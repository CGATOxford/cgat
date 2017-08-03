"""
gtf2fasta.py - annotate genomic bases from a gene set
=====================================================

:Tags: Genomics Genesets Sequences GTF FASTA Transformation


Purpose
-------
This script can be used for a quick-and-dirty annotation of variants
in a genome. It is most appropriately used in exploratory analyses
of the effect of variants/alleles.

For a better prediction of variant effects in coding sequences,
see :doc:`snp2counts` and :doc:`gtf2alleles`.

If you wish to convert gtf intervals into fasta sequences, use gff2fasta.py.

This script takes a :term:`gtf` formatted file from ENSEMBL and
annotates each base in the genome according to its *function*. The
script multiplexes both strands with lower- case characters referring
to the forward strand and upper-case characters referring to the
reverse strand.

The codes and their meaning are:

+-----+----------------------------------------------------------------------+
|code | description                                                          |
+-----+----------------------------------------------------------------------+
|a    | first codon position within a complete codon                         |
+-----+----------------------------------------------------------------------+
|b    | second codon position within a complete codon                        |
+-----+----------------------------------------------------------------------+
|c    | third codon position within a complete codon                         |
+-----+----------------------------------------------------------------------+
|d    | coding base, but in multiple frames or strands                       |
+-----+----------------------------------------------------------------------+
|e    | non-coding base in exon                                              |
+-----+----------------------------------------------------------------------+
|f    | frame-shifted base                                                   |
+-----+----------------------------------------------------------------------+
|g    | intergenic base                                                      |
+-----+----------------------------------------------------------------------+
|i    | intronic base                                                        |
+-----+----------------------------------------------------------------------+
|l    | base in other RNA                                                    |
+-----+----------------------------------------------------------------------+
|m    | base in miRNA                                                        |
+-----+----------------------------------------------------------------------+
|n    | base in snRNA                                                        |
+-----+----------------------------------------------------------------------+
|o    | base in snoRNA                                                       |
+-----+----------------------------------------------------------------------+
|r    | base in rRNA (both genomic and mitochondrial)                        |
+-----+----------------------------------------------------------------------+
|p    | base in pseudogene (including transcribed, unprocessed and processed)|
+-----+----------------------------------------------------------------------+
|q    | base in retrotransposon                                              |
+-----+----------------------------------------------------------------------+
|s    | base within a splice signal (GT/AG)                                  |
+-----+----------------------------------------------------------------------+
|t    | base in tRNA (both genomic and mitochondrial)                        |
+-----+----------------------------------------------------------------------+
|u    | base in 5' UTR                                                       |
+-----+----------------------------------------------------------------------+
|v    | base in 3' UTR                                                       |
+-----+----------------------------------------------------------------------+
|x    | ambiguous base with multiple functions.                              |
+-----+----------------------------------------------------------------------+
|y    | unknown base                                                         |
+-----+----------------------------------------------------------------------+



Output files
++++++++++++

The annotated genome is output on stdout.

The script creates the following additional output files:

counts
   Counts for each annotations

junctions
   Splice junctions. This is a tab separated table linking residues that are
   joined via features. The coordinates are forward/reverse coordinates.

   The columns are:

   contig
      the contig
   strand
      direction of linkage
   end
      last base of exon in direction of strand
   start
      first base of exon in direction of strand
   frame
      frame base at second coordinate (for coding sequences)
    
Known problems
--------------

The stop-codon is part of the UTR. This has the following effects:

   * On the mitochondrial chromosome, the stop-codon might be used for
     ncRNA transcripts and thus the base is recorded as ambiguous.

   * On the mitochondrial chromosome, alternative transcripts might
     read through a stop-codon (RNA editing). The codon itself will be
     recorded as ambiguous.

Usage
-----

For example::

   zcat hg19.gtf.gz | python gtf2fasta.py --genome-file=hg19 > hg19.annotated

Type::

   python gtf2fasta.py --help

for command line help.

Command line options
--------------------

``--genome-file``
    required option. filename for genome fasta file

``--ignore-missing``
    transcripts on contigs not in the genome file will be ignored

``--min-intron-length``
    intronic bases in introns less than specified length
    will be marked "unknown"

"""

import sys
import collections
import array

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Genomics as Genomics
import CGAT.Intervals as Intervals

MAP_ENSEMBL = {'miRNA': 'm',
               'misc_RNA': 'l',
               'pseudogene': 'p',
               'transcribed_pseudogene': 'p',
               'unprocessed_pseudogene': 'p',
               'processed_pseudogene': 'p',
               'retrotransposed': 'q',
               'rRNA': 'r',
               'snRNA': 'n',
               'snoRNA': 'o',
               'Mt_rRNA': 'r',
               'Mt_tRNA': 't'}

DEFAULT_CODE = "g"
AMBIGUOUS_CODE = "x"
CODING_CODE = "d"
ALL_CODES = DEFAULT_CODE + AMBIGUOUS_CODE + \
    "abcdefilmnorpqstuvy" + "abcdefilmnorpqstuvy".upper()
CODING_CODES = "abcdABCD"
NONCODING_CODES = "iuvIUV"
INTRON_CODES = "iI"
UTR_CODES = "uvUV"


def setCode(annotation, pos, code):
    """set *pos* to *code* in annotation.

    This method performs conflict resolution in the following cases:

    1. Introns and UTRs do not cause a conflict if they overlap with other
    features

    1. Coding bases take precedence over intronic and UTR bases.
    2. UTR base take precedence over intronic bases.
    3. Coding bases in different frames/strands are set to CODING_CODE
    4. Non-coding features of the same type but different strand are permitted

    All other conflicts are marked as ambiguous bases.
    """
    c = annotation[pos]
    if c == DEFAULT_CODE or c in NONCODING_CODES:
        annotation[pos] = code
    elif c == AMBIGUOUS_CODE:
        return
    elif code in NONCODING_CODES:
        # only set introns/UTR if no other code is present
        return
    elif c == code:
        return
    elif code in CODING_CODES and c in CODING_CODES:
        # ambiguous frame/strand in coding sequence
        annotation[pos] = CODING_CODE
    elif code in NONCODING_CODES and c in CODING_CODES:
        # permit alternative transcripts
        return
    elif c not in CODING_CODES and code not in CODING_CODES and c.upper() == code.upper():
        # permit features of the same type on different strands to overlap (for
        # example, tRNAs)
        return
    else:
        annotation[pos] = AMBIGUOUS_CODE
        E.warn("ambiguous position %i: %s - %s" % (pos, c, code))


def addSegments(annotation, intervals, is_positive, code):
    """add intervals."""
    if not intervals:
        return
    if not is_positive:
        code = code.upper()

    for start, end in intervals:
        for x in range(start, end):
            setCode(annotation, x, code)


def addIntrons(annotation, intervals, is_positive, max_frameshift_length):
    """add introns for intervals.

    Intervals need to be sorted in incremental order.
    """
    if not intervals:
        return
    intervals.sort()
    last = intervals[0][1]
    code_i, code_f, code_s = "i", "f", "s"

    if not is_positive:
        code_i = code_i.upper()
        code_f = code_f.upper()
        code_s = code_s.upper()

    for start, end in intervals[1:]:
        d = start - last
        if d < max_frameshift_length:
            code = code_f
        else:
            code = code_i
            # add splice sites
            setCode(annotation, last, code_s)
            setCode(annotation, last + 1, code_s)
            setCode(annotation, start - 2, code_s)
            setCode(annotation, start - 1, code_s)
            last += 2
            start -= 2

        for x in range(last, start):
            setCode(annotation, x, code)

        last = end


def addCDS(annotation, gtfs, is_positive):
    """add coding sequence to contig.

    Also adds the splice sites.
    """

    if not gtfs:
        return

    if is_positive:
        chars = "abc"
    else:
        chars = "ABC"

    last = None
    for cds in gtfs:

        c = int(cds.frame)
        if c != 0:
            c = 3 - c

        if is_positive:
            r = range(cds.start, cds.end)
        else:
            r = range(cds.end - 1, cds.start - 1, -1)

        for x in r:
            code = chars[c]
            c += 1
            if c == len(chars):
                c = 0
            setCode(annotation, x, code)


def outputCounts(outfile, annotations):
    """output table into outfile with annotations."""

    total_counts = collections.defaultdict(int)

    outfile.write("contig\ttotal\t%s\n" % ("\t".join(ALL_CODES)))

    for k in sorted(annotations.keys()):
        counts = collections.defaultdict(int)
        for x in annotations[k]:
            counts[x] += 1
        outfile.write("\t".join((k,
                                 str(len(annotations[k])),
                                 "\t".join([str(counts[x]) for x in ALL_CODES]))) + "\n")

        for k, v in counts.items():
            total_counts[k] += v

    outfile.write("\t".join(("total",
                             str(sum([len(x) for x in list(annotations.values())])),
                             "\t".join([str(total_counts[x]) for x in ALL_CODES]))) + "\n")


def annotateGenome(iterator, fasta, options, default_code=DEFAULT_CODE):
    """annotate a genome given by the indexed *fasta* file and 
    an iterator over gtf annotations.
    """

    annotations = {}
    contig_sizes = fasta.getContigSizes(with_synonyms=False)
    E.info("allocating memory for %i contigs and %i bytes" %
           (len(contig_sizes), sum(contig_sizes.values()) * array.array("B").itemsize))
    # AString.AString( "a").itemsize ))

    for contig, size in list(contig_sizes.items()):
        E.debug("allocating %s: %i bases" % (contig, size))
        # annotations[contig] = AString.AString( default_code * size )
        # annotations[contig] = array.array("", default_code * size)
        # Go to list for py3 compatibility, patch
        annotations[contig] = [default_code] * size

    E.info("allocated memory for %i contigs" % len(fasta))

    counter = E.Counter()

    # output splice junctions
    outfile_junctions = E.openOutputFile("junctions")
    outfile_junctions.write(
        "contig\tstrand\tpos1\tpos2\tframe\tgene_id\ttranscript_id\n")
    for gtfs in iterator:

        counter.input += 1

        if counter.input % options.report_step == 0:
            E.info("iteration %i" % counter.input)

        try:
            contig = fasta.getToken(gtfs[0].contig)
        except KeyError as msg:
            E.warn("contig %s not found - annotation ignored" % gtfs[0].contig)
            counter.skipped_contig += 1
            continue

        lcontig = fasta.getLength(contig)

        # make sure that exons are sorted by coordinate
        gtfs.sort(key=lambda x: x.start)

        is_positive = Genomics.IsPositiveStrand(gtfs[0].strand)
        source = gtfs[0].source

        # process non-coding data
        if source in MAP_ENSEMBL:
            code = MAP_ENSEMBL[source]

            intervals = [(x.start, x.end) for x in gtfs]
            addSegments(annotations[contig],
                        intervals,
                        is_positive,
                        code)

        elif source == "protein_coding":

            # collect exons for utr
            exons = [(x.start, x.end) for x in gtfs if x.feature == "exon"]
            cds = [(x.start, x.end) for x in gtfs if x.feature == "CDS"]
            if len(cds) == 0:
                counter.skipped_transcripts += 1
                E.warn("protein-coding transcript %s without CDS - skipped" %
                       gtfs[0].transcript_id)
                continue

            exons = Intervals.truncate(exons, cds)
            start, end = cds[0][0], cds[-1][1]

            UTR5 = [x for x in exons if x[1] < start]
            UTR3 = [x for x in exons if x[0] >= end]

            if not is_positive:
                UTR5, UTR3 = UTR3, UTR5
                splice_code = "S"
            else:
                splice_code = "s"

            addSegments(annotations[contig],
                        UTR5,
                        is_positive,
                        "u")

            addIntrons(annotations[contig],
                       UTR5,
                       is_positive,
                       options.max_frameshift_length)

            addSegments(annotations[contig],
                        UTR3,
                        is_positive,
                        "v")

            addIntrons(annotations[contig],
                       UTR3,
                       is_positive,
                       options.max_frameshift_length)

            # output CDS according to frame
            addCDS(annotations[contig],
                   [x for x in gtfs if x.feature == "CDS"],
                   is_positive)

            # add introns between CDS
            addIntrons(annotations[contig],
                       cds,
                       is_positive,
                       options.max_frameshift_length)

            # output splice junctions
            cds = [x for x in gtfs if x.feature == "CDS"]

            # apply corrections for 1-past end coordinates
            # to point between residues within CDS
            if is_positive:
                ender = lambda x: x.end - 1
                starter = lambda x: x.start
                out_positive = "+"
            else:
                ender = lambda x: lcontig - x.start - 1
                starter = lambda x: lcontig - x.end
                out_positive = "-"
                cds.reverse()

            end = ender(cds[0])
            for c in cds[1:]:
                start = starter(c)
                outfile_junctions.write("%s\t%s\t%i\t%i\t%s\t%s\t%s\n" %
                                        (contig,
                                         out_positive,
                                         end,
                                         start,
                                         c.frame,
                                         c.gene_id,
                                         c.transcript_id,
                                         ))
                end = ender(c)

    E.info("finished reading genes: %s" % str(counter))

    outfile_junctions.close()

    E.info("started counting")
    outfile = E.openOutputFile("counts")
    outputCounts(outfile, annotations)
    outfile.close()

    E.info("started output")
    for k in sorted(annotations.keys()):
        # options.stdout.write(">%s\n%s\n" % (k, annotations[k].tostring()))
        options.stdout.write(">%s\n%s\n" % (k, "".join(annotations[k])))


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id: gtf2fasta.py 2861 2010-02-23 17:36:32Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default].")

    parser.add_option("-i", "--ignore-missing", dest="ignore_missing", action="store_true",
                      help="Ignore transcripts on contigs that are not in the genome-file [default=%default].")

    parser.add_option("--min-intron-length", dest="min_intron_length", type="int",
                      help="minimum intron length. If the distance between two consecutive exons is smaller, the region will be marked 'unknown' [default=%default].")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("full", ),
                      help="method to apply [default=%default].")

    parser.set_defaults(
        genome_file=None,
        flank=1000,
        max_frameshift_length=4,
        min_intron_length=30,
        ignore_missing=False,
        restrict_source=None,
        method="full",
        report_step=1000,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    if not options.genome_file:
        raise ValueError("an indexed genome is required.")

    fasta = IndexedFasta.IndexedFasta(options.genome_file)

    iterator = GTF.transcript_iterator(GTF.iterator(options.stdin))

    annotateGenome(iterator, fasta, options)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

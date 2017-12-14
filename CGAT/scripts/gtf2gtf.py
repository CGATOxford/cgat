'''gtf2gtf.py - manipulate transcript models
============================================

:Tags: Genomics Genesets GTF Manipulation

Purpose
-------

This script reads a gene set in :term:`gtf` format from stdin, applies some
transformation, and outputs a new gene set in :term:`gtf` format to stdout.
The transformation is chosen by the ``--method`` command line option.

Transformations available for use in this script can broadly be
classified into four categories:

1. sorting gene sets
2. manipulating gene models
3. filtering gene sets
4. setting/resetting fields within a gtf file

Further options for working with gtf files are available in gff2gff.py,
which can be run with the specification --is-gtf


Sorting gene sets
+++++++++++++++++

``sort``

   Sorts entries in gtf file by one or more fields

   +-----------------+---------------------------------------+
   | option          | order in which fields are sorted      |
   +-----------------|---------------------------------------+
   | gene            | gene_id, contig, start                |
   +-----------------+---------------------------------------+
   | gene+transcript | gene_id, transcript_id, contig, start |
   +-----------------+---------------------------------------+
   | contig+gene     | contig, gene_id, transcript_id, start |
   +-----------------+---------------------------------------+
   | transcript      | transcript_id, contig, start          |
   +-----------------+---------------------------------------+
   | position        | contig, start                         |
   +-----------------+---------------------------------------+
   | position+gene   | contig( gene_id, start )              |
   +-----------------+---------------------------------------+
   | gene+position   | gene_id, contig, start                |
   +-----------------+---------------------------------------+
   | gene+exon       | gene_id, exon_id                      |
   +-----------------+---------------------------------------+

   N.B. position+gene sorts by gene_id, start, then subsequently sorts
   flattened gene lists by contig, start


Manipulating gene-models
++++++++++++++++++++++++

Options that can be used to alter the features represented in a :term:`gtf`
file. Only one method can be specified at once.

Input gtfs need to be sorted so that features for a gene or transcript
appear consecutively within the file. This can be achevied using
``--method=sort``.


``genes-to-unique-chunks```
    Divide the complete length of a gene up into chunks that represent
    ranges of bases that are all present in the same set of transcripts.
    E.g. for two overlapping exons an entry will be output representing
    the overlap and a seperate entry each for the sequences only present
    in one. Ranges which are between the first TSS and last TTS but not
    present in any transcript (i.e. merged introns) are also output.
    Useful for DEXSeq like splicing analysis

``find-retained-introns``
    Finds intervals within a transcript that represent retained-introns,
    here a retained intron is considered to be an intron in one transcript
    that is entirely contianed within the exon of another. The retained
    intron will be assigned to the transcript with the containing exon. Where
    multiple, overlapping introns are contained within a single exon of a
    transcript, the union of the introns will be output. Thus when considering
    an indevidual transcript, outputs will be non-overlapping. However,
    overlapping, or even identical feature can be output if they belong to
    different transcripts.

``merge-exons``
    Merges overlapping exons for all transcripts of a gene, outputting
    the merged exons. Can be used in conjunction with
    ``merge-exons-distance`` to set the minimum distance that may
    appear between two exons before they are merged.If
    ``--mark-utr`` is set, the UTR regions will be output separately.

``merge-transcripts``
    Merges all transcripts of a gene. Outputs contains a single interval that
    spans the original gene (both introns and exons). If ``--with-utr`` is
    set, the output interval will also contain UTR.

``merge-genes``

    Merges genes that have overlapping exons, outputting a single
    gene_id and transcript_id for all exons of overlapping genes. The
    input needs te sorted by transcript " (Does not merge intervals on
    different strands).

``join-exons``
    Joins together all exons of a transcript, outputting a single
    interval that spans the original transcript (both introns and
    exons). Input needs to be sorted by transcript.

``intersect-transcripts``
    Finds regions representing the intersection of all transcripts of a gene.
    Output will contain intervals spanning only those bases covered by all
    transcripts. If ``--with-utr`` is set, the UTR will also be included in the
    intersect. This method only uses ``exon`` or ``CDS`` features.

``merge-introns``
    Outputs a single interval that spans the region between the start
    of the first intron and the end of last intron. Single exons genes
    will not be output. The input needs to be sorted by gene

``exons2introns``
    Merges overlapping introns for all transcripts of a gene,
    outputting the merged introns. Use ``--intron-min-length`` to
    ignore merged introns below a specified length. Use
    ``--intron-border`` to specify a number of residues to remove at
    either end of output introns (residues are removed prior to
    filtering on size when used in conjunction with
    ``--intron-min-length``).

``transcripts2genes``
    Cluster transcripts into genes by exon overlap ignoring any
    gene_ids in the :term:`gtf` file. May be used in conjunction with
    ``reset-strand``

The option ``permit-duplicates`` may be specified in order to
allow gene-ids to be duplicated within the input :term:`gtf` file
(i.e. for the same gene-id to appear non-consecutively within the
input file). However, this option currently only works for
``merge-exons``, ``merge-transcripts``, ``merge-introns``, and
``intersect-transcripts``. It DOES NOT work for ``merge-genes``,
``join-exons``, or ``exons-file2introns``.

Filtering gene sets
+++++++++++++++++++

Options that can be used to filter :term:`gtf` files. For further
detail see command line options.

Input gtfs need to be sorted so that features for a gene or transcript
appear consecutively within the file. This can be achevied using
``--method=sort --sort-order``.

``filter``
    When filtering on the basis of 'gene-id' or 'transcript-id' a
    filename containing ids to be removed may provided using
    ``--map-tsv-file``. Alternatively, a random subsample of
    genes/transcripts may be retained using
    ``--sam-fileple-size``. Use ``--min-exons-length`` in conjunction
    with ``--sam-fileple-size`` to specify a minimum length for
    genes/transcripts to be retained. Use ``--ignore-strand`` to set
    strand to '.' in output.

    Other filter options include longest-gene, longest-transcript,
    or representative-transcript.

    When filtering on the basis of gene-id, transcript-id or longest-gene,
    ``--invert-filter`` may be used to invert the selection.

``remove-overlapping``
    Given a second :term:`gff` formatted file (``--file-gff``) removes
    any features overlapping. Any transcripts that intersect intervals
    in the supplied file are removed.  (Does not account for strand.)

``remove-duplicates``
    Remove duplicate features from :term:`gtf` file. The type of
    feature to be removed is set by the option ``-duplicate-feature``.
    Setting ``--duplicate-feature`` to 'gene', 'transcript', or
    'coordinates' will remove any interval for which non-consecutive
    occurrances of specified term appear in input :term:`gtf` file.
    Setting to 'ucsc', will remove any interval for which
    transcript-id contains '_dup'.


Setting fields
++++++++++++++

Options for altering fields within :term:`gtf`.

``rename-genes``
    With a mapping file is provided using ``--map-tsv-file``, renames
    the gene_id to the one supplied. Outputs a :term:`gtf` file with
    field renamed. Any entry in input :term:`gtf` not appearing in
    mapping file is discarded.

``rename-transcripts``
    as ``rename-genes``, but renames the transcript_id.

``add-protein-id``
    Takes a map of transcript_id to protein_id from the a tsv file
    (see option ``--map-tsv-file``) and appends the protein_id
    provided to the attributes field.  Any entry with a transcript_id
    not appearing in the tsv file is discarded.

``renumber-genes``
    Renumber genes from 1 using the pattern provided in
    ``--pattern-identifier``.

``renumber-transcripts``
    Renumber transcripts from 1 using the pattern provided in
    ``--pattern-identifier``.

``unset-genes``
    Renumber genes from 1 using the pattern provided in
    ``--pattern-identifier``. Transcripts with the same gene-id in the
    input :term:`gtf` file will have different gene-ids in the output
    :term:`gtf` file.

``set-transcript-to-gene``
    Will set the transcript-id to the gene-id for each feature.

``set-gene-to-transcript``
    Will set the gene-id to the transcript-id for each each feature.

``set-protein-to-transcript``
    Will append transcript_id to attributes field as 'protein_id'

``set-score-to-distance``
    Will reset the score field (field 6) of each feature in input
    :term:`gtf` to be the distance from transcription start site to
    the start of the feature.  (Assumes input file is sorted by
    transcript-id)

``set-gene_biotype-to-source``
    Sets the ``gene_biotype`` attribute from the source column. Will only set
    if biotype attribute is not present in the current record.

``rename-duplicates``
    Rename duplicate gene_ids and transcript_ids by addition of
    numerical suffix

``set-source-to-transcript_biotype``
    Sets the source attribute to the ``transcript_biotype``
    attribute. Will only set if ``transcript_biotype`` attribute is
    present in the current record.

Usage
-----

The following example sorts the input gene set by gene
(``method=sort``) so that it can be used as input for
``method=intersect-transcripts`` that outputs genomic the genomic
regions within a gene that is covered by all transcripts in a gene.
Finally, the resultant transcripts are renamed with the pattern
"MERGED_%i"::

    cgat gtf2gtf
            --method=sort
            --sort-order=gene \
    | cgat gtf2gtf
               --method=intersect-transcripts
               --with-utr
    | cgat gtf2gtf
               --method=renumber-transcripts
               --pattern-identifier=MERGED_%i

Type::

    cgat gtf2gtf --help

for command line options.

Command line Options
--------------------

'''
import sys
import re
import random
import collections
import itertools

import CGAT.GTF as GTF
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.Intervals as Intervals
import CGAT.IOTools as IOTools


def find_retained_introns(gene):
    '''Given a bundle of transcripts, find intervals matching retained
    introns. A retained intron is defined as an interval from an exon/intron
    boundary to the next where both boundaries are in the same exon of another
    transcript'''

    intron_intervals = [GTF.toIntronIntervals(transcript)
                        for transcript in gene]
    intron_intervals = list(set(
        itertools.chain.from_iterable(intron_intervals)))
    intron_intervals.sort()

    for transcript in gene:
        exons = iter(sorted(GTF.asRanges(transcript)))
        introns = iter(intron_intervals)
        retained_introns = []
        try:
            intron = next(introns)
            exon = next(exons)
            while True:

                if exon[1] < intron[0]:

                    exon = next(exons)
                    continue

                if intron[0] >= exon[0] and intron[1] <= exon[1]:
                    E.debug("exon %s of transcript %s contains intron %s" %
                            (exon, transcript[0].transcript_id, intron))
                    retained_introns.append(intron)
                intron = next(introns)
        except StopIteration:
            pass

        retained_introns = Intervals.combine(retained_introns)

        for intron in retained_introns:
            entry = GTF.Entry()
            entry = entry.copy(transcript[0])
            entry.start = intron[0]
            entry.end = intron[1]
            yield entry


def gene_to_blocks(gene):
    '''Given a bundle of all exons in a gene, create a seperate exon
    for each unqiue part of a exon, as well as one for introns. '''

    exons = [(e.start, e.end)
             for e in gene if e.feature == "exon"]

    exons = list(set(sum(exons, ())))
    exons.sort()

    entry = GTF.Entry()
    entry = entry.copy(gene[0])
    entry.transcript_id = "merged"
    entry.feature = "exon"
    entry.source = "merged"

    for i in range(len(exons)-1):
        entry.start = exons[i]
        entry.end = exons[i+1]
        entry.attributes["exon_id"] = str(i + 1)
        yield entry


def main(argv=None):

    if not argv:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--merge-exons-distance",
                      dest="merge_exons_distance",
                      type="int",
                      help="distance in nucleotides between "
                      "exons to be merged [%default].")

    parser.add_option("--pattern-identifier", dest="pattern", type="string",
                      help="pattern to use for renaming genes/transcripts. "
                      "The pattern should contain a %i, for example "
                      "--pattern-identifier=ENSG%010i [%default].")

    parser.add_option("--sort-order",
                      dest="sort_order",
                      type="choice",
                      choices=("gene",
                               "gene+transcript",
                               "transcript",
                               "position",
                               "contig+gene",
                               "position+gene",
                               "gene+position",
                               "gene+exon"),
                      help="sort input data [%default].")

    parser.add_option("--mark-utr",
                      dest="mark_utr",
                      action="store_true",
                      help="mark utr for method --merge-exons. "
                      "[%default].")

    parser.add_option(
        "--without-utr",
        dest="with_utr",
        action="store_false",
        help="exclude UTR in methods --merge-exons, merge-transcripts "
        "and intersect-transripts. Setting this option will remove "
        "non-coding transcripts. "
        "[%default].")

    parser.add_option(
        "--filter-method", dest="filter_method",
        type="choice",
        choices=("gene",
                 "transcript",
                 "longest-gene",
                 "longest-transcript",
                 "representative-transcript",
                 "proteincoding",
                 "lincrna"),
        help="Filter method to apply. Available filters are: "
        "'gene': filter by gene_id given in ``--map-tsv-file``, "
        "'transcript': filter by transcript_id given in ``--map-tsv-file``, "
        "'longest-gene': output the longest gene for overlapping genes ,"
        "'longest-transcript': output the longest transcript per gene,"
        "'representative-transcript': output the representative transcript "
        "per gene. The representative transcript is the transcript "
        "that shares most exons with other transcripts in a gene. "
        "The input needs to be sorted by gene. "
        "'proteincoding': only output protein coding features. "
        "'lincrna': only output lincRNA features. "
        "[%default].")

    parser.add_option("-a", "--map-tsv-file", dest="filename_filter",
                      type="string",
                      metavar="tsv",
                      help="filename of ids to map/filter [%default].")

    parser.add_option(
        "--gff-file", dest="filename_gff", type="string",
        metavar="GFF",
        help="second filename of features (see --remove-overlapping) "
        "[%default]")

    parser.add_option("--invert-filter",
                      dest="invert_filter",
                      action="store_true",
                      help="when using --filter, invert selection "
                      "(like grep -v). "
                      "[%default].")

    parser.add_option("--sample-size", dest="sample_size", type="int",
                      help="extract a random sample of size # if the option "
                      "'--method=filter --filter-method' is set "
                      "[%default].")

    parser.add_option(
        "--intron-min-length",
        dest="intron_min_length", type="int",
        help="minimum length for introns (for --exons-file2introns) "
        "[%default].")

    parser.add_option("--min-exons-length",
                      dest="min_exons_length",
                      type="int",
                      help="minimum length for gene (sum of exons) "
                      "(--sam-fileple-size) [%default].")

    parser.add_option(
        "--intron-border",
        dest="intron_border",
        type="int",
        help="number of residues to exclude at intron at either end "
        "(--exons-file2introns) [%default].")

    parser.add_option("--ignore-strand",
                      dest="ignore_strand",
                      action="store_true",
                      help="remove strandedness of features (set to '.') when "
                      "using ``transcripts2genes`` or ``filter``"
                      "[%default].")

    parser.add_option("--permit-duplicates", dest="strict",
                      action="store_false",
                      help="permit duplicate genes. "
                      "[%default]")

    parser.add_option(
        "--duplicate-feature",
        dest="duplicate_feature",
        type="choice",
        choices=("gene", "transcript", "both", "ucsc", "coordinates"),
        help="remove duplicates by gene/transcript. "
        "If ``ucsc`` is chosen, transcripts ending on _dup# are "
        "removed. This is necessary to remove duplicate entries "
        "that are next to each other in the sort order "
        "[%default]")

    parser.add_option("--use-gene-id", dest="use_geneid", action="store_true",
                      help="when merging transcripts, exons or introns, use "
                      "the parent gene_id as the transcript id.")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      action="append",
                      choices=(
                          "add-protein-id",
                          "exons2introns",
                          "filter",
                          "find-retained-introns",
                          "genes-to-unique-chunks",
                          "intersect-transcripts",
                          "join-exons",
                          "merge-exons",
                          "merge-transcripts",
                          "merge-genes",
                          "merge-introns",
                          "remove-overlapping",
                          "remove-duplicates",
                          "rename-genes",
                          "rename-transcripts",
                          "rename-duplicates",
                          "renumber-genes",
                          "renumber-transcripts",
                          "set-transcript-to-gene",
                          "set-gene-to-transcript",
                          "set-protein-to-transcript",
                          "set-score-to-distance",
                          "set-gene_biotype-to-source",
                          "set-source-to-transcript_biotype",
                          "sort",
                          "transcript2genes",
                          "unset-genes"),
                      help="Method to apply [%default]."
                      "Please only select one.")

    parser.set_defaults(
        sort_order="gene",
        filter_method="gene",
        pattern="%i",
        merge_exons_distance=0,
        filename_filter=None,
        intron_border=None,
        intron_min_length=None,
        sample_size=0,
        min_exons_length=0,
        ignore_strand=False,
        mark_utr=False,
        with_utr=True,
        invert_filter=False,
        duplicate_feature=None,
        strict=True,
        method=None,
        use_geneid=False,
    )

    (options, args) = E.Start(parser, argv=argv)

    ninput, noutput, nfeatures, ndiscarded = 0, 0, 0, 0

    if options.method is None:
        raise ValueError("please specify a --method")

    if len(options.method) > 1:
        raise ValueError("multiple --method arguements specified")
    else:
        options.method = options.method[0]

    if options.method == "set-transcript-to-gene":

        for gff in GTF.iterator(options.stdin):

            ninput += 1

            gff.transcript_id = gff.gene_id
            options.stdout.write("%s\n" % str(gff))

            noutput += 1
            nfeatures += 1

    elif options.method == "set-gene_biotype-to-source":

        for gff in GTF.iterator(options.stdin):

            ninput += 1

            if "gene_biotype" not in gff.attributes:
                gff.gene_biotype = gff.source

            options.stdout.write("%s\n" % str(gff))

            noutput += 1
            nfeatures += 1

    elif options.method == "set-source-to-transcript_biotype":

        for gff in GTF.iterator(options.stdin):

            ninput += 1

            try:
                gff.source = gff.transcript_biotype
            except AttributeError:
                pass

            options.stdout.write("%s\n" % str(gff))

            noutput += 1
            nfeatures += 1

    elif options.method == "remove-duplicates":

        counts = collections.defaultdict(int)

        if options.duplicate_feature == "ucsc":
            store = []
            remove = set()
            f = lambda x: x[0].transcript_id

            gffs = GTF.transcript_iterator(
                GTF.iterator(options.stdin), strict=False)
            outf = lambda x: "\n".join([str(y) for y in x])

            for entry in gffs:
                ninput += 1
                store.append(entry)
                id = f(entry)
                if "_dup" in id:
                    remove.add(re.sub("_dup\d+", "", id))
                    remove.add(id)

            for entry in store:
                id = f(entry)
                if id not in remove:
                    options.stdout.write(outf(entry) + "\n")
                    noutput += 1
                else:
                    ndiscarded += 1
                    E.info("discarded duplicates for %s" % (id))
        else:

            if options.duplicate_feature == "gene":
                gffs = GTF.gene_iterator(
                    GTF.iterator(options.stdin), strict=False)
                f = lambda x: x[0][0].gene_id
                outf = lambda x: "\n".join(
                    ["\n".join([str(y) for y in xx]) for xx in x])
            elif options.duplicate_feature == "transcript":
                gffs = GTF.transcript_iterator(
                    GTF.iterator(options.stdin), strict=False)
                f = lambda x: x[0].transcript_id
                outf = lambda x: "\n".join([str(y) for y in x])
            elif options.duplicate_feature == "coordinates":
                gffs = GTF.chunk_iterator(GTF.iterator(options.stdin))
                f = lambda x: x[0].contig + "_" + \
                    str(x[0].start) + "-" + str(x[0].end)
                outf = lambda x: "\n".join([str(y) for y in x])

            store = []

            for entry in gffs:
                ninput += 1
                store.append(entry)
                id = f(entry)
                counts[id] += 1

            # Assumes GTF file sorted by contig then start
            last_id = ""
            if options.duplicate_feature == "coordinates":
                for entry in store:
                    id = f(entry)
                    if id == last_id:
                        ndiscarded += 1
                        E.info("discarded duplicates for %s: %i" %
                               (id, counts[id]))
                    else:
                        options.stdout.write(outf(entry) + "\n")
                        noutput += 1
                    last_id = id

            else:
                for entry in store:
                    id = f(entry)
                    if counts[id] == 1:
                        options.stdout.write(outf(entry) + "\n")
                        noutput += 1
                    else:
                        ndiscarded += 1
                        E.info("discarded duplicates for %s: %i" %
                               (id, counts[id]))

    elif "sort" == options.method:

        for gff in GTF.iterator_sorted(GTF.iterator(options.stdin),
                                       sort_order=options.sort_order):
            ninput += 1
            options.stdout.write("%s\n" % str(gff))
            noutput += 1
            nfeatures += 1

    elif "set-gene-to-transcript" == options.method:

        for gff in GTF.iterator(options.stdin):

            ninput += 1

            gff.gene_id = gff.transcript_id
            options.stdout.write("%s\n" % str(gff))

            noutput += 1
            nfeatures += 1

    elif "set-protein-to-transcript" == options.method:

        for gff in GTF.iterator(options.stdin):
            ninput += 1
            gff.protein_id = gff.transcript_id
            options.stdout.write("%s\n" % str(gff))
            noutput += 1
            nfeatures += 1

    elif "add-protein-id" == options.method:

        transcript2protein = IOTools.readMap(
            IOTools.openFile(options.filename_filter, "r"))

        missing = set()
        for gff in GTF.iterator(options.stdin):
            ninput += 1
            if gff.transcript_id not in transcript2protein:
                if gff.transcript_id not in missing:
                    E.debug(
                        ("removing transcript '%s' due to "
                         "missing protein id") % gff.transcript_id)
                    missing.add(gff.transcript_id)
                ndiscarded += 1
                continue

            gff.protein_id = transcript2protein[gff.transcript_id]
            options.stdout.write("%s\n" % str(gff))
            noutput += 1
            nfeatures += 1

        E.info("transcripts removed due to missing protein ids: %i" %
               len(missing))

    elif "join-exons" == options.method:

        for exons in GTF.transcript_iterator(GTF.iterator(options.stdin)):
            ninput += 1
            strand = Genomics.convertStrand(exons[0].strand)
            contig = exons[0].contig
            transid = exons[0].transcript_id
            geneid = exons[0].gene_id
            biotype = exons[0].source
            all_start, all_end = min([x.start for x in exons]), max(
                [x.end for x in exons])
            y = GTF.Entry()
            y.contig = contig
            y.source = biotype
            y.feature = "transcript"
            y.start = all_start
            y.end = all_end
            y.strand = strand
            y.transcript_id = transid
            y.gene_id = geneid
            options.stdout.write("%s\n" % str(y))

    elif "merge-genes" == options.method:
        # merges overlapping genes
        #
        gffs = GTF.iterator_sorted_chunks(
            GTF.flat_gene_iterator(GTF.iterator(options.stdin)),
            sort_by="contig-strand-start")

        def iterate_chunks(gff_chunks):

            last = next(gff_chunks)
            to_join = [last]

            for gffs in gff_chunks:
                d = gffs[0].start - last[-1].end

                if gffs[0].contig == last[0].contig and \
                   gffs[0].strand == last[0].strand:
                    assert gffs[0].start >= last[0].start, \
                        ("input file should be sorted by contig, strand "
                         "and position: d=%i:\nlast=\n%s\nthis=\n%s\n") % \
                        (d,
                         "\n".join([str(x) for x in last]),
                         "\n".join([str(x) for x in gffs]))

                if gffs[0].contig != last[0].contig or \
                        gffs[0].strand != last[0].strand or \
                        d > 0:
                    yield to_join
                    to_join = []

                last = gffs
                to_join.append(gffs)

            yield to_join
            raise StopIteration

        for chunks in iterate_chunks(gffs):
            ninput += 1
            if len(chunks) > 1:
                gene_id = "merged_%s" % chunks[0][0].gene_id
                transcript_id = "merged_%s" % chunks[0][0].transcript_id
                info = ",".join([x[0].gene_id for x in chunks])
            else:
                gene_id = chunks[0][0].gene_id
                transcript_id = chunks[0][0].transcript_id
                info = None

            intervals = []
            for c in chunks:
                intervals += [(x.start, x.end) for x in c]

            intervals = Intervals.combine(intervals)
            # take single strand
            strand = chunks[0][0].strand

            for start, end in intervals:
                y = GTF.Entry()
                y.fromGTF(chunks[0][0], gene_id, transcript_id)
                y.start = start
                y.end = end
                y.strand = strand

                if info:
                    y.addAttribute("merged", info)
                options.stdout.write("%s\n" % str(y))
                nfeatures += 1

            noutput += 1

    elif options.method == "renumber-genes":

        map_old2new = {}
        for gtf in GTF.iterator(options.stdin):
            ninput += 1
            if gtf.gene_id not in map_old2new:
                map_old2new[gtf.gene_id] = options.pattern % (
                    len(map_old2new) + 1)
            gtf.gene_id = map_old2new[gtf.gene_id]
            options.stdout.write("%s\n" % str(gtf))
            noutput += 1

    elif options.method == "unset-genes":

        map_old2new = {}
        for gtf in GTF.iterator(options.stdin):
            ninput += 1
            key = gtf.transcript_id
            if key not in map_old2new:
                map_old2new[key] = options.pattern % (len(map_old2new) + 1)
            gtf.gene_id = map_old2new[key]
            options.stdout.write("%s\n" % str(gtf))
            noutput += 1

    elif options.method == "renumber-transcripts":

        map_old2new = {}
        for gtf in GTF.iterator(options.stdin):
            ninput += 1
            key = (gtf.gene_id, gtf.transcript_id)
            if key not in map_old2new:
                map_old2new[key] = options.pattern % (
                    len(map_old2new) + 1)
            gtf.transcript_id = map_old2new[key]
            options.stdout.write("%s\n" % str(gtf))
            noutput += 1

    elif options.method == "transcripts2genes":

        transcripts = set()
        genes = set()
        ignore_strand = options.ignore_strand
        for gtfs in GTF.iterator_transcripts2genes(
                GTF.iterator(options.stdin)):

            ninput += 1
            for gtf in gtfs:
                if ignore_strand:
                    gtf.strand = "."
                options.stdout.write("%s\n" % str(gtf))
                transcripts.add(gtf.transcript_id)
                genes.add(gtf.gene_id)
                nfeatures += 1
            noutput += 1

        E.info("transcripts2genes: transcripts=%i, genes=%i" %
               (len(transcripts), len(genes)))

    elif options.method in ("rename-genes", "rename-transcripts"):

        map_old2new = IOTools.readMap(
            IOTools.openFile(options.filename_filter, "r"))

        if options.method == "rename-transcripts":
            is_gene_id = False
        elif options.method == "rename-genes":
            is_gene_id = True

        for gff in GTF.iterator(options.stdin):
            ninput += 1

            if is_gene_id:
                if gff.gene_id in map_old2new:
                    gff.gene_id = map_old2new[gff.gene_id]
                else:
                    E.debug("removing missing gene_id %s" % gff.gene_id)
                    ndiscarded += 1
                    continue

            else:
                if gff.transcript_id in map_old2new:
                    gff.transcript_id = map_old2new[gff.transcript_id]
                else:
                    E.debug("removing missing transcript_id %s" %
                            gff.transcript_id)
                    ndiscarded += 1
                    continue

            noutput += 1
            options.stdout.write("%s\n" % str(gff))

    elif options.method == "filter":

        keep_genes = set()
        if options.filter_method == "longest-gene":
            iterator = GTF.flat_gene_iterator(GTF.iterator(options.stdin))
            coords = []
            gffs = []
            for gff in iterator:
                gff.sort(key=lambda x: x.start)
                coords.append((gff[0].contig,
                               min([x.start for x in gff]),
                               max([x.end for x in gff]),
                               gff[0].gene_id))
                gffs.append(gff)
            coords.sort()

            last_contig = None
            max_end = 0
            longest_gene_id = None
            longest_length = None

            for contig, start, end, gene_id in coords:
                ninput += 1
                if contig != last_contig or start >= max_end:
                    if longest_gene_id:
                        keep_genes.add(longest_gene_id)
                    longest_gene_id = gene_id
                    longest_length = end - start
                    max_end = end
                else:
                    if end - start > longest_length:
                        longest_length, longest_gene_id = end - start, gene_id
                last_contig = contig
                max_end = max(max_end, end)

            keep_genes.add(longest_gene_id)
            invert = options.invert_filter
            for gff in gffs:
                keep = gff[0].gene_id in keep_genes

                if (keep and not invert) or (not keep and invert):
                    noutput += 1
                    for g in gff:
                        nfeatures += 1
                        options.stdout.write("%s\n" % g)
                else:
                    ndiscarded += 1
        elif options.filter_method in ("longest-transcript",
                                       "representative-transcript"):

            iterator = GTF.gene_iterator(GTF.iterator(options.stdin))

            def selectLongestTranscript(gene):
                r = []
                for transcript in gene:
                    transcript.sort(key=lambda x: x.start)
                    length = transcript[-1].end - transcript[0].start
                    r.append((length, transcript))
                r.sort()
                return r[-1][1]

            def selectRepresentativeTranscript(gene):
                '''select a representative transcript.

                The representative transcript represent the largest number
                of exons over all transcripts.
                '''
                all_exons = []
                for transcript in gene:
                    all_exons.extend([(x.start, x.end)
                                      for x in transcript
                                      if x.feature == "exon"])
                exon_counts = {}
                for key, exons in itertools.groupby(all_exons):
                    exon_counts[key] = len(list(exons))
                transcript_counts = []
                for transcript in gene:
                    count = sum([exon_counts[(x.start, x.end)]
                                 for x in transcript if x.feature == "exon"])
                    # add transcript id to sort to provide a stable
                    # segmentation.
                    transcript_counts.append((count,
                                              transcript[0].transcript_id,
                                              transcript))
                transcript_counts.sort()
                return transcript_counts[-1][-1]

            if options.filter_method == "longest-transcript":
                _select = selectLongestTranscript
            elif options.filter_method == "representative-transcript":
                _select = selectRepresentativeTranscript

            for gene in iterator:
                ninput += 1
                # sort in order to make reproducible which
                # gene is chosen.
                transcript = _select(sorted(gene))
                noutput += 1
                for g in transcript:
                    nfeatures += 1
                    options.stdout.write("%s\n" % g)

        elif options.filter_method in ("gene", "transcript"):

            if options.filename_filter:

                ids = IOTools.readList(
                    IOTools.openFile(options.filename_filter, "r"))
                E.info("read %i ids" % len(ids))

                ids = set(ids)
                by_gene = options.filter_method == "gene"
                by_transcript = options.filter_method == "transcript"
                invert = options.invert_filter

                ignore_strand = options.ignore_strand
                for gff in GTF.iterator(options.stdin):

                    ninput += 1

                    keep = False
                    if by_gene:
                        keep = gff.gene_id in ids
                    if by_transcript:
                        keep = gff.transcript_id in ids
                    if (invert and keep) or (not invert and not keep):
                        continue

                    if ignore_strand:
                        gff.strand = "."

                    options.stdout.write("%s\n" % str(gff))
                    nfeatures += 1
                    noutput += 1

            elif options.sample_size:

                if options.filter_method == "gene":
                    iterator = GTF.flat_gene_iterator(
                        GTF.iterator(options.stdin))
                elif options.filter_method == "transcript":
                    iterator = GTF.transcript_iterator(
                        GTF.iterator(options.stdin))
                if options.min_exons_length:
                    iterator = GTF.iterator_min_feature_length(
                        iterator,
                        min_length=options.min_exons_length,
                        feature="exon")

                data = [x for x in iterator]
                ninput = len(data)
                if len(data) > options.sample_size:
                    data = random.sample(data, options.sample_size)

                for d in data:
                    noutput += 1
                    for dd in d:
                        nfeatures += 1
                        options.stdout.write(str(dd) + "\n")

            else:
                assert False, "please supply either a filename "
                "with ids to filter with (--map-tsv-file) or a sample-size."

        elif options.filter_method in ("proteincoding", "lincrna",
                                       "processed-pseudogene"):
            # extract entries by transcript/gene biotype.
            # This filter uses a test on the source field (ENSEMBL pre v78)
            # a regular expression on the attributes (ENSEMBL >= v78).
            tag = {"proteincoding": "protein_coding",
                   "processed-pseudogene": "processed_pseudogene",
                   "lincrna": "lincRNA"}[options.filter_method]
            rx = re.compile('"%s"' % tag)
            if not options.invert_filter:
                f = lambda x: x.source == tag or rx.search(x.attributes)
            else:
                f = lambda x: x.source != tag and not rx.search(x.attributes)

            for gff in GTF.iterator(options.stdin):
                ninput += 1
                if f(gff):
                    options.stdout.write(str(gff) + "\n")
                    noutput += 1
                else:
                    ndiscarded += 1

    elif options.method == "exons2introns":

        for gffs in GTF.flat_gene_iterator(GTF.iterator(options.stdin)):

            ninput += 1

            cds_ranges = GTF.asRanges(gffs, "CDS")
            exon_ranges = GTF.asRanges(gffs, "exon")
            input_ranges = Intervals.combine(cds_ranges + exon_ranges)

            if len(input_ranges) > 1:
                last = input_ranges[0][1]
                output_ranges = []
                for start, end in input_ranges[1:]:
                    output_ranges.append((last, start))
                    last = end

                if options.intron_border:
                    b = options.intron_border
                    output_ranges = [(x[0] + b, x[1] - b)
                                     for x in output_ranges]

                if options.intron_min_length:
                    l = options.intron_min_length
                    output_ranges = [
                        x for x in output_ranges if x[1] - x[0] > l]

                for start, end in output_ranges:

                    entry = GTF.Entry()
                    entry.copy(gffs[0])
                    entry.clearAttributes()
                    entry.transcript_id = "merged"
                    entry.feature = "intron"
                    entry.start = start
                    entry.end = end
                    options.stdout.write("%s\n" % str(entry))
                    nfeatures += 1
                noutput += 1
            else:
                ndiscarded += 1

    elif options.method == "set-score-to-distance":

        for gffs in GTF.transcript_iterator(GTF.iterator(options.stdin)):
            ninput += 1
            strand = Genomics.convertStrand(gffs[0].strand)
            all_start, all_end = min([x.start for x in gffs]), max(
                [x.end for x in gffs])

            if strand != ".":
                t = 0
                if strand == "-":
                    gffs.reverse()
                for gff in gffs:
                    gff.score = t
                    t += gff.end - gff.start

                if strand == "-":
                    gffs.reverse()
            for gff in gffs:
                options.stdout.write("%s\n" % str(gff))
                nfeatures += 1
            noutput += 1

    elif options.method == "remove-overlapping":

        index = GTF.readAndIndex(
            GTF.iterator(IOTools.openFile(options.filename_gff, "r")))

        for gffs in GTF.transcript_iterator(GTF.iterator(options.stdin)):
            ninput += 1
            found = False
            for e in gffs:
                if index.contains(e.contig, e.start, e.end):
                    found = True
                    break

            if found:
                ndiscarded += 1
            else:
                noutput += 1
                for e in gffs:
                    nfeatures += 1
                    options.stdout.write("%s\n" % str(e))

    elif options.method == "intersect-transcripts":

        for gffs in GTF.gene_iterator(GTF.iterator(options.stdin),
                                      strict=options.strict):

            ninput += 1
            r = []
            for g in gffs:
                if options.with_utr:
                    ranges = GTF.asRanges(g, "exon")
                else:
                    ranges = GTF.asRanges(g, "CDS")
                r.append(ranges)

            result = r[0]
            for x in r[1:]:
                result = Intervals.intersect(result, x)

            entry = GTF.Entry()
            entry.copy(gffs[0][0])
            entry.clearAttributes()
            entry.transcript_id = "merged"
            entry.feature = "exon"
            for start, end in result:
                entry.start = start
                entry.end = end
                options.stdout.write("%s\n" % str(entry))
                nfeatures += 1

            noutput += 1

    elif "rename-duplicates" == options.method:
        # note: this will only rename entries with "CDS" in feature column

        assert options.duplicate_feature in ["gene", "transcript", "both"],\
            ("for renaming duplicates, --duplicate-feature must be set to one "
             "of 'gene', transcript' or 'both'")

        gene_ids = list()
        transcript_ids = list()
        gtfs = list()

        for gtf in GTF.iterator(options.stdin):
            gtfs.append(gtf)
            if gtf.feature == "CDS":
                gene_ids.append(gtf.gene_id)
                transcript_ids.append(gtf.transcript_id)

        dup_gene = [item for item in set(gene_ids) if gene_ids.count(item) > 1]
        dup_transcript = [item for item in set(transcript_ids)
                          if transcript_ids.count(item) > 1]

        E.info("Number of duplicated gene_ids: %i" % len(dup_gene))
        E.info("Number of duplicated transcript_ids: %i" % len(dup_transcript))

        gene_dict = dict(list(zip(dup_gene, ([0] * len(dup_gene)))))
        transcript_dict = dict(list(zip(dup_transcript,
                                        ([0] * len(dup_transcript)))))

        for gtf in gtfs:
            if gtf.feature == "CDS":
                if options.duplicate_feature in ["both", "gene"]:
                    if gtf.gene_id in dup_gene:
                        gene_dict[gtf.gene_id] = gene_dict[gtf.gene_id] + 1
                        # TS. patch until pysam.ctabixproxies.pyx bugfixed
                        gtf.attributes = gtf.attributes.strip()
                        gtf.gene_id = str(gtf.gene_id) + "." + str(gene_dict[gtf.gene_id])

                if options.duplicate_feature in ["both", "transcript"]:
                    if gtf.transcript_id in dup_transcript:
                        transcript_dict[gtf.transcript_id] = \
                            transcript_dict[gtf.transcript_id] + 1
                        # TS. patch until pysam.ctabixproxies.pyx bugfixed
                        gtf.attributes = gtf.attributes.strip()
                        gtf.transcript_id = str(gtf.transcript_id) + "." + str(transcript_dict[gtf.transcript_id])

            options.stdout.write("%s\n" % gtf)

    elif options.method in ("merge-exons",
                            "merge-introns",
                            "merge-transcripts"):
        for gffs in GTF.flat_gene_iterator(
                GTF.iterator(options.stdin),
                strict=options.strict):
            ninput += 1

            cds_ranges = GTF.asRanges(gffs, "CDS")
            exon_ranges = GTF.asRanges(gffs, "exon")
            # sanity checks
            strands = set([x.strand for x in gffs])
            contigs = set([x.contig for x in gffs])
            if len(strands) > 1:
                E.warn('Skipping gene %s on multiple strands' % gffs[0].gene_id)
                break
                raise ValueError(
                    "can not merge gene '%s' on multiple strands: %s" % (
                        gffs[0].gene_id, str(strands)))

            if len(contigs) > 1:
                raise ValueError(
                    "can not merge gene '%s' on multiple contigs: %s" % (
                        gffs[0].gene_id, str(contigs)))

            strand = Genomics.convertStrand(gffs[0].strand)
            utr_ranges = []

            if cds_ranges and options.mark_utr:
                cds_start, cds_end = cds_ranges[0][0], cds_ranges[-1][1]
                midpoint = (cds_end - cds_start) / 2 + cds_start

                utr_ranges = []
                for start, end in Intervals.truncate(exon_ranges, cds_ranges):
                    if end - start > 3:
                        if strand == ".":
                            feature = "UTR"
                        elif strand == "+":
                            if start < midpoint:
                                feature = "UTR5"
                            else:
                                feature = "UTR3"
                        elif strand == "-":
                            if start < midpoint:
                                feature = "UTR3"
                            else:
                                feature = "UTR5"
                        utr_ranges.append((feature, start, end))

            try:
                biotypes = [x["gene_biotype"] for x in gffs]
                biotype = ":".join(set(biotypes))
            except (KeyError, AttributeError):
                biotype = None

            def output_ranges(ranges, gffs, biotype=None,
                              use_geneid=False):
                result = []
                for feature, start, end in ranges:
                    entry = GTF.Entry()
                    entry.copy(gffs[0])
                    entry.clearAttributes()
                    entry.feature = feature
                    if use_geneid:
                        entry.transcript_id = entry.gene_id
                    else:
                        entry.transcript_id = "merged"
                    if biotype:
                        entry.addAttribute("gene_biotype", biotype)
                    entry.start = start
                    entry.end = end
                    result.append(entry)
                return result

            result = []

            if options.method == "merge-exons":
                if options.with_utr:
                    if options.mark_utr:
                        result.extend(output_ranges(utr_ranges, gffs, biotype,
                                                    options.use_geneid))
                        r = [("CDS", x, y) for x, y in
                             Intervals.combineAtDistance(
                                 cds_ranges, options.merge_exons_distance)]
                    else:
                        r = [("exon", x, y) for x, y in
                             Intervals.combineAtDistance(
                                 exon_ranges, options.merge_exons_distance)]
                else:
                    r = [("CDS", x, y) for x, y in
                         Intervals.combineAtDistance(
                             cds_ranges, options.merge_exons_distance)]

            elif options.method == "merge-transcripts":

                if options.with_utr:
                    r = [("exon", exon_ranges[0][0],
                          exon_ranges[-1][1])]
                elif cds_ranges:
                    r = [("exon", cds_ranges[0][0],
                          cds_ranges[-1][1])]
                else:
                    ndiscarded += 1
                    continue

            elif options.method == "merge-introns":

                if len(exon_ranges) >= 2:
                    r = [("exon",
                          exon_ranges[0][1],
                          exon_ranges[-1][0])]
                else:
                    ndiscarded += 1
                    continue

            result.extend(output_ranges(r, gffs, biotype, options.use_geneid))

            result.sort(key=lambda x: x.start)

            for x in result:
                options.stdout.write("%s\n" % str(x))
                nfeatures += 1
            noutput += 1

    elif options.method == "find-retained-introns":

        for gene in GTF.gene_iterator(GTF.iterator(options.stdin)):
            ninput += 1
            found_any = False
            for intron in find_retained_introns(gene):
                found_any = True
                options.stdout.write("%s\n" % str(intron))
                nfeatures += 1
            if found_any:
                noutput += 1

    elif options.method == "genes-to-unique-chunks":

        for gene in GTF.flat_gene_iterator(GTF.iterator(options.stdin)):
            ninput += 1
            for exon in gene_to_blocks(gene):
                options.stdout.write("%s\n" % str(exon))
                nfeatures += 1
            noutput += 1

    else:
        raise ValueError("unknown method '%s'" % options.method)

    E.info("ninput=%i, noutput=%i, nfeatures=%i, ndiscarded=%i" %
           (ninput, noutput, nfeatures, ndiscarded))
    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

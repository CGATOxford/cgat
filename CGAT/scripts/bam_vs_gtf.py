'''bam_vs_gtf.py - compare bam file against gene set
=================================================

:Tags: Genomics NGS Genesets BAM GTF Summary

Purpose
-------

Compare RNASeq reads in a BAM file and compares it against reference exons to quantify exon overrun / underrun.

Documentation
-------------

This script is for validation purposes:
   * Exon overrun should be minimal - reads should not extend beyond known exons.
   * Spliced reads should link known exons.
   
Please note:
   * For unspliced reads, any bases extending beyond exon boundaries are counted. 
   * For spliced reads, both parts of the reads are examined for their overlap.
        As a consequence, counts are doubled for spliced reads.
   * The script requires a list of non-overlapping exons as input.
   * For read counts to be correct the NH (number of hits) flag needs to be set correctly.

Usage
-----

Example::

   # Preview the BAM file using Samtools view
   samtools view tests/bam_vs_gtf.py/small.bam | head
   # Pipe input bam to script and specify gtf file as argument
   cat tests/bam_vs_gtf.py/small.bam | cgat bam_vs_gtf.py --gtf-file=tests/bam_vs_gtf.py/hg19.chr19.gtf.gz

+--------------------+---------+   
|category            |counts   |
+====================+=========+
|spliced_bothoverlap |0        |
+--------------------+---------+
|unspliced_overlap   |0        |
+--------------------+---------+
|unspliced_nooverrun |0        |
+--------------------+---------+
|unspliced	     |207      |
+--------------------+---------+
|unspliced_nooverlap |207      |
+--------------------+---------+
|spliced_overrun     |0        |
+--------------------+---------+
|spliced_halfoverlap |0        |
+--------------------+---------+
|spliced_exact	     |0        |
+--------------------+---------+
|spliced_inexact     |0        |
+--------------------+---------+
|unspliced_overrun   |0        |
+--------------------+---------+
|spliced	     |18       |
+--------------------+---------+
|spliced_underrun    |0        |
+--------------------+---------+
|mapped	             |225      |
+--------------------+---------+
|unmapped	     |0        |
+--------------------+---------+
|input	             |225      |
+--------------------+---------+
|spliced_nooverlap   |18       |
+--------------------+---------+
|spliced_ignored     |0        |
+--------------------+---------+

Type::

   python bam_vs_gtf.py --help

for command line help.

Command line options
--------------------

filename-exons / filename-gtf: a gtf formatted file containing the
genomic coordinates of a set of non-overlapping exons, such as from a
reference genome annotation database (Ensembl, UCSC etc.).

'''

import sys
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import CGAT.GTF as GTF


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-e", "--exons-file", "--gtf-file",
        dest="filename_exons", type="string", metavar="gtf",
        help="gtf formatted file with non-overlapping exon "
        "locations (required). [%default]")

    parser.set_defaults(
        filename_exons=None,
        read_length=200,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    exons = GTF.readAndIndex(
        GTF.iterator(IOTools.openFile(options.filename_exons)))

    pysam_in = pysam.AlignmentFile("-", "rb")

    nspliced = 0
    nspliced_ignored = 0
    nspliced_nooverlap = 0
    nspliced_halfoverlap = 0
    nspliced_bothoverlap = 0
    nspliced_overrun = [0] * 2 * (options.read_length + 10)
    nspliced_exact = 0
    nspliced_inexact = 0
    nunspliced = 0
    nunspliced_overlap = 0
    nunspliced_ignored = 0
    nunspliced_nooverlap = 0
    nunspliced_overrun = [0] * (options.read_length + 10)
    overrun_offset = options.read_length + 10
    ninput = 0
    nunmapped = 0

    c = E.Counter()

    def _splice_overrun(start, end, overlap):
        '''return splicesite over/underrun.

        positive values: overrun
        negative values: underrun
        0: no over/underrun
        '''

        exon_start = min([x[0] for x in overlap])
        exon_end = max([x[1] for x in overlap])

        if start <= exon_start and end > exon_start:
            # overrun at start or match
            r = exon_start - start
        elif start < exon_end and end >= exon_end:
            # overrun at end or match
            r = end - exon_end
        else:
            # underrun - distance to closest exon boundary
            r = -min(start - exon_start, exon_end - end)

        return r

    for read in pysam_in:
        ninput += 1
        if read.is_unmapped:
            nunmapped += 1
            continue

        # check for BAM_CREF_SKIP code in cigar string
        cigar = read.cigar
        is_spliced = 3 in [x[0] for x in cigar]

        contig = pysam_in.getrname(read.tid)
        start = read.pos
        end = read.aend
        if is_spliced:
            # count both ends
            nspliced += 1

            if len(cigar) != 3:
                nspliced_ignored += 1
                continue

            start5, end5 = start, start + cigar[0][1]
            start3, end3 = end - cigar[2][1], end
            try:
                overlap3 = list(exons.get(contig, start3, end3))
                overlap5 = list(exons.get(contig, start5, end5))
            except KeyError:
                overlap3 = overlap5 = []

            ovl3 = len(overlap3)
            ovl5 = len(overlap5)
            o3 = o5 = None
            if not ovl3 and not ovl5:
                nspliced_nooverlap += 1
            elif ovl3 and not ovl5:
                nspliced_halfoverlap += 1
                o3 = _splice_overrun(start3, end3, overlap3)
            elif ovl5 and not ovl3:
                nspliced_halfoverlap += 1
                o5 = _splice_overrun(start5, end5, overlap5)
            else:
                # both overlap
                nspliced_bothoverlap += 1
                o3 = _splice_overrun(start3, end3, overlap3)
                o5 = _splice_overrun(start5, end5, overlap5)

            if o3 is not None:
                if o3 == 0:
                    nspliced_exact += 1
                else:
                    nspliced_inexact += 1
                nspliced_overrun[max(0, overrun_offset + o3)] += 1
            if o5 is not None:
                if o5 == 0:
                    nspliced_exact += 1
                else:
                    nspliced_inexact += 1
                nspliced_overrun[max(0, overrun_offset + o5)] += 1
        else:
            nunspliced += 1
            try:
                overlap = list(exons.get(contig, start, end))
            except KeyError:
                overlap = []

            if len(overlap) == 0:
                nunspliced_nooverlap += 1
            elif len(overlap) >= 1:
                nunspliced_overlap += 1
                # multiple overlap - merge exons (usually: small introns)
                exon_start = min([x[0] for x in overlap])
                exon_end = max([x[1] for x in overlap])
                ostart = max(0, exon_start - start)
                oend = max(0, end - exon_end)
                o = min(end, exon_end) - max(start, exon_start)
                overrun = ostart + oend
                nunspliced_overrun[overrun] += 1

    # output histograms
    outfile = E.openOutputFile("overrun")
    outfile.write(
        "bases\tunspliced_overrun_counts\tspliced_overrun_counts\tspliced_underrun_counts\n")
    _nspliced_overrun = nspliced_overrun[overrun_offset:]
    _nspliced_underrun = nspliced_overrun[:overrun_offset + 1]
    _nspliced_underrun.reverse()
    for x, v in enumerate(zip(nunspliced_overrun, _nspliced_overrun, _nspliced_underrun)):
        outfile.write("%i\t%s\n" % (x, "\t".join(map(str, v))))
    outfile.close()

    # output summary
    # convert to counter
    c.input = ninput
    c.unmapped = nunmapped
    c.mapped = ninput - nunmapped

    c.unspliced = nunspliced
    c.unspliced_nooverlap = nunspliced_nooverlap
    c.unspliced_nooverrun = nunspliced_overrun[0]
    c.unspliced_overlap = nunspliced_overlap
    c.unspliced_overrun = sum(nunspliced_overrun[1:])

    c.spliced = nspliced
    c.spliced_nooverlap = nspliced_nooverlap
    c.spliced_halfoverlap = nspliced_halfoverlap
    c.spliced_bothoverlap = nspliced_bothoverlap
    c.spliced_exact = nspliced_exact
    c.spliced_inexact = nspliced_inexact
    c.spliced_ignored = nspliced_ignored
    c.spliced_underrun = sum(_nspliced_underrun[1:])
    c.spliced_overrun = sum(_nspliced_overrun[1:])

    outfile = options.stdout
    outfile.write("category\tcounts\n")
    for k, v in sorted(c.items()):
        outfile.write("%s\t%i\n" % (k, v))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

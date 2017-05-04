'''gtfs2tsv.py - compare two genesets
==================================

:Tags: Python

Purpose
-------

This script compares two genesets (required) in :term:`gtf`-formatted
files and output lists of shared and unique genes.

It outputs the results of the comparison into various sections. The
sections are split into separate output files whose names are
determined by the ``--output-filename-pattern`` option. The sections
are:

``genes_ovl``
   Table with overlapping genes

``genes_total``
   Summary statistic of overlapping genes

``genes_uniq1``
   List of genes unique in set 1

``genes_uniq2``
   List of genes unique in set 2

Options
-------

``--output-filename-pattern``
   This option defines how the output filenames are determined for the
   sections described in the :term:`Purpose` section above.


Usage
-----

Example::

   head a.gtf::

     19 processed_transcript exon 66346 66509 . - . gene_id "ENSG00000225373";
     transcript_id "ENST00000592209"; exon_number "1"; gene_name "AC008993.5";
     gene_biotype "pseudogene"; transcript_name "AC008993.5-002";
     exon_id "ENSE00001701708";

     19 processed_transcript exon 60521 60747 . - . gene_id "ENSG00000225373";
     transcript_id "ENST00000592209"; exon_number "2"; gene_name "AC008993.5";
     gene_biotype "pseudogene"; transcript_name "AC008993.5-002";
     exon_id "ENSE00002735807";

     19 processed_transcript exon 60105 60162 . - . gene_id "ENSG00000225373";
     transcript_id "ENST00000592209"; exon_number "3"; gene_name "AC008993.5";
     gene_biotype "pseudogene"; transcript_name "AC008993.5-002";
     exon_id "ENSE00002846866";

   head b.gtf::

     19 transcribed_processed_pseudogene exon 66320 66492 . - .
     gene_id "ENSG00000225373"; transcript_id "ENST00000587045"; exon_number "1";
     gene_name "AC008993.5"; gene_biotype "pseudogene";
     transcript_name "AC008993.5-001"; exon_id "ENSE00002739353";

     19 lincRNA exon 68403 69146 . + . gene_id "ENSG00000267111";
     transcript_id "ENST00000589495"; exon_number "1"; gene_name "AC008993.2";
     gene_biotype "lincRNA"; transcript_name "AC008993.2-001";
     exon_id "ENSE00002777656";

     19 lincRNA exon 71161 71646 . + . gene_id "ENSG00000267588";
     transcript_id "ENST00000590978"; exon_number "1"; gene_name "MIR1302-2";
     gene_biotype "lincRNA"; transcript_name "MIR1302-2-001";
     exon_id "ENSE00002870487";

   python gtfs2tsv.py a.gtf b.gtf > out.tsv

   head out.tsv::

     contigs source feature start end score strand frame gene_id transcript_id attributes
     19 processed_transcript exon 66345 66509 . - . ENSG00000225373 ENST00000592209 exon_number "1";
     gene_name "AC008993.5"; gene_biotype "pseudogene"; transcript_name "AC008993.5-002";
     exon_id "ENSE00001701708"
     19 processed_transcript exon 60520 60747 . - . ENSG00000225373 ENST00000592209 exon_number "2";
     gene_name "AC008993.5"; gene_biotype "pseudogene"; transcript_name "AC008993.5-002";
     exon_id "ENSE00002735807"

Type::

   python gtfs2tsv.py --help

for command line help.

Command line options
--------------------

'''
import sys
import os
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import bx.intervals.intersection


def GetNextLine(infile):
    for line in infile:
        if line[0] == "#":
            continue
        return line
    return None


class Counts:

    mPercentFormat = "%.2f"

    def __init__(self, add_percent):

        self.nleft, self.nright, self.noverlap = 0, 0, 0
        self.nunique_left = 0
        self.nunique_right = 0
        self.nidentical = 0
        self.nhalf = 0
        self.nsplit_left = 0
        self.nsplit_right = 0
        self.mAddPercent = add_percent

    def __add__(self, other):
        self.nleft += other.nleft
        self.nright += other.nright
        self.noverlap += other.noverlap
        self.nunique_left += other.nunique_left
        self.nunique_right += other.nunique_right
        self.nidentical += other.nidentical
        self.nhalf += other.nhalf
        self.nsplit_left += other.nsplit_left
        self.nsplit_right += other.nsplit_right
        return self

    def getHeader(self):
        h = "total_left\ttotal_right\tnoverlap\tnidentical\tnhalf\tunique_left\tunique_right\tsplit_left\tsplit_right"
        if self.mAddPercent:
            h += "\t" + self.getHeaderPercent()
        return h

    def __str__(self):
        h = "\t".join(map(str,
                          (self.nleft, self.nright,
                              self.noverlap, self.nidentical, self.nhalf,
                              self.nunique_left, self.nunique_right,
                              self.nsplit_left, self.nsplit_right)))
        if self.mAddPercent:
            h += "\t" + self.asPercent()

        return h

    def getHeaderPercent(self):
        return "\t".join(["pl%s\tpr%s" % (x, x) for x in ("overlap", "identical", "half", "unique", "split")])

    def asPercent(self):
        return "\t".join([self.mPercentFormat % (100.0 * x) for x in (
            float(self.noverlap) / self.nleft,
            float(self.noverlap) / self.nright,
            float(self.nidentical) / self.nleft,
            float(self.nidentical) / self.nright,
            float(self.nhalf) / self.nleft,
            float(self.nhalf) / self.nright,
            float(self.nunique_left) / self.nleft,
            float(self.nunique_right) / self.nright,
            float(self.nsplit_left) / self.nleft,
            float(self.nsplit_right) / self.nright)])


def getFile(options, section):

    if options.output_filename_pattern:
        outfile = IOTools.openFile(
            options.output_filename_pattern % section, "w")
        E.info("output for section '%s' goes to file %s" %
               (section, options.output_filename_pattern % section))
    else:
        outfile = options.stdout
        outfile.write("## section: %s\n" % section)
    return outfile


def writeDiff(outfile, symbol, genes):

    for gene in sorted(genes):
        for exon in sorted(gene):
            outfile.write("%s\t%s\n" % (symbol, str(exon)))


def main(argv=None):

    if not argv:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option(
        "-e", "--output-equivalent", dest="write_equivalent",
        action="store_true",
        help="write equivalent entries [default=%default].")

    parser.add_option(
        "-f", "--output-full", dest="write_full",
        action="store_true",
        help="write full gff entries [default=%default].")

    parser.add_option("-p", "--add-percent", dest="add_percent",
                      action="store_true",
                      help="add percentage columns [default=%default].")

    parser.add_option("-s", "--ignore-strand", dest="ignore_strand",
                      action="store_true",
                      help="ignore strand information [default=%default].")

    parser.set_defaults(
        write_equivalent=False,
        write_full=False,
        add_percent=False,
        ignore_strand=False,
        as_gtf=False,
    )

    (options, args) = E.Start(parser, argv, add_output_options=True)

    if len(args) != 2:
        raise ValueError("two arguments required")

    input_filename1, input_filename2 = args

    # duplicated features cause a problem. Make sure
    # features are non-overlapping by running
    # gff_combine.py on GFF files first.

    E.info("reading data started")

    idx, genes2 = {}, set()
    for e in GTF.readFromFile(IOTools.openFile(input_filename2, "r")):
        genes2.add(e.gene_id)
        if e.contig not in idx:
            idx[e.contig] = bx.intervals.intersection.Intersecter()
        idx[e.contig].add_interval(
            bx.intervals.Interval(e.start, e.end, value=e))

    overlaps_genes = []

    E.info("reading data finished: %i contigs" % len(idx))

    # outfile_diff and outfile_overlap not implemented
    # outfile_diff = getFile( options, "diff" )
    # outfile_overlap = getFile( options, "overlap" )
    overlapping_genes = set()

    genes1 = set()

    # iterate over exons
    with IOTools.openFile(input_filename1, "r") as infile:
        for this in GTF.iterator(infile):

            genes1.add(this.gene_id)

            try:
                intervals = idx[this.contig].find(this.start, this.end)
            except KeyError:
                continue

            others = [x.value for x in intervals]
            for other in others:
                overlapping_genes.add((this.gene_id, other.gene_id))

            # check for identical/half-identical matches
            output = None
            for other in others:
                if this.start == other.start and this.end == other.end:
                    output, symbol = other, "="
                    break
            else:
                for other in others:
                    if this.start == other.start or this.end == other.end:
                        output, symbol = other, "|"
                        break
                else:
                    symbol = "~"

    # if outfile_diff != options.stdout: outfile_diff.close()
    # if outfile_overlap != options.stdout: outfile_overlap.close()

    outfile = None
    ##################################################################
    ##################################################################
    ##################################################################
    # print gene based information
    ##################################################################
    if overlapping_genes:
        outfile = getFile(options, "genes_ovl")
        outfile.write("gene_id1\tgene_id2\n")
        for a, b in sorted(overlapping_genes):
            outfile.write("%s\t%s\n" % (a, b))
        if outfile != options.stdout:
            outfile.close()

        outfile_total = getFile(options, "genes_total")
        outfile_total.write(
            "set\tngenes\tnoverlapping\tpoverlapping\tnunique\tpunique\n")

        outfile = getFile(options, "genes_uniq1")
        b = set([x[0] for x in overlapping_genes])
        d = genes1.difference(b)
        outfile.write("gene_id1\n")
        outfile.write("\n".join(sorted(d)) + "\n")
        if outfile != options.stdout:
            outfile.close()
        outfile_total.write("%s\t%i\t%i\t%5.2f\t%i\t%5.2f\n" % (
            os.path.basename(input_filename1), len(
                genes1), len(b), 100.0 * len(b) / len(a),
            len(d), 100.0 * len(d) / len(genes1)))

        outfile = getFile(options, "genes_uniq2")
        b = set([x[1] for x in overlapping_genes])
        d = genes2.difference(b)
        outfile.write("gene_id2\n")
        outfile.write("\n".join(sorted(d)) + "\n")
        if outfile != options.stdout:
            outfile.close()

        outfile_total.write("%s\t%i\t%i\t%5.2f\t%i\t%5.2f\n" % (
            os.path.basename(input_filename2), len(
                genes2), len(b), 100.0 * len(b) / len(a),
            len(d), 100.0 * len(d) / len(genes2)))
        if outfile_total != options.stdout:
            outfile_total.close()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

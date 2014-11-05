##########################################################################
#   Gene prediction pipeline
#
#   $Id: diff_gff.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''diff_gff.py - compare two gff files
======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script compares the intervals in two :term:`gff` formatted files.

.. note::

    This script only works with non-overlapping features.

The script outputs the following symbols::

    >: an entry unique to the second file
    <: an entry unique to the first file

    if --output-equivalent is set, then it also outputs overlapping entries:

    ~: partially overlapping entries
    /: an entry that is split in the first file (split_right)
    \: an entry that is split in the second file (split_left)
    =: identical entries
    |: half-identical entries, only one boundary is matching

Usage
-----

Example::

   python diff_gff.py a.gtf.gz b.gtf.gz

Type::

   python diff_gff.py --help

for command line help.

Documentation
-------------

Code
----
'''
import sys
import re
import os

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools


def GetNextLine(infile):
    for line in infile:
        if line[0] == "#":
            continue
        return line
    return None


class Counts:

    mPercentFormat = "%5.2f"

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
        return "\t".join(map(lambda x: "pl%s\tpr%s" % (x, x), ("overlap", "identical", "half", "unique", "split")))

    def asPercent(self):

        def toPercent(a, b):
            if b > 0:
                return self.mPercentFormat % (100.0 * float(a) / b)
            else:
                return "na"

        return "\t".join(map(lambda x: toPercent(x[0], x[1]),
                             ((self.noverlap, self.nleft),
                              (self.noverlap, self.nright),
                              (self.nidentical, self.nleft),
                              (self.nidentical, self.nright),
                              (self.nhalf, self.nleft),
                              (self.nhalf, self.nright),
                              (self.nunique_left, self.nleft),
                              (self.nunique_right, self.nright),
                              (self.nsplit_left, self.nleft),
                              (self.nsplit_right, self.nright))))


def getFile(options, section):

    if options.output_filename_pattern:
        outfile = open(options.output_filename_pattern % section, "w")
        E.info("output for section '%s' goes to file %s" %
               (section, options.output_filename_pattern % section))
    else:
        outfile = options.stdout
        outfile.write("## section: %s\n" % section)
    return outfile


def _cmp(this, other):
    '''compare to gtf entries'''
    return cmp((this.contig, this.strand, this.start),
               (other.contig, other.strand, other.start))


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: diff_gff.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-e", "--output-equivalent", dest="write_equivalent",
                      help="write equivalent entries [default=%default].", action="store_true")

    parser.add_option("-f", "--output-full", dest="write_full",
                      help="write full gff entries [default=%default].", action="store_true")

    parser.add_option("-o", "--format=", dest="format",
                      help="output format [flat|multi-line] [default=%default]")

    parser.add_option("-p", "--add-percent", dest="add_percent", action="store_true",
                      help="add percentage columns [default=%default].")

    parser.add_option("-a", "--as-gtf", "--is-gtf", dest="as_gtf", action="store_true",
                      help="input is in gtf format. Output on overlapping genes will be output [default=%default].")

    parser.add_option("-s", "--ignore-strand", dest="ignore_strand", action="store_true",
                      help="ignore strand information [default=%default].")

    parser.set_defaults(
        write_equivalent=False,
        write_full=False,
        format="flat",
        add_percent=False,
        ignore_strand=False,
        as_gtf=False,
    )

    (options, args) = E.Start(parser, add_output_options=True)

    if len(args) != 2:
        raise ValueError("two arguments required")

    input_filename1, input_filename2 = args

    # duplicated features cause a problem. Make sure
    # features are non-overlapping by running
    # gff_combine.py on GFF files first.

    E.info("reading data")

    if options.as_gtf:
        gff1 = GTF.readFromFile(IOTools.openFile(input_filename1, "r"))
        gff2 = GTF.readFromFile(IOTools.openFile(input_filename2, "r"))
        overlaps_genes = []
    else:
        gff1 = GTF.readFromFile(IOTools.openFile(input_filename1, "r"))
        gff2 = GTF.readFromFile(IOTools.openFile(input_filename2, "r"))

    E.info("reading data finished: %i, %i" % (len(gff1), len(gff2)))

    # removing everything but exons
    gff1 = [x for x in gff1 if x.feature == "exon"]
    gff2 = [x for x in gff2 if x.feature == "exon"]

    E.info("after keeping only 'exons': %i, %i" % (len(gff1), len(gff2)))

    if options.ignore_strand:
        for e in gff1:
            e.strand = "."
        for e in gff2:
            e.strand = "."

    E.info("sorting exons")

    gff1.sort(key=lambda x: (x.contig, x.strand, x.start, x.end))
    gff2.sort(key=lambda x: (x.contig, x.strand, x.start, x.end))

    E.info("sorting exons finished")

    subtotals = []
    subtotal = Counts(add_percent=options.add_percent)

    outfile_diff = getFile(options, "diff")
    outfile_overlap = getFile(options, "overlap")

    if options.as_gtf:
        overlapping_genes = []
    else:
        overlapping_genes = None

    i1, i2 = 0, 0
    n1 = len(gff1)
    n2 = len(gff2)
    first_entry2, first_entry1 = None, None

    while i1 < n1 and i2 < n2:

        entry1 = gff1[i1]
        entry2 = gff2[i2]

        E.debug("1: i1=%i n1=%i entry1=%s" % (i1, n1, str(entry1)))
        E.debug("2: i2=%i n2=%i entry2=%s" % (i2, n2, str(entry2)))

        # when chromosome/strand have changed in both (and are the same), print
        # summary info:
        if first_entry1:

            if (first_entry1.contig != entry1.contig or
                first_entry1.strand != entry1.strand) and \
                (first_entry2.contig != entry2.contig or
                 first_entry2.strand != entry2.strand) and \
                    entry1.contig == entry2.contig and \
                    entry1.strand == entry2.strand:

                subtotals.append(
                    (first_entry1.contig, first_entry1.strand, subtotal))
                subtotal = Counts(add_percent=options.add_percent)
                first_entry1 = entry1
                first_entry2 = entry2

        else:
            first_entry1 = entry1
            first_entry2 = entry2

        output_1, output_2 = None, None

        if GTF.Overlap(entry1, entry2):

            # collect multiple matches
            last_l = True
            while GTF.Overlap(entry1, entry2):

                if overlapping_genes is not None:
                    overlapping_genes.append((entry1.gene_id, entry2.gene_id))

                write_last = True
                subtotal.noverlap += 1
                if entry1.start == entry2.start and entry1.end == entry2.end:
                    symbol = "="
                    subtotal.nidentical += 1
                elif entry1.start == entry2.start or entry1.end == entry2.end:
                    symbol = "|"
                    subtotal.nhalf += 1
                else:
                    symbol = "~"

                output_1 = entry1
                output_2 = entry2

                if entry1.end < entry2.end:
                    i1 += 1
                    subtotal.nleft += 1
                    last_l = True

                    if i1 >= n1:
                        i2 += 1
                        break

                    entry1 = gff1[i1]
                    if GTF.Overlap(entry1, entry2):
                        symbol = "/"
                        # outfile.write( "# split right\n" )
                        subtotal.nsplit_right += 1

                else:
                    i2 += 1
                    subtotal.nright += 1
                    last_l = False

                    if i2 >= n2:
                        i1 += 1
                        break

                    entry2 = gff2[i2]
                    if GTF.Overlap(entry1, entry2):
                        symbol = "\\"
                        # outfile.write("# split left\n")
                        subtotal.nsplit_left += 1

                # output at the end, so that symbol is known
                if options.write_equivalent:
                    if options.format == "flat":
                        outfile_overlap.write(
                            "%s\t%s\t%s\n" % (symbol, str(output_1), str(output_2)))
                    elif options.format == "multi-line":
                        outfile_overlap.write(
                            "%s\t%s\n\t%s\n" % (symbol, str(output_1), str(output_2)))

                write_last = False

            if write_last and output_1 and output_2 and options.write_equivalent:
                if options.format == "flat":
                    outfile_overlap.write(
                        "%s\t%s\t%s\n" % (symbol, str(output_1), str(output_2)))
                elif options.format == "multi-line":
                    outfile_overlap.write(
                        "%s\t%s\n\t%s\n" % (symbol, str(output_1), str(output_2)))

            # if last advance was left, go right, and vice versa
            if last_l:
                i2 += 1
                subtotal.nright += 1
            else:
                i1 += 1
                subtotal.nleft += 1

        elif _cmp(entry1, entry2) < 0:
            outfile_diff.write("<\t%s\n" % str(entry1))
            subtotal.nunique_left += 1
            i1 += 1
            subtotal.nleft += 1

        elif _cmp(entry1, entry2) > 0:
            outfile_diff.write(">\t%s\n" % str(entry2))
            subtotal.nunique_right += 1
            i2 += 1
            subtotal.nright += 1

    while i1 < n1:
        outfile_diff.write("<\t%s\n" % str(entry1))
        subtotal.nunique_left += 1
        i1 += 1
        if i1 >= n1:
            break
        entry1 = gff1[i1]
        subtotal.nleft += 1

    while i2 < n2:
        outfile_diff.write(">\t%s\n" % str(entry2))
        subtotal.nunique_right += 1
        i2 += 1
        if i2 >= n2:
            break
        entry2 = gff2[i2]
        subtotal.nright += 1

    subtotals.append((entry1.contig, entry1.strand, subtotal))

    if outfile_diff != options.stdout:
        outfile_diff.close()
    if outfile_overlap != options.stdout:
        outfile_overlap.close()

    ##################################################################
    ##################################################################
    ##################################################################
    # print gene based information
    ##################################################################
    if overlapping_genes:
        outfile = getFile(options, "genes_ovl")
        s = set(overlapping_genes)
        outfile.write("gene_id1\tgene_id2\n")
        for a, b in s:
            outfile.write("%s\t%s\n" % (a, b))
        if outfile != options.stdout:
            outfile.close()

        outfile_total = getFile(options, "genes_total")
        outfile_total.write(
            "set\tngenes\tnoverlapping\tpoverlapping\tnunique\tpunique\n")

        outfile = getFile(options, "genes_uniq1")
        a = set([x.gene_id for x in gff1])
        b = set([x[0] for x in s])
        d = a.difference(b)
        outfile.write("gene_id1\n")
        outfile.write("\n".join(d) + "\n")
        if outfile != options.stdout:
            outfile.close()
        outfile_total.write("%s\t%i\t%i\t%5.2f\t%i\t%5.2f\n" % (
            os.path.basename(input_filename1), len(
                a), len(b), 100.0 * len(b) / len(a),
            len(d), 100.0 * len(d) / len(a)))

        outfile = getFile(options, "genes_uniq2")
        a = set([x.gene_id for x in gff2])
        b = set([x[1] for x in s])
        d = a.difference(b)
        outfile.write("gene_id2\n")
        outfile.write("\n".join(d) + "\n")
        if outfile != options.stdout:
            outfile.close()

        outfile_total.write("%s\t%i\t%i\t%5.2f\t%i\t%5.2f\n" % (
            os.path.basename(input_filename2), len(
                a), len(b), 100.0 * len(b) / len(a),
            len(d), 100.0 * len(d) / len(a)))
        if outfile_total != options.stdout:
            outfile_total.close()

    ##################################################################
    ##################################################################
    ##################################################################
    # print totals
    ##################################################################
    outfile = getFile(options, "total")
    outfile.write("chr\tstrand\t%s\n" %
                  Counts(add_percent=options.add_percent).getHeader())

    total = Counts(add_percent=options.add_percent)
    for x in subtotals:
        outfile.write("\t".join((x[0], x[1], str(x[2]))) + "\n")
        total += x[2]

    outfile.write("\t".join(("all", "all", str(total))) + "\n")

    if outfile != options.stdout:
        outfile.close()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

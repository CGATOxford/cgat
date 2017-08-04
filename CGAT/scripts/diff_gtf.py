'''diff_gtf.py - compute overlap between multiple gtf files
===========================================================

:Tags: Genomics Intervals Genesets GTF Comparison

Purpose
-------

This script compares multiple set of gtf files. It computes
the overlap between bases, exons and genes between each pair
of gtf files.

If results from a previous run are present, existing
pairs are not re-computed but simply echoed.

The output is a tab-separated table with counts for each pair
of files being compared. The fields are:

+--------------+----------------------------------+
|*Column*      |*Content*                         |
+--------------+----------------------------------+
|set           |Name of the set                   |
+--------------+----------------------------------+
|ngenes_total  |number of genes in set            |
+--------------+----------------------------------+
|ngenes_ovl    |number of genes overlapping       |
+--------------+----------------------------------+
|ngenes_unique |number of unique genes            |
+--------------+----------------------------------+
|nexons_total  |number of exons in set            |
+--------------+----------------------------------+
|nexons_ovl    |number of exons overlapping       |
+--------------+----------------------------------+
|nexons_unique |number of unique exons            |
+--------------+----------------------------------+
|nbases_total  |number of bases in gene set       |
+--------------+----------------------------------+
|nbases_ovl    |number of bases overlapping       |
+--------------+----------------------------------+
|nbases_unique |number of unique bases            |
+--------------+----------------------------------+

Each of these fields will appear twice, once for each of the pair of files.
Hence ngenes_unique1 will be the number of genes in set1 that have no exons
that overlap with any exons in set2, and vice versa for ngenes_unique2. And
on for each field in the table above. This makes a total of 9*2=18 fields
containing counts, each starting with an n.

A further set of 18 fields each start with a ``p`` and are the corresponding
percentage values.

Options
-------

-s, --ignore-strand
    Ignore strand infomation so that bases overlap even if exons/genes are
    on different strands

-u, --update=FILENAME
    Read in previous results from FILENAME and only output comparisons that
    are missing.

-p, --pattern-identifier=PATTERN
    Provide a regular expression pattern for converting a filename into a
    set name for the output. The regular expression should capture at least
    one group. That group will be used to identify that file in the output
    table (see examples)


Examples
--------

For example if we have two gtf_files that look like::

   first_set_of_genes.gtf:
   1	protein_coding	exon	1	10	.	+	.	gene_id "1"; transcript_id "1"
   1	protein_coding	exon	20	30	.	+	.	gene_id "1"; transcript_id "1"

   second_set_of_genes.gtf:
   1	protein_coding	exon	25	35	.	+	.	gene_id "1"; transcript_id "1"
   2	protein_coding	exon	100	200	.	+	.	gene_id "2"; transcript_id "3"

Then the command::

   python diff_gtf.py *.gtf --pattern-identifier='(.+)_of_genes.gtf' > out.tsv

would produce an output file that has a single row with set1 being "second_set"
and set2 being "first_set" (these are extracted using that --pattern-identifier
option). It will report that set1 contains 2 genes and set2 1 gene. That for
each set one of these genes overlaps with the other set. For set1 it will
report that 1 gene is unique and that no genes are unique for set2 and so on
for exons and bases.

If we want to add a third file to the comparison,
"third_set_of_genes.gtf", we don't need to redo the comparison between
first_set_of_genes.gtf and second_set_of_genes.gtf::

   python diff_gtf.py --update=out.tsv *.gtf.gz > new.tsv

This will output a table with a row for third_set vs second_set and
third_set vs second_set, along with the comparison of first_set and
second_set that will simply be copied from the previous results. It is
important to include all files on the command line that are to be
output. Any comparisons in ``out.tsv`` that correspond to files that
are not given on the command line will not be output.

Usage
-----

::
   cgat diff_gtf.py GTF GTF [GTF [GTF [...]]] [OPTIONS]
   cgat diff_gtf GTF1 --update=OUTFILE [OPTIONS]

where GTF is a gtf or compressed gtf formated file and OUTFILE is the results
from a previous run.  At least two must be provided unless --update is present.

Type::

   python diff_gtf.py --help

for command line help.

Command line options
--------------------

'''

import sys
import re
import numpy

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.NCL as NCL


class Counter:

    mPercentFormat = "%5.2f"

    def __init__(self):
        pass

    def getHeader(self):
        h = []
        for a in ("genes", "exons", "bases"):
            for b in ("total", "ovl", "unique"):
                for c in ("1", "2"):
                    h.append("n" + a + "_" + b + c)
        for a in ("genes", "exons", "bases"):
            for b in ("ovl", "unique"):
                for c in ("1", "2"):
                    h.append("p" + a + "_" + b + c)

        return "\t".join(h)

    @E.cachedmethod
    def buildIndex(self, filename):
        """read and index."""

        idx = {}
        infile = IOTools.openFile(filename, "r")
        for e in GTF.readFromFile(infile):
            if e.contig not in idx:
                idx[e.contig] = NCL.NCLSimple()
            idx[e.contig].add(e.start, e.end)
        infile.close()
        return idx

    def _count(self, filename, idx):

        overlapping_genes = set()
        genes = set()
        # iterate over exons
        infile = IOTools.openFile(filename, "r")
        it = GTF.iterator(infile)

        nexons, nexons_overlapping = 0, 0
        nbases, nbases_overlapping = 0, 0
        for this in it:
            nexons += 1
            nbases += this.end - this.start
            genes.add(this.gene_id)

            try:
                intervals = list(idx[this.contig].find(this.start, this.end))
            except KeyError:
                continue

            if len(intervals) == 0:
                continue

            overlapping_genes.add(this.gene_id)
            nexons_overlapping += 1
            start, end = this.start, this.end
            counts = numpy.zeros(end - start, numpy.int)
            for other_start, other_end, other_value in intervals:
                for x in range(max(start, other_start) - start, min(end, other_end) - start):
                    counts[x] += 1
            nbases_overlapping += sum([1 for x in counts if x > 0])

        infile.close()

        return len(genes), len(overlapping_genes), nexons, nexons_overlapping, nbases, nbases_overlapping

    def count(self, filename1, filename2):
        """count overlap between two gtf files."""

        E.info("counting started for %s versus %s" % (filename1, filename2))

        idx2 = self.buildIndex(filename2)

        (self.mGenes1, self.mGenesOverlapping1,
         self.mExons1, self.mExonsOverlapping1,
         self.mBases1, self.mBasesOverlapping1 ) = \
            self._count(filename1, idx2)

        self.mGenesUnique1 = self.mGenes1 - self.mGenesOverlapping1
        self.mExonsUnique1 = self.mExons1 - self.mExonsOverlapping1
        self.mBasesUnique1 = self.mBases1 - self.mBasesOverlapping1

        idx1 = self.buildIndex(filename1)

        (self.mGenes2, self.mGenesOverlapping2,
         self.mExons2, self.mExonsOverlapping2,
         self.mBases2, self.mBasesOverlapping2 ) = \
            self._count(filename2, idx1)

        self.mGenesUnique2 = self.mGenes2 - self.mGenesOverlapping2
        self.mExonsUnique2 = self.mExons2 - self.mExonsOverlapping2
        self.mBasesUnique2 = self.mBases2 - self.mBasesOverlapping2

    def __str__(self):

        return "\t".join(map(str, (
            self.mGenes1, self.mGenes2,
            self.mGenesOverlapping1, self.mGenesOverlapping2,
            self.mGenesUnique1, self.mGenesUnique2,
            self.mExons1, self.mExons2,
            self.mExonsOverlapping1, self.mExonsOverlapping2,
            self.mExonsUnique1, self.mExonsUnique2,
            self.mBases1, self.mBases2,
            self.mBasesOverlapping1, self.mBasesOverlapping2,
            self.mBasesUnique1, self.mBasesUnique2 ) ) ) + "\t" +\
            "\t".join([IOTools.prettyPercent(*x) for x in (
                (self.mGenesOverlapping1, self.mGenes1),
                (self.mGenesOverlapping2, self.mGenes2),
                (self.mGenesUnique1, self.mGenes1),
                (self.mGenesUnique2, self.mGenes2),
                (self.mExonsOverlapping1, self.mExons1),
                (self.mExonsOverlapping2, self.mExons2),
                (self.mExonsUnique1, self.mExons1),
                (self.mExonsUnique2, self.mExons2),
                (self.mBasesOverlapping1, self.mBases1),
                (self.mBasesOverlapping2, self.mBases2),
                (self.mBasesUnique1, self.mBases1),
                (self.mBasesUnique2, self.mBases2))])


class CounterGenes(Counter):

    """output only genes."""

    mSeparator = ";"

    def __init__(self, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs)

    def getHeader(self):
        h = ["ngenes1", "ngenes2", "novl1", "novl2", "nuniq1",
             "nuniq2", "ovl1", "ovl2", "uniq1", "uniq2"]
        return "\t".join(h)

    def _count(self, filename, idx):

        overlapping_genes = set()
        genes = set()
        # iterate over exons
        infile = IOTools.openFile(filename, "r")
        it = GTF.iterator(infile)

        for this in it:
            genes.add(this.gene_id)

            try:
                intervals = idx[this.contig].find(this.start, this.end)
            except KeyError:
                continue

            if len(intervals) == 0:
                continue
	
            overlapping_genes.add(this.gene_id)

        infile.close()

        return genes, overlapping_genes

    def count(self, filename1, filename2):
        """count overlap between two gtf files."""

        E.info("counting started for %s versus %s" % (filename1, filename2))

        idx2 = self.buildIndex(filename2)

        (self.mGenes1, self.mGenesOverlapping1) = self._count(filename1, idx2)

        idx1 = self.buildIndex(filename1)
        (self.mGenes2, self.mGenesOverlapping2) = self._count(filename2, idx1)

    def __str__(self):

        uniq1 = self.mGenes1.difference(self.mGenesOverlapping)
        uniq2 = self.mGenes2.difference(self.mGenesOverlapping)

        return "\t".join(map(str, (
            len(self.mGenes1),
            len(self.mGenes2),
            len(self.mGenesOverlapping1),
            len(self.mGenesOverlapping2),
            len(uniq1),
            len(uniq2),
            self.mSeparator.join(self.mGenesOverlapping1),
            self.mSeparator.join(self.mGenesOverlapping2),
            self.mSeparator.join(uniq1),
            self.mSeparator.join(uniq2))))


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-s", "--ignore-strand", dest="ignore_strand",
                      action="store_true",
                      help="ignore strand information [default=%default].")

    parser.add_option(
        "-u", "--update", dest="filename_update", type="string",
        help="if filename is given, previous results will be read"
        "from there and only changed sets will be computed "
        "[default=%default].")

    parser.add_option(
        "-p", "--pattern-identifier", dest="pattern_id", type="string",
        help="pattern to convert a filename to an id"
        "[default=%default].")

    parser.add_option(
        "-g", "--output-only-genes", dest="output_only_genes",
        action="store_true",
        help="only output gene stats (includes gene lists)"
        " [default=%default].")

    parser.set_defaults(
        ignore_strand=False,
        filename_update=None,
        pattern_id="(.*).gtf",
        output_only_genes=False,
    )

    (options, args) = E.Start(parser)

    if len(args) < 2:
        print(USAGE)
        raise ValueError("at least two arguments required")

    if options.filename_update:
        infile = IOTools.openFile(options.filename_update, "r")
        previous_results = {}
        for line in infile:
            if line.startswith("#"):
                continue
            if line.startswith("set1"):
                continue
            data = line[:-1].split("\t")
            set1, set2 = data[0], data[1]

            if set1 not in previous_results:
                previous_results[set1] = {}
            if set2 not in previous_results:
                previous_results[set2] = {}

            previous_results[set1][set2] = "\t".join(data[2:])
            rev = [(data[x + 1], data[x]) for x in range(2, len(data), 2)]
            previous_results[set2][set1] = "\t".join(IOTools.flatten(rev))
    else:
        previous_results = {}

    if options.output_only_genes:
        counter = CounterGenes()
    else:
        counter = Counter()

    options.stdout.write("set1\tset2\t%s\n" % counter.getHeader())

    pattern_id = re.compile(options.pattern_id)

    def getTitle(x):
        try:
            return pattern_id.search(x).groups()[0]
        except AttributeError:
            return x

    ncomputed, nupdated = 0, 0
    for x in range(len(args)):
        title1 = getTitle(args[x])
        for y in range(0, x):
            title2 = getTitle(args[y])
            if previous_results:
                try:
                    prev = previous_results[title1][title2]
                except KeyError:
                    pass
                else:
                    options.stdout.write(
                        "%s\t%s\t%s\n" % ((title1, title2, prev)))
                    nupdated += 1
                    continue

            counter.count(args[x], args[y])
            options.stdout.write(
                "%s\t%s\t%s\n" % ((title1, title2, str(counter))))
            ncomputed += 1

    E.info("nupdated=%i, ncomputed=%i" % (nupdated, ncomputed))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

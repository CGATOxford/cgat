"""
diff_bed.py - count differences between several bed files
=========================================================

:Tags: Genomics Intervals BED Comparison

Purpose
-------

Compute overlap statistics between multiple bed files. For each pairwise
comparison, this script outputs the number of intervals (exons) and
bases overlapping.

Using the ``--update`` option, a table can be incrementally updated with
additional comparisons.

The strand of intervals is ignored in comparisons.

+--------------+----------------------------------+
|*Column*      |*Content*                         |
+--------------+----------------------------------+
|set           |Name of the set                   |
+--------------+----------------------------------+
|nexons_total  |number of intervals in set        |
+--------------+----------------------------------+
|nexons_ovl    |number of intervals overlapping   |
+--------------+----------------------------------+
|nexons_unique |number of unique intervals        |
+--------------+----------------------------------+
|nbases_total  |number of bases in gene set       |
+--------------+----------------------------------+
|nbases_ovl    |number of bases overlapping       |
+--------------+----------------------------------+
|nbases_unique |number of unique bases            |
+--------------+----------------------------------+

Usage
-----

For example::

   python diff_bed.py *.bed.gz > out.tsv

To update results from a previous run, type::

   python diff_bed.py --update=out.tsv *.bed.gz > new.tsv

Type::

   python diff_bed.py --help

for command line help.

Command line options
--------------------

"""

import sys
import re
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Bed as Bed
import numpy


class Counter:

    mPercentFormat = "%5.2f"

    def __init__(self):
        pass

    def getHeader(self):
        h = []
        for a in ("exons", "bases"):
            for b in ("total", "ovl", "unique"):
                for c in ("1", "2"):
                    h.append("n" + a + "_" + b + c)
        for a in ("exons", "bases"):
            for b in ("ovl", "unique"):
                for c in ("1", "2"):
                    h.append("p" + a + "_" + b + c)

        return "\t".join(h)

    @E.cachedmethod
    def buildIndex(self, filename):
        return Bed.readAndIndex(IOTools.openFile(filename, "r"))

    def _count(self, filename, idx):
        '''count filename against idx.'''

        overlapping_genes = set()
        genes = set()

        # iterate over exons
        infile = IOTools.openFile(filename, "r")
        it = Bed.bed_iterator(infile)

        nexons, nexons_overlapping = 0, 0
        nbases, nbases_overlapping = 0, 0
        for this in it:
            nexons += 1
            nbases += this.end - this.start

            try:
                intervals = list(
                    idx[this.contig].find(max(0, this.start), this.end))
            except KeyError:
                continue
            except Exception as msg:
                raise Exception(
                    "error while processing %s, msg=%s" % (filename, msg))
            if len(intervals) == 0:
                continue

            nexons_overlapping += 1
            start, end = this.start, this.end
            counts = numpy.zeros(end - start, numpy.int)
            for other_start, other_end, other_value in intervals:
                for x in range(max(start, other_start) - start, min(end, other_end) - start):
                    counts[x] += 1
            nbases_overlapping += sum([1 for x in counts if x > 0])

        infile.close()

        return nexons, nexons_overlapping, nbases, nbases_overlapping

    def count(self, filename1, filename2):
        """count overlap between two bed files."""

        E.info("counting started for %s versus %s" % (filename1, filename2))

        idx2 = self.buildIndex(filename2)

        (self.mExons1, self.mExonsOverlapping1,
         self.mBases1, self.mBasesOverlapping1 ) = \
            self._count(filename1, idx2)

        self.mExonsUnique1 = self.mExons1 - self.mExonsOverlapping1
        self.mBasesUnique1 = self.mBases1 - self.mBasesOverlapping1

        idx1 = self.buildIndex(filename1)

        (self.mExons2, self.mExonsOverlapping2,
         self.mBases2, self.mBasesOverlapping2 ) = \
            self._count(filename2, idx1)

        self.mExonsUnique2 = self.mExons2 - self.mExonsOverlapping2
        self.mBasesUnique2 = self.mBases2 - self.mBasesOverlapping2

    def __str__(self):

        return "\t".join(map(str, (
            self.mExons1, self.mExons2,
            self.mExonsOverlapping1, self.mExonsOverlapping2,
            self.mExonsUnique1, self.mExonsUnique2,
            self.mBases1, self.mBases2,
            self.mBasesOverlapping1, self.mBasesOverlapping2,
            self.mBasesUnique1, self.mBasesUnique2 ) ) ) + "\t" +\
            "\t".join([IOTools.prettyPercent(*x) for x in (
                (self.mExonsOverlapping1, self.mExons1),
                (self.mExonsOverlapping2, self.mExons2),
                (self.mExonsUnique1, self.mExons1),
                (self.mExonsUnique2, self.mExons2),
                (self.mBasesOverlapping1, self.mBases1),
                (self.mBasesOverlapping2, self.mBases2),
                (self.mBasesUnique1, self.mBases1),
                (self.mBasesUnique2, self.mBases2))])


class CounterTracks(Counter):

    def __init__(self, filename):
        self.mIndices = Bed.readAndIndex(IOTools.openFile(filename, "r"),
                                         per_track=True)

    def getTracks(self):
        return sorted(self.mIndices.keys())

    def _countIndices(self, idx_in, idx):
        '''count filename against idx.'''

        overlapping_genes = set()
        genes = set()

        # iterate over exons

        nexons, nexons_overlapping = 0, 0
        nbases, nbases_overlapping = 0, 0
        for contig, ix in idx_in.items():

            # note: add a findall function to ncl
            for start, end, value in ix.find(0, 1000000000):
                nexons += 1
                nbases += end - start

                try:
                    intervals = list(idx[contig].find(start, end))
                except KeyError:
                    continue

                if len(intervals) == 0:
                    continue

                nexons_overlapping += 1
                counts = numpy.zeros(end - start, numpy.int)
                for other_start, other_end, other_value in intervals:
                    for x in range(max(start, other_start) - start, min(end, other_end) - start):
                        counts[x] += 1
                nbases_overlapping += sum([1 for x in counts if x > 0])

        return nexons, nexons_overlapping, nbases, nbases_overlapping

    def count(self, filename, track):
        """count overlap between two gtf files."""

        E.info("counting started for %s versus %s" % (filename, track))

        (self.mExons1, self.mExonsOverlapping1,
         self.mBases1, self.mBasesOverlapping1 ) = \
            self._count(filename, self.mIndices[track])

        self.mExonsUnique1 = self.mExons1 - self.mExonsOverlapping1
        self.mBasesUnique1 = self.mBases1 - self.mBasesOverlapping1

        idx = self.buildIndex(filename)

        # count index against index
        (self.mExons2, self.mExonsOverlapping2,
         self.mBases2, self.mBasesOverlapping2 ) = \
            self._countIndices(self.mIndices[track], idx)

        self.mExonsUnique2 = self.mExons2 - self.mExonsOverlapping2
        self.mBasesUnique2 = self.mBases2 - self.mBasesOverlapping2


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id: diff_bed.py 2866 2010-03-03 10:18:49Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-u", "--update", dest="filename_update", type="string",
                      help="if filename is given, previous results will be read from there and only changed sets will be computed [default=%default].")

    parser.add_option("-p", "--pattern-identifier", dest="pattern_id", type="string",
                      help="pattern to convert a filename to an id [default=%default].")

    parser.add_option("-t", "--tracks", dest="tracks", action="store_true",
                      help="compare files against all tracks in the first file [default=%default]")

    parser.set_defaults(
        filename_update=None,
        pattern_id="(.*).bed",
        tracks=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) < 2:
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

    pattern_id = re.compile(options.pattern_id)

    def getTitle(x):
        try:
            return pattern_id.search(x).groups()[0]
        except AttributeError:
            return x

    ncomputed, nupdated = 0, 0

    if options.tracks:
        counter = CounterTracks(args[0])
        options.stdout.write("set1\tset2\t%s\n" % counter.getHeader())
        for filename in args[1:]:
            title1 = getTitle(filename)
            for title2 in counter.getTracks():

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

                counter.count(filename, title2)
                options.stdout.write(
                    "%s\t%s\t%s\n" % ((title1, title2, str(counter))))
                ncomputed += 1
    else:
        counter = Counter()
        options.stdout.write("set1\tset2\t%s\n" % counter.getHeader())

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

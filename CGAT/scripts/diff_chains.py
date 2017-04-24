"""diff_chains.py - compare to chain formatted files
=================================================

:Tags: Genomics GenomeAlignment CHAIN Comparison

Purpose
-------

Compare two genomic alignment files and calculate statistics from the
comparison.

Documentation
--------------

Operates on two `chain
<https://genome.ucsc.edu/goldenPath/help/chain.html>`_ formatted
files.

Outputs a table with the following columns:

+-----------+-------------------------------+
|*Column*   |*Content*                      |
+-----------+-------------------------------+
|contig1    |contig name                    |
+-----------+-------------------------------+
|contig2    |contig name                    |
+-----------+-------------------------------+
|strand     |strand                         |
+-----------+-------------------------------+
|mapped1    |mapped residues                |
+-----------+-------------------------------+
|identical1 |identically mapped residues    |
+-----------+-------------------------------+
|different1 |differently mapped residues    |
+-----------+-------------------------------+
|unique1    |residues mapped only from set1 |
+-----------+-------------------------------+
|pmapped1   |percentage of mapped residues  |
+-----------+-------------------------------+
|pidentical1|percentage of identically      |
|           |mapped residues                |
+-----------+-------------------------------+
|pdifferent1|percentage of differently      |
|           |mapped residues                |
+-----------+-------------------------------+

Similar columns exist for data set 2

Usage
-----

Example::

   cgat diff_chains.py hg19ToMm10v1.chain.over.gz hg19ToMm10v2.chain.over.gz

This will compare the locations that regions within the genome hg19
map to between two different mappings to the genome mm10.

Type::

   python diff_chains.py --help

for command line help.

Command line options
--------------------

"""

import sys
import collections

import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import alignlib_lite


def chain_iterator(infile):
    lines = []
    for line in infile:

        if line.startswith("#"):
            continue
        if line.strip() == "":
            continue
        if line.startswith("chain"):
            if lines:
                yield lines
            lines = []
        lines.append(line)

    yield lines


def validateChain(infile):
    '''validate a chain file.

    No overlapping target coordinates.
    '''

    pairs_t2q = collections.defaultdict(alignlib_lite.py_makeAlignmentBlocks)
    pairs_q2t = collections.defaultdict(alignlib_lite.py_makeAlignmentBlocks)

    for lines in chain_iterator(infile):

        (_,
         _,
         target_contig,
         target_length,
         target_strand,
         target_start,
         target_end,
         query_contig,
         query_length,
         query_strand,
         query_start,
         query_end,
         alignment_id) = lines[0][:-1].split()

        (query_start,
         query_end,
         query_length,
         target_start,
         target_end,
         target_length) = \
            [int(x) for x in
             (query_start,
              query_end,
              query_length,
              target_start,
              target_end,
              target_length)]

        # target_strand is always positive
        assert(target_strand == "+")

        map_target2query = pairs_t2q[target_contig]
        map_query2target = pairs_q2t[query_contig]

        qstart, tstart = query_start, target_start

        for line in lines[1:-1]:
            size, dt, dq = [int(x) for x in line[:-1].split()]
            map_target2query.addDiagonal(tstart,
                                         tstart + size,
                                         0)
            map_query2target.addDiagonal(qstart,
                                         qstart + size,
                                         0)

            qstart += size + dq
            tstart += size + dt

        size = int(lines[-1][:-1])
        map_target2query.addDiagonal(tstart,
                                     tstart + size,
                                     0)
        map_query2target.addDiagonal(qstart,
                                     qstart + size,
                                     0)

    try:
        x = map_query2target.mapRowToCol(0)
    except RuntimeError:
        E.info("query is not unique - this is ok.")

    try:
        x = map_target2query.mapRowToCol(0)
    except RuntimeError:
        E.info("target is not unique")
        return False

    return True


def buildPairs(infile):
    '''read a chain file.

    build target2query alignments.
    The target is always on the positive strand.
    '''
    pairs = collections.defaultdict(alignlib_lite.py_makeAlignmentBlocks)

    def chain_iterator(infile):
        lines = []
        for line in infile:

            if line.startswith("#"):
                continue
            if line.strip() == "":
                continue
            if line.startswith("chain"):
                if lines:
                    yield lines
                lines = []
            lines.append(line)

        yield lines

    for lines in chain_iterator(infile):

        (_,
         _,
         target_contig,
         target_length,
         target_strand,
         target_start,
         target_end,
         query_contig,
         query_length,
         query_strand,
         query_start,
         query_end,
         alignment_id) = lines[0][:-1].split()

        (query_start,
         query_end,
         query_length,
         target_start,
         target_end,
         target_length) = \
            [int(x) for x in
             (query_start,
              query_end,
              query_length,
              target_start,
              target_end,
              target_length)]

        # target_strand is always positive
        assert(target_strand == "+")

        map_target2query = pairs[(target_contig, query_contig, query_strand)]

        qstart, tstart = query_start, target_start

        for line in lines[1:-1]:
            size, dt, dq = [int(x) for x in line[:-1].split()]
            map_target2query.addDiagonal(tstart,
                                         tstart + size,
                                         qstart - tstart)
            qstart += size + dq
            tstart += size + dt

        size = int(lines[-1][:-1])

        map_target2query.addDiagonal(tstart,
                                     tstart + size,
                                     qstart - tstart)

    return pairs

DiffResult = collections.namedtuple(
    "DiffResult", "total same different unique")


def compareChains(pairs1, pairs2):
    '''compare chains in pairs1 versus those in pairs2'''

    result = {}
    for key1, chain1 in pairs1.items():
        E.debug("comparing %s" % str(key1))

        ntotal = chain1.getNumAligned()

        if key1 not in pairs2:
            result[key1] = DiffResult._make((ntotal, 0, 0, ntotal))
            continue

        chain2 = pairs2[key1]
        nsame = alignlib_lite.py_getAlignmentIdentity(
            chain1, chain2, alignlib_lite.py_RR)
        noverlap = alignlib_lite.py_getAlignmentOverlap(
            chain1, chain2, alignlib_lite.py_RR)
        ndifferent = noverlap - nsame
        nunique = ntotal - noverlap

        result[key1] = DiffResult._make((ntotal, nsame, ndifferent, nunique))

    return result


def outputMismatches(pairs1, pairs2,
                     output_mismatches=False,
                     output_unique=False,
                     output_matches=False):
    '''output mismatches.

    This is a very slow operation.
    '''

    outfile = sys.stdout

    for key1, chain1 in pairs1.items():
        E.debug("comparing %s" % str(key1))

        if key1 not in pairs2:
            outfile.write("%s\t%s\t%s\t%i\t%i\t-1\n" %
                          (key1 + (chain1.getRowFrom(), chain2.getRowTo())))
            continue

        chain2 = pairs2[key1]
        reg_start = chain1.getRowFrom()
        for pos in range(chain1.getRowFrom(), chain2.getRowTo()):
            x = chain1.mapRowToCol(pos)
            y = chain2.mapRowToCol(pos)

            if x == y:
                if output_matches:
                    outfile.write("%s\t%s\t%s\t%i\t%i\t%i\t%s\n" %
                                  (key1 + (pos, x, y, "=")))
                continue

            mismatch = x != -1 and y != -1

            if mismatch:
                if output_mismatches:
                    outfile.write("%s\t%s\t%s\t%i\t%i\t%i\t%i\n" %
                                  (key1 + (pos, x, y, x - y)))
            else:
                if output_unique:
                    outfile.write("%s\t%s\t%s\t%i\t%i\t%i\t%s\n" %
                                  (key1 + (pos, x, y, "-")))


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: chain2psl.py 2899 2010-04-13 14:37:37Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--output-mismatches", dest="output_mismatches", action="store_true",
                      help="output mismatches [%default]")

    parser.add_option("-a", "--output-matches", dest="output_matches", action="store_true",
                      help="output matches [%default]")

    parser.add_option("-u", "--output-unique", dest="output_unique", action="store_true",
                      help="output unique positions [%default]")

    parser.add_option("-r", "--restrict", dest="restrict", type="string",
                      help="restrict analysis to a chromosome pair (chr1:chr1:+) [%default]")

    parser.set_defaults(
        output_mismatches=False,
        output_unique=False,
        restrict=None
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) != 2:
        raise ValueError("expected two chain files")

    filename_chain1, filename_chain2 = args

    E.info("validating chain 1")
    if not validateChain(IOTools.openFile(filename_chain1)):
        E.warn("validation failed - exiting")
        return 1

    E.info("validating chain 2")
    if not validateChain(IOTools.openFile(filename_chain2)):
        E.warn("validation failed - exiting")
        return 1

    E.info("building pairs for %s" % filename_chain1)
    pairs1 = buildPairs(IOTools.openFile(filename_chain1))
    E.info("read %i pairs" % len(pairs1))

    E.info("building pairs for %s" % filename_chain2)
    pairs2 = buildPairs(IOTools.openFile(filename_chain2))
    E.info("read %i pairs" % len(pairs2))

    if options.restrict:
        restrict = tuple(options.restrict.split(":"))
        pairs1 = {restrict: pairs1[restrict]}
        pairs2 = {restrict: pairs2[restrict]}

    E.info("comparing 1 -> 2")
    comparison1 = compareChains(pairs1, pairs2)
    E.info("comparing 2 -> 1")
    comparison2 = compareChains(pairs2, pairs1)

    all_keys = sorted(list(set(list(comparison1.keys()) + list(comparison2.keys()))))

    outfile = options.stdout
    headers = ("mapped", "identical", "different", "unique")
    outfile.write("contig1\tcontig2\tstrand\t%s\t%s\t%s\t%s\n" %
                  (
                      "\t".join(["%s1" % x for x in headers]),
                      "\t".join(["p%s1" % x for x in headers]),
                      "\t".join(["%s2" % x for x in headers]),
                      "\t".join(["p%s2" % x for x in headers])))

    totals = E.Counter()

    for key in all_keys:
        outfile.write("%s\t%s\t%s" % key)

        if key in comparison1:
            c = comparison1[key]
            outfile.write("\t%i\t%i\t%i\t%i\t" %
                          (c.total, c.same, c.different, c.unique))
            outfile.write(
                "\t".join([IOTools.prettyPercent(x, c.total) for x in c]))

            totals.total1 += c.total
            totals.same1 += c.same
            totals.different1 += c.different
            totals.unique1 += c.unique
        else:
            outfile.write("\t%i\t%i\t%i\t%i\t" % (0, 0, 0, 0))
            outfile.write("\t%i\t%i\t%i\t%i" % (0, 0, 0, 0))

        if key in comparison2:
            c = comparison2[key]
            outfile.write("\t%i\t%i\t%i\t%i\t" %
                          (c.total, c.same, c.different, c.unique))
            outfile.write(
                "\t".join([IOTools.prettyPercent(x, c.total) for x in c]))

            totals.same2 += c.same
            totals.total2 += c.total
            totals.different2 += c.different
            totals.unique2 += c.unique
        else:
            outfile.write("\t%i\t%i\t%i\t%i\t" % (0, 0, 0, 0))
            outfile.write("\t%i\t%i\t%i\t%i" % (0, 0, 0, 0))

        outfile.write("\n")

    outfile.write("total\ttotal\t.\t")
    outfile.write("\t".join(map(str, (totals.total1,
                                      totals.same1,
                                      totals.different1,
                                      totals.unique1,
                                      IOTools.prettyPercent(
                                          totals.total1, totals.total1),
                                      IOTools.prettyPercent(
                                          totals.same1, totals.total1),
                                      IOTools.prettyPercent(
                                          totals.different1, totals.total1),
                                      IOTools.prettyPercent(
                                          totals.unique1, totals.total1),
                                      totals.total2,
                                      totals.same2,
                                      totals.different2,
                                      totals.unique2,
                                      IOTools.prettyPercent(
                                          totals.total2, totals.total2),
                                      IOTools.prettyPercent(
                                          totals.same2, totals.total2),
                                      IOTools.prettyPercent(
                                          totals.different2, totals.total2),
                                      IOTools.prettyPercent(
                                          totals.unique2, totals.total2),
                                      ))) + "\n")

    # output mismapped residues
    if options.output_mismatches or options.output_unique:
        outputMismatches(pairs1, pairs2,
                         output_mismatches=options.output_mismatches,
                         output_unique=options.output_unique,
                         output_matches=options.output_matches,
                         )

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

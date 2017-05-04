'''
diamond2counts.py - count alignments to reference
===================================================

:Tags: Python

Purpose
-------

Count the number of alignments to each reference in outfmt6.

Counts are based on various options specified by --method.

best       This will take the best alignment as judged by the highest
           bitscore.




TODO::
Add additional options

Usage
-----

Example::

   python diamond2counts.py

Type::

   python diamond2counts.py --help

for command line help.

Command line options
--------------------

'''

import sys

import CGAT.Experiment as E
from CGAT.Diamond import *
import collections
import CGAT.IOTools as IOTools


def readCogMap(cog_map):
    '''
    return a dictionary mapping gene to cog
    '''
    gene2cog = {}
    for line in IOTools.openFile(cog_map):
        data = line[:-1].split("\t")
        gene2cog[data[0]] = data[1]
    return gene2cog


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("best", None),
                      help="method for determing what to count")

    parser.add_option("--sum-cog", dest="sum_cog", action="store_true",
                      help="sum counts over functions (COGs) in --cog-map")

    parser.add_option("--evaluate-cog",
                      dest="evaluate_cog",
                      action="store_true",
                      help="""output the percent of
                              alignments for each read = best hit""")

    parser.add_option("--cog-map", dest="cog_map", type="string",
                      help="file with gene to cog map")

    parser.add_option("-n", "--nsamples", dest="nsamples", type="int",
                      help="""number of queries to evaluate-
                              will take the first n in the file""")

    parser.set_defaults(method=None,
                        sum_cog=False,
                        evaluate_cog=False,
                        cog_map=None,
                        nsamples=10000)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.evaluate_cog:
        assert options.cog_map, """must specify an annotation
                                   mapping gene to function (COG)"""
        assert not options.method, """evaluation performed
                                      in the absence of counting"""

        E.info("reading gene to function (COG) map %s" % options.cog_map)
        gene2cog = readCogMap(options.cog_map)
        E.info("loaded gene to function (COG) map")

        E.info("retrieving alignment data")
        options.stdout.write("query\tpbest\tnalignments\n")

        c = 0
        for alignments in query_iterator(alignment_iterator(options.stdin)):
            c += 1
            scores = []
            if c <= options.nsamples:
                for alignment in alignments:
                    scores.append(alignment.score)
                    best = max(scores)
                    best_alignments = [
                        x for x in alignments if x.score == best]
                    if len(best_alignments) > 1:
                        best_alignments = random.sample(best_alignments, 1)
                    best_alignment = best_alignments[0]
                    best_cog = gene2cog[best_alignment.ref]
                pbest = float(len(
                    [gene2cog[x.ref]
                     for x in alignments
                     if gene2cog[x.ref] == best_cog])) / len(alignments) * 100
                nalignments = len(alignments)
                options.stdout.write(
                    "\t".join(map(
                        str, [alignments[0].qid,
                              pbest,
                              nalignments])) + "\n"
                )
            else:
                break
        return

    # container for counts
    counts = collections.defaultdict(int)
    E.info("counting alignments")
    assert options.method, "required option --method"
    if options.method == "best":
        if options.sum_cog:
            E.warn("""summing over functions (COGS)
                      will remove genes with no annotations
                      and those with multiple COG assignments""")
            assert options.cog_map, """a mapping between gene and
                                       function (COG) is required"""

            E.info("""reading gene to function (COG) mapping from %s"""
                   % options.cog_map)
            gene2cog = readCogMap(options.cog_map)
            E.info("loaded gene to function (COG) mapping")

            E.info("summing functional assignments")
            query_it = query_iterator(alignment_iterator(options.stdin))
            for best in best_alignment_iterator(query_it):
                cog = gene2cog[best.ref]
                # removing uassigned or multiple assignments
                if cog == "unknown" or cog.find(";") != -1:
                    continue
                counts[cog] += 1
        else:
            E.info("counting best alignments")
            query_it = query_iterator(alignment_iterator(options.stdin))
            for best in best_alignment_iterator(query_it):
                counts[best.ref] += 1
        E.info("finished counting")

        E.info("writing results")
        options.stdout.write("ref\tcount\n")
        for ref, count in sorted(counts.items()):
            options.stdout.write("\t".join([ref, str(count)]) + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

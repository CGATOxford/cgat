'''
psl2fasta.py - output alignment pairs from psl formatted alignments
===================================================================

:Tags: Python

Purpose
-------

analyze sequence pairs from a psl formatted table.

The sequences are assumed to be nucleotide sequences.

Usage
-----

Example::

   python psl2fasta.py --help

Type::

   python psl2fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys

import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Blat as Blat
import CGAT.Genomics as Genomics

import alignlib_lite



def getAlignmentFull(m, q, t, options):
    """print alignment with gaps in both query and target."""
    a = alignlib_lite.py_AlignmentFormatExplicit(
        m, alignlib_lite.py_makeSequence(q), alignlib_lite.py_makeSequence(t))
    return a.mRowAlignment, a.mColAlignment



def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--query-psl-file", dest="filename_query", type="string",
                      help="fasta filename with queries.")

    parser.add_option("--target-psl-file", dest="filename_target", type="string",
                      help="fasta filename with target.")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=(
                          "full", "pileup-query", "pileup-target", "gapless"),
                      help="method to use for constructing the alignment [%default].")

    parser.add_option("--forward-query", dest="forward_query", action="store_true",
                      help="reverse-complement sequences such that query is always on forward strand [%default]")

    parser.add_option("--target-prefix", dest="target_prefix", type="string",
                      help="prefix to use for target [%default].")

    parser.add_option("--query-prefix", dest="query_prefix", type="string",
                      help="prefix to use for query [%default].")

    parser.add_option("--id", dest="id", type="choice",
                      choices=("numeric", "query"),
                      help="choose type of identifier to use [%default]")

    parser.set_defaults(
        filename_query=None,
        filename_target=None,
        method="full",
        output_format_id="%06i",
        target_prefix="",
        query_prefix="",
        forward_query=False,
    )

    (options, args) = E.Start(parser)

    if options.filename_query:
        query = IndexedFasta.IndexedFasta(options.filename_query)

    if options.filename_target:
        target = IndexedFasta.IndexedFasta(options.filename_target)

    if options.method == "full":
        getAlignment = getAlignmentFull

    id = 0
    for match in Blat.iterator(options.stdin):
        if options.loglevel >= 2:
            options.stdout.write("# %s\n" % str(match))

        m = match.getMapQuery2Target()
        m.moveAlignment(-min(match.mQueryBlockStarts), -
                        min(match.mSbjctBlockStarts))
        q = query.getSequence(
            match.mQueryId, match.strand, match.mQueryFrom, match.mQueryTo)
        t = target.getSequence(
            match.mSbjctId, "+", match.mSbjctFrom, match.mSbjctTo)
        query_ali, sbjct_ali = getAlignment(m, q, t, options)

        if match.strand == "-" and options.forward_query:
            query_ali = Genomics.complement(query_ali)
            sbjct_ali = Genomics.complement(sbjct_ali)

        options.stdout.write(">%s%s:%s/%i-%i\n%s\n>%s%s:%s%s/%i-%i\n%s\n" %
                             (options.query_prefix,
                              options.output_format_id % id,
                              match.mQueryId, match.mQueryFrom, match.mQueryTo,
                              query_ali,
                              options.target_prefix,
                              options.output_format_id % id,
                              match.mSbjctId, match.strand,
                              match.mSbjctFrom, match.mSbjctTo,
                              sbjct_ali))
        id += 1

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

"""
bed2psl.py - convert a bed file to a psl file
=============================================

:Tags: Python

Purpose
-------

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Command line options
--------------------

"""

import sys
import CGAT.Experiment as E

import CGAT.Blat as Blat
import CGAT.Bed as Bed
import CGAT.IndexedFasta as IndexedFasta


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: bed2psl.py 2899 2010-04-13 14:37:37Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-q", "--query", dest="query", type="string",
                      help="sequence to use for query [default=%default].")

    parser.add_option("-t", "--target", dest="target", type="string",
                      help="sequence to use for target [default=%default].")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome.")

    parser.set_defaults(
        genome_file=None,
        query=None,
        target=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # do sth
    ninput, nskipped, noutput = 0, 0, 0

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)
    else:
        fasta = None

    psl = Blat.Match()

    for bed in Bed.iterator(options.stdin):

        ninput += 1

        start, end = bed.start, bed.end

        if "blockSizes" in bed:
            psl.mQueryId = bed["name"]
            blocksizes = [int(x) for x in bed["blockSizes"].split(",")[:-1]]
            sbjctblockstarts = [
                int(x) + start for x in bed["blockStarts"].split(",")[:-1]]
            strand = bed["strand"]
        else:
            psl.mQueryId = "%i" % ninput
            blocksizes = [end - start]
            sbjctblockstarts = [start, ]

            strand = "+"

        psl.mSbjctId = bed.contig
        psl.mSbjctFrom, psl.mSbjctTo = start, end
        psl.mQueryFrom, psl.mQueryTo = 0, end - start

        psl.mBlockSizes = blocksizes
        psl.mNBlocks = len(blocksizes)
        psl.strand = strand
        q, qp = [], 0
        for x in blocksizes:
            q.append(qp)
            qp += x

        psl.mQueryBlockStarts = q
        psl.mSbjctBlockStarts = sbjctblockstarts
        psl.mQueryLength = sum(psl.mBlockSizes)
        if fasta:
            psl.mSbjctLength = fasta.getLength(bed.contig)

        options.stdout.write("%s\n" % str(psl))
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

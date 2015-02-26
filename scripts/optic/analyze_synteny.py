##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
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
'''
optic/analyze_synteny.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python optic/analyze_synteny.py --help

Type::

   python optic/analyze_synteny.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E
import CGAT.Synteny as Synteny
import scipy
import scipy.stats
import CGAT.GTF as GTF


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/analyze_synteny.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("--method", dest="methods", type="string",
                      help="methods to apply [off-diagonal].")
    parser.add_option("--max-synteny-distance", dest="max_synteny_distance", type="int",
                      help="maximal synteny distance in genes for local rearrangments/deletions.")
    parser.add_option("--max-look-ahead", dest="max_look_ahead", type="int",
                      help="maximal distance in genes to continue a block.")
    parser.add_option("--filter-junk", dest="filter_junk", action="store_true",
                      help="do not take into account orthologs on junk chromosomes.")

    parser.set_defaults(
        separator="|",
        methods="counts",
        max_synteny_distance=5,
        max_look_ahead=1,
        filter_junk=False,
        output_pattern_gff="%s.gff",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if len(args) != 2:
        raise "please supply two files with synteny information."

    options.methods = options.methods.split(",")

    options.filename1, options.filename2 = args

    orthologs1 = Synteny.ReadOrthologs(open(options.filename1, "r"))
    orthologs2 = Synteny.ReadOrthologs(open(options.filename2, "r"))

    if options.loglevel >= 1:
        print "# read orthologs: %i/%i" % (len(orthologs1), len(orthologs2))

    # make sure, that there are no orthologs without a corresponding
    # entry in the other set.
    Synteny.CleanOrthologs(orthologs1, orthologs2,
                           filter_junk=options.filter_junk)

    if options.loglevel >= 1:
        print "# cleaned orthologs: %i/%i" % (len(orthologs1), len(orthologs2))

    # get sorted list of contigs and list of orthologs per contig
    sorted_contigs1, sorted_contigs2, contigs1, contigs2 = Synteny.SortAssignOrthologs(
        orthologs1, orthologs2)

    # set rank for orthologs within a contig
    Synteny.SetRankToPosition(orthologs1, contigs1, relative=False)
    Synteny.SetRankToPosition(orthologs2, contigs2, relative=False)

    # go through analysis methods
    for method in options.methods:

        outfile = options.stdout

        if method == "counts":

            outfile.write("contig1\tcontig2\tngenes\n")

            counts = Synteny.GetContigCounts(orthologs1, orthologs2)

            values = counts.items()

            values.sort(lambda x, y: cmp(x[1], y[1]))
            values.reverse()

            for key, value in values:
                outfile.write("%s\t%s\t%i\n" % (key[0], key[1], value))

        elif method == "contents":

            contents = Synteny.GetContigContents(orthologs1, orthologs2)

            values = contents.items()

            values.sort(lambda x, y: cmp(len(x[1]), len(y[1])))
            values.reverse()

            for key, value in values:
                outfile.write("%s\t%s\t%s\n" %
                              (key[0], key[1], str(list(set(value)))))

        elif method in ("blocks", "blocks-gff"):

            # compute synteny blocks
            blocks = Synteny.GetSyntenyBlocks(orthologs1, orthologs2,
                                              max_synteny_distance=options.max_synteny_distance,
                                              max_look_ahead=options.max_look_ahead,
                                              loglevel=options.loglevel)

            if method == "blocks":
                outfile.write(Synteny.Block().getHeader() + "\n")

                for b in range(len(blocks)):
                    block = blocks[b]
                    outfile.write(str(block) + "\n")

            elif method == "blocks-gff":

                def writeGFF(blocks, first, filename):

                    outfile.write("writing gff entries to %s\n" % filename)

                    outfile_gff = open(filename, "w")

                    entry = GTF.Entry()
                    entry.source = "gpipe"
                    entry.feature = "synteny"

                    for b in range(len(blocks)):
                        block = blocks[b]
                        if first:
                            entry.name = block.contig1
                            entry.start = block.mFrom1
                            entry.end = block.mTo1
                        else:
                            entry.name = block.contig2
                            entry.start = block.mFrom2
                            entry.end = block.mTo2

                        entry.info = "Block=%i" % block.mBlockId

                        outfile_gff.write(str(entry) + "\n")

                    outfile_gff.close()

                writeGFF(blocks, True, options.output_pattern_gff % "species1")
                writeGFF(blocks, False, options.output_pattern_gff %
                         "species2")

        elif method == "members":

            # compute synteny blocks
            blocks = Synteny.GetSyntenyBlocks(orthologs1, orthologs2,
                                              max_synteny_distance=options.max_synteny_distance,
                                              max_look_ahead=options.max_look_ahead,
                                              loglevel=options.loglevel)

            outfile.write(Synteny.Block().getHeader() + "\n")

            for block in blocks:
                outfile.write(">" + str(block.mBlockId) + "\n")
                outfile.write("\n".join(map(str, block.mMembers1)) + "\n")
                outfile.write("\n".join(map(str, block.mMembers2)) + "\n")

        elif method == "breakers":
            # output information about synteny breakers
            blocks = Synteny.GetSyntenyBlocks(orthologs1, orthologs2,
                                              max_synteny_distance=options.max_synteny_distance,
                                              max_look_ahead=options.max_look_ahead,
                                              loglevel=options.loglevel)

            ntotal_same, ntotal_different, nassigned = 0, 0, 0
            total_block_sizes = []
            total_breakers_ortholog_ids = {}

            for block in blocks:

                nassigned += len(block.mMembers2)

                if block.mBreakers2:
                    nsame = 0
                    ndifferent = 0
                    block_sizes = []
                    block.mBreakers2.sort(lambda x, y: cmp(x.mRank, y.mRank))
                    last_rank = None
                    n = 0

                    for b in block.mBreakers2:
                        if b.contig == block.contig2:
                            nsame += 1
                        else:
                            ndifferent += 1

                        if last_rank and b.mRank - last_rank > 1:
                            block_sizes.append(n)
                            n = 0

                        n += 1
                        last_rank = b.mRank

                        if b.mOrthologId not in total_breakers_ortholog_ids:
                            total_breakers_ortholog_ids[b.mOrthologId] = 0

                        total_breakers_ortholog_ids[b.mOrthologId] += 1

                    block_sizes.append(n)

                    ntotal_same += nsame
                    ntotal_different += ndifferent
                    outfile.write(">%i %s:%i..%i %s:%i..%i\n" % (
                        block.mBlockId, block.contig1, block.mFrom1, block.mTo1, block.contig2, block.mFrom2, block.mTo2))
                    outfile.write("stats\t%i\t%s\t%i\t%i\t%i\t%s\n" % (
                        block.mBlockId, block.contig2, nsame, ndifferent, len(block_sizes), ";".join(map(str, block_sizes))))
                    outfile.write("\n".join(map(str, block.mBreakers2)) + "\n")

                    total_block_sizes += block_sizes

            r = range(1, max(total_block_sizes) + 1)
            h = scipy.stats.histogram2(total_block_sizes, r)

            outfile.write(
                "# Histogram of block sizes of synteny breakers in genome2.\n")
            outfile.write("bin\tcounts\n")
            for bin, val in zip(r, h):
                outfile.write("%i\t%i\n" % (bin, val))

            ntotal = ntotal_same + ntotal_different

            outfile.write("""# Gene summary of synteny breakers (counts in genes).
# nsyn:         total number of genes assigned to synteny blocks in genome2
# ntotal:       total number of synteny breakers in genome2.
# ptotal:       percentage of synteny breakers.
# nsame:        number of synteny breakers in genome2 on the same contig.
# psame:        percentage of synteny breakers in genome2 on the same contig.
# ndiff:        number of synteny breakers in genome2 on a different contig.
# pdiff:        percentage of synteny breakers in genome2 on a different contig.
# rsd:          ratio of same to different contigs
""")

            outfile.write(
                "gene_total\tnsyn\tntotal\tptotal\tnsame\tpsame\tndiff\tpdiff\trsd\n")
            outfile.write("gene_total\t%i\t%i\t%5.2f\t%i\t%5.2f\t%i\t%5.2f\t%5.2f\n" % (nassigned,
                                                                                        ntotal, float(
                                                                                            ntotal) * 100 / nassigned,
                                                                                        ntotal_same, float(
                                                                                            ntotal_same) * 100 / nassigned,
                                                                                        ntotal_different, float(
                                                                                            ntotal_different) * 100 / nassigned,
                                                                                        float(ntotal_same) / float(ntotal_different)))

            outfile.write("""# Ortholog summary of synteny breakers (counts in ortholog clusters).
# nsyn:         total number of orthologs assigned to synteny blocks in genome2
# ntotal:       total number of synteny breakers in genome2.
# ptotal:       percentage of synteny breakers.
""")

            nassigned = len(orthologs2)
            ntotal = len(total_breakers_ortholog_ids)
            outfile.write("ortholog_total\tnsyn\tntotal\tptotal\n")
            outfile.write("ortholog_total\t%i\t%i\t%5.2f\n" % (nassigned,
                                                               ntotal, float(ntotal) * 100 / nassigned))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

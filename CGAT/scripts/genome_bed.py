"""=================================================================
genome_bed.py - Create a bed file tiling a genome from a fai file
=================================================================

:Tags: Python

This program takes an indexed genome and creates windows of a certain
size.

It also takes two input parameters: the window/tile size (bases) and
the shift size.  By default the shift size is equal to the window
size.  The default window size is 1000.

Usage
=====

   python genome_bed -g <genome.fai> -o <output.bed> -w window size -s shift size


Command line options
====================

"""


import CGAT.Experiment as E
import sys
import CGAT.IndexedFasta as IndexedFasta


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option(
        "-g", "--genome-file", dest="genome_file", type="string",
        help="filename with Samtools indexed genome [default=%default].")
    parser.add_option(
        "-w", "--window-size", dest="window", type="int",
        help="Window size for tiling [default=%default].")
    parser.add_option(
        "-s", "--shift-size", dest="shift", type="int",
        help="Window shift for tiling [default=%default].")

    parser.set_defaults(
        genome_file=None,
        window=1000,
        shift=1000,
        output_file=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    ninput, nunchanged, nchanged = 0, 0, 0

    # Open input file
    E.info("Opening input file: %s" % options.genome_file)
    fasta = IndexedFasta.IndexedFasta(options.genome_file)
    contigs = fasta.getContigSizes(with_synonyms=False)

    # Open output file
    bed = options.stdout
    shift = options.shift

    # Loop over input files and convert to soft clipped
    nwindows = 0
    ncontigs = 0
    for contig, stop in contigs.items():
        ncontigs += 1
        i = 0
        while (i < stop):
            j = min(i + options.window, stop)
            #bed.write( """%(contig)s\t%(i)i\t%(j)i\t%(contig)s:%(i)i..%(j)i\n""" % locals() )
            bed.write( """%(contig)s\t%(i)i\t%(j)i\n""" % locals() )
            nwindows += 1
            i += shift

    # Report statistics
    E.info("ncontigs=%i, nwindows=%i" % (ncontigs, nwindows))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

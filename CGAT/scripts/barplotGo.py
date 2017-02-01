'''
barplotGo.py
=============================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes the output from GO.py (section = results) and produces a barplot of the 
enrichment results. The plotting is performed using Rpy.

Usage
-----

Example::

   python barplotGo.py --help

Type::

   python barplotGo.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
from rpy2.robjects import r as R


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-i", "--infile", dest="infile", type="string",
                      help="supply infile name")
    parser.add_option("-o", "--outfile", dest="outfile", type="string",
                      help="supply outfile name")
    parser.add_option("-t", "--height", dest="height", type="int",
                      help="height of plotting window")
    parser.add_option("-w", "--width", dest="width", type="int",
                      help="number of bars")
    parser.add_option("-n", "--number", dest="number", type="int",
                      help="width of plotting window")
    parser.add_option("-c", "--colour", dest="colour", type="choice", choices=("red",
                                                                               "purple", "blue", "grey", "black", "pink", "green"), help="colour choices")
    parser.add_option(
        "--side", dest="side", action="store_true", help="set to draw sideways barplot")
    parser.add_option("-m", "--plot-margins", dest="margins",
                      type="string", help="set margins - follows R convention for mai")
    parser.add_option("--cexn", dest="cex_names", type="float",
                      help="set the size of the names on the barplot")
    parser.add_option("--cexa", dest="cex_axis", type="float",
                      help="set the size of the axis labels on the barplot")

    parser.set_defaults(
        height=500, width=500, number=5, colour="blue", margins="1,1,1,1", cex_names=1, cex_axis=1)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # don't plot more than 20
    if options.number > 20:
        sys.exit(0)

    inf = options.infile
    outf = options.outfile

    R('''data <- read.table("%s", sep = "\t", header = T, stringsAsFactors = F)''' %
      inf)
    R('''pdf("%s", height = %i, width = %i)''' %
      (outf, options.height, options.width))

    R('''ratio <- as.matrix(data$ratio[c(1:%i)])''' % options.number)
    R('''par(mai = c(%s))''' % options.margins)
    if options.side:
        horiz = "horiz = T"
    else:
        horiz = "horiz = F"
    R('''barplot(ratio
                     , beside = T
                     , names = data$description[c(1:%i)]
                     , col = "%s"
                     , las = 2
                     , %s
                     , cex.names = %f
                     , cex.axis = %f)''' % (options.number, options.colour, horiz, options.cex_names, options.cex_axis))
    R["dev.off"]()

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

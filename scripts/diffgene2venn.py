'''
diffgene2venn.py - template for CGAT scripts
====================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------
A plotting script for venn diagrams for overlapping differentially
expressed genes.  Uses the R VennDiagram library.

Usage
-----
Input are results of differential analysis, plots are written in png format to
a user supplied directory.  Can be use to plot up to 5 overlapping sets

Example::

   python diffgene2venn.py --alpha=0.01
                           --file-list=file1.txt,file2.txt,file3.txt
                           --output-directory=images.dir

Type::

   python diffgene2venn.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.Timeseries as TS


def main(argv=None):
    """script main.

parses command line options in sys.argv, unless *argv* is given.
"""

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--alpha", dest="alpha", type="string",
                      help="false positive rate for differentially"
                      " expressed genes")

    parser.add_option("--file-list", dest="infiles", type="string",
                      help="comma separated list of input files")

    parser.add_option("--output-directory", dest="out_dir", type="string",
                      help="output directory for png images")

    # add common options (-h/--help, ...) and parse command line

    (options, args) = E.Start(parser, argv=argv)

    infiles = options.infiles.split(",")
    TS.genSigGenes(file_list=infiles,
                   alpha=float(options.alpha),
                   out_dir=options.out_dir)

    # Write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

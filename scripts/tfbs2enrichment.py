'''
tfbs2enrichment.py - template for CGAT scripts
====================================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. Runs transcription factor motif enrichment using the
Transfac(R) database.  Signficance testing is currently
carried out by Fisher's Exact test based on a background
gene set matched for %CG content.

Currently only for use in pipeline_transfacmatch.py, this is
not a stand-alone script.

Usage
-----

.. Example use case

Example::

   python tfbs2enrichment.py

Type::

   python tfbs2enrichment.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGATPipelines.PipelineTransfacMatch as PipelineTFM


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--foreground", dest="foreground", type="string",
                      help="foreground file")

    parser.add_option("--background", dest="background", type="string",
                      help="matched background file")

    parser.add_option("--database", dest="database", type="string",
                      help="PATH to sqlite database containing Transfac Match"
                      " output.")

    parser.add_option("--match-table", dest="match_table", type="string",
                      help="tablename containing Transfac Match output")

    parser.add_option("--outfile", dest="outfile", type="string",
                      help="output file name")

    parser.add_option("--direction", dest="direction", type="string",
                      help="direction to test for significance in Fisher's "
                      "exact test. Default=2-tailed")

    parser.add_option("--geneset-header", dest="genesets", type="string",
                      help="geneset header")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    PipelineTFM.testSignificanceOfMatrices(options.background,
                                           options.foreground,
                                           options.database,
                                           options.match_table,
                                           options.outfile,
                                           options.genesets,
                                           options.direction)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

'''
distance2merge.py
=============================================

:Tags: Python

Purpose
-------

merge distance matrices from split jobs into a single distance
matrix.  Merge in the correct order based on the filename sub.
Only used for parallelising distance matrix calculations, not
designed as a stand-alone script.  See usage in pipeline_timeseries.py
for details

Usage
-----

Example::

   python distance2merge.py

Type::

   python distance2merge.py --help

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

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--outfile", dest="outfile", type="string",
                      help="output filename")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infiles = argv[-1]

    files_list = infiles.split(",")

    if not options.outfile:
        outfile = options.stdout
    else:
        outfile = options.outfile

    TS.mergeFiles(file_list=files_list,
                  outfile=outfile)

    # write footer and output benchmark information
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

'''
table2bed.py
=============================================

:Tags: Python

Purpose
-------

Not generic. Takes a tab delimited table and converts to bed file

Usage
-----

Example::

   python table2bed.py --help

Type::

   python table2bed.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--table", dest="table", type="string",
                      help="supply input table name")
    parser.add_option("-o", "--outfile", dest="outfile", type="string",
                      help="supply output file name")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    outfile = IOTools.openFile(options.outfile, "w")
    for line in IOTools.openFile(options.table).readlines():
        contig = line.split("\t")[0].split('"')[1]
        outfile.write(
            "\t".join((contig, line.split("\t")[1], line.split("\t")[2])))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

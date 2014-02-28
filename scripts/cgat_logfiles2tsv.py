'''
cgat_logfiles2tsv.py - create summary from logfiles
===================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts are collection of logfiles and collates
summary information.

This script uses the ``# job finished`` tag that
is added by scripts using the module :mod:`Experiment`.

Usage
-----

Example::

   python cgat_logfiles2tsv.py --help

Type::

   python cgat_logfiles2tsv.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import tempfile
import subprocess
import optparse
import gzip
import glob

from types import *

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Logfile as Logfile


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: cgat_logfiles2tsv.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-g", "--glob", dest="glob_pattern", type="string",
                      help="glob pattern to use for collecting files [%default].")

    parser.add_option("-f", "--file-pattern", dest="file_pattern", type="string",
                      help="only check files matching this pattern [%default].")

    parser.add_option("-m", "--mode", dest="mode", type="choice",
                      choices=("file", "node"),
                      help="analysis mode [%default].")

    parser.add_option("-r", "--recursive", action="store_true",
                      help="recursively look for logfiles from current directory [%default].")

    parser.set_defaults(
        truncate_sites_list=0,
        glob_pattern="*.log",
        mode="file",
        recursive=False,
    )

    (options, args) = E.Start(parser)

    if args:
        filenames = args
    elif options.glob_pattern:
        filenames = glob.glob(options.glob_pattern)

    if len(filenames) == 0:
        raise "no files to analyse"

    if options.mode == "file":
        totals = Logfile.LogFileData()

        options.stdout.write("file\t%s\n" % totals.getHeader())

        for filename in filenames:
            if filename == "-":
                infile = sys.stdin
            elif filename[-3:] == ".gz":
                infile = gzip.open(filename, "r")
            else:
                infile = open(filename, "r")

            subtotals = Logfile.LogFileData()
            for line in infile:
                subtotals.add(line)

            infile.close()

            options.stdout.write("%s\t%s\n" % (filename, str(subtotals)))
            totals += subtotals

        options.stdout.write("%s\t%s\n" % ("total", str(totals)))

    elif options.mode == "node":

        chunks_per_node = {}

        rx_node = re.compile("# job started at .* \d+ on (\S+)")

        for filename in filenames:
            if filename == "-":
                infile = sys.stdin
            elif filename[-3:] == ".gz":
                infile = gzip.open(filename, "r")
            else:
                infile = open(filename, "r")

            data = Logfile.LogFileDataLines()

            for line in infile:

                if rx_node.match(line):
                    node_id = rx_node.match(line).groups()[0]
                    data = Logfile.LogFileDataLines()
                    if node_id not in chunks_per_node:
                        chunks_per_node[node_id] = []
                    chunks_per_node[node_id].append(data)
                    continue

                data.add(line)

        options.stdout.write("node\t%s\n" % data.getHeader())
        total = Logfile.LogFileDataLines()

        for node, data in sorted(chunks_per_node.items()):
            subtotal = Logfile.LogFileDataLines()
            for d in data:
                # options.stdout.write( "%s\t%s\n" % (node, str(d) ) )
                subtotal += d

            options.stdout.write("%s\t%s\n" % (node, str(subtotal)))

            total += subtotal

        options.stdout.write("%s\t%s\n" % ("total", str(total)))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

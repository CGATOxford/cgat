'''
mali_phylip2fasta.py - 
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

   python mali_phylip2fasta.py --help

Type::

   python mali_phylip2fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys
import string
import re
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    line = sys.stdin.readline()
    num_lines, width = map(int, re.search("(\d+)\s+(\d+)", line[:-1]).groups())

    alignment = []

    for x in range(num_lines):
        line = sys.stdin.readline()
        id, ali = re.search("^(\S+)\s+(.+)", line[:-1]).groups()
        alignment.append((id, [ali]))

    while 1:
        line = sys.stdin.readline()
        if not line:
            break
        for x in range(num_lines):
            line = sys.stdin.readline()
            alignment[x][1].append(line[:-1])

    for x in range(num_lines):
        print ">%s\n%s" % (alignment[x][0],
                           re.sub("\s", "",
                                  string.join(alignment[x][1], "")))

    # write footer and output benchmark information.
    E.Stop()

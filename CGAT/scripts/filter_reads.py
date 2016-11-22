'''
filter_reads.py -
=============================================

:Author: ???
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import CGAT.Experiment as E


def readread(s):
    return [s.readline(), s.readline(), s.readline(), s.readline()]


def writeread(sread, s):
    for i in range(4):
        s.write(sread[i])


def filter_pattern(query, s, o):
    sread = readread(s)
    while (sread[0]):
        if (query.search(sread[1])):
            writeread(sread, o)
        sread = readread(s)


def filter_low_complexity(s, o):
    sread = readread(s)
    biased = 0
    low_complexity = 0
    total_reads = 0
    remaining_reads = 0
    while (sread[0]):
        total_reads += 1
        my_read = Seq(sread[1], generic_dna)
        a = my_read.count("A")
        c = my_read.count("C")
        t = my_read.count("T")
        g = my_read.count("G")
        seq_len = len(my_read)
        count_list = [a, c, t, g]
        if (count_list.count(0) < 2):
            if (max(a, c, t, g) / seq_len < 0.9):
                writeread(sread, o)
                remaining_reads += 1
            else:
                biased += 1
        else:
            low_complexity += 1
        sread = readread(s)
    removed = biased + low_complexity
    sys.stderr.write("Total reads processed: %s\\n" % total_reads)
    sys.stderr.write(r"Low complexity reads removed: %s\n" % low_complexity)
    sys.stderr.write(r"Biased reads removed: %s\n" % biased)
    sys.stderr.write(r"Total reads removed: %s\n" % removed)
    sys.stderr.write(r"Total reads remaining: %s\n" % remaining_reads)


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

    filter_low_complexity(sys.stdin, sys.stdout)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

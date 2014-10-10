'''
correlate_fasta_identifier.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Given two :term:`fasta` formatted files, substitute
identifiers in stream with those given in file with filename.

Usage
-----

Example::

   python correlate_fasta_identifier.py --help

Type::

   python correlate_fasta_identifier.py --help

for command line help.

Command line options
--------------------

'''
import sys
import string
import getopt

param_loglevel = 0

param_long_options = ["verbose=", "help", "version"]
param_short_options = "v:h"


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(
            sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print globals()["__doc__"], msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print globals()["__doc__"]
            sys.exit(0)

    if len(args) != 1:
        print "please supply filename with replacement identifiers."
        print globals()["__doc__"]
        sys.exit(1)

    identifiers = []
    infile = open(args[0], "r")
    for line in infile:
        if line[0] == ">":
            identifiers.append(line[1:string.find(" ", line)])
    infile.close()

    x = 0
    for line in sys.stdin:
        if line[0] == ">":
            if x >= len(identifiers):
                raise "different number of sequences."
            line = ">" + identifiers[x] + "\n"
            x += 1
        print line[:-1]

if __name__ == "__main__":
    sys.exit(main(sys.argv))

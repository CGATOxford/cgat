##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
gpipe/translate_forward2backward.py - 
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

   python gpipe/translate_forward2backward.py --help

Type::

   python gpipe/translate_forward2backward.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import getopt
import tempfile
import time
import popen2

USAGE = """python %s [OPTIONS] 

Version: $Id: gpipe/translate_forward2backward.py 18 2005-08-09 15:32:24Z andreas $

Wrapper for running gene predictions.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
""" % sys.argv[0]

param_long_options = ["verbose=", "help",
                      "bracket-increment=", "query-border=",
                      "border-refinement=",
                      "exit-identical", "min-score=", "method=",
                      "recursive", "refinement", "probe", "incremental",
                      "exons=", "mask-probe", "format=",
                      "probe-options=", "version"]

param_short_options = "v:hi:b:em:procx:af:"

param_columns = (1, 2, 3, 4)

param_filename_contigs = "contig_sizes"


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
        print USAGE, msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o in ("-b", "--query-border"):
            param_query_border = int(a)

    contig_sizes = {}

    infile = open(param_filename_contigs, "r")
    for line in infile:
        if line[0] == "#":
            continue

        sbjct_token, size, offset = line[:-1].split("\t")
        contig_sizes[sbjct_token] = int(size)

    for line in sys.stdin:
        if line[0] == "#":
            continue

        data = line[:-1].split("\t")
        sbjct_token, sbjct_strand, sbjct_from, sbjct_to = (
            data[param_columns[0]],
            data[param_columns[1]],
            data[param_columns[2]],
            data[param_columns[3]])

        sbjct_from, sbjct_to = int(sbjct_from), int(sbjct_to)

        if sbjct_strand == "-" or sbjct_strand == "-1":
            if contig_sizes.has_key(sbjct_token):
                size = contig_sizes[sbjct_token]
                sbjct_from, sbjct_to = size - sbjct_to, size - sbjct_from

            data[param_columns[2]] = sbjct_from
            data[param_columns[3]] = sbjct_to

        print string.join(map(str, data), "\t")


if __name__ == "__main__":
    sys.exit(main(sys.argv))

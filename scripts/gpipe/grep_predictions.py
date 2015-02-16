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
gpipe/grep_predictions.py - 
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

   python gpipe/grep_predictions.py --help

Type::

   python gpipe/grep_predictions.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import getopt
import CGAT.PredictionParser as PredictionParser

USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: gpipe/grep_predictions.py 18 2005-08-09 15:32:24Z andreas $

Grep predictions from a predictions file.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-f, --file=                     filename with patterns
-k, --keys=                     filename with keys
-o, --format=                   format [predictions|matches]
""" % sys.argv[0]

param_long_options = [
    "verbose=", "help", "file=", "invert-match", "keys=", "format=", "version"]
param_short_options = "v:hf:k:"

param_keep = 1

param_filename = None
param_filename_keys = None
param_format = "predictions"


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
        elif o in ("-f", "--file"):
            param_filename = a
        elif o in ("-k", "--keys"):
            param_filename_keys = a
        elif o in ("-v", "--invert-match"):
            param_keep = 0
        elif o in ("-o", "--format"):
            param_format = a

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    patterns = {}
    if param_filename:
        infile = open(param_filename, "r")
        for line in infile:
            if line[0] == "#":
                continue
            data = line[:-1].split("\t")
            patterns["%s_%s_%s_%s_%s" % tuple(data)] = 1
        infile.close()

    if param_filename_keys:
        infile = open(param_filename_keys, "r")
        for line in infile:
            if line[0] == "#":
                continue
            patterns[line.strip()] = 1
        infile.close()

    if param_format == "predictions":

        p = PredictionParser.PredictionParserEntry()
        for line in sys.stdin:
            if line[0] == "#":
                continue

            p.Read(line)

            key = "%s_%s_%s_%i_%i" % (p.mQueryToken, p.mSbjctToken,
                                      p.mSbjctStrand,
                                      p.mSbjctGenomeFrom, p.mSbjctGenomeTo)

            if patterns.has_key(key):
                keep = param_keep
            else:
                keep = not param_keep

            if keep:
                print line[:-1]

    elif param_format == "matches":

        keep = False
        for line in sys.stdin:
            if re.match("# START:", line):
                key = re.search("# START: key=(\S+)", line).groups()[0]
                if key in patterns:
                    keep = True
                else:
                    keep = False

            if keep:
                print line[:-1]


if __name__ == "__main__":
    sys.exit(main(sys.argv))

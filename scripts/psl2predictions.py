'''
psl2predictions.py - 
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

   python psl2predictions.py --help

Type::

   python psl2predictions.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import getopt
import CGAT.Experiment as E
import CGAT.PredictionParser as PredictionParser

USAGE = """python %s [OPTIONS] < psl > predictions

Convert BLAT output to predictions.

Version: $Id: psl2predictions.py 14 2005-08-09 15:24:07Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-t, --trans                     input is translated DNA
""" % sys.argv[0]


param_long_options = ["verbose=", "help", "trans", "version"]
param_short_options = "v:ht"

param_trans = None


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
        elif o in ("-t", "--trans"):
            param_trans = 1

    print E.GetHeader()
    print E.GetParams()

    if param_trans:
        parser = PredictionParser.PredictionParserBlatTrans()
    else:
        parser = PredictionParser.PredictionParserBlatCDNA()

    nmatches = 1
    for line in sys.stdin:
        if line[0] == "#":
            continue
        if not re.match("^[0-9]", line):
            continue

        try:
            entries = parser.Parse((line,))
        except PredictionParser.AlignmentError, e:
            print "# %s" % str(e)
            print "#", line[:-1]
            sys.exit(1)

        for entry in entries:
            entry.mPredictionId = nmatches
            nmatches += 1

        print str(entries)

    print E.GetFooter()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

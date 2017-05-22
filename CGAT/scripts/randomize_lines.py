'''
randomize_lines.py - randomize lines from stdin
===============================================

:Tags: Python

Purpose
-------

This script reads lines from stdin and outputs them
in randomized order.

Usage
-----

Example::

   cgat randomize-lines < in.lines > out.lines

Command line options
--------------------

'''

import sys
import random
import CGAT.Experiment as E


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-k", "--keep-header", dest="keep_header", type="int",
                      help="randomize, but keep header in place [%default]")

    parser.set_defaults(keep_header=0)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    inf = options.stdin
    outf = options.stdout
    c = E.Counter()
    for x in range(options.keep_header):
        c.header += 1
        outf.write(inf.readline())

    lines = inf.readlines()
    c.lines_input = len(lines)
    random.shuffle(lines)
    for line in lines:
        outf.write(line)
    c.lines_output = len(lines)

    E.info(c)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

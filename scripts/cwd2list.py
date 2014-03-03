'''
cwd2list.py
====================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

create a flat file that represents the current state of the current working directory.
Useful for when data is deleted due to space contraints. If files need to be recreated
then the difference between the current present and previous states can be assessed.

Usage
-----

Example::

   python cwd2list.py 

Type::

   python cwd2list.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import collections
import time
import datetime

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

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    dir2files = {}
    for root, directory, files in os.walk("."):
        dir2files[root] = files

    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H:%M:%S')
    outf = open("CWD_%s" % st, "w")
    outf.write("##state of cwd on %s\n\n" % st)
    for directory, files in dir2files.iteritems():
        for file in files:
            path = os.path.join(directory, file)
            outf.write(path + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

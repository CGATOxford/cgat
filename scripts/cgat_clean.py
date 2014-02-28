'''
cgat_clean.py - remove incomplete files
=======================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script looks at files matching a certain pattern and will remove
incomplete files. File completeness is determined by the file itself
or an associated log-files.

Usage
-----

Example::

   python cgat_clean.py *.tsv

Type::

   python cgat_clean.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys

import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-n", "--dry-run", dest="dry_run", action="store_true",
                      help="dry run, do not delete any files [%default]")

    parser.set_defaults(dry_run=False)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    filenames = args

    c = E.Counter()
    for filename in filenames:
        c.checked += 1
        if IOTools.isComplete(filename):
            c.complete += 1
            continue

        c.incomplete += 1
        E.info('deleting %s' % filename)
        if options.dry_run:
            continue
        os.unlink(filename)
        c.deleted += 1

    E.info(c)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

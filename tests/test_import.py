'''test_import - test importing all modules and pipelines
=========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script attempts to import all the python libraries and
pipeline scripts in the CGAT code collection.

Importing a script/module is a pre-requisite for building
documentation with sphinx. A script/module that can not be imported
will fail within sphinx.

Usage
-----

Example::

   python test_import.py

Type::

   python cgat_import_extensions.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import glob
import traceback
import imp

import CGAT.Experiment as E

# DIRECTORIES to examine for python modulesl/scripts
DIRECTORIES = ('scripts',
               'scripts/optic',
               'scripts/gpipe',
               'CGAT',
               'CGATPipelines')

def main(argv = None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id",
                            usage=globals()["__doc__"])

    parser.set_defaults(
    )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start(parser, argv=argv)

    totals = E.Counter()

    for directory in DIRECTORIES:

        sys.path.insert(0, os.path.abspath(directory))

        files = glob.glob(os.path.join(directory, "*.py"))
        files.sort()

        c = E.Counter()

        for f in files:
            # if f != "pipeline_ancestral_repeats.py" : continue
            E.debug( "importing %s" % f )
            c.input += 1
            prefix, suffix = os.path.splitext(f)

            dirname, basename = os.path.split(prefix)

            if os.path.exists( prefix + ".pyc"):
                os.remove( prefix + ".pyc" )

            success = False
            try:
                __import__(basename, globals(), locals() )
                c.success += 1
                success = True
                options.stdout.write("PASS %s\n" % basename )
                options.stdout.flush()
            except ImportError, msg:
                c.import_fail += 1
                options.stdout.write( "FAIL %s\n%s\n" % (basename, msg) )
                options.stdout.flush()
                traceback.print_exc()
            except Exception, msg:
                c.other_fail += 1
                options.stdout.write( "FAIL %s\n%s\n" % (basename, msg))
                options.stdout.flush()
                traceback.print_exc()

            c.output += 1
            
        E.info('%s: %s' % (directory, c))
        totals += c

    E.info('%s: %s' % ('totals', c))

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

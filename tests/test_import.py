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

def check_import(directory):

    files = glob.glob(os.path.join(directory, "*.py"))
    files.sort()
        
    c = E.Counter()

    #import IndexedGenome
    #files = ('IndexedGenome.py',)

    sys.stderr.write(str(sys.path) + "\n")
    for filename in files:
        sys.stderr.write("importing %s\n" % filename)
        c.input += 1
        prefix, suffix = os.path.splitext(filename)

        dirname, basename = os.path.split(prefix)

        if os.path.exists(prefix + ".pyc"):
            os.remove(prefix + ".pyc")

        pyxfile = os.path.join(dirname, "_") + basename + "x"
        # ignore script with pyximport for now, something does not work
        if os.path.exists(pyxfile):
            c.skipped_pyx += 1
            continue

        success = False
        try:
            imp.load_source(basename, os.path.join(directory,filename))
            c.success += 1
            success = True
            sys.stderr.write("PASS %s\n" % basename)
            sys.stderr.flush()
        except ImportError, msg:
            c.import_fail += 1
            sys.stderr.write("FAIL %s\n%s\n" % (basename, msg))
            sys.stderr.flush()
            traceback.print_exc()
        except Exception, msg:
            c.other_fail += 1
            sys.stderr.write("FAIL %s\n%s\n" % (basename, msg))
            sys.stderr.flush()
            traceback.print_exc()
            
        c.output += 1

    sys.stderr.write('%s: %s\n' % (os.path.dirname(directory), c))


def test_imports():
    '''test importing'''

    # add scripts directories to path.  imp.load_source does not import
    # modules that are in the same directory as the module being
    # loaded from source. This causes a problem in the scripts
    # directory where some scripts import function from other
    # scripts.

    for directory in ('scripts',):
        sys.path.insert(0, os.path.abspath(directory))

    for directory in DIRECTORIES:
        check_import.description = directory
        yield(check_import, os.path.abspath(directory))

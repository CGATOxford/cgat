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

This script is best run within nosetests::

   nosetests tests/test_import.py


'''

import os
import glob
import traceback
import imp

from nose.tools import ok_

import CGAT.Experiment as E

# DIRECTORIES to examine for python modulesl/scripts
DIRECTORIES = ('scripts',
               'scripts/optic',
               'scripts/gpipe',
               'CGAT',
               'CGATPipelines')


# Scripts to exclude as they fail imports.
EXCLUDE = (
    # The following fail because of pybedtools
    # compilation fails. Reason why it triggers
    # recompilation or why it fails is unknown
    # (it seems using C compiler for C++ code).
    'pipeline_intervals',
    'bam2transcriptContribution',
    'beds2counts',
    'fasta2bed',
    # The following fail because of pyximport
    # problems
    'bed2table',)


def check_import(directory):

    files = glob.glob(os.path.join(directory, "*.py"))
    files.sort()

    c = E.Counter()

    outfile = open('test_import.log', 'a')
    outfile.write('testing %s' % directory)

    for filename in files:

        c.input += 1
        # outfile.write('testing %s\n' % filename)

        prefix, suffix = os.path.splitext(filename)

        dirname, basename = os.path.split(prefix)

        if basename in EXCLUDE:
            c.skipped += 1
            continue

        if os.path.exists(prefix + ".pyc"):
            os.remove(prefix + ".pyc")

        pyxfile = os.path.join(dirname, "_") + basename + "x"
        # ignore script with pyximport for now, something does not work
        if os.path.exists(pyxfile):
            c.skipped_pyx += 1
            continue

        try:
            imp.load_source(basename, os.path.join(directory, filename))
            c.success += 1
            outfile.write("PASS %s\n" % basename)
            outfile.flush()
        except ImportError, msg:
            c.fail += 1
            c.import_fail += 1
            outfile.write("FAIL %s\n%s\n" % (basename, msg))
            outfile.flush()
            traceback.print_exc(file=outfile)
        except Exception, msg:
            c.fail += 1
            c.other_fail += 1
            outfile.write("FAIL %s\n%s\n" % (basename, msg))
            outfile.flush()
            traceback.print_exc(file=outfile)

        c.output += 1

    outfile.write('completed %s: %s\n' % (os.path.dirname(directory), c))
    outfile.close()

    ok_(c.fail == 0, '%i scripts/modules could not be imported' % c.fail)


def test_imports():
    '''test importing

    Relative imports will cause a failure because
    imp.load_source does not import modules that are in the same
    directory as the module being loaded from source.
    '''

    for directory in DIRECTORIES:
        check_import.description = directory
        yield(check_import, os.path.abspath(directory))



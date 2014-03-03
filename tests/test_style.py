'''test_import - test importing all modules and pipelines
=========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script runs pep8 on all scripts in the CGAT
code collection.

This script is best run within nosetests::

   nosetests tests/test_style.py

'''
import pep8
import glob
import os
from nose.tools import ok_

# DIRECTORIES to examine for python modules/scripts
EXPRESSIONS = (
    ('tests', 'tests/*.py'),
    ('scripts', 'scripts/*.py'),
    ('optic', 'scripts/optic/*.py'),
    ('gpipe', 'scripts/gpipe/*.py'),
    ('CGAT', 'CGAT/*.py'),
    ('CGATPipelines', 'CGATPipelines/*.py'),
    ('trackers', 'CGATPipelines/pipeline_docs/*/trackers/*.py'))

# Counters to ignore in the pep8 BaseReport
IGNORE = set(('E101',
              'E201',
              'E202',
              'E501',
              'E502',
              'W191',
              'W291',
              'W293',
              'W391',
              'W601',
              'W602',
              'files',
              'directories',
              'physical lines',
              'logical lines',))


def check_style(filename):
    '''check style of filename.
    '''

    p = pep8.StyleGuide(quiet=True)
    report = p.check_files([filename])

    # count errors/warning excluding
    # those to ignore
    take = [y for x, y in report.counters.items() if x not in IGNORE]
    total = sum(take)
    ok_(total == 0, 'pep8 style violations')


def test_style():
    '''test style of scripts
    '''

    for label, expression in EXPRESSIONS:

        files = glob.glob(expression)
        files.sort()

        for f in files:
            if os.path.isdir(f):
                continue
            check_style.description = os.path.abspath(f)
            yield(check_style, os.path.abspath(f))

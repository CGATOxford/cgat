'''test_commandline - test coding style confirmation of CGAT code
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script test the command line usage of all scripts in the
CGAT code collection.

This script is best run within nosetests::

   nosetests tests/test_commandline.py


.. note::

   Make sure to run::

       python setup.py develop

   Before running these tests.

'''
import glob
import os
import imp
import sys
from nose.tools import ok_
import CGAT.Experiment as E

# handle to original E.Start function
ORIGINAL_START = None

# Parser object collected from child script
PARSER = None

# DIRECTORIES to examine for python modules/scripts
EXPRESSIONS = (
    ('scripts', 'scripts/*.py'),)
# ('optic', 'scripts/optic/*.py'),
# ('gpipe', 'scripts/gpipe/*.py'))

EXCLUDE = ("__init__.py",)


class DummyError(Exception):
    pass


def LocalStart(parser, *args, **kwargs):
    '''stub for E.Start - set return_parser argument to true'''
    global PARSER
    PARSER = ORIGINAL_START(parser,
                            return_parser=True,
                            **kwargs
                            )
    raise DummyError()


def loadScript(script_name):

    # call other script
    prefix, suffix = os.path.splitext(script_name)

    dirname = os.path.dirname(script_name)
    basename = os.path.basename(script_name)[:-3]

    if os.path.exists(prefix + ".pyc"):
        os.remove(prefix + ".pyc")

    try:
        module = imp.load_source(basename, script_name)
    except ImportError, msg:
        E.warn('could not import %s - skipped: %s' % (basename, msg))
        return

    return module


def check_options(script_name):
    '''import script and get command line options.

    Test command line options for conformity.
    '''

    module = loadScript(script_name)

    E.Start = LocalStart

    try:
        module.main(argv=["--help"])
    except AttributeError:
        ok_(False, "no main method in %s" % script_name)
    except DummyError:
        pass

    for option in PARSER.option_list:
        # ignore options added by optparse
        if option.dest is None:
            continue

        optstring = option.get_opt_string()
        if optstring.startswith("--"):
            optstring = optstring[2:]

        if "file" in optstring:
            components = optstring.split("-")
            ok_(len(components) == 3,
                "option %s has not three components" % optstring)


def test_cmdline():
    '''test style of scripts
    '''
    # start script in order to build the command line parser
    global ORIGINAL_START
    ORIGINAL_START = E.Start

    x = 0
    for label, expression in EXPRESSIONS:

        files = glob.glob(expression)
        files.sort()

        for f in files:
            if os.path.isdir(f):
                continue
            if os.path.basename(f) in EXCLUDE:
                continue

            check_options.description = os.path.abspath(f)
            yield(check_options, os.path.abspath(f))
            x += 1
            if x > 10:
                break

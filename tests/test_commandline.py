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

   nosetests tests/test_commandline.py --nocapture


.. note::

   Make sure to run::

       python setup.py develop

   Before running these tests.

'''
import glob
import os
import importlib
import yaml
import re
import sys

from nose.tools import ok_
import CGAT.Experiment as E
import CGAT.IOTools as IOTools

# handle to original E.Start function
ORIGINAL_START = None

# Parser object collected from child script
PARSER = None

# DIRECTORIES to examine for python modules/scripts
EXPRESSIONS = (
    ('scripts', 'scripts/*.py'),)
# ('optic', 'scripts/optic/*.py'),
# ('gpipe', 'scripts/gpipe/*.py'))

EXCLUDE = ("__init__.py",
           "version.py",
           "cgat.py",
           "fasta2bed.py",   # fails because of pybedtools rebuild
           )

# Filename with the black/white list of options.
# The file is a tab-separated with the first column
# an option name and the second field a marker.
# Possible markers are:
# ok = whitelist - this option is ok.
# 'bad', 'rename', '?', '' - this option is not ok.
FILENAME_OPTIONLIST = "tests/option_list.tsv"


class DummyError(Exception):
    pass


def filterFiles(files):
    '''filter list of files according to filters set in
    configuration file tests/_test_commandline.yaml'''

    if os.path.exists("tests/_test_commandline.yaml"):
        config = yaml.load(open("tests/_test_commandline.yaml"))
        if config is not None:
            if "restrict" in config and config["restrict"]:
                values = config["restrict"]
                if "manifest" in values:
                    # take scripts defined in the MANIFEST.in file
                    scriptdirs = [x for x in open("MANIFEST.in")
                                  if x.startswith("include scripts") and
                                  x.endswith(".py\n")]
                    take = set([re.sub("include\s*", "",
                                       x[:-1]) for x in scriptdirs])
                    files = [x for x in files if x in take]

                if "regex" in values:
                    rx = re.compile(values["regex"])
                    files = filter(rx.search, files)
    return files


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

    dirname = os.path.relpath(os.path.dirname(script_name))
    basename = os.path.basename(script_name)[:-3]

    if os.path.exists(prefix + ".pyc"):
        os.remove(prefix + ".pyc")

    modulename = ".".join((re.sub("/", ".", dirname), basename))
    try:
        module = importlib.import_module(modulename)
    except ImportError, msg:
        sys.stderr.write('could not import %s - skipped: %s\n' %
                         (modulename, msg))
        module = None

    return module, modulename


def check_option(option, script_name, map_option2action):
    '''import script and get command line options.

    Test command line options for conformity.
    '''
    if option in map_option2action:
        ok_(option in map_option2action,
            'option %s:%s unknown')
        ok_(map_option2action[option] == "ok",
            'option %s:%s wrong: action="%s"' %
            (script_name, option, map_option2action[option]))


def failTest(msg):
    '''create test that fails with *msg*.'''
    ok_(False, msg)


def test_cmdline():
    '''test style of scripts
    '''
    # start script in order to build the command line parser
    global ORIGINAL_START
    ORIGINAL_START = E.Start

    # read the first two columns
    map_option2action = IOTools.readMap(
        IOTools.openFile(FILENAME_OPTIONLIST),
        columns=(0, 1),
        has_header=True)

    files = []
    for label, expression in EXPRESSIONS:
        f = glob.glob(expression)
        files.extend(sorted(f))

    files = filterFiles(files)

    # make sure to use the current working directory as
    # primary lookup.
    sys.path.insert(0, ".")

    # files = [
    #    'scripts/check_db.py',
    #    'scripts/cgat_build_report_page.py']

    for f in files:
        if os.path.isdir(f):
            continue
        if os.path.basename(f) in EXCLUDE:
            continue

        script_name = os.path.abspath(f)
        pyxfile = (os.path.join(os.path.dirname(f), "_") +
                   os.path.basename(f) + "x")

        failTest.description = script_name
        # check if script contains getopt
        with IOTools.openFile(script_name) as inf:
            if "getopt" in inf.read():
                yield (failTest,
                       "script uses getopt directly: %s" % script_name)
                continue

        module, modulename = loadScript(script_name)
        if module is None:
            yield (failTest,
                   "module could not be imported: %s\n" % script_name)
            continue
        E.Start = LocalStart

        try:
            module.main(argv=["--help"])
        except AttributeError:
            yield (failTest,
                   "no main method in %s\n" % script_name)
            ok_(False, "no main method in %s" % script_name)
        except SystemExit:
            yield (failTest,
                   "script does not use E.Start() %s\n" % script_name)
        except DummyError:
            pass

        for option in PARSER.option_list:
            # ignore options added by optparse
            if option.dest is None:
                continue

            optstring = option.get_opt_string()
            if optstring.startswith("--"):
                optstring = optstring[2:]

            check_option.description = script_name + ":" + optstring

            yield(check_option, optstring, os.path.abspath(f),
                  map_option2action)

        # clear up
        del sys.modules[modulename]

        # scripts with pyximport need special handling.
        #
        # Multiple imports of pyximport seems to create
        # some confusion - here, clear up sys.meta_path after
        # each script
        if os.path.exists(pyxfile):
            sys.meta_path = []


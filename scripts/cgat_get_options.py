'''
cgat_get_options.py - build a sorted list of all options used in scripts
========================================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Go through all scripts in the CGAT code collection and collect
options used in the scripts.

This script expects to be executed at the root of the
CGAT code repository.


Usage
-----

.. Example use case

Example::

   python cgat_get_options.py

Type::

   python cgat_get_options.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import glob
import imp
import collections
import pandas
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

# scripts to exclude from collection
EXCLUDE = ("__init__.py",
           "cgat.py",)


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


def collectOptionsFromScript(script_name):
    '''collect options used in script *script_name*.'''

    # call other script
    prefix, suffix = os.path.splitext(script_name)

    dirname = os.path.dirname(script_name)
    basename = os.path.basename(script_name)[:-3]

    if os.path.exists(prefix + ".pyc"):
        os.remove(prefix + ".pyc")

    # check if script contains getopt
    with IOTools.openFile(script_name) as inf:
        if "getopt" in inf.read():
            E.warn("script %s uses getopt directly" % script_name)
            return []

    try:
        module = imp.load_source(basename, script_name)
    except ImportError, msg:
        E.warn('could not import %s - skipped: %s' % (basename, msg))
        return []

    E.Start = LocalStart

    try:
        module.main(argv=["--help"])
    except AttributeError:
        E.warn("no main method in %s" % script_name)
        return []
    except SystemExit:
        E.warn("script exits - possibly does not use E.Start()")
        return []
    except DummyError:
        pass

    result = []
    for option in PARSER.option_list:
        # ignore options added by optparse
        if option.dest is None:
            continue

        optstring = option.get_opt_string()
        if optstring.startswith("--"):
            optstring = optstring[2:]
        result.append(optstring)

    return result


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "--inplace", dest="inplace", action="store_true",
        help="update option list in place. New options will"
        "be added to the list given by --options-tsv-file. "
        "Options will only be added, not removed [%default]")

    parser.add_option(
        "--options-tsv-file", dest="tsv_file", type="string",
        help="existing table with options. Will be updated if "
        "--in-place is set [default]")

    parser.set_defaults(
        inplace=False,
        tsv_file=None)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    old_options = None
    if options.tsv_file:
        if not os.path.exists(options.tsv_file):
            raise OSError(
                "filename %s not found, see --options-tsv-file" %
                options.tsv_file)
        old_options = pandas.read_csv(
            IOTools.openFile(options.tsv_file),
            sep="\t",
            index_col=0,
        )
        old_options = old_options.fillna("")

    global ORIGINAL_START
    ORIGINAL_START = E.Start

    all_options = collections.defaultdict(list)

    for label, expression in EXPRESSIONS:

        files = glob.glob(expression)
        files.sort()

        for f in files:

            E.debug("processing %s" % f)
            if os.path.isdir(f):
                continue
            if os.path.basename(f) in EXCLUDE:
                continue
            collected_options = collectOptionsFromScript(os.path.abspath(f))
            for o in collected_options:
                all_options[o].append(f)

    # add old options
    for x in old_options.index:
        if x not in all_options:
            all_options[x].append("--")

    if options.inplace:
        outfile = IOTools.openFile(options.tsv_file, "w")
        E.info("updating file '%s'" % options.tsv_file)
    else:
        outfile = options.stdout

    outfile.write("option\taction\tcomment\talternative\tfiles\n")
    for o, v in sorted(all_options.items()):
        try:
            action, comment, alternative, ff = old_options.xs(o)

        except KeyError:
            action, comment, alternative, ff = "", "", "", ""

        if comment == "nan":
            comment = ""
        if alternative == "nan":
            alternative = ""

        outfile.write("\t".join((map(
            str, (o, action, comment, alternative, ",".join(v))))) + "\n")

    if outfile != options.stdout:
        outfile.close()

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

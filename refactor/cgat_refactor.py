##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''cgat_refactor.py - refactor CGAT Code
=====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Overview
--------

This script automates some of the refactoring tasks within the CGAT
code base. It currently has two capabilities - renaming scripts and
changing command line options.

Renaming scripts
++++++++++++++++

Renaming command line options
+++++++++++++++++++++++++++++

All command line options in the CGAT code collection are listed in
the file :file:`tests/option_list.tsv`. This list is basically a
white-list of permitted options and a black-list of invalid options.
The black-list contains instructions for changing an invalid option
into a valid option. The script :file:`cgat_refactor.py` can apply
these rules recursively in a collection of files. For example, to update
all options in CGAT code collection, type::

   python refactor/cgat_refactor.py --options=tests/option_list.tsv \
            --suffix=.py --suffix=.rst --suffix=.yaml --suffix=.ipynb \
            --directories=CGAT,CGATPipelines,scripts,makefiles,doc,tests,recipes

This statement will update all files ending in ``.py``, ``.rst`` or ``.yaml``
in the listed directories.

Note that the script only changes long options and implements a simple
pattern match where the option must be preceded by two hyphens and
constitutes a full word. It might get confused by non-standard option
usage.

Usage
-----

The script expects to be executed in the root of the CGAT code repository.

Type::

   python cgat_refactor.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import glob
import pandas

import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def checkUnique(l):
    '''check if elements is list l are unique.'''
    # check for unique mapping
    values = list(sorted(l))
    unique = set(values)
    if len(values) != len(unique):
        raise ValueError(
            "non-unique option mappings")
    else:
        E.info("option list is unique")


def checkOverlap(a, b):
    '''check if a and b not overlap.'''
    aa = set(list(a))
    bb = set(list(b))

    if len(aa.intersection(bb)) != 0:
        raise ValueError(
            "option lists are not separate: shared=%s" %
            ",".join(aa.intersection(bb)))
    else:
        E.info("no overlap between option lists")


def updateFiles(dirs, map_old2new, counter,
                suffixes,
                regex_restrict=None,
                dry_run=False):
    '''iterate through all files in dirs and
    replace patterns is map_old2new'''

    if regex_restrict:
        rx = re.compile(regex_restrict)
    else:
        rx = None

    for d in dirs:
        for root, dirs, files in os.walk(d):
            for f in files:
                _, ext = os.path.splitext(f)
                if rx and not rx.search(f):
                    continue
                if ext not in suffixes:
                    continue

                counter.files_examined += 1
                fn = os.path.join(root, f)
                with IOTools.openFile(fn, "r") as inf:
                    old_data = inf.read()

                changed = False
                for old_name, new_name in map_old2new.items():
                    # only replace at word boundaries
                    old_name += """(['`\s"=])"""
                    new_name += r"\1"
                    new_data = re.sub(old_name, new_name, old_data)
                    if old_data != new_data:
                        changed = True
                        E.info("changed: %s : %s to %s" %
                               (fn, old_name, new_name))
                    old_data = new_data

                if changed:
                    counter.files_changed += 1

                    if not dry_run:
                        with IOTools.openFile(fn, "w") as outf:
                            outf.write(new_data)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--scripts", dest="rename_scripts", type="string",
                      help="rename scripts")

    parser.add_option("--options", dest="rename_options", type="string",
                      help="rename command line options")

    parser.add_option("--split-prefix", dest="split_prefix",
                      type="string",
                      help="move scripts with prefix to subdirectory")

    parser.add_option("--suffix", dest="suffixes", action="append",
                      type="string",
                      help="file suffixes to use.")

    parser.add_option("-n", "--dry-run", dest="dry_run",
                      action="store_true",
                      help="dry run, do not implement any changes")

    parser.add_option("--restrict-regex", dest="regex_restrict", type="string",
                      help="regular expression to restrict refactoring to")

    parser.add_option("-d", "--directories", dest="dirs", action="append",
                      type="string",
                      help="directories to change files in [%defaul]")

    parser.set_defaults(
        rename_scripts=None,
        rename_options=None,
        split_prefix=None,
        scriptsdir="scripts",
        dirs=[],
        suffixes=[],
        dry_run=False,
        regex_restrict=None,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(options.suffixes) == 0:
        raise ValueError("please supply --suffixes")

    if len(options.dirs) == 0:
        raise ValueError("please supply --directories")

    scriptsdir = options.scriptsdir

    counter = E.Counter()

    map_old2new = {}

    if options.rename_scripts or options.split_prefix:

        if options.rename:
            with IOTools.openFile(options.rename_scripts, "r") as inf:
                for line in inf:
                    if line.startswith("#"):
                        continue
                    if line.startswith("old"):
                        continue
                    try:
                        old, new = line[:-1].split("\t")
                    except ValueError:
                        continue
                    if not os.path.exists(os.path.join(scriptsdir, old)):
                        E.warn("%s does not exist - no renaming" % old)
                        continue
                    map_old2new[old] = new

        elif options.split_prefix:
            if not os.path.exists(os.path.join(scriptsdir,
                                               options.split_prefix)):
                E.warn("destination %s does not exist - no renaming" %
                       options.split_prefix)
                return

            scripts = glob.glob("%s/%s_*.py" % (scriptsdir,
                                                options.split_prefix))
            if len(scripts) == 0:
                E.info("nothing to change")
                return

            for script in scripts:
                scriptname = os.path.basename(script)
                newname = scriptname[len(options.split_prefix) + 1:]
                map_old2new[scriptname] = "%s/%s" % (options.split_prefix,
                                                     newname)

        if len(map_old2new) == 0:
            E.info("nothing to change")
            return

        for old, new in map_old2new.items():
            statement = "git mv %(scriptsdir)s/%(old)s %(scriptsdir)s/%(new)s" % locals()
            counter.renamed += 1

            if options.dry_run:
                E.info(statement)
            else:
                E.run(statement)

        updateFiles(options.dirs,
                    map_old2new, counter,
                    suffixes=options.suffixes,
                    dry_run=options.dry_run)

    elif options.rename_options:
        # read refactoring guides
        table = pandas.read_csv(
            IOTools.openFile(options.rename_options),
            sep="\t")

        # select all options that need to renamed
        selected = table[table.action == "rename"]

        # check if all are unique
        checkUnique(selected["option"])
        checkOverlap(selected["option"],
                     selected["alternative"])

        # build map adding "--" prefix
        map_old2new = dict(zip(
            ["--%s" % x for x in selected["option"]],
            ["--%s" % x for x in selected["alternative"]]))

        updateFiles(options.dirs, map_old2new, counter,
                    suffixes=options.suffixes,
                    regex_restrict=options.regex_restrict,
                    dry_run=options.dry_run)

    E.info(str(counter))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

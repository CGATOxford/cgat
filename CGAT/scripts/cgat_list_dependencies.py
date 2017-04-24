'''cgat_list_dependencies.py - collect all external dependencies
=============================================================

:Author:
:Tags: Python

Purpose
-------

This script iterates over all scripts and modules in the CGAT
script collection and checks/list the dependencies.

Usage
-----

The script expects to be executed in the root directory of the CGAT
repository.

Example::

   python scripts/cgat_list_dependencies.py

Type::

   python cgat_list_dependencies.py --help

for command line help.

Command line options
--------------------

'''

import sys
import glob
import os
import collections
import numpy
import CGAT.Experiment as E
import CGAT.Requirements as Requirements


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.set_defaults(
        directories=["scripts", "CGAT", "CGAT/scripts", "CGATPipelines"])

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    # collect all dependencies
    dependencies = {}
    for directory in options.directories:
        files = glob.glob(os.path.join(directory, "*.py"))
        E.debug("processing %i files in %s" % (len(files), directory))

        for f in files:
            dependencies[f] = Requirements.checkRequirementsFromFile(f)

    # group by tool
    tools = collections.defaultdict(list)
    for key, deps in list(dependencies.items()):
        for r in deps:
            tools[r.tool].append((r, key))

    options.stdout.write(
        "tool\trequired\tinstalled\tis_required\tlocations\n")
    for tool, deps in sorted(tools.items()):
        options.stdout.write("\t".join((
            tool,
            ",".join(numpy.unique(
                [x[0].operation + x[0].required_version for x in deps])),
            ",".join(numpy.unique(
                [x[0].installed_version for x in deps])),
            str(not numpy.all([x[0].optional for x in deps])),
            ",".join(numpy.unique(
                [x[1] for x in deps])))) + "\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

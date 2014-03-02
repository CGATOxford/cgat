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
'''
setup_test.py - setup a testing stub for a python script
========================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script sets up a testing stub for a python script
in the CGAT code collection.

The script will:

1. create a testing directory in the current working directory
2. create a stub test.yaml file in the testing directory

Usage
-----

Example::

   python setup_test.py ../scripts/bam2bam.py

Type::

   python setup_test.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E

# common tests to perform
# basic_test: call script, ask for version
YAML_TEMPLATE = '''
version:
    stdin: null
    outputs: [stdout]
    references: []
    options: --version
'''


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("script", "module"),
                      help="type of tests to create [%default].")

    parser.set_defaults(method="script")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) == 0:
        raise ValueError(
            "setup_test.py requires one or more command line arguments")

    targetdir = os.path.dirname(__file__)

    counter = E.Counter()

    for arg in args:
        counter.input += 1
        script_dirname, basename = os.path.split(arg)

        dirname = os.path.join(targetdir, basename)

        if os.path.exists(dirname):
            E.warn("%s already exists - skipping" % basename)
            counter.skipped += 1
            continue

        os.mkdir(dirname)

        with open(os.path.join(dirname, "tests.yaml"), "w") as outf:
            outf.write(YAML_TEMPLATE)

        counter.created += 1

    E.info("%s" % str(counter))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

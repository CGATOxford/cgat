'''
bigwig2hilbert.py - template for CGAT scripts
====================================================

:Author:
:Tags: Python

Purpose
-------

.. Generate hilbert curves for each contig in a bigwig file

Options
-------

.. Options for the script, with detail of how they are combined
.. to provide intended functionality

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import CGAT.Experiment as E
from rpy2.robjects import r as R


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

    parser.add_option("--images-dir", dest="images_dir", type="string",
                      help="directory to save hilbert curves image files to")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    infile = argv[-1]
    pref = infile.split("/")[1].split(".")[0]
    spec = infile.split("/")[1].split(".")[1].split("-")[1]
    header = "%s-%s" % (pref, spec)
    image_dir = options.images_dir

    if os.path.exists(image_dir):
        pass
    else:
        os.mkdir(image_dir)
    
    # set path for R scripts to source
    lib_dir = os.path.dirname(__file__)
    root_dir = os.path.dirname(lib_dir)
    r_dir = os.path.join(root_dir, "R")

    # test R scripts directory - fail if not present
    assert r_dir

    R('''suppressPackageStartupMessages(library(rtracklayer))''')
    R('''data.rle <- rtracklayer::import.bw(con="%(infile)s", '''
      '''as="Rle")''' % locals())
    R('''source("%(r_dir)s/wiggle2hilbert.R")''' % locals())
    R('''wiggle2Hilbert(wiggleRle=data.rle, '''
      '''image.dir="%(image_dir)s", datName="%(header)s")''' % locals())

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

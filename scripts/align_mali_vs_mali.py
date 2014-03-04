'''
align_mali_vs_mali.py - align two multiple alignments
=====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script does a profile vs profile alignment of two
multiple alignments.

This script is only parameterized for the alignment of
protein sequences, but could be extended to the alignment
of DNA sequences as well.

Usage
-----

Example::

   python align_mali_vs_mali.py mali1.fasta mali2.fasta > mali_out.fasta

Type::

   python align_mali_vs_mali.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.Mali as Mali
import CGAT.IOTools as IOTools
import alignlib_lite


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-o", "--gop", dest="gop", type="float",
                      help="gap opening penalty [default=%default].")

    parser.add_option("-e", "--gep", dest="gep", type="float",
                      help="gap extension penalty [default=%default].")

    parser.add_option("-m", "--mode", dest="mode", type="choice",
                      choices=("global", "local"),
                      help="alignment mode, global=nw, local=sw [default=%default].")

    parser.set_defaults(
        gop=-12.0,
        gep=-2.0,
        format="fasta",
        mode="local",
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if len(args) != 2:
        raise ValueError(
            "please supply two multiple alignments in FASTA format.")

    mali1 = Mali.Mali()
    mali2 = Mali.Mali()

    E.info("read 2 multiple alignments")

    mali1.readFromFile(IOTools.openFile(args[0], "r"), format=options.format)
    mali2.readFromFile(IOTools.openFile(args[1], "r"), format=options.format)

    cmali1 = Mali.convertMali2Alignlib(mali1)
    cmali2 = Mali.convertMali2Alignlib(mali2)

    if options.mode == "local":
        mode = alignlib_lite.py_ALIGNMENT_LOCAL
    elif options.mode == "global":
        mode = alignlib_lite.py_ALIGNMENT_GLOBAL

    alignator = alignlib_lite.py_makeAlignatorDPFull(mode,
                                                     options.gop, options.gep)

    alignlib_lite.py_setDefaultEncoder(
        alignlib_lite.py_getEncoder(alignlib_lite.py_Protein20))
    alignlib_lite.py_setDefaultLogOddor(
        alignlib_lite.py_makeLogOddorDirichlet(0.3))
    alignlib_lite.py_setDefaultRegularizor(
        alignlib_lite.py_makeRegularizorDirichletPrecomputed())

    cprofile1 = alignlib_lite.py_makeProfile(cmali1)
    cprofile2 = alignlib_lite.py_makeProfile(cmali2)

    result = alignlib_lite.py_makeAlignmentVector()

    alignator.align(result, cprofile1, cprofile2)

    E.debug("result=\n%s" % alignlib_lite.py_AlignmentFormatEmissions(result))

    cmali1.add(cmali2, result)

    outmali = Mali.convertAlignlib2Mali(cmali1,
                                        identifiers=mali1.getIdentifiers() + mali2.getIdentifiers())

    outmali.writeToFile(options.stdout, format=options.format)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

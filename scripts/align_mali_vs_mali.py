################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
#################################################################################
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
import CGAT.Mali as Mali
import CGAT.IOTools as IOTools
import alignlib

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-o", "--gop", dest="gop", type="float",
                      help="gap opening penalty [default=%default]."  )

    parser.add_option("-e", "--gep", dest="gep", type="float",
                      help="gap extension penalty [default=%default]."  )

    parser.add_option("-m", "--mode", dest="mode", type="choice",
                      choices = ("global", "local" ),
                      help="alignment mode, global=nw, local=sw [default=%default]."  )

    parser.set_defaults(
        gop = -12.0,
        gep = -2.0,
        format= "fasta",
        mode = "local",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 2: raise ValueError("please supply two multiple alignments in FASTA format.")

    mali1 = Mali.Mali()
    mali2 = Mali.Mali()

    E.info( "read 2 multiple alignments" )

    mali1.readFromFile( IOTools.openFile( args[0], "r" ), format=options.format )
    mali2.readFromFile( IOTools.openFile( args[1], "r" ), format=options.format )

    cmali1 = Mali.convertMali2Alignlib( mali1 )
    cmali2 = Mali.convertMali2Alignlib( mali2 )

    if options.mode == "local":
        mode = alignlib.ALIGNMENT_LOCAL
    elif options.mode == "global":
        mode = alignlib.ALIGNMENT_GLOBAL
        
    alignator = alignlib.makeAlignatorDPFull( mode,
                                              options.gop, options.gep )

    alignlib.setDefaultEncoder( alignlib.getEncoder( alignlib.Protein20) )
    alignlib.setDefaultLogOddor( alignlib.makeLogOddorDirichlet( 0.3 ) )
    alignlib.setDefaultRegularizor( alignlib.makeRegularizorDirichletPrecomputed() )

    cprofile1 = alignlib.makeProfile( cmali1 )
    cprofile2 = alignlib.makeProfile( cmali2 )

    result = alignlib.makeAlignmentVector()

    alignator.align( result, cprofile1, cprofile2 )

    E.debug( "result=\n%s" % alignlib.AlignmentFormatEmissions( result) )

    cmali1.add( cmali2, result )

    outmali = Mali.convertAlignlib2Mali( cmali1,
                                         identifiers = mali1.getIdentifiers() + mali2.getIdentifiers() )
    
    outmali.writeToFile( options.stdout, format=options.format)

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

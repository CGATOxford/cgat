################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
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
WrapperHmmer.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile, subprocess, optparse, shutil

from types import *

"""Wrapper for Hmmer
"""

# raise NotImplementedError("incomplete")

import Experiment
import TreeTools, MatrixTools

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

    def __str__(self):
        return str(self.message)

class ParsingError(Error):
    """Exception raised for errors while parsing

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message, line):
        self.message = message + " at line " + line

class UsageError(Error):
    """Exception raised for errors while starting

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class HmmerSearchResult:
    def __init__(self):
        self.mNexus = None
        self.mOutfile = None

class Hmmer:

    def __init__( self ):
        
        self.mMapInput2Phylip = {}
        self.mMapPhylip2Input = {}
        
        self.mTempdir = None
        
        self.mInputTree = None
        self.mInputMatrix = None
        
        self.mOutputTree = "-"
        
        self.mOptions = ""

        self.mOutputStdout = None
        self.mOutputStderr = None

        self.mPruneTree = False

    def __del__(self):

        if self.mTempdir and self.mTempdir != "tmp":
            os.system( "rm -rf %s" % self.mTempdir )

    def setProgram( self, program ):
        self.mProgram = program

    def setOptions( self, options ):
        self.mOptions = options
        
    def setTree( self, tree ):
        self.mInputTree = tree

    def setMatrix( self, matrix):
        self.mInputMatrix = matrix

    def setPruneTree( self, flag = True ):
        self.mPruneTree = flag

    def updateMaps( self, taxa ):
        """update maps."""
        for taxon in taxa:
            if taxon not in self.mMapInput2Phylip:
                key = "t_%05i" % (len(self.mMapInput2Phylip) )
                self.mMapInput2Phylip[taxon] = key
                self.mMapPhylip2Input[key] = taxon


    def __reset(self):
        """clear internal state variables."""
        self.mMapInput2Phylip = {}
        self.mMapPhylip2Input = {}
        
        self.mTempdir = None
        

    def parseSearch( self, lines ):
        pass


if __name__ == "__main__":

    
    phylip = Phylip()

    parser = E.OptionParser( version = "%prog version: $Id: WrapperHmmer.py 2807 2009-10-22 13:42:10Z andreas $" )

    parser.add_option("-t", "--filename-input-tree", dest="filename_input_tree", type="string",
                      help="filename with tree information."  )
    parser.add_option("-T", "--filename-output-tree", dest="filename_output_tree", type="string",
                      help="output filename with tree information."  )
    parser.add_option("-p", "--program", dest="program", type="string",
                      help="program to use."  )
    parser.add_option("-o", "--options", dest="options", type="string",
                      help="input options."  )

    parser.set_defaults(
        filename_input_tree = None,
        filename_output_tree = None,
        program = None,
        options = "",
        )

    (options, args) = Experiment.Start( parser )

    if options.filename_input_tree != "-":
        phylip.setTree( open(options.filename_input_tree, "r").readlines() )
    elif options.filename_input_tree == "-":
        phylip.setTree( sys.stdin.readlines() )

    phylip.setOptions( options.options )
        
    phylip.setProgram( options.program )

    phylip.run()

    Experiment.Stop()
    

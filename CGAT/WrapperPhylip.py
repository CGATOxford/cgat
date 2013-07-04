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
WrapperPhylip.py - 
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

"""Wrapper for CodeML
"""

import Experiment
import TreeTools, MatlabTools

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

class PhylipResult:
    def __init__(self):
        self.mNexus = None
        self.mOutfile = None
        self.mMatrix = None

    def parseContrasts( self, infile ):
        """parse results from contrasts run."""
        ## remove empty lines and strip
        
        lines = filter( lambda x: x, map( lambda x: x.strip(), infile.readlines()))

        sections = filter( lambda x: lines[x+1][:3] == "---", range(0,len(lines)-1))
        sections.append(len(lines))
        
        def readTable( lines ):
            data = []
            for line in lines:
                data.append( map( float, re.split( "\s+", line[:-1] ) ) )
            return data

        for s in range(len(sections)-1):
            start = sections[s]
            end = sections[s+1]

            header = lines[start]
            if header[:len("Contrasts")] == "Contrasts":
                self.mContrasts = readTable( lines[start+2:end] )
            elif header[:len("Covariance")] == "Covariance":
                self.mCovariance = readTable( lines[start+2:end] )
            elif header[:len("Correlations")] == "Correlations":
                self.mCorrelations = readTable( lines[start+2:end] )
            elif header[:len("Regressions")] == "Regressions":
                self.mRegressions = readTable( lines[start+2:end] )
                
class Phylip:

    def __init__( self ):
        
        self.mMapInput2Phylip = {}
        self.mMapPhylip2Input = {}
        
        self.mTempdir = None
        
        self.mInputTree = None
        self.mInputTrees = None
        self.mInputData = None
        self.mInputMatrix = None
        self.mInputMali = None
        
        self.mOutputTree = "-"
        
        self.mOptions = ""

        self.mOutputStdout = None
        self.mOutputStderr = None

        self.mPruneTree = False

        self.mLogLevel = 0

    def __del__(self):

        if self.mTempdir and self.mTempdir != "tmp" and self.mLogLevel == 0:
            os.system( "rm -rf %s" % self.mTempdir )

    def setLogLevel( self, loglevel ):
        self.mLogLevel = loglevel 

    def setProgram( self, program ):
        self.mProgram = program

    def setOptions( self, options ):
        self.mOptions = options
        
    def setTree( self, tree ):
        self.mInputTree = tree

    def setTrees( self, trees ):
        self.mInputTrees = trees

    def setData( self, data ):
        """set input data. This is a matrix with row headers."""
        self.mInputData = data

    def setMatrix( self, matrix):
        self.mInputMatrix = matrix

    def setMali( self, mali):
        self.mInputMali = mali

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
        
    def prepareRun(self):

        self.__reset()

        self.mTempdir = tempfile.mkdtemp()
        # self.mTempdir = "tmp"
        if not os.path.exists( self.mTempdir ):
            os.mkdir(self.mTempdir)

        if self.mInputMatrix and self.mInputData:
            raise ValueError("please specify either input matrix or input data, but not both.")
        
        ## prepare input matrix. Should already be in phylip like
        ## format, but long identifiers are shortened and tabs are
        ## replaced by spaces.
        if self.mInputMatrix:

            outfile = open(self.mTempdir + "/infile", "w")
            
            identifiers = map( lambda x: re.split("\s+", x[:-1])[0], self.mInputMatrix[1:])
            self.updateMaps( identifiers )

            outfile.write(self.mInputMatrix[0])
            for line in self.mInputMatrix[1:]:
                data = re.split("\s+", line[:-1])
                new_line = self.mMapInput2Phylip[data[0]] + "       " + "  ".join(data[1:])
                outfile.write(new_line + "\n" )
                
            outfile.close()

            if self.mLogLevel >= 1:
                print "# written input matrix with %i taxa to %s" % (len(identifiers), self.mTempdir + "/infile" )
                os.system("cat %s" %  self.mTempdir + "/infile" )

        elif self.mInputData:
            
            outfile = open(self.mTempdir + "/infile", "w")
            outfile.write( "%i %i\n" % (len(self.mInputData), len(self.mInputData[0]) -1 ) )
            identifiers = map( lambda x: x[0], self.mInputData) 
            self.updateMaps( identifiers)
            
            for x in range(len(identifiers)):
                outfile.write( "%-10s %s\n" % (self.mMapInput2Phylip[identifiers[x]], " ".join(self.mInputData[x][1:]) ) )

            outfile.close()
            
            if self.mLogLevel >= 1:
                print "# written input matrix with %i taxa to %s" % (len(identifiers), self.mTempdir + "/infile" )
                os.system("cat %s" %  self.mTempdir + "/infile" )
                
        ## prepare input tree or trees
        self.mNInputTrees = 0
        if self.mInputTree or self.mInputTrees:

            outfile = open(self.mTempdir + "/intree", "w")

            if self.mInputTree and self.mInputTrees:
                raise UsageError("please supply either one or mupltiple trees, but not both." )

            if self.mInputTree:
                trees = [ self.mInputTree ]
            else:
                trees = self.mInputTrees

            for tree in trees:
                if self.mPruneTree:
                    taxa = self.mMapInput2Phylip.keys()
                    TreeTools.PruneTree( tree, taxa )

                taxa = TreeTools.GetTaxa( tree )
                self.updateMaps( taxa )
                TreeTools.MapTaxa( tree, self.mMapInput2Phylip )

                ## check if taxa are unique
                taxa = tree.get_taxa()
                staxa = set()
                
                skip = False
                for t in taxa:
                    if t in staxa:
                        if self.mLogLevel >= 1:
                            print "# skipping tree %s because of duplicate taxa." % (tree.name)
                        skip = True
                    staxa.add(t)
                    
                if skip:
                    continue
                
                outfile.write( TreeTools.Tree2Newick( tree ) + "\n" )
                self.mNInputTrees += 1
                
                if self.mLogLevel >= 1:
                    print "# written input tree with %i taxa to %s" % (len(TreeTools.GetTaxa( tree )), self.mTempdir + "/intree" )
                    print "#", TreeTools.Tree2Newick( tree ) 

            outfile.close()

        ## prepare input multiple alignment
        if self.mInputMali:

            if self.mInputMatrix:
                raise "both mali and matrix supplied - infile conflict."
            
            outfile = open(self.mTempdir + "/infile", "w" )

            identifiers = self.mInputMali.getIdentifiers()
            self.updateMaps( identifiers )
            self.mInputMali.mapIdentifiers( self.mMapInput2Phylip )
            self.mInputMali.writeToFile( outfile, format="phylip" )

            outfile.close()

            if self.mLogLevel >= 1:
                print "# written input multiple alignments with %i taxa and with %i to %s" %\
                      (self.mInputMali.getLength(), self.mInputMali.getWidth(), self.mTempdir + "/intree" )

    def run( self ):

        self.prepareRun()

        if not self.mProgram:
            raise UsageError( "no program specified." )

        s = subprocess.Popen( "%s" % (self.mProgram),
                              shell = True,
                              stdin = subprocess.PIPE,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = self.mTempdir,
                              close_fds = True)                              

        (out, err) = s.communicate( "\n".join(self.mOptions) + "\n")

        if s.returncode != 0:
            raise UsageError, "Error in running phylip.\n%s\n%s\nTemporary directory was %s" % (out, err, self.mTempdir)

        ## Parse output files that might have been created:
        result = PhylipResult()

        ## parse tree file
        if os.path.exists( "%s/outtree" % self.mTempdir):
            
            nexus = TreeTools.Newick2Nexus( open("%s/outtree" % self.mTempdir, "r") )
            for tree in nexus.trees:
                TreeTools.MapTaxa( tree, self.mMapPhylip2Input )
            result.mNexus = nexus
            if self.mLogLevel >= 1:
                print "# received tree with %i taxa" % (len(TreeTools.GetTaxa(nexus.trees[0])))
                
        elif os.path.exists( "%s/outfile" % self.mTempdir ):
            
            if self.mProgram in ("dnadist", "protdist" ):
                infile = open( "%s/outfile" % self.mTempdir, "r")
                result.mMatrix, row_headers, col_headers = MatlabTools.readMatrix( infile, format="phylip" )
                result.mRowHeaders = []
                for x in row_headers: result.mRowHeaders.append( self.mMapPhylip2Input[ x ] )
                result.mColHeaders = result.mRowHeaders
            elif self.mProgram == "contrast":

                infile = open( "%s/outfile" % self.mTempdir, "r" )
                result.parseContrasts( infile )
                infile.close()

        else:
            raise "other return types not implemented"

        if self.mLogLevel >= 2:
            print out

        if self.mLogLevel == 0:
            shutil.rmtree( self.mTempdir )

        return result

if __name__ == "__main__":
    
    phylip = Phylip()

    parser = E.OptionParser( version = "%prog version: $Id: WrapperPhylip.py 2784 2009-09-10 11:41:14Z andreas $" )

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
    

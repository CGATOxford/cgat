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
WrapperNJTree.py - 
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

"""Wrapper for NJTree
"""

import Experiment as E
import Mali
import IOTools
import TreeTools

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

class Result:
    
    mMaxLineLength = 100
    
    def __init__(self):
        self.mLog = None
        self.mErr = None
        
    def __str__(self):
        return TreeTools.Tree2Newick(self.mTree)
        
class IdentifierParser:

    def __init__(self, map_species2sp = None):
        self.mMapSpecies2Sp = map_species2sp

class IdentifierParserGPipe(IdentifierParser):

    mSeparator = "|"
    
    def __init__(self, *args, **kwargs):
        IdentifierParser.__init__(self, *args, **kwargs)

    def __call__(self, id ):
        data = id.split( self.mSeparator )

        species = data[0]
        
        if len(data) == 2:
            gene = data[1]
            transcript = None
        elif len(data) >= 3:
            gene = data[2]
            transcript = data[1]
            
        if self.mMapSpecies2Sp:
            species = self.mMapSpecies2Sp[species]

        return species, gene, transcript

class NJTree:
    """
    """
    mExecutable = "treebest"

    def __init__( self, identifier_parser ):

        self.mFilenameSpeciesTree = None
        self.mFilenameMali = None
        self.mFilenameTree = None
        
        self.mIdentifierParser = identifier_parser

        self.mDefaults = { "positive_only": "0",
                           "initial_kappa": "2.0",
                           "initial_omega": "0.1" }

        self.mOptions = {
            "positive_only" : {  "0": "output positve and neutral sites",
                                 "1": "only output positive sites"
                                 },
            "initial_kappa"   : { "#" : "initial value for kappa",
                                  },
            "initial_omega"  : { "#" : "initial value for omega",
                                 },
            }

        self.mLog = sys.stdout
        self.mErr = sys.stderr

    ##------------------------------------------------------------------------
    def SetLog( self, log ):
        """set log stream."""
        self.mLog = log

    ##------------------------------------------------------------------------
    def SetErr( self, err ):
        """set error stream."""
        self.mErr = err
        
    ##------------------------------------------------------------------------    
    def GetOptions( self ):
        """return options in pretty format"""
        result = ["# Options for NJTree" ]
        for var in self.mOptions.keys():
            result.append("# %-40s: %s" % (var, self.mOptions[var]))
        return string.join(result, "\n")


    ##------------------------------------------------------------------------
    def WriteAlignment( self, mali ):
        """write alignment in Phylip format."""

        mali.mapIdentifiers( self.mMapOld2New )
        outfile = open( self.mTempdir + "/" + self.mFilenameMali, "w" )
        mali.writeToFile( outfile, format="plain-fasta" )
        outfile.close()

    ##------------------------------------------------------------------------
    def BuildMap( self, mali ):
        """build a map of identifiers that conform to the NJ standard."""

        self.mMapOld2New = {}
        self.mMapNew2Old = {}
        for old_id in mali.getIdentifiers():

            species, gene, transcript = self.mIdentifierParser( old_id )

            if transcript:
                new_id = "%s_%s GENEid=%s" % ( transcripts, species, gene )
            else:
                new_id = "%s_%s" % ( gene, species )                

            self.mMapOld2New[old_id] = new_id
            self.mMapNew2Old[new_id] = old_id

    ##------------------------------------------------------------------------
    def SetSpeciesTree( self, filename_tree ):
        """filename with species tree."""
        self.mFilenameSpeciesTree = filename_tree
    ##------------------------------------------------------------------------            
    def WriteTree( self, tree ):
        """write tree to file.
        """

        nexus = TreeTools.Newick2Nexus( tree )
        t = nexus.trees[0]
        TreeTools.MapTaxa( t, self.mMapOld2New )
        
        outfile = open( self.mTempdir + "/" + self.mFilenameTree, "w" )
        outfile.write("%i 1\n" % self.mNumSequences )
        outfile.write("%s\n" % TreeTools.Tree2Newick(t))
        outfile.close()

    ##------------------------------------------------------------------------        
    def Run( self, alignment, tree = None,
             dump = 0,
             test = False,
             options = {} ):

        self.mTempdir = tempfile.mkdtemp()
        self.mFilenameMali = "input"
        self.mFilenameOutput = "output"
        
        self.BuildMap( alignment )
        self.WriteAlignment( alignment )

        if test:
            print "# temporary directory is %s" % self.mTempdir
            
        if tree:
            self.mFilenameTree = "tree"
            
            ## check what kind of tree is given.
            if type(tree) == StringType:
                t = tree.strip()
                if not (t[0] == "(" and t[-1] in ");"):
                    tree = "".join(open( tree, "r").readlines())
                self.WriteTree( tree )

        extra_options = []

        if self.mFilenameSpeciesTree:
            extra_options.append( "-f %s" % os.path.abspath(self.mFilenameSpeciesTree ))
            
        statement = " ".join( ( self.mExecutable,
                                "best",
                                "-o %s" % self.mFilenameOutput,
                                " ".join(extra_options),
                                self.mFilenameMali) )

        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = self.mTempdir,
                              close_fds = True)                              

        (out, err) = s.communicate()
        
        if s.returncode != 0:
            raise UsageError, "Error in running %s \n%s\n%s\nTemporary directory in %s" % (statement, err, out, self.mTempdir)

        if dump:
            self.mLog.write( "# stdout output of %s:\n%s\n######################################" % (self.mExecutable, out) )
            self.mLog.write( "# stderr output of %s:\n%s\n######################################" % (self.mExecutable, err) )

        lines = open( "%s/%s" % (self.mTempdir, self.mFilenameOutput), "r").readlines()

        if len(lines) == 0:
            raise UsageError, "Empty result from %s \n%s\n%s\nTemporary directory in %s" % (statement, err, out, self.mTempdir)            

        if dump:
            self.mLog.write( "# result output of %s:\n%s\n######################################\n" % (self.mExecutable, "".join(lines)) )

        if not test:
            shutil.rmtree( self.mTempdir )

        return self.parseOutput( lines, out, err )

    def parseOutput( self, lines, out, err ):

        lines = re.sub("\s", "", "".join(lines))
        lines = re.sub("\[[^\]]+\]", "", lines)
        
        t = TreeTools.Newick2Nexus( "".join(lines) )

        result = Result()
        t = t.trees[0]
        
        TreeTools.MapTaxa( t, self.mMapNew2Old )
        
        result.mTree = t

        result.mLog = out
        result.mErr = err
        
        return result

if __name__ == "__main__":
    
    parser = E.OptionParser( version = "%prog version: $Id: WrapperNJTree.py 2781 2009-09-10 11:33:14Z andreas $" )

    parser.add_option( "-m", "--map", dest="filename_map", type="string",
                       help="filename with mapping of species ids to swissprot species ids." )
    
    parser.add_option( "-t", "--tree", dest="filename_tree", type="string",
                       help="filename of trees." )

    parser.add_option("-a", "--filename-alignment", dest="filename_alignment", type="string",
                      help="filename with aligned codon sequences."  )

    parser.add_option( "--dump", dest="dump", action="store_true",
                      help="dump output."  )

    parser.set_defaults(
        separator="|",        
        dump = False,
        filename_map = None,
        filename_alignment = "-",
        filename_tree = None,
        )

    (options, args) = E.Start( parser )

    if options.filename_map:
        map_species2sp = IOTools.ReadMap( open(options.filename_map, "r"))

    E.debug( "species map: %s" % str(map_species2sp) )

    identifier_parser = IdentifierParserGPipe( map_species2sp = map_species2sp )

    njtree = NJTree( identifier_parser = identifier_parser )

    njtree.SetLog( options.stdlog )
    njtree.SetErr( options.stderr )    

    if options.filename_tree:
        njtree.SetSpeciesTree( options.filename_tree )
        
    mali = Mali.Mali()
    if options.filename_alignment == "-":
        infile = sys.stdin
    else:
        infile = open(options.filename_alignment, "r")
        
    mali.readFromFile( infile, format = "fasta")

    if mali.getLength() == 1:
        if options.loglevel >= 1:
            options.stdlog.write("# Warning: single gene tree\n" )
        options.stdout.write( "(%s:1);\n" % tuple(mali.getIdentifiers()) )
    elif mali.getLength() == 2:
        if options.loglevel >= 1:
            options.stdlog.write("# Warning: two gene tree\n" )
        options.stdout.write( "(%s:1,%s:1);\n" % tuple(mali.getIdentifiers()) )
    else:
        result = njtree.Run( mali,
                             dump = options.dump )
        
        options.stdout.write( str(result) + "\n" )

        if result.mLog:
            options.stdlog.write( str(result.mLog) + "\n" )

        if result.mErr:
            options.stderr.write( str(result.mErr) + "\n" )

    E.Stop()

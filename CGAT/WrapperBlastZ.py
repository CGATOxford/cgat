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
WrapperBlastZ.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile, optparse, subprocess

"""Wrapper for BlastZ.
"""

import Genomics
import Experiment as E

class BlastZError(Exception):
    pass

class BlastZ:

    # options for sensitive alignment of short (<10Mb) sequences
    # according to blastz manual.
    mOptions = "C=2 B=0 T=0 W=6 K=2200"
    mExecutable = "blastz"
    mStderr = sys.stderr
    
    def __init__( self, options=None):

        if options:
            self.mOptions = options

    def CreateTemporaryFiles( self ):
        """create temporary files."""
        self.mTempDirectory = tempfile.mkdtemp()
        
        self.mFilenameTempSeq1 = self.mTempDirectory + "/seq1"
        self.mFilenameTempSeq2 = self.mTempDirectory + "/seq2"

    def DeleteTemporaryFiles( self ):
        """clean up."""
        os.remove( self.mFilenameTempSeq1 )
        os.remove( self.mFilenameTempSeq2 )
        os.rmdir( self.mTempDirectory )

    def SetStderr( self, file = None):
        """set file for dumping stderr."""
        self.mStderr= file
        
    def ParseResult( self, out, result ):
        """parse BLASTZ output.

        TODO: check for subopt alignments - does blastz do these?
        """
        result.clear()

        keep = False

        self.mReverseComplement = False

        def block_parser( out ):
            keep = False
            lines = []
            for line in out.split("\n"):
                if not line: continue

                if line[0] == "}":
                    yield id, lines
                    
                if re.match( "[shaxmd] {", line):
                    id = line[0]
                    keep = True
                    lines = []
                    continue
                
                if keep:
                    lines.append( line.strip() )
                    
        for id, lines in block_parser(out):

            if id == "a":
                for line in lines:
                    d = line.split(" ")
                    if d[0] == "l":
                        row_from, col_from, row_to, col_to = map(int, d[1:5])
                        result.addDiagonal( row_from-1, row_to, col_from - row_from )
            
            elif id == "s":
                d = lines[1].split(" ")
                self.mReverseComplement = d[3] == "1"

    def isReverseComplement( self ):
        return self.mReverseComplement

    def Align( self, seq1, seq2, result ):
        """align two sequences with BlastZ."""
        
        self.CreateTemporaryFiles()

        f = open(self.mFilenameTempSeq1, "w")
        f.write(">seq1\n%s" % seq1 )
        f.close()

        f = open(self.mFilenameTempSeq2, "w")
        f.write(">seq2\n%s" % seq2 )
        f.close()
        
        statement = string.join( ( self.mExecutable,
                                   self.mFilenameTempSeq1,
                                   self.mFilenameTempSeq2,                                   
                                   self.mOptions),
                                 " ")

        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = self.mTempDirectory,
                              close_fds = True)                              

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise BlastZError, "Error in running BLASTZ \n%s" % err

        self.ParseResult( out, result )

        self.mStderr.write( err )

        self.DeleteTemporaryFiles()
     
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: WrapperBlastZ.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-i", "--input-file-seq1", dest="input_filename_seq1", type="string",
                      help="input filename for sequence 1.",
                      metavar="FILE" )

    parser.add_option("-j", "--input-file-seq2", dest="input_filename_seq2", type="string",
                      help="input filename for sequence 2.",
                      metavar="FILE" )

    parser.add_option("-o", "--options", dest="options", type="string",
                      help="BlastZ options." )

    parser.set_defaults( \
        input_filename_seq1 = None,
        input_filename_seq2 = None,
        options = "B=0 C=2")
    
    (options, args) = E.Start( parser ) 
    
    wrapper = BlastZ( options.options )

    import alignlib_lite
    seqs1 = Genomics.ReadPeptideSequences( open(options.input_filename_seq1, "r") )
    seqs2 = Genomics.ReadPeptideSequences( open(options.input_filename_seq2, "r") )
    seq1 = seqs1[seqs1.keys()[0]]
    seq2 = seqs2[seqs2.keys()[0]]    
    result = alignlib_lite.py_makeAlignmentVector()
    wrapper.Align( seq1, seq2, result) 

    print str( alignlib_lite.py_AlignmentFormatExplicit( result,
                                                 alignlib_lite.py_makeSequence( seq1 ),
                                                 alignlib_lite.py_makeSequence( seq2 ) ) )
    
    E.Stop()
        
            
        
                                 
        
        
        
        
    
        

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
WrapperBl2Seq.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile, subprocess, optparse

"""Wrapper for adaptive codon bias program
"""

import Experiment
import FastaIterator

class Bl2SeqError(Exception):
    pass
            
class Bl2Seq:

    mOptions = ""
    mExecutable = "bl2seq"
    mStderr = sys.stderr
    
    def __init__( self, options=""):
        
        self.mOptions = options

    def CreateTemporaryFiles( self ):
        """create temporary files."""
        self.mTempDirectory = tempfile.mkdtemp()
        
        self.mFilenameTempInput = self.mTempDirectory + "/input"
        self.mFilenameTempOutput = self.mTempDirectory + "/output"

    def DeleteTemporaryFiles( self ):
        """clean up."""
        os.remove( self.mFilenameTempInput )
        os.remove( self.mFilenameTempOutput )
        os.rmdir( self.mTempDirectory )

    def SetStderr( self, file = None):
        """set file for dumping stderr."""
        self.mStderr= file
        
    def WriteOutput( self, lines, filename_output = None):
        """write output to file.

        If file is not given, lines are written to stdout.
        """
        
        if filename_output:
            outfile = open(filename_output, "w")
        else:
            outfile = sys.stdout
            
        outfile.write( string.join( lines, "") )

        if filename_output:
            outfile.close()

    def ParseResult( self, trace_file = None, information_file = None ):

        result = AdaptiveCAIResult()

        result.Read( trace_file, information_file )
        
        return result
    
    def RunOnFile( self, infile, outfile, errfile ):

        self.CreateTemporaryFiles()

        statement = string.join( ( self.mExecutable,
                                   self.mFilenameTempInput,
                                   self.mFilenameTempOutput ),
                                 " ")

        
        i = FastaIterator.FastaIterator( infile )

        outfile.write( "GENE\tBl2Seq\n")
        
        while 1:
            f = i.next()
            if f == None: break

            file = open(self.mFilenameTempInput, "w")
            file.write( ">%s\n%s" % (f.title, f.sequence) )
            file.close()

            
            s = subprocess.Popen( statement,
                                  shell = True,
                                  stdout = subprocess.PIPE,
                                  stderr = subprocess.PIPE,
                                  cwd = self.mTempDirectory,
                                  close_fds = True)                              

            (out, err) = s.communicate()

            if s.returncode != 0:
                raise Bl2SeqError, "Error in calculating Bl2Seq\n%s" % err

            d = open(self.mFilenameTempOutput).readlines()[2][:-1]
            enc = d.split( " " )[2]
            
            outfile.write( (string.join((f.title, enc), "\t") )  + "\n")
            
            errfile.write( err )

        self.DeleteTemporaryFiles()
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: WrapperBl2Seq.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-f", "--input-file", dest="input_filename", type="string",
                      help="input filename. If '-', stdin is used [default=%default].",
                      metavar="FILE" )

    parser.add_option("-o", "--output-file", dest="output_filename", type="string",
                      help="output filename for codon usage. If '-', output is stdout [default=%default].",
                      metavar="FILE" )

    parser.add_option("-e", "--error-file", dest="error_filename", type="string",
                      help="output filename for error messages. If '-', output is stderr [default=%default].",
                      metavar="FILE" )


    parser.set_defaults( \
        input_filename = "-",
        output_filename = "-",
        error_filename = "/dev/null",
        )
    
    (options, args) = Experiment.Start( parser ) 

    wrapper = Bl2Seq( )

    if options.input_filename == "-":
        file_stdin = sys.stdin
    else:
        file_stdin = open( options.input_filename, "r" )

    if options.output_filename:
        if options.output_filename == "-":
            file_stdout = sys.stdout
        else:
            file_stdout = open( options.output_filename, "w")

    if options.error_filename:
        if options.error_filename == "-":
            file_stderr = sys.stderr
        else:
            file_stderr = open( options.error_filename, "w")
        
    wrapper.RunOnFile( file_stdin, file_stdout, file_stderr )

    if file_stdin and file_stdin != sys.stdin:
        file_stdin.close()

    if file_stdout and file_stdout != sys.stdout:
        file_stdout.close()

    if file_stderr and file_stderr != sys.stderr:
        file_stderr.close()

    Experiment.Stop()
        

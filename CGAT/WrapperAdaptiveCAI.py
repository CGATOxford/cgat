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
WrapperAdaptiveCAI.py - 
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

import Experiment as E

class CAIError(Exception):
    pass

class AdaptiveCAIResult:

    def __init__(self):

        self.mIterations = 0
        self.mCodonUsageSphereStart = []        
        self.mCodonUsageSpheres = []
        self.mSetSizes = []
        self.mNumSequences = 0
        self.mWeights = []
        self.mFinalWeights = {}
        self.mGeneInfo = {}
        self.mCodonUsages = []

    def GetDominantThreshold( self, fraction = 0.1 ):
        """get CAI threshold for dominant set."""

        cai_values = map( lambda x: x["CAICLASS"], self.mGeneInfo.values())

        cai_values.sort()
        cai_values.reverse()

        return cai_values[int(len(cai_values) * fraction)]

    def Read( self, trace_file = None, gene_file = None, codon_file = None):
        
        """parse data from CAI output files."""
        
        if trace_file:
            iteration = None
            final_weights = None
            
            for line in trace_file:
                if line[:3] == "Got": continue

                if re.match( "CodonUsageIt", line ):
                    data = re.split("\s", line.strip())[2:]
                    self.mCodonUsageSpheres.append( map(float, data) )
                elif re.match( "CodonUsage sphere", line ):
                    data = re.split("\s", line.strip())[2:]
                    self.mCodonUsageSphereStart = map(float, data)
                elif re.match("First \d+ genes", line):
                    size = int(line.split(" ")[1])
                    self.mSetSizes.append( size )
                elif re.match("Number of points = \d+", line):
                    self.mNumSequences = int(line.split(" ")[4])
                elif re.match("Iteration \d", line):
                    iteration = int(line.split(" ")[1])
                    self.mWeights.append( {} )
                elif re.match("AdditionalInfoFile", line):
                    continue
                elif re.match("Looking", line):
                    continue
                elif re.match("Loaded", line):
                    continue
                elif re.match("Loading", line):
                   continue
                elif re.match("Gene", line):
                    continue
                elif re.match("External W-table", line):
                    ## ignore this table.
                    break
                elif re.match("GTFFinal", line):
                    final_weights = True
                    continue
                else:
                    data = re.split("\s", line.strip())
                    if iteration == None: continue
                    if final_weights:
                        for d in data:
                            codon, f = d.split(":")
                            self.mFinalWeights[codon.upper()] = float(f)
                    else:
                        for d in data:
                            codon, f = d.split(":")
                            self.mWeights[iteration][codon.upper()] = float(f)

        if gene_file:
            lines = gene_file.readlines()
            headers = []
            num_columns, num_rows = map(int, lines[0][:-1].split(" "))
            
            for x in range( 1, num_columns+1):
                header, type = lines[x][:-1].split( " ")
                header = header.upper()
                if header == "GENENAME":
                    genename_index = len(headers)
                headers.append((header, type))

            if num_columns+1+num_rows != len(lines):
                raise ValueError, "parsing error: not enough lines in %s" % gene_file

            for x in range( num_columns+1,num_columns+num_rows+1):
                data = re.split( "\s+", lines[x].strip())
                v = {}
                
                if len(data) != len(headers):
                    raise ValueError, "not enough fields in line %s" % lines[x]
                
                for c in range( len(headers) ):
                    if c == genename_index:
                        genename = data[c].replace( '"', '')
                        continue
                    
                    header,type = headers[c]

                    if type == "FLOAT":
                        v[header] = float(data[c])
                    elif type == "STRING":
                        v[header] = data[c]

                self.mGeneInfo[genename] = v

        if codon_file:

            for line in codon_file:

                data = re.split("\s+", line.strip())

                aa, codon = data[:2]

                if not self.mCodonUsages:
                    for x in range(len(data)-2):
                        self.mCodonUsages.append({})

                iteration = 0
                for v in data[2:]:

                    self.mCodonUsages[iteration][codon.upper()] = float(v)
                    iteration += 1

                
                
            
class AdaptiveCAI:

    mOptions = ""
    mIterations = 10
    mExecutable = "caijava"
    mTrace  = sys.stdout
    mStderr = sys.stderr
    mExtra  = sys.stdout
    mCodons = sys.stdout
    
    def __init__( self, iterations=10, options=""):
        
        self.mOptions = options
        self.mIterations = iterations

    def CreateTemporaryFiles( self ):
        """create temporary files."""
        self.mTempDirectory = tempfile.mkdtemp()
        
        self.mFilenameTempExtra = self.mTempDirectory + "/extra"
        self.mFilenameTempInput = self.mTempDirectory + "/input"

        ## note: this file is created automatically by caijava
        self.mFilenameTempCodons = self.mTempDirectory + "/codonusage"

        self.mFilenameTempWeights = self.mTempDirectory + "/weights"        

    def DeleteTemporaryFiles( self ):
        """clean up."""
        os.remove( self.mFilenameTempExtra )
        os.remove( self.mFilenameTempInput )
        os.remove( self.mFilenameTempCodons )
        if os.path.exists( self.mFilenameTempWeights ):
            os.remove( self.mFilenameTempWeights )                
        os.rmdir( self.mTempDirectory )

    def SetCodons( self, file = None):
        """set file for the information file."""
        self.mCodons = file
        
    def SetExtra( self, file = None):
        """set file for the information file."""
        self.mExtra = file
        
    def SetTrace( self, file = None):
        """set file for dumping stdout.

        If file is None, stdout is parsed but not dumped.
        """
        self.mTrace = file

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
    
    def RunOnFile( self, infile, infile_weights = None ):

        self.CreateTemporaryFiles()

        f = open(self.mFilenameTempInput, "w")
        lines = filter(lambda x: x[0] != "#", infile.readlines() )
        for line in lines:
            f.write(line)
        
        f.close()

        statement = string.join( ( self.mExecutable,
                                   self.mFilenameTempInput,
                                   "-i %i" % self.mIterations,
                                   "-s -g",
                                   "-f %s" % self.mFilenameTempExtra),
                                 " ")

        if infile_weights:
            f = open(self.mFilenameTempWeights, "w")
            for line in infile_weights: f.write(line)
            f.close()
            statement += " -ew %s" % self.mFilenameTempWeights

        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = self.mTempDirectory,
                              close_fds = True)                              

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise CAIError, "Error in calculating CAI\n%s" % err
            
        self.mTrace.write( out )
        self.mStderr.write( err )
        self.mExtra.write( string.join(open(self.mFilenameTempExtra, "r").readlines(), ""))
        self.mCodons.write( string.join(open(self.mFilenameTempCodons, "r").readlines(), ""))

        self.DeleteTemporaryFiles()

    
    
     
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: WrapperAdaptiveCAI.py 756 2006-09-20 16:38:02Z andreas $")

    parser.add_option("-f", "--input-file", dest="input_filename", type="string",
                      help="input filename. If not given, stdin is used.",
                      metavar="FILE" )

    parser.add_option("-w", "--input-file-weights", dest="input_filename_weights", type="string",
                      help="input filename with weights. If not given, no CAIEXT values will be there.",
                      metavar="FILE" )

    parser.add_option ("-i", "--iterations", dest="iterations", type="int",
                       help="number of iterations." )

    parser.add_option("-o", "--output-file-trace", dest="output_filename_trace", type="string",
                      help="output filename for cai. If not given, no file is produced.",
                      metavar="FILE" )

    parser.add_option("-e", "--output-file-extra", dest="output_filename_extra", type="string",
                      help="output filename for extra information from cai. If not given, no file is produced.",
                      metavar="FILE" )

    parser.add_option("-c", "--output-file-codons", dest="output_filename_codons", type="string",
                      help="output filename for codon usage information.",
                      metavar="FILE" )


    parser.set_defaults( \
        input_filename = "-",
        input_filename_weights = None,
        output_filename_trace = None,
        output_filename_extra = None,
        output_filename_codons = None,                
        iterations = 10),

    (options, args) = E.Start( parser ) 

    wrapper = AdaptiveCAI( iterations=options.iterations )

    if options.input_filename == "-":
        file_stdin = sys.stdin
    else:
        file_stdin = open( options.input_filename, "r" )

    file_trace, file_extra, file_codons, file_weights = None, None, None, None
    
    if options.output_filename_trace:
        if options.output_filename_trace == "-":
            file_trace = sys.stdout
        else:
            file_trace = open( options.output_filename_trace, "w")
        wrapper.SetTrace( file_trace )
        
    if options.output_filename_extra:
        if options.output_filename_extra == "-":
            file_extra = sys.stdout
        else:
            file_extra = open( options.output_filename_extra, "w")
        wrapper.SetExtra( file_extra )

    if options.output_filename_codons:
        if options.output_filename_codons == "-":
            file_codons = sys.stdout
        else:
            file_codons = open( options.output_filename_codons, "w")
        wrapper.SetCodons( file_codons )

    if options.input_filename_weights:
        if options.input_filename_weights == "-":
            file_weights = sys.stdin
        else:
            file_weights = open( options.input_filename_weights, "r")

    wrapper.RunOnFile( file_stdin, infile_weights = file_weights)

    if file_codons and file_codons != sys.stdout:
        file_codons.close()
            
    if file_extra and file_extra != sys.stdout:
        file_extra.close()
        
    if file_trace and file_trace != sys.stdout:
        file_trace.close()

    if file_stdin and file_stdin != sys.stdin:
        file_stdin.close()

    if file_weights and file_weights != sys.stdin:
        file_weights.close()

    E.Stop()
        

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
WrapperExonerate.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile, subprocess, optparse

"""Wrapper for Exonerate
"""

import Experiment as E
import PredictionParser
import Genomics

class ExonerateError(Exception):
    pass
            
class Exonerate:

    mDefaultOptions = "-m p2g --showvulgar FALSE --showsugar FALSE --showcigar FALSE"
    mOptions = ""
    mOutputOptions = '--showalignment FALSE ' + \
                     ' --ryo "diy\t%S\t%ql\t%r\t%pi\t%ps\t%V\n" --showtargetgff FALSE --showquerygff FALSE'
    
    mExecutable = "exonerate"
    mStdout = sys.stdout
    mStderr = sys.stderr
    
    def __init__( self, options="", output_options=[]):
        
        self.mOptions = options
        self.mDoParse = True

        if output_options:
            o = []
            if "gff" in output_options:
                o.append( "--showtargetgff TRUE --showquerygff TRUE")
            else:
                o.append( "--showtargetgff FALSE --showquerygff FALSE")
                    
            if "alignment" in output_options:
                o.append( "--showalignment TRUE")
            else:
                o.append( "--showalignment FALSE")

            if "parsed" in output_options:
                o.append( '--ryo "diy\t%S\t%ql\t%r\t%pi\t%ps\t%V\n"')
                self.mDoParse = True
            else:
                self.mDoParse = False
                
            self.mOutputOptions = " ".join(o)
        
    def CreateTemporaryFiles( self ):
        """create temporary files."""
        self.mTempDirectory = tempfile.mkdtemp()
        self.mFilenameTempPeptide = self.mTempDirectory + "/peptide.fasta"
        self.mFilenameTempGenome = self.mTempDirectory + "/genome.fasta"        
        
    def DeleteTemporaryFiles( self ):
        """clean up."""
        os.remove( self.mFilenameTempPeptide )
        os.remove( self.mFilenameTempGenome )
        os.rmdir( self.mTempDirectory )

    def SetStderr( self, file = None):
        """set file for dumping stderr."""
        self.mStderr= file

    def SetStdout( self, file = None):
        """set file for dumping stderr."""
        self.mStdout= file
        
    def Run( self,
             peptide_sequence,
             genomic_sequence ):

        self.CreateTemporaryFiles()

        outfile = open(self.mFilenameTempPeptide, "w")
        outfile.write(">%s\n%s\n" % ("peptide", peptide_sequence))
        outfile.close()

        outfile = open(self.mFilenameTempGenome, "w")
        outfile.write(">%s\n%s\n" % ("genome", genomic_sequence))
        outfile.close()
        
        statement = string.join( map(str, (
            self.mExecutable,
            self.mDefaultOptions,
            self.mOptions,
            self.mOutputOptions,
            "--query",
            self.mFilenameTempPeptide,
            "--target",
            self.mFilenameTempGenome,
            )), " " )

        if self.mLogLevel >= 3:
            print "# statement: %s" % statement

        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = self.mTempDirectory,
                              close_fds = True)                              

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise ExonerateError, "Error in calculating Exonerate\n%s" % err
            
        self.mStdout.write( out )
        self.mStderr.write( err )

        self.DeleteTemporaryFiles()

        parser = PredictionParser.PredictionParserExonerate()
        if self.mDoParse:
            return parser.Parse( out, peptide_sequence, genomic_sequence )
        

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: WrapperExonerate.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-p", "--peptide", dest="input_filename_peptide", type="string",
                      help="input filename with peptide sequence.",
                      metavar="FILE" )

    parser.add_option("-g", "--genome", dest="input_filename_genome", type="string",
                      help="input filename with genome sequence.",
                      metavar="FILE" )

    parser.add_option("-R", "--range-genome", dest="range_genome", type="string",
                      help="range on genome sequence.")

    parser.add_option("-P", "--range-peptide", dest="range_peptide", type="string",
                      help="range on peptide sequence.")

    parser.add_option("-d", "--output", dest="output_options", action="append",
                      help="output options [gff|alignment]" )

    parser.add_option("-o", "--options", dest="options", type="string",
                      help="options to exonerate" )

    parser.set_defaults( \
        input_filename_peptide = None,
        input_filename_genome = None,
        output_filename_stdout = "-",
        output_filename_stderr = "-",        
        options = "",
        range_genome = None,
        range_peptide = None,
        id_peptides = None,
        id_genomes = None,
        output_options = []
        )

    (options, args) = E.Start( parser ) 

    if options.range_genome: options.range_genome = map(int, options.range_genome.split(","))
    if options.range_peptide: options.range_peptide = map(int, options.range_peptide.split(","))    

    wrapper = Exonerate( options=options.options, output_options=options.output_options )
    wrapper.mLogLevel = options.loglevel

    if options.loglevel >= 2:
        print "# reading peptide sequence."
    peptide_sequences = Genomics.ReadPeptideSequences( open(options.input_filename_peptide, "r") )
    
    if options.loglevel >= 2:
        print "# reading genome sequence."
    genome_sequences = Genomics.ReadGenomicSequences( open(options.input_filename_genome, "r"), do_reverse = 0 )    

    if not options.id_peptides:
        options.id_peptides= peptide_sequences.keys()
    if not options.id_genomes:
        options.id_genomes= genome_sequences.keys()
        
    for x in options.id_peptides:
        ps = peptide_sequences[x]
        if options.range_peptide:
            ps = ps[options.range_peptide[0]:options.range_peptide[1]]
        for y in options.id_genomes:
            gs = genome_sequences[y]
            if options.range_genome:
                gs = gs[options.range_genome[0]:options.range_genome[1]]

            if options.loglevel >= 2:
                print "# aligning peptide of length %i: %s to genome of length %i: %s" % (len(ps), x, len(gs), y)
            result = wrapper.Run( ps, gs )

            print result
            
    E.Stop()
        


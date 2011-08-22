################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
WrapperMEDIPS.py - wrap MEDIPS methylation analysis
===================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

perform Medip-Seq data analysis using the MEDIPS R package.

Note that the UCSC genome file has to have been downloaded 
previously in Bioconductor::

   source("http://bioconductor.org/biocLite.R")
   biocLite("BSgenome.Hsapiens.UCSC.hg19")

Usage
-----

Documentation
-------------

Code
----

'''

import os, sys, re, optparse, tempfile, shutil, subprocess

import Experiment as E

## for zinba
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError

def bamToMEDIPS( infile, outfile ):
    '''convert bam to medips format

    contig, start, end, strand

    Start is 1-based.
    .'''

    statement = '''bamToBed -i %(infile)s | awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\n", $1,$2+1,$3,$6)}' > %(outfile)s''' % locals()

    E.debug( "executing statement '%s'" % statement )

    retcode = subprocess.call(  statement,
                                cwd = os.getcwd(), 
                                shell = True )
    if retcode < 0:
        raise OSError( "Child was terminated by signal %i: \n%s\n" % (-retcode, statement ))

    return outfile

def compress( infile ):
    '''gzip infile'''

    statement = "gzip %(infile)s" % locals() 

    E.debug( "executing statement '%s'" % statement )

    retcode = subprocess.call(  statement,
                                cwd = os.getcwd(), 
                                shell = True )
    if retcode < 0:
        raise OSError( "Child was terminated by signal %i: \n%s\n" % (-retcode, statement ))

    return outfile

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-f", "--input-format", dest="input_format", type="choice",
                      choices = ("bed", "bam"),
                      help="input file format [default=%default]."  )
    
    parser.add_option("-g", "--genome", dest="genome", type="string",
                      help="UCSC genome identifier [default=%default]."  )

    parser.add_option("-e", "--extension", dest="extension", type="int",
                      help="extension size [default=%default]."  )

    parser.add_option("-b", "--bin-size", dest="bin_size", type="int",
                      help="bin size of genome vector [default=%default]."  )

    parser.add_option("-l", "--fragment-length", dest="fragment_length", type="int",
                      help="bin size of genome vector [default=%default]."  )

    parser.add_option("-s", "--saturation-iterations", dest="saturation_iterations", type="int",
                      help = "iterations for saturation analysis [default=%default]."  )
    
    parser.set_defaults(
        input_format = "bam",
        genome = "hg19",
        extension = 400,
        bin_size = 50,
        saturation_iterations = 10,
        fragment_length = 700,
        )


    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    if len(args) != 1:
        raise ValueError("please specify a filename with sample data")

    filename_sample = args[0]
    
    # load MEDIPS
    R.library( 'MEDIPS' )
    genome_file = 'BSgenome.Hsapiens.UCSC.%s' % options.genome 
    R.library( genome_file )

    #tmpdir = tempfile.mkdtemp( )
    tmpdir = "/tmp/tmps6hP4h"

    E.debug( "temporary files are in %s" % tmpdir )

    if options.input_format == "bam" and 0:
        E.info( "converting bam files to bed" )
        filename_sample = bamToMEDIPS( filename_sample, os.path.join( tmpdir, "sample.bed" ) )

    E.info( "loading data" )
    R('''CONTROL.SET = MEDIPS.readAlignedSequences(
                       BSgenome = "%(genome_file)s", 
                       file = "%(filename_sample)s" ) ''' % locals() )

    E.info( "computing genome vector" )
    R('''CONTROL.SET = MEDIPS.genomeVector(data = CONTROL.SET, bin_size = 50, extend=400 )''')
    
    if options.saturation_analysis:
        E.info( "saturation analysis" )
        R('''sr.control = MEDIPS.saturationAnalysis(data = CONTROL.SET, bin_size = 50, extend = 400, no_iterations = 10, no_random_iterations = 1)''')

        R.png( Experiment.getOutputFile( "calibration.png" ) )
        R('''MEDIPS.plotSaturation(sr.control)''')

    E.info( "computing CpG positions" )
    R('''CONTROL.SET = MEDIPS.getPositions(data = CONTROL.SET, pattern = "CG")''' )

    if options.coverage_analysis:
        E.info( "CpG coverage analysis" )
        R('''cr.control = MEDIPS.coverageAnalysis(data = CONTROL.SET, extend = 400, no_iterations = 10)''')
        R.png( Experiment.getOutputFile( "cpg_coverage.png" ) )
        MEDIPS.plotCoverage(cr.control)
        R('''er.control = MEDIPS.CpGenrich(data = CONTROL.SET)''')

    E.info( "compute coupling vector" )
    R('''CONTROL.SET = MEDIPS.couplingVector(data = CONTROL.SET, fragmentLength = 700, func = "count")''')
    
    E.info( "plotting calibration" )
    R.png( Experiment.getOutputFile( "calibration.png" ) )
    R('''MEDIPS.plotCalibrationPlot(data = CONTROL.SET, linearFit = T)''')
    
    outputfile = Experiment.getOutputFile( "rpm.wig" )
    R('''MEDIPS.exportWIG(file = %(outputfile)s, data = CONTROL.SET, raw = T, descr = "rpm")''' % locals())
    compress( outputfile )
    
    outputfile = Experiment.getOutputFile( "rms.wig" )
    R('''MEDIPS.exportWIG(file = %(outputfile)s, data = CONTROL.SET, raw = F, descr = "rms")''' % locals())
    compress( outputfile )

    shutil.rmtree( tmpdir )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

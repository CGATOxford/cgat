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

import os, sys, re, optparse, tempfile, shutil, subprocess, tempfile

import Experiment as E
import IOTools
import IndexedFasta

## for zinba
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError


def bamToMEDIPS( infile, outfile ):
    '''convert bam to medips format

    contig, start, end, strand

    Start is 1-based.
    '''

    statement = '''bamToBed -i %(infile)s | awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\n", $1,$2+1,$3,$6)}' > %(outfile)s''' % locals()

    E.debug( "executing statement '%s'" % statement )

    E.run( statement )

    return outfile

def bedToMEDIPS( infile, outfile ):
    '''convert bam to medips format

    contig, start, end, strand

    Start is 1-based.
    '''

    if infile.endswith( ".gz" ): cat = "zcat"
    else: cat = "cat"

    statement = '''%(cat)s %(infile)s | awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\n", $1,$2+1,$3,$6)}' > %(outfile)s''' % locals()

    E.run( statement )

    return outfile

def compress( infile ):
    '''gzip infile'''

    statement = "gzip -f %(infile)s" % locals() 

    E.debug( "executing statement '%s'" % statement )

    return E.run( statement )

def bigwig( infile, contig_sizes ):
    '''convert infile to bigwig file'''

    if infile.endswith( ".wig"):
        outfile = infile[:-4] + ".bigwig"
    else:
        outfile = infile + ".bigwig"
        
    tmp, filename_sizes = tempfile.mkstemp() 

    os.write( tmp, "\n".join( [ "\t".join(map(str,x)) for x in contig_sizes.iteritems() ] ) )
    os.close( tmp )

    statement = "wigToBigWig -clip %(infile)s %(filename_sizes)s %(outfile)s " % locals() 

    E.debug( "executing statement '%s'" % statement )

    if E.run( statement ):
        os.unlink( infile )

    os.unlink( filename_sizes )

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-f", "--input-format", dest="input_format", type="choice",
                      choices = ("bed", "bam"),
                      help="input file format [default=%default]."  )
    
    parser.add_option("-u", "--ucsc-genome", dest="ucsc_genome", type="string",
                      help="UCSC genome identifier [default=%default]."  )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option("-e", "--extension", dest="extension", type="int",
                      help="extension size [default=%default]."  )

    parser.add_option("-b", "--bin-size", dest="bin_size", type="int",
                      help="bin size of genome vector [default=%default]."  )

    parser.add_option("-l", "--fragment-length", dest="fragment_length", type="int",
                      help="bin size of genome vector [default=%default]."  )

    parser.add_option("-s", "--saturation-iterations", dest="saturation_iterations", type="int",
                      help = "iterations for saturation analysis [default=%default]."  )
    
    parser.add_option( "-t", "--toolset", dest="toolset", type="choice", action="append",
                       choices = ("saturation", "coverage", "rms", "rpm", "all"),
                       help = "actions to perform [default=%default]." )
    
    parser.add_option( "-w", "--bigwig", dest="bigwig", action = "store_true",
                       help = "store wig files as bigwig files - requires a genome file [default=%default]" )

    parser.set_defaults(
        input_format = "bam",
        ucsc_genome = "hg19",
        genome_file = None,
        extension = 400,
        bin_size = 50,
        saturation_iterations = 10,
        fragment_length = 700,
        toolset = [],
        bigwig = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    if len(args) != 1:
        raise ValueError("please specify a filename with sample data")

    if options.bigwig and not options.genome_file:
        raise ValueError("please provide a genome file when outputting bigwig")

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
        contig_sizes = fasta.getContigSizes()
        
    filename_sample = args[0]

    if len(options.toolset) == 0: options.toolset = ["all"]

    do_all = "all" in options.toolset
    
    # load MEDIPS
    R.library( 'MEDIPS' )
    genome_file = 'BSgenome.Hsapiens.UCSC.%s' % options.ucsc_genome 
    R.library( genome_file )
    
    tmpdir = tempfile.mkdtemp( )

    E.debug( "temporary files are in %s" % tmpdir )

    bin_size = options.bin_size
    extension = options.extension
    fragment_length = options.fragment_length
    saturation_iterations = options.saturation_iterations

    if options.input_format == "bam":
        E.info( "converting bam files" )
        filename_sample = bamToMEDIPS( filename_sample, os.path.join( tmpdir, "sample.medips" ) )
    elif options.input_format == "bed":
        E.info( "converting bed files" )
        filename_sample = bedToMEDIPS( filename_sample, os.path.join( tmpdir, "sample.medips" ) )

    E.info( "loading data" )
    R('''CONTROL.SET = MEDIPS.readAlignedSequences(
                       BSgenome = "%(genome_file)s", 
                       file = "%(filename_sample)s" ) ''' % locals() )
    slotnames = ( ( "extend", "extend", "%i"),
                  ( "distFunction", "distance_function", "%s"),
                  ( "slope", "slope", "%f"),
                  ( "fragmentLength", "fragment_length", "%i" ),
                  ( "bin_size", "bin_size", "%i"),
                  ( "seq_pattern", "pattern", "%s" ),
                  ( "number_regions", "nregions", "%i"),
                  ( "number_pattern", "npatterns", "%i" ),
                  ( "cali_chr", "calibration_contig", "%s"),
                  ( "genome_name", "genome", "%s") )


    E.info( "computing genome vector" )
    R('''CONTROL.SET = MEDIPS.genomeVector(data = CONTROL.SET, 
                       bin_size = %(bin_size)i, 
                       extend=%(extension)i )''' % locals())

    E.info( "computing CpG positions" )
    R('''CONTROL.SET = MEDIPS.getPositions(data = CONTROL.SET, pattern = "CG")''' )

    E.info( "compute coupling vector" )
    R('''CONTROL.SET = MEDIPS.couplingVector(data = CONTROL.SET, 
                       fragmentLength = %(fragment_length)i, 
                       func = "count")''' % locals() )
    
    E.info( "compute calibration curve" )
    R('''CONTROL.SET = MEDIPS.calibrationCurve(data = CONTROL.SET)''')

    E.info( "normalizing" )
    R('''CONTROL.SET = MEDIPS.normalize(data = CONTROL.SET)''')

    outfile = IOTools.openFile( E.getOutputFile( "summary.tsv.gz" ), "w" )
    outfile.write( "category\tvalue\n" )

    if "saturation" in options.toolset or do_all:
        E.info( "saturation analysis" )
        R('''sr.control = MEDIPS.saturationAnalysis(data = CONTROL.SET, 
                            bin_size = %(bin_size)i, 
                            extend = %(extension)i, 
                            no_iterations = %(saturation_iterations)i, 
                            no_random_iterations = 1)''' % locals() )

        R.png( E.getOutputFile( "saturation.png" ) )
        R('''MEDIPS.plotSaturation(sr.control)''')
        R('''dev.off()''')

        R('''write.csv( sr.control$estimation, file ='%s' )'''% E.getOutputFile( "saturation_estimation.csv" ) )
        outfile.write( "estimated_correlation\t%f\n" % R('''sr.control$maxEstCor''')[1] )
        outfile.write( "true_correlation\t%f\n" % R('''sr.control$maxTruCor''')[1] )

    if "coverage" in options.toolset or do_all:
        E.info( "CpG coverage analysis" )
        R('''cr.control = MEDIPS.coverageAnalysis(data = CONTROL.SET, 
                                extend = %(extension)i, 
                                no_iterations = 10)''' % locals())

        R.png( E.getOutputFile( "cpg_coverage.png" ) )
        R('''MEDIPS.plotCoverage(cr.control)''')
        R('''dev.off()''')

        # three rows
        R('''write.csv( cr.control$coveredPos, file ='%s' )'''% E.getOutputFile( "saturation_coveredpos.csv" ) )
        # coverage threshold
        # number of CpG covered
        # percentage of CpG covered

        R('''write.csv( cr.control$matrix, file ='%s' )'''% E.getOutputFile( "saturation_matrix.csv" ) )

        # R('''er.control = MEDIPS.CpGenrich(data = CONTROL.SET)''')

    if "calibration" in options.toolset or do_all:
        E.info( "plotting calibration" )
        R.png( E.getOutputFile( "calibration.png" ) )
        R('''MEDIPS.plotCalibrationPlot(data = CONTROL.SET, linearFit = T, xrange=250)''')
        R('''dev.off()''')

    
    for slotname, label, pattern in slotnames:
        value = tuple(R('''CONTROL.SET@%s''' % slotname ))
        if len(value) == 0: continue
        outfile.write( "%s\t%s\n" % (label, pattern % tuple(R('''CONTROL.SET@%s''' % slotname ))[0] ) )
        
    outfile.close()
        
    if "rpm" in options.toolset or do_all:
        outputfile = E.getOutputFile( "rpm.wig" )
        R('''MEDIPS.exportWIG(file = '%(outputfile)s', data = CONTROL.SET, raw = T, descr = "rpm")''' % locals())
        if options.bigwig:
            bigwig( outputfile, contig_sizes )
        else:
            compress( outputfile )
    
    if "rms" in options.toolset or do_all:
        outputfile = E.getOutputFile( "rms.wig" )
        R('''MEDIPS.exportWIG(file = '%(outputfile)s', data = CONTROL.SET, raw = F, descr = "rms")''' % locals())
        if options.bigwig:
            bigwig( outputfile, contig_sizes )
        else:
            compress( outputfile )

    shutil.rmtree( tmpdir )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

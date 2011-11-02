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
WrapperZinba.py - wrap zinba peak caller
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Runs the zinba peak caller.

Note that for zinba to work, a mappability file needs to have
been created.

The output bed files contain the P-value as their score field.

Usage
-----

Documentation
-------------

Code
----

'''

import os, sys, re, optparse, tempfile, shutil, subprocess
import collections

import Experiment as E

## for zinba
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError

def bamToBed( infile, outfile ):
    '''convert bam to bed with bedtools.'''

    statement = "bamToBed -i %(infile)s > %(outfile)s" % locals()

    E.debug( "executing statement '%s'" % statement )

    retcode = subprocess.call(  statement,
                                cwd = os.getcwd(), 
                                shell = True )
    if retcode < 0:
        raise OSError( "Child was terminated by signal %i: \n%s\n" % (-retcode, statement ))

    return outfile

ZinbaPeak = collections.namedtuple( "ZinbaPeak", "contig unrefined_start unrefined_end strand posterior summit height refined_start refined_end median fdr" )

def iteratePeaks( infile ):
    '''iterate of zinba peaks in infile.'''
    
    for line in infile:

        if line.startswith("#"): continue
        if line.startswith("PEAKID\tChrom"): continue
        # skip empty lines
        if line.startswith("\n"): continue

        data = line[:-1].split("\t")

        if len(data) != 12:
            raise ValueError( "could not parse line %s" % line )

        # I assume these are 1-based coordinates
        data[2] = max(int(data[2]) - 1, 0)
        # end
        data[3] = int(data[3])
        # posterior
        data[5] = float(data[5])
        # summit
        data[6] = max(int(data[6]) - 1, 0)
        # height
        data[7] = int(data[7])
        # refined_start
        data[8] = max(int(data[8]) - 1, 0)
        # end
        data[9] = int(data[9])
        # median
        data[10] = int(data[10])
        # qvalue
        data[11] = float(data[11])

        yield ZinbaPeak._make( data[1:] )

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
    
    parser.add_option("-s", "--fragment-size", dest="fragment_size", type="int",
                      help="fragment size [default=%default]."  )

    parser.add_option("-m", "--mappability-dir", dest="mappability_dir", type="string",
                      help="mappability_dir [default=%default]."  )

    parser.add_option("-b", "--bit-filename", dest="bit_filename", type="string",
                      help="2bit genome filename [default=%default]."  )

    parser.add_option("-c", "--control-filename", dest="control_filename", type="string",
                      help="filename of input/control data in bed format [default=%default]."  )

    parser.add_option("-i", "--index-dir", dest="index_dir", type="string",
                      help="index directory [default=%default]."  )

    parser.add_option("-t", "--threads", dest="threads", type="int",
                      help="number of threads to use [default=%default]."  )

    parser.add_option("-q", "--fdr-threshold", dest="fdr_threshold", type="float",
                      help="fdr threshold [default=%default]."  )

    parser.add_option("-a", "--alignability-threshold", dest="alignability_threshold", type="int",
                      help="alignability threshold [default=%default]."  )

    parser.set_defaults(
        input_format = "bed",
        fragment_size = 200,
        mappability_dir = None,
        threads = 1,
        alignability_threshold = 1,
        bit_filename = None,
        fdr_threshold = 0.05,
        )


    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 2:
        raise ValueError("please specify a filename with sample data and an output file")

    filename_sample, filename_output = args[0], args[1]
    filename_control = options.control_filename
    
    # load Zinba
    R.library( 'zinba' )

    tmpdir = tempfile.mkdtemp( )

    E.debug( "temporary files are in %s" % tmpdir )

    if options.input_format == "bam":
        E.info( "converting bam files to bed" )
        filename_sample = bamToBed( filename_sample, os.path.join( tmpdir, "sample.bed" ) )
        if filename_control:
            filename_control = bamToBed( filename_control, os.path.join( tmpdir, "control.bed" ) )

    fragment_size = options.fragment_size
    threads = options.threads
    bit_filename = options.bit_filename
    mappability_dir = options.mappability_dir
    fdr_threshold = options.fdr_threshold

    E.info( "computing counts" )
    R( '''basealigncount( inputfile='%(filename_sample)s', 
                          outputfile='%(tmpdir)s/basecount', 
                          extension=%(fragment_size)i, 
                          filetype='bed', 
                          twoBitFile='%(bit_filename)s' )
                          '''  % locals() )
    
    E.info( "computing peaks" )

    R( '''zinba( refinepeaks=1, 
                 seq='%(filename_sample)s', 
                 input='%(filename_control)s', 
                 filetype='bed',  
                 align='%(mappability_dir)s',
                 twoBit='%(bit_filename)s', 
                 outfile='%(filename_output)s', 
                 extension=%(fragment_size)s, 
                 basecountfile='%(tmpdir)s/basecount', 
                 numProc=%(threads)i, 
                 threshold=%(fdr_threshold)f, 
                 broad=FALSE, 
                 printFullOut=0, 
                 interaction=FALSE, 
                 mode='peaks', 
                 FDR=TRUE) '''  % locals() )
    
    shutil.rmtree( tmpdir )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


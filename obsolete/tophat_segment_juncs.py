#! /bin/env python
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
tophat_segment_juncs - speeding up tophat's segment_juncs
=====================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

A bottle neck in the tophat pipeline is :file:`segment_juncs`. This script
serves as a drop-in. It splits the input files by chromosome/contig
and calls :file:`segment_juncs` for contig separately and in parallel.

Other than renaming the original executable :file:`segment_juncs` to
:file:`segment_juncs.original` and linking to this script as
:file:`segment_juncs` no modification of the tophat pipeline
is required.

This script uses 8 threads.

.. note:: 
   The original :file:`segment_juncs` is still required and should
   be renamed to :file:`segment_juncs.original` and need to reside
   within the user`s path.

Usage
-----

No usage - is called within tophat.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse
import subprocess
import signal
import glob
import argparse

import multiprocessing.pool

import CGAT.Experiment as E
import CGAT.IOTools as IOTools

# If true, the original segment_juncs will be called without splitting
# the input data.
DISABLE = False

def subprocess_setup():
    # Python installs a SIGPIPE handler by default, which causes
    # gzip or other de/compression pipes to complain about "stdout: Broken pipe"
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)

def runCommand( args ):

    cmd, logfile = args
    log = open(logfile, "w" )
    log.write("# %s" % " ".join(cmd))
    retcode = subprocess.call( cmd,
                               preexec_fn=subprocess_setup,
                               stderr=log)

    if retcode != 0:
        die( "Error: segment-based junction search failed with err ="+str(retcode))

    log.close()

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if DISABLE:
        print "# tophat_segment_juncs.py disabled"
        argv[0] = "segment_juncs.original"
        runCommand( argv , "segment_juncs.log" )
        return 0

    E.Start( no_parsing = True )

    # collect arguments
    parser = argparse.ArgumentParser(description='Process tophat options.')
    parser.add_argument('-p', '--num-threads', metavar='N', type=int, dest='nthreads',
                         help='number of threads')
    parser.add_argument('--version', action='version', version='%(prog)s')
    options, args = parser.parse_known_args( argv[1:] )

    E.info( "parallelizing segment juncs with %i threads" % options.nthreads )
    
    x = argv.index("--ium-reads") + 1
    
    all_options = argv[1:x]

    (input_missing_reads, input_genome, 
     output_junctions, 
     output_insertions, output_deletions,
     input_left_all_reads,
     input_left_all_map,
     input_left_segments_maps ) = argv[x:x + 8]

    input_left_segments_maps = input_left_segments_maps.split(",")

    if len(argv) > x + 8:
        ( input_right_all_reads,
          input_right_all_map,
          input_right_segments_maps ) = argv[x+8:x+11]
        input_right_segments_maps = input_right_segments_maps.split(",")
    else:
        input_right_all_reads = ""
        input_right_all_map = ""
        input_right_segments_maps = []

    keys = set()
    
    # some filenames might appear multiple times
    files_to_split = set([input_left_all_map, \
                              input_right_all_map ] +\
                             input_left_segments_maps +\
                             input_right_segments_maps )

    E.info( "splitting %i files" % len(files_to_split))

    ## split all map files by chromosome
    for filename in files_to_split:
        if filename == "": continue
        E.info("splitting %s" % filename )
        base, ext = os.path.splitext( filename )

        f = glob.glob( "%s.input.*%s" % (filename, ext) )
        if f:
            E.info("files already exist - skipping" )
            keys.update( [ re.match("%s.input.(\S+)%s" % (filename,ext), x ).groups()[0] for x in f ] )
            continue
        
        infile = IOTools.openFile( filename )

        outfiles = IOTools.FilePool( filename + ".input.%s" + ext )

        for line in infile:
            key = line.split("\t")[2]
            keys.add( key )
            outfiles.write( key, line )

        outfiles.close()

    # keys = set( ["chr1", "chr2", "chr3", "chr4", "chr5",
    #              "chr6", "chr7", "chr8", "chr9", "chr10",
    #              "chr11", "chr12", "chr13", "chr14", "chr15",
    #              "chr16", "chr17", "chr18", "chr19", "chr20",
    #              "chr21", "chr22", "chrX", "chrY", "chrM" ] )

    E.info( "working on %i contigs: %s" % (len(keys), list(keys)))

    pool = multiprocessing.pool.ThreadPool( options.nthreads )
    #pool = threadpool.ThreadPool( THREADS )

    tmpdir = os.path.dirname( input_left_all_reads )
    logdir = os.path.join( tmpdir[:-len("tmp")], "logs" )

    if not os.path.exists(logdir):
        raise IOError( "can not find logdir %s" % logdir )

    args = []
    for key in keys:

        def modout( old, key ):
            if not old:return ""
            _, ext = os.path.splitext( old )
            return old + ".output.%s%s" % (key, ext)

        def modin( old, key ):
            if not old:return ""
            _, ext = os.path.splitext( old )
            return old + ".input.%s%s" % (key,ext)

        def modgenome( old, key ):
            dirname, filename = os.path.split(old)
            genome, ext = os.path.splitext( filename )
            if genome.lower().endswith("_cs"): genome = genome[:-3]
            new = os.path.join( dirname, genome + ".perchrom", key + ext )
            if not os.path.exists(new):
                raise ValueError( "can not find chromoseme file %s" % new )
            return new

        cmd = ["segment_juncs"] +\
            all_options +\
            [input_missing_reads,  \
                 modgenome(input_genome,key), \
                 modout(output_junctions,key),\
                 modout(output_insertions,key),\
                 modout(output_deletions,key),\
                 input_left_all_reads,\
                 modin( input_left_all_map, key ),\
                 ",".join( [ modin( x, key ) for x in input_left_segments_maps ] ),\
                 input_right_all_reads,\
                 modin( input_right_all_map, key ),\
                 ",".join( [ modin( x, key ) for x in input_right_segments_maps ] ) ]


        logfile = os.path.join(logdir, "segment_juncs_%s.log" % key )
        args.append( (cmd,logfile) )

    E.info( "submitting %i jobs" % len(keys) )

    pool.map( runCommand, args, chunksize = 1 )
    pool.close()
    pool.join()

    E.info("all jobs finished successfully" )

    E.info("merging results")
    ## merge results
    for filename in (output_junctions, output_insertions, output_deletions):
        outfile = open(filename, "w")
        for inf in glob.glob( filename + ".output.*" ):
            infile = open( inf, "r" )
            outfile.write( infile.read() )
            infile.close()
        outfile.close()
        
    E.info("results merged")

    ## cleaning up is done automatically by tophat
    E.info("cleaning up" )
    for f in glob.glob( os.path.join( tmpdir, "*.output.*") ) +\
            glob.glob( os.path.join( tmpdir, "*.input.*") ):
        os.remove(f)

    ## write footer and output benchmark information.
    E.Stop()
    
if __name__ == "__main__":
    sys.exit( main( sys.argv) )


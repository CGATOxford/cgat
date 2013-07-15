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
gtf2reads.py - sample reads from genes
======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script converts a :term:`gtf` formatted file into a :term:`fasta` formatted file of 
reads. The reads are sampled from the gene structures defined in the file.

Reads are sampled from a normal distribution with (read_length_mean,
read_length_stddev).

There are two modeling schemas:
   1 uniform expression: the sampling aims at a given
     statistical coverage. Each position in a transcript
     is equally likely to be a starting point of a sequencing
     read. The sampling procedure::

         1. estimate lvsequence = lsequence * gauss( coverage_mean, coverage_stddev)
         Sample until sum of residues in reads > lvsequence

            1. positions of read: uniformly from (1/lsequence)
            2. length of read: drawn from N(<read length>,stddev(read_length))
            3. strand of read: uniformly from (1/2)
            4. reject read if less than minimum size (min_read_length)

   Note that the effective read length and coverage will lie below the input
   values, as read lengths are constricted by the transcript length.

   2 sample according to expression profile
      TO BE IMPLEMENTED

For cross-species comparison, a mutation rate can be applied
from a normal distribution with (mutation_mean, mutation_stddev).


Usage
-----

Example::

   python gtf2reads.py --help

Type::

   python gtf2reads.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import optparse
import types
import random
import CGAT.GTF as GTF
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Genomics as Genomics
import CGAT.SequencePairProperties as SequencePairProperties
import CGAT.Iterators as Iterators

import numpy
import numpy.random
import random

def getMutatedSequence( sequence, divergence ):
    """sample number of events from a Poisson distribution.

    This is a very hacky, simple mutator that does not take into account
    multiple substitutions per site and/or a substitution matrix.
    Only use for small divergences.
    """

    lsequence = len(sequence)
    nmutate = numpy.random.poisson( float(lsequence) * divergence )
    sequence = list(sequence.upper())

    for pos in random.sample(xrange(lsequence), nmutate):
        c = sequence[pos]
        x = c
        while x == c: x = random.choice( "ACGT" )
        sequence[pos] = x

    return "".join(sequence)

##------------------------------------------------------------
def main():

    parser = E.OptionParser( version = "%prog version: $Id: gtf2reads.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option( "-p", "--output-filename-pattern", dest="output_filename_pattern", type="string" ,
                       help="OUTPUT filename with histogram information on aggregate coverages [%default].")

    parser.add_option( "--read-length-mean", dest="read_length_mean", type="float",
                      help="simulation parameter [default=%default]."  )

    parser.add_option( "--read-length-std", dest="read_length_stddev", type="float",
                      help="simulation parameter [default=%default]."  )

    parser.add_option( "--coverage-mean", dest="coverage_mean", type="float",
                      help="simulation parameter [default=%default]."  )

    parser.add_option( "--coverage-std", dest="coverage_stddev", type="float",
                      help="simulation parameter [default=%default]."  )

    parser.add_option( "--ds-mean", dest="ds_mean", type="float",
                      help="simulation parameter [default=%default]."  )

    parser.add_option( "--ds-std", dest="ds_stddev", type="float",
                      help="simulation parameter [default=%default]."  )

    parser.add_option( "--error-mean", dest="error_mean", type="float",
                      help="simulation parameter [default=%default]."  )

    parser.add_option( "--error-std", dest="error_stddev", type="float",
                      help="simulation parameter [default=%default]."  )

    parser.add_option( "--min-read-length", dest="min_read_length", type="int",
                      help="minimum read length [default=%default]."  )

    parser.add_option( "--sample-size", dest="sample_size", type="int",
                      help="randomly sample from selected transcripts [default=%default]."  )

    parser.add_option( "--test", dest="test", type="int",
                      help="test with # first entries [default=%default]."  )

    parser.add_option( "--mode", dest="mode", type="choice",
                       choices=("genes", "transcripts" ),
                       help="use genes or transcripts [default=%default]."  )



    parser.set_defaults(
        genome_file = None,
        read_length_mean = 200.0,
        read_length_stddev = 20.0,
        coverage_mean = 2.0,
        coverage_stddev = 1.0,
        ds_mean = None,
        ds_stddev = None,
        error_mean = None,
        error_stddev = None,
        min_read_length = 50,
        test = None,
        mode = "transcripts",
        output_filename_pattern = None,
        output_format_id = "%010i",
        sample_size = 0,
        )

    (options, args) = E.Start( parser )
    
    assert options.genome_file, "please supply an indexed genome." 

    if options.output_filename_pattern:
        outfile_stats = open( options.output_filename_pattern % "stats", "w" )
        outfile_stats.write( "id\tlen\tnreads\tlen_mean\tlen_std\tcov_mean\tcov_std\n" )
        outfile_map = open( options.output_filename_pattern % "map", "w" )
        outfile_map.write( "id\ttranscript\n" )
    else:
        outfile_stats = None
        outfile_map = None

    genome = IndexedFasta.IndexedFasta( options.genome_file )

    ninput, noutput, nskipped = 0, 0, 0

    total_counts, total_read_lengths, total_len = [], [], 0
    total_pids = []
    total_error_pids = []

    if options.mode == "transcripts":
        iterator = GTF.transcript_iterator( GTF.iterator_filtered( GTF.iterator(options.stdin), feature="exon" ))
        getId = lambda x: x.transcript_id
    elif options.mode == "genes":
        iterator = GTF.flat_gene_iterator( GTF.iterator_filtered( GTF.iterator(options.stdin), feature="exon" ))
        getId = lambda x: x.gene_id

    if options.sample_size:
        iterator = Iterators.sample( iterator ) 
        
    if options.ds_mean:
        do_mutate = True
        pid_calc = SequencePairProperties.SequencePairPropertiesPID()
    else:
        do_mutate = False

    if options.error_mean:
        do_error = True
        pid_calc = SequencePairProperties.SequencePairPropertiesPID()
    else:
        do_error = False

    for gtfs in iterator:

        id = getId(gtfs[0])

        try:
            sequence = GTF.toSequence( gtfs, genome )
        except KeyError, msg:
            if options.loglevel >= 2:
                options.stdlog.write( "# skipping %s: %s\n" % (id, msg))
            nskipped += 1
            continue

        lsequence = len(sequence)
        
        if lsequence <= options.min_read_length * 2:
            if options.loglevel >= 2:
                options.stdlog.write( "# skipping %s - sequence is too short: %i\n" % (id, lsequence) )
            nskipped += 1
            continue

        ninput += 1

        if do_mutate:
            new_sequence = getMutatedSequence( sequence, options.ds_mean )
            pid_calc.loadPair( sequence, new_sequence)
            pid = pid_calc.mPID 
            total_pids.append( pid )
            sequence = new_sequence
        else:
            pid = 100.0

        if options.loglevel >= 2:
            options.stdlog.write("# processing %s - len=%i\n" % (id, lsequence) )
            options.stdlog.flush()

        total_len += lsequence
        lvsequence = lsequence * random.gauss( options.coverage_mean, options.coverage_stddev )
        
        covered = 0
        counts = numpy.zeros( lsequence )
        nreads = 0

        error_pids, read_lengths = [], []
        
        while covered < lvsequence:
            
            read_length = int(random.gauss( options.read_length_mean, options.read_length_stddev ))
            positive = random.randint(0,1)

            if positive:
                start = random.randint(0,lsequence)
                end = min( lsequence, start + read_length )
            else:
                end = random.randint(0,lsequence)
                start = max(0, end - read_length)
                
            read_length = end - start
            if read_length < options.min_read_length:
                continue
            
            segment = sequence[start:end]

            if not positive:
                segment = Genomics.complement( segment )

            noutput += 1

            if do_error:
                new_segment = getMutatedSequence( segment, options.error_mean )
                pid_calc.loadPair( segment, new_segment)
                pid = pid_calc.mPID 
                error_pids.append( pid )
                segment = new_segment
            else:
                pid = 100.0
                
            options.stdout.write( ">%s\n%s\n" % (options.output_format_id % noutput, segment ) )
            
            if outfile_map:
                outfile_map.write( "%s\t%s\n" % (id, options.output_format_id % noutput) )

            for x in range(start,end):
                counts[x] += 1

            nreads += 1

            covered += read_length
            read_lengths.append(read_length)
            
        if options.loglevel >= 2:
            options.stdout.write( "# transcript %s: len=%i, nreads=%i, len_mean=%.2f, len_std=%.2f, cov_mean=%.2f, cov_stddev=%.2f\n" % (id, 
                                                                                                                                         lsequence,
                                                                                                                                         nreads, 
                                                                                                                                         numpy.mean( read_lengths ),
                                                                                                                                         numpy.std( read_lengths ),
                                                                                                                                         numpy.mean( counts ),
                                                                                                                                         numpy.std( counts) ))
        if outfile_stats:
            outfile_stats.write( "%s\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\n" % (id, 
                                                                           lsequence,
                                                                           nreads, 
                                                                           numpy.mean( read_lengths ),
                                                                           numpy.std( read_lengths ),
                                                                           numpy.mean( counts ),
                                                                           numpy.std( counts) ))


        total_counts += list(counts)
        total_read_lengths += read_lengths
        total_error_pids += error_pids

        if options.test and ninput >= options.test:
            break

        if options.sample_size and ninput >= options.sample_size:
            break

    if options.loglevel >= 1:
        output = ["len=%i, nreads=%i" % ( total_len,
                                          noutput ) ]
        output.append( "len_mean=%.2f, len_std=%.2f, cov_mean=%.2f, cov_stddev=%.2f" % ( 
                numpy.mean( total_read_lengths ),
                numpy.std( total_read_lengths ),
                numpy.mean( total_counts ),
                numpy.std( total_counts) ))

        no_uncovered = [ x for x in total_counts if x > 0 ]

        output.append( "cov0_mean=%.2f, cov0_stddev=%.2f" % (numpy.mean( no_uncovered ),
                                                             numpy.std( no_uncovered ) ) )

        if do_mutate:
            output.append( "pid_mean=%.2f, pid_std=%.2f" % (numpy.mean( total_pids ), numpy.std( total_pids ) ))

        if do_error:
            output.append( "pid_error_mean=%.2f, pid_error_std=%.2f" % (numpy.mean( total_error_pids ), numpy.std( total_error_pids ) ))
            
            
        options.stdlog.write( "# effective: %s\n" % ", ".join( output ))

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped ) )

    E.Stop()

if __name__ == '__main__':
    main()

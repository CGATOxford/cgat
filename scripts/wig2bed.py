####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: bed2bed.py 2861 2010-02-23 17:36:32Z andreas $
##
##
####
####
'''
wig2bed.py - convert densities to intervals
===========================================

Purpose
-------

define intervals based on densities within a bigwig file.

The script currently implements the following methods (``--method``):

threshold
    output windows that contain values above a certain
    threshold.

std-above-mean
    output windows that are a certain number of standard deviations 
    above the mean. 

multiple-of-mean
    output windows that are a certain times above the mean.

Usage
-----

Bigwig files need to be supplied by the --bigwig-file options.

For example::

    python wig2bed.py --threshold=10 --method=threshold --genome-file=mm10 --bigwig-file=in.bw > out.bed

Command line options
--------------------

'''   

import sys
import re
import string
import optparse
import time
import os
import itertools
import tempfile
import subprocess
import shutil
import collections

import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools

import bx
import bx.bbi.bigwig_file

def applyThreshold( infile, fasta, threshold, max_distance = 0 ):
    '''apply threshold to a wig file writing a
    bed-formatted file as output.'''
    
    c = E.Counter()

    for contig, size in fasta.getContigSizes( with_synonyms = False ).items():
        c.contigs += 1

        E.debug( "processing %s" % contig )
        
        iterator = infile.get( contig, 0, size )
        if iterator == None: 
            c.skipped += 1
            E.warn( "skipping %s" % contig )
            continue

        last_start, last_end = -1, 0

        for start, end, value in iterator:
            d = start - last_end
            if (d > 0 or value < threshold): 
                if last_start >= 0:
                    yield contig, last_start, last_end
                    c.intervals += 1
                last_start = -1
            elif last_start < 0 and value >= threshold:
                last_start = start

            last_end = end

        if last_start >= 0:
            yield contig, last_start, end
            c.intervals += 1

        c.output += 1

    E.info( str(c) )

def getBigwigSummary( bigwig_file ):
    '''return summary of bigwig contents.

    This method uses the bigWigInfo UCSC utility
    '''
    
    results = E.run( "bigWigInfo %(bigwig_file)s" % locals(), return_stdout = True )

    data = [ x.split( ":") for x in results.split( "\n" ) if x != "" ]
    fields = [ x[0] for x in data]
    Results = collections.namedtuple( "BigwigInfo", fields )
    def conv( v ):
        return IOTools.str2val( re.sub( ",", "", v.strip()) )

    results = Results( *[ conv(x[1]) for x in data] )
    return results

def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id: bed2bed.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-m", "--method", dest="methods", type="choice", action="append",
                       choices=("threshold", "stddev-above-mean", "multiple-of-mean"),
                       help="method to apply [default=%default]"  )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("-t", "--threshold", dest="threshold", type="float",
                      help="threshold to apply [default=%default]"  )

    parser.add_option("-i", "--bigwig-file", dest="bigwig_file", type="string", metavar="bigwig",
                      help="filename with bigwig information [default=%default]."  )

    parser.add_option("-b", "--bigwig", dest="bigwig", action = "store_true",
                      help="input is bigwig [default=%default]."  )

    parser.set_defaults( methods = [],
                         genome_file = None,
                         threshold = 10,
                         bigwig = False,
                         max_distance = 0)
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    if options.bigwig_file:
        bigwig_file = bx.bbi.bigwig_file.BigWigFile( open(options.bigwig_file ) )
    else:
        bigwig_file = None

    if options.genome_file:
        genome_fasta = IndexedFasta.IndexedFasta( options.genome_file )
        contigs = genome_fasta.getContigSizes()

    for method in options.methods:
        if method ==  "threshold":
            if not contigs: raise ValueError("please supply contig sizes" )
            if not bigwig_file: raise NotImplementedError( "threshold not implemented for wig files")
            processor = applyThreshold( bigwig_file, 
                                        genome_fasta, 
                                        threshold = options.threshold,
                                        max_distance = options.max_distance )
        elif method ==  "stddev-above-mean":
            if not contigs: raise ValueError("please supply contig sizes" )
            if not bigwig_file: raise NotImplementedError( "threshold not implemented for wig files")
            summary = getBigwigSummary( options.bigwig_file )
            threshold = summary.mean + options.threshold * summary.std
            E.info( "applying threshold %f: mean=%f, std=%f" % (threshold, summary.mean, summary.std) )
            processor = applyThreshold( bigwig_file, 
                                        genome_fasta, 
                                        threshold = threshold,
                                        max_distance = options.max_distance )

        elif method ==  "multiple-of-mean":
            if not contigs: raise ValueError("please supply contig sizes" )
            if not bigwig_file: raise NotImplementedError( "threshold not implemented for wig files")
            summary = getBigwigSummary( options.bigwig_file )
            threshold = summary.mean * options.threshold 
            E.info( "applying threshold %f: mean=%f, std=%f" % (threshold, summary.mean, summary.std) )
            processor = applyThreshold( bigwig_file, 
                                        genome_fasta, 
                                        threshold = threshold,
                                        max_distance = options.max_distance )


    outfile = options.stdout
    
    outfile.write( "".join( [ "%s\t%i\t%i\n" % x for x in processor ] ) )
    outfile.write( "\n" )

    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv ) )

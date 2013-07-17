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
wig2wig.py - manipulate wiggle files
====================================

Purpose
-------

manipulate wig-formatted files.

The script currently implements the following methods:

1. sanitize-genome: remove all empty intervals and intervals on unknown contigs. 
   Intervals extending beyond a contig a truncated. 
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

import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools

def sanitizeGenome( infile, outfile, contigs ):
    """truncate bed intervals that extend beyond contigs.

    removes empty intervals (start == end).

    throws an error if start > end.
    """
    
    ninput, noutput = 0, 0
    ntruncated_contig, nskipped_contig, nskipped_empty = 0, 0, 0

    for line in infile:
        if line.startswith( "track"):
            outfile.write(line)
            continue
        elif line.startswith( "variableStep" ):
            contig = re.search( "chrom=(\S+)", line ).groups()[0]
            span = int(re.search( "span=(\d+)", line ).groups()[0])
            
            if contig not in contigs:
                nskipped_contig += 1
                skip_contig = True
            else:
                skip_contig = False
                outfile.write(line)
                contig_size = contigs[contig]

            continue

        if skip_contig: continue
    
        start, val = map(int, line[:-1].split("\t") )
        
        if start + span > contig_size:
            ntruncated_contig += 1
            continue
        
        outfile.write( line )

    E.info( "ninput=%i, noutput=%i, nskipped_contig=%i, ntruncated=%i, nskipped_empty=%i" % \
                (ninput, noutput, nskipped_contig, ntruncated_contig, nskipped_empty) )


def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id: bed2bed.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage = globals()["__doc__"] )


    parser.add_option( "-m", "--method", dest="methods", type="choice", action="append",
                       choices=("merge", "filter-genome", "bins", "block", "sanitize-genome", "shift", "extend" ),
                       help="method to apply [default=%default]"  )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.set_defaults( methods = [],
                         merge_distance = 0,
                         binning_method = "equal-bases",
                         genome_file = None,
                         bam_file = None,
                         num_bins = 5,
                         merge_min_intervals = 1,
                         bin_edges = None,
                         offset = 10000,
                         test = None,
                         extend_distance=1000)
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    contigs = None

    if options.genome_file:
        genome_fasta = IndexedFasta.IndexedFasta( options.genome_file )
        contigs = genome_fasta.getContigSizes()

    for method in options.methods:
        if method ==  "filter-genome":
            if not contigs: raise ValueError("please supply contig sizes" )
            processor = filterGenome( processor, contigs )
        elif method == "sanitize-genome":
            if not contigs: raise ValueError("please supply contig sizes" )
            processor = sanitizeGenome( options.stdin, options.stdout, contigs )

    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv ) )

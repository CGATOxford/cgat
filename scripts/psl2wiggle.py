####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: psl2wiggle.py 2834 2009-11-24 16:11:23Z andreas $
##
##
####
####
'''
psl2wiggle.py - convert from psl to wiggle
==========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script converts from a :term:`psl` formatted
file to a :term:`wiggle` formatted file by stacking 
alignments on the target on top of each other.

This script uses the UCSC tools for bigwig/bigbed output.

Usage
-----

Example::

   python <script_name>.py --help

Type::

   python <script_name>.py --help

for command line help.

Documentation
-------------

Code
----

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
import CGAT.Stats as Stats
import CGAT.Blat as Blat
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools
import numpy

def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id: psl2wiggle.py 2834 2009-11-24 16:11:23Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option("-b", "--output-filename", dest="output_filename", type="string",
                      help="filename for output [default=%default]" )

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("bedgraph", "wiggle", "bigbed", "bigwig"),
                      help="output format [default=%default]" )

    parser.set_defaults( genome_file = None,
                         typecode = numpy.int16,
                         output_filename = None,
                         output_format = "wiggle",
                         test = None )

    (options, args) = E.Start( parser, add_pipe_options = True )

    typecode = options.typecode

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file )
        counts = {}
        contig_sizes = fasta.getContigSizes( with_synonyms = False )
        E.info("allocating memory for %i contigs and %i bytes" % (len(contig_sizes), sum(contig_sizes.values()) * typecode().itemsize ))
        for contig,size in contig_sizes.items():
            E.debug( "allocating %s: %i bases" % (contig, size ) )
            counts[contig] = numpy.zeros( size, typecode )

        E.info("allocated memory for %i contigs" % len(fasta))

    else:
        fasta = None
        contig_sizes = {}

    if options.output_format in ("bigwig", "bigbed"):
        
        if not options.genome_file:
            raise ValueError("please supply genome file for bigwig/bigbed computation.")

        if not options.output_filename:
            raise ValueError("please output file for bigwig/bigbed computation.")

        if options.output_format == "bigwig":
            executable_name = "wigToBigWig"
        elif options.output_format == "bigbed":
            executable_name = "bedToBigBed"
        else:
            raise ValueError("unknown output format `%s`" % options.output_format)

        executable = IOTools.which( executable_name )

        if not executable:
            raise OSError( "could not find %s in path." % executable_name )

        tmpdir = tempfile.mkdtemp()
        E.debug( "temporary files are in %s" % tmpdir)

        tmpfile_wig = os.path.join( tmpdir, "wig" )
        tmpfile_sizes = os.path.join( tmpdir, "sizes" )

        # write contig sizes
        outfile_size = open( tmpfile_sizes, "w")
        for contig, size in contig_sizes.items():
            outfile_size.write("%s\t%s\n" % (contig, size) )
        outfile_size.close()    
        
        outfile = open( tmpfile_wig, "w" )

    else:
        outfile = options.stdout

    iterator = Blat.BlatIterator( sys.stdin )

    ninput, ncontigs, nskipped = 0, 0, 0

    E.info( "started counting" )

    while 1:

        if options.test and ninput >= options.test:
            break

        match = iterator.next()
        
        if match == None: break
        
        ninput += 1

        contig = match.mSbjctId

        for start, length in zip( match.mSbjctBlockStarts, match.mBlockSizes):
            counts[contig][start:start+length] += 1

    E.info( "finished counting" )

    if options.output_format in ("wig", "bigwig"):
        E.info( "starting wig output" )

        for contig, vals in counts.items():

            E.debug("output for %s" % contig )
            for val, iter in itertools.groupby( enumerate( vals ), lambda x: x[1] ):
                l = list(iter)
                start,end = l[0][0],l[-1][0]
                val = vals[start]
                if val > 0:
                    outfile.write("variableStep chrom=%s span=%i\n" % (contig, end-start+1))
                    outfile.write("%i\t%i\n" % (start, val) )

            ncontigs += 1
    elif options.output_format in ("bedgraph", "bigbed"):
        
        E.info( "starting bedgraph output" )

        for contig, vals in counts.items():
            E.debug("output for %s" % contig )
            for val, iter in itertools.groupby( enumerate( vals ), lambda x: x[1] ):
                l = list(iter)
                start,end = l[0][0],l[-1][0]
                val = vals[start]
                if val > 0:
                    outfile.write("%s\t%i\t%i\t%i\n" % (contig, start, end+1,val))
            
            ncontigs += 1

    E.info( "finished output" )

    if options.output_format in ("bigwig", "bigbed"):
        outfile.close()

        E.info( "starting bigwig conversion" )
        try:
            retcode = subprocess.call( " ".join( (executable,
                                                 tmpfile_wig,
                                                 tmpfile_sizes,
                                                 os.path.abspath( options.output_filename )), ),
                                       shell=True)
            if retcode < 0:
                warn( "wigToBigWig terminated with signal: %i" % -retcode)
                return -retcode
        except OSError, msg:
            warn( "Error while executing bigwig: %s" % e)
            return 1

        shutil.rmtree( tmpdir )
        
        E.info( "finished bigwig conversion" )

    E.info( "ninput=%i, ncontigs=%i, nskipped=%i\n" % (ninput, ncontigs, nskipped) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )

'''
fastq2fastq.py - manipulate fastq files
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS Sequences FASTQ Manipulation

Purpose
-------

This script performs manipulations on :term:`fastq` formatted
files. For example it can be used to change the quality score format
or sample a subset of reads.

The script predominantly is used for manipulation of single fastq
files. However, for some of its functionality it will take paired data
using the --pair and --outfile-pair options. This applies to the
--sample and --sort options.

Usage
-----

Example::

   cat in.fastq.1 | python fastq2fastq.py --sample 0.5 --pair in.fastq.2 --outfile-pair out.fastq.2  > out.fastq.1

In this example we randomly sample 50% of reads from paired data provided in two files.

Type::

   python fastq2fastq.py --help

for command line help.


Command line options
--------------------

'''

import os
import sys
import re
import optparse
import math
import random
import itertools
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGAT.Fastq as Fastq

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-f", "--change-format", dest="change_format", type="choice",
                      choices = ('sanger', 'solexa', 'phred64', 'integer' ),
                      help="guess quality score format and set quality scores to format [default=%default]."  )

    parser.add_option( "--guess-format", dest="guess_format", type="choice",
                      choices = ('sanger', 'solexa', 'phred64', 'integer' ),
                      help="quality score format to assume if ambiguous [default=%default]."  )

    parser.add_option( "--sample", dest="sample", type="float",
                       help="sample a proportion of reads [default=%default]."  )

    parser.add_option( "--pair", dest="pair", type="string",
                       help="if data is paired, filename with second pair. "
                       "Implemented for sampling [default=%default]."  )

    parser.add_option( "--outfile-pair", dest="outfile_pair", type="string",
                       help="if data is paired, filename for second pair. "
                       "Implemented for sampling [default=%default]."  )

    parser.add_option( "--uniq", dest="uniq", action="store_true",
                       help="remove duplicate reads (by name) [default=%default]."  )

    parser.add_option( "--apply", dest="apply", type="string",
                       help="apply a filter to fastq file (taking only reads in filename) [default=%default]."  )

    parser.add_option( "--trim3", dest="trim3", type="int",
                       help="trim # bases from 3' end [default=%default]."  )

    parser.add_option( "--sort", dest="sort", action="store_true",
                       help="sort fastq by sequence id [default=%default]."  )

    parser.add_option( "--seed", dest="seed", type="int",
                       help="seed for random number generator [default=%default]."  )


    parser.set_defaults(
        change_format = None,
        guess_format = None,
        sample = None,
        trim3 = None,
        pair = None,
        apply = None,
        uniq = False,
        outfile_pair = None,
        sort = None,
        seed = None )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    c = E.Counter()

    if options.change_format:
        for record in Fastq.iterate_convert( options.stdin, 
                                             format = options.change_format,
                                             guess = options.guess_format ):
            c.input += 1
            options.stdout.write( "%s\n" % record )
            c.output += 1
                                                                   
    elif options.sample:
        sample_threshold = min( 1.0, options.sample)

        random.seed( options.seed )
        
        if options.pair:
            if not options.outfile_pair:
                raise ValueError( "please specify output filename for second pair (--outfile-pair)")

            outfile1 = options.stdout
            outfile2 = IOTools.openFile( options.outfile_pair, "w" )
            
            for record1, record2 in itertools.izip( Fastq.iterate( options.stdin ), Fastq.iterate( IOTools.openFile( options.pair) ) ):
                c.input += 1
                if random.random() <= sample_threshold:
                    c.output += 1
                    outfile1.write( "%s\n" % record1 )
                    outfile2.write( "%s\n" % record2 )

        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            if random.random() <= sample_threshold:
                c.output += 1
                options.stdout.write( "%s\n" % record )

    elif options.apply:
        ids = set(IOTools.readList( IOTools.openFile( options.apply ) ))
        
        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            if re.sub(" .*", "", record.identifier).strip() in ids:
                c.output += 1
                options.stdout.write( "%s\n" % record )

    elif options.trim3:
        trim3 = options.trim3
        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            record.trim( trim3 )
            options.stdout.write( "%s\n" % record )
            c.output += 1
            
    elif options.uniq:
        keys = set()
        for record in Fastq.iterate( options.stdin ):
            c.input += 1
            if record.identifier in keys: continue
            else: keys.add( record.identifier )
            options.stdout.write( "%s\n" % record )
            c.output += 1
        
    # Need to change this to incorporate both pairs
    elif options.sort:
        if not options.pair:
            # This is quicker for a single fastq file
            statement = "paste - - - - | sort -k1,1 -t ' ' | tr '\t' '\n'"
            os.system(statement)
        else:
            if not options.outfile_pair:
                raise ValueError( "please specify output filename for second pair (--outfile-pair)")
            E.warn("consider sorting individual fastq files - this is memory intensive")
            entries1 = {}
            entries2 = {}
            for record1, record2 in itertools.izip( Fastq.iterate( options.stdin ), Fastq.iterate( IOTools.openFile( options.pair) ) ):
                entries1[record1.identifier[:-2]] = (record1.seq, record1.quals)
                entries2[record2.identifier[:-2]] = (record2.seq, record2.quals)
            
            outfile1 = options.stdout
            outfile2 = IOTools.openFile(options.outfile_pair, "w")
            assert len(set(entries1.keys()).intersection(set(entries2.keys()))) == len(entries1), """paired files do not contain the same reads
                                                                                                     need to reconcile files"""
            for entry in sorted(entries1):
                outfile1.write("@%s/1\n%s\n+\n%s\n" % (entry, entries1[entry][0], entries1[entry][1]))
                outfile2.write("@%s/2\n%s\n+\n%s\n" % (entry, entries2[entry][0], entries2[entry][1]))

    ## write footer and output benchmark information.
    E.info( "%s" % str(c) )
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

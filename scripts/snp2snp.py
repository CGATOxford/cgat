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
snp2snp.py - manipulate lists of SNPs
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------


Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse
import glob

import pysam

import CGAT.Experiment as E

def validateSNPs( options, fastafile ):
    '''read SNPs in pileup format.'''

    callers = []

    for filename in options.filename_references:
        for f in glob.glob( filename ):
            callers.append( pysam.SNPCaller( pysam.Samfile( f, "rb" ), fastafile ) )

    if len(callers) == 0: 
        E.warning( "no transcript data available" )
    else:
        E.info( "validating against %i reference files" % (len(callers)))

    c = E.Counter()

    outf = options.stdout

    nfiles = len(callers)

    outf.write( "\t".join( \
            ("contig", "pos", "reference", "genotype", "consensus_quality", "genotype_quality", "mapping_quality", "coverage" ) ) )
    outf.write( "\tstatus\tcalls\tfiltered\tmin_coverage\tmax_coverage\tmin_quality\tmax_quality\t" )
    outf.write( "\t".join( ("genotypes", "consensus_qualities", "genotype_qualities", "mapping_qualities", "coverages" ) ) )
    outf.write( "\n" )

    min_coverage = options.min_coverage

    for line in options.stdin:

        c.input += 1

        if c.input % options.report_step == 0:
            E.info("iteration %i: %s" % (c.input, str(c)))

        data = line[:-1].split()

        contig, pos, reference, genotype, consensus_quality, genotype_quality, mapping_quality, read_depth = data[:8]
        pos, consensus_quality, genotype_quality, mapping_quality, read_depth = \
            map( int, (pos, consensus_quality, genotype_quality, mapping_quality, read_depth))

        if reference == "*":
            # todo: treat indels
            c.skipped += 1
            continue

        # go to 0-based coordinates
        pos -= 1
        calls = []
        for x, caller in enumerate(callers):
            try:
                calls.append( caller.call( contig, pos ) )
            except ValueError, msg:
                continue
        outf.write( "\t".join( map(str, 
                                   (contig, pos, reference, genotype, 
                                    consensus_quality, genotype_quality, mapping_quality, read_depth ) )))

        gt = [ x.genotype for x in calls ]
        cq = [ x.consensus_quality for x in calls]
        sq = [ x.snp_quality for x in calls]
        mq = [ x.mapping_quality for x in calls ]
        co = [ x.coverage for x in calls ]

        # filter calls by quality
        filtered_calls = [ x for x in calls if x.coverage >= min_coverage ] 

        nidentical = len( [x for x in filtered_calls if x.genotype == genotype ] )
        nreference = len( [x for x in filtered_calls if x.genotype == reference ] )
        ncalls = len(filtered_calls)

        if ncalls == 0:
            # no data
            status = "?"
        elif nidentical == ncalls:
            # ok
            status = "O"
        elif nidentical:
            # ambiguous
            status = "A"
        elif nreference == ncalls:
            # wildtype
            status = "W"
        else:
            # conflict
            status = "C"

        c[status] += 1

        outf.write( "\t%s\t%i\t%i" % ( status, len(calls), len(filtered_calls) ) )

        if len(calls) > 1:
            outf.write( "\t%i\t%i\t%i\t%i" % (
                    min( co ), max(co),
                    min( sq ), max(sq) ))
        else:
            outf.write( "\t" + "\t".join( (
                        "".join( map(str,co) ), "".join(map(str,co)),
                        "".join( map(str,sq) ), "".join(map(str,sq))) ) )

        outf.write( "\t%s" % ",".join( map(str, gt ) ) )
        outf.write( "\t%s" % ",".join( map(str, cq ) ) )
        outf.write( "\t%s" % ",".join( map(str, sq ) ) )
        outf.write( "\t%s" % ",".join( map(str, mq ) ) )
        outf.write( "\t%s" % ",".join( map(str, co ) ) )
        outf.write("\n")

        c.output += 1

    return c

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices = ("validate", ),
                      help="methods to apply [default=%default]."  )

    parser.add_option("-g", "--filename-genome", dest="filename_genome", type="string",
                      help="filename of faidx indexed genome [default=%default]."  )

    parser.add_option("-r", "--filename-reference", dest="filename_references", action="append",
                      help="reference filenames for validation step [default=%default]."  )

    parser.add_option( "--min-coverage", dest="min_coverage", type="int",
                      help="minimum coverage [default=%default]."  )

    parser.set_defaults(
        min_coverage = 0,
        method = None,
        filename_genome = None,
        filename_references = [],
        report_step = 10000,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if options.method == None:
        raise ValueError("please supply a method")

    if options.filename_genome != None:
        fastafile = pysam.Fastafile( options.filename_genome )    
    else:
        fastafile = None

    if options.method == "validate":
        counter = validateSNPs( options, fastafile )

    E.info( "%s" % str(counter) )
    
    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: bams2correlation.py 2861 2010-02-23 17:36:32Z andreas $
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
"""
bams2correlation.py
===================

:Author: Andreas Heger
:Release: $Id: bams2correlation.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Compare per base coverage between BAM files.

Usage
-----

Example::

   python bams2correlation.py in1.bam in2.bam -g reference.fasta

This command generates a tab delimited output with columns chromosome, base coordinate, number of overlapping reads
in in1.bam, and number of overlapping reads in in2.bam.

Type::

   python <script_name>.py --help

for command line help.

Documentation
-------------

This tools allows users to obtain per base coverage over the genome for one or more BAM files.  The output includes
all bases in the supplied reference fasta except those with no coverage in the input BAM or BAMs.

At present the --interval or -i option has not been implemented.  

Code
----

""" 

import os
import sys
import re
import optparse
import collections

import CGAT.Experiment as E
import pysam
import CGAT.IndexedFasta as IndexedFasta

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: bams2correlation.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-i", "--intervals", dest="filename_intervals", type="string",
                      help="filename with intervals to use [default=%default]."  )
    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.set_defaults(
        filename_intervals = None,
        genome_file = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) < 1: raise ValueError("please supply at least two BAM files.")

    samfiles = []
    for f in args:
        samfiles.append( pysam.Samfile( f, "rb" ) )
        
    if options.filename_intervals:
        raise NotImplementedError( "It is not yet possible to specify intervals of interest.  Repeat command without intervals option." )
    
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )

    options.stdout.write( "contig\tpos\t%s\n" % "\t".join(args))

    ninput, nskipped, noutput = 0, 0, 0
    for contig in fasta.getContigs():
        
        missing_contig = False

        positions = {}

        # lazy way: use dictionary
        for x, f in enumerate(samfiles): 
            try:
                i = f.pileup( contig )
            except ValueError:
                missing_contig = True
                break

            for v in i: 
                vp = v.pos
                if vp in positions:
                    positions[vp].append( v.n )
                else:
                    positions[vp] = [0] * x + [v.n]

            # fill with 0 those not touched in this file
            for p in positions.keys():
                if len(positions[p]) <= x:
                    positions[p].append( 0 )

        if missing_contig:
            nskipped += 1
            continue

        noutput += 1
        for pos in sorted(positions.keys()):
            vals = positions[pos]
            options.stdout.write( "%s\t%i\t%s\n" % (contig, pos, 
                                                    "\t".join(map(str, vals))))

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

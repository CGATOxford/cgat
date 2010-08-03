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
gff2psl.py - convert from gff to psl
====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts converts from a :term:`gff` formatted
file to a :term:`psl` formatted file.

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

import sys, re, string, optparse, time, os, glob

import Experiment
import IndexedFasta
import Blat
import Genomics
import GFF, GTF
import alignlib
import Intervals

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: gff2psl.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option( "--is-gtf", dest="is_gtf", action="store_true",
                      help="input is gtf."  )

    parser.add_option( "--no-header", dest="with_header", action="store_false",
                      help="do not output BLAT header [default=%default]."  )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("--input-filename-queries", dest="input_filename_queries", type="string",
                      help="fasta filename with queries [default=%default]."  )

    parser.add_option("--allow-duplicates", dest="allow_duplicates", action="store_true",
                      help="""permit duplicate entries. Adjacent exons of a transcript will still be merged [default=%default]."""  )

    parser.set_defaults( is_gtf = False,
                         genome_file = None,
                         with_header = True,
                         allow_duplicates = False,
                         test = None )
    
    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    if options.genome_file:
        genome_fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        genome_fasta = None

    if options.input_filename_queries:
        queries_fasta =  IndexedFasta.IndexedFasta( options.input_filename_queries )
    else:
        queries_fasta = None

    ninput, noutput, nskipped = 0, 0, 0

    if options.is_gtf:
        iterator = GTF.transcript_iterator( GTF.iterator_filtered( GTF.iterator( sys.stdin ), feature="exon" ), 
                                            strict = not options.allow_duplicates )
    else:
        iterator = GFF.joined_iterator( GFF.iterator(sys.stdin) )

    if options.with_header:
        options.stdout.write( Blat.Match().getHeader() + "\n" )

    for gffs in iterator:
        
        if options.test and ninput >= options.test:
            break

        ninput += 1

        result = alignlib.makeAlignmentBlocks()

        xstart = 0

        intervals = Intervals.combine( [ (gff.start, gff.end) for gff in gffs ] )

        for start, end in intervals:
            xend = xstart + end - start
            
            result.addDiagonal( xstart, xend,
                                start - xstart )
            xstart = xend

        entry = Blat.Match()
        entry.mQueryId = gff.transcript_id
        entry.mSbjctId = gff.contig
        entry.strand  = gff.strand

        if genome_fasta:
            if entry.mSbjctId in genome_fasta:
                entry.mSbjctLength = genome_fasta.getLength( entry.mSbjctId )
            else:
                entry.mSbjctLength = result.getColTo()

        if queries_fasta:
            if entry.mQueryId in queries_fasta:
                entry.mQueryLength = queries_fasta.getLength( entry.mQueryId )
        else:
            entry.mQueryLength = result.getRowTo()

        entry.fromMap( result )
        
        options.stdout.write (str(entry) + "\n" )
        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped) )

    Experiment.Stop()

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
bed2table.py - annotate intervals
=================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

compute sequence properties for genes of given by a gtf file and output them in 
tabular format.

Usage
-----

Example::

   python gtf2table.py --help

Type::

   python gtf2table.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os, sys, string, re, optparse, math, time, tempfile, subprocess, types, bisect, array, collections
import GFF, GTF, Bed
import Experiment as E
import IndexedFasta
import Stats
import SequenceProperties
import Genomics
import Intervals

import bx.intervals.io
import bx.intervals.intersection
import alignlib
import numpy
import IndexedGenome
import pysam

class Counter( object ):
    
    def __init__(self, *args, **kwargs ):
        pass

class CounterOverlap( Counter ):
    '''count overlap for each interval in tracks.'''

    def __init__(self, filename, *args, **kwargs ):
        
        assert filename != None, "please supply filename for CounterOverlap"

        Counter.__init__(self, *args, **kwargs )

        self.filename = filename

        E.info( "reading intervals from %s" % self.filename )

        self.index = Bed.readAndIndex( open( self.filename, "r"),
                                       per_track = True )
        
        E.info( "read intervals for %s tracks" % len(self.index) )

        self.tracks = self.index.keys()
        self.headers = []
        for track in self.tracks:
            self.headers.extend( ["%s_nover" % track, "%s_bases" % track] )
            
    
    def update( self, bed ):
        '''update internal counts.'''

        results = []
        for track in self.tracks:
            try:
                overlaps = [ (x[0],x[1]) for x in self.index[track][bed.contig].find( bed.start, bed.end ) ]
            except KeyError:
                overlaps = []

            results.append( (len(overlaps), 
                             Intervals.calculateOverlap( [(bed.start, bed.end),],
                                                         Intervals.combine( overlaps ) ) ) )

        self.data = results

    def __str__(self):
        '''output overlap of interval in *bed*'''
                        
        r = []
        for track, result in zip(self.tracks, self.data):
            r.append( "\t".join( (str(result[0]), str(result[1])) ) )

        return "\t".join(r)

##------------------------------------------------------------
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: gtf2table.py 2888 2010-04-07 08:48:36Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option("-q", "--quality-file", dest="quality_file", type="string",
                      help="filename with genomic base quality information [default=%default]."  )

    parser.add_option("-b", "--bam-file", dest="bam_files", type="string",
                      help="filename with read mapping information. Multiple files can be submitted in a comma-separated list [default=%default]."  )

    parser.add_option("-f", "--filename-bed", dest="filename_bed", type="string",
                      help="filename with extra bed file (for counter: overlap) [default=%default]."  )

    parser.add_option( "--filename-format", dest="filename_format", type="choice",
                       choices=("bed", "gff", "gtf" ),
                       help="format of secondary stream [default=%default]."  )

    parser.add_option( "--gff-source", dest="gff_sources", type="string", action="append",
                      help="restrict to this source in extra gff file (for counter: overlap) [default=%default]."  )

    parser.add_option( "--gff-feature", dest="gff_features", type="string", action="append",
                      help="restrict to this feature in extra gff file (for counter: overlap) [default=%default]."  )

    parser.add_option("-r", "--reporter", dest="reporter", type="choice",
                      choices=("genes", "transcripts" ),
                      help="report results for 'genes' or 'transcripts' [default=%default]."  )

    parser.add_option("-s", "--section", dest="sections", type="choice", action="append",
                      choices=("exons", "introns" ),
                      help="select range on which counters will operate [default=%default]."  )

    parser.add_option("-c", "--counter", dest="counters", type="choice", action="append",
                      choices=( "overlap", ),
                      help="select counters to apply [default=%default]."  )

    parser.add_option( "--add-gtf-source", dest="add_gtf_source", action="store_true",
                      help="add gtf field of source to output [default=%default]."  )

    parser.add_option( "--proximal-distance", dest="proximal_distance", type="int",
                      help="distance to be considered proximal to an interval [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        reporter = "genes",
        with_values = True,
        sections = [],
        counters = [],
        filename_gff = None,
        filename_format = "gtf",
        gff_features = [],
        gff_sources = [],
        add_gtf_source = False,
        proximal_distance = 10000, 
        bam_files = None,
        )

    (options, args) = E.Start( parser )

    # get files
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    if options.quality_file:
        quality = IndexedFasta.IndexedFasta( options.quality_file )
        quality.setTranslator( IndexedFasta.TranslatorBytes() )
    else:
        quality = None

    counters = []

    for c in options.counters:
        if c == "overlap":
            counters.append( CounterOverlap( filename = options.filename_bed,
                                             fasta=fasta,
                                             options = options) )
            

    options.stdout.write( "\t".join( [ "contig", "start", "end", "name" ] ) )
    for counter in counters: 
        options.stdout.write("\t%s" % "\t".join( counter.headers ) )
    options.stdout.write("\n")

    for bed in Bed.iterator(options.stdin):
        options.stdout.write( "\t".join( (bed.contig, 
                                          str(bed.start), 
                                          str(bed.end), 
                                          bed.mFields[0]) ) )
        for counter in counters: 
            counter.update(bed)
            options.stdout.write("\t%s" % str(counter) )
        options.stdout.write("\n")

    E.Stop()

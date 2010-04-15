################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: bed2graph.py 2861 2010-02-23 17:36:32Z andreas $
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

:Author: Andreas Heger
:Release: $Id: bed2graph.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

compute the overlap graph between two bed files.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os, sys, re, optparse

import Experiment as E
import IOTools
import Bed
import ncl
import numpy

class Counter:

    mPercentFormat = "%5.2f"

    def __init__(self):
        pass
        
    def getHeader( self ):
        h = []
        for a in ( "exons", "bases" ):
            for b in ("total", "ovl", "unique"):
                for c in ("1","2"):
                    h.append( "n" + a + "_" + b + c )
        for a in ( "exons", "bases" ):
            for b in ("ovl", "unique"):
                for c in ("1","2"):
                    h.append( "p" + a + "_" + b + c )

        return "\t".join(h)

    @E.cachedmethod
    def buildIndex( self, filename ):
        return Bed.readAndIndex( open(filename, "r"), simple = True )

    def _count( self, filename, idx ):
        '''count filename against idx.'''

        overlapping_genes = set()
        genes = set()

        # iterate over exons
        infile = open( filename, "r" )
        it = Bed.bed_iterator( infile )

        nexons, nexons_overlapping = 0, 0
        nbases, nbases_overlapping = 0, 0
        for this in it:
            nexons += 1
            nbases += this.end - this.start 
            
            try:
                intervals = list(idx[this.contig].find( this.start, this.end ))
            except KeyError:
                continue
            
            if len(intervals) == 0: 
                continue

            nexons_overlapping += 1
            start, end = this.start, this.end 
            counts = numpy.zeros( end - start, numpy.int )
            for other_start,other_end,other_value in intervals:
                for x in range( max(start, other_start ) - start, min(end, other_end ) - start ): 
                    counts[x] += 1
            nbases_overlapping += sum( [1 for x in counts if x > 0 ] )

        infile.close()

        return nexons, nexons_overlapping, nbases, nbases_overlapping

    def count( self, filename1, filename2 ):
        """count overlap between two gtf files."""

        E.info( "counting started for %s versus %s" % (filename1, filename2))

        idx2 = self.buildIndex( filename2 )

        (self.mExons1, self.mExonsOverlapping1,
         self.mBases1, self.mBasesOverlapping1 ) = \
         self._count( filename1, idx2 )

        self.mExonsUnique1 = self.mExons1 - self.mExonsOverlapping1
        self.mBasesUnique1 = self.mBases1 - self.mBasesOverlapping1

        idx1 = self.buildIndex( filename1 )

        (self.mExons2, self.mExonsOverlapping2,
         self.mBases2, self.mBasesOverlapping2 ) = \
         self._count( filename2, idx1 ) 

        self.mExonsUnique2 = self.mExons2 - self.mExonsOverlapping2
        self.mBasesUnique2 = self.mBases2 - self.mBasesOverlapping2

    def __str__(self):

        return "\t".join( map( str, (
                        self.mExons1, self.mExons2,
                        self.mExonsOverlapping1, self.mExonsOverlapping2,
                        self.mExonsUnique1, self.mExonsUnique2,
                        self.mBases1, self.mBases2,
                        self.mBasesOverlapping1, self.mBasesOverlapping2,
                        self.mBasesUnique1, self.mBasesUnique2 ) ) ) + "\t" +\
              "\t".join( map( lambda x: IOTools.prettyPercent( *x), (
                    (self.mExonsOverlapping1, self.mExons1), 
                    (self.mExonsOverlapping2, self.mExons2),
                    (self.mExonsUnique1, self.mExons1), 
                    (self.mExonsUnique2, self.mExons2),
                    (self.mBasesOverlapping1, self.mBases1), 
                    (self.mBasesOverlapping2, self.mBases2),
                    (self.mBasesUnique1, self.mBases1), 
                    (self.mBasesUnique2, self.mBases2)  ) ) )

class CounterTracks(Counter):

    def __init__(self, filename ):
        self.mIndices = Bed.readAndIndex( open(filename, "r"), 
                                          simple = True,
                                          per_track = True )
        
        
    def getTracks( self ):
        return sorted(self.mIndices.keys())

    def _countIndices( self, idx_in, idx ):
        '''count filename against idx.'''

        overlapping_genes = set()
        genes = set()

        # iterate over exons

        nexons, nexons_overlapping = 0, 0
        nbases, nbases_overlapping = 0, 0
        for contig, ix in idx_in.iteritems():

            # note: add a findall function to ncl
            for start, end, value in ix.find(0, 1000000000):
                nexons += 1
                nbases += end - start
            
                try:
                    intervals = list(idx[contig].find( start, end ))
                except KeyError:
                    continue
            
                if len(intervals) == 0: 
                    continue

                nexons_overlapping += 1
                counts = numpy.zeros( end - start, numpy.int )
                for other_start,other_end,other_value in intervals:
                    for x in range( max(start, other_start ) - start, min(end, other_end ) - start ): 
                        counts[x] += 1
                nbases_overlapping += sum( [1 for x in counts if x > 0 ] )

        return nexons, nexons_overlapping, nbases, nbases_overlapping

    def count( self, filename, track ):
        """count overlap between two gtf files."""

        E.info( "counting started for %s versus %s" % (filename, track))

        (self.mExons1, self.mExonsOverlapping1,
         self.mBases1, self.mBasesOverlapping1 ) = \
         self._count( filename, self.mIndices[track] )

        self.mExonsUnique1 = self.mExons1 - self.mExonsOverlapping1
        self.mBasesUnique1 = self.mBases1 - self.mBasesOverlapping1

        idx = self.buildIndex( filename )

        # count index against index
        (self.mExons2, self.mExonsOverlapping2,
         self.mBases2, self.mBasesOverlapping2 ) = \
         self._countIndices( self.mIndices[track], idx )

        self.mExonsUnique2 = self.mExons2 - self.mExonsOverlapping2
        self.mBasesUnique2 = self.mBases2 - self.mBasesOverlapping2

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: bed2graph.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-s", "--ignore-strand", dest="ignore_strand", action="store_true",
                      help="ignore strand information [default=%default]." )

    parser.add_option("-o", "--output", dest="output", type="choice",
                      choices = ("full", "name"),
                      help="output either ``full`` overlapping entries, only the ``name``s. [default=%default]." )

    parser.add_option("-p", "--pattern-id", dest="pattern_id", type="string",
                      help="pattern to convert a filename to an id [default=%default]." )

    parser.add_option("-t", "--tracks", dest="tracks", action="store_true",
                      help="compare files against all tracks in the first file [default=%default]" )
    
    parser.set_defaults(
        ignore_strand = False,
        filename_update = None,
        pattern_id = "(.*).bed",
        output = "full",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 2:
        raise ValueError( "two arguments required" )

    if args[0] == "-":
        infile1 = options.stdin
    else:
        infile1 = open( args[0], "r")

    infile2 = open(args[1], "r")        
    
    idx = Bed.readAndIndex( infile2, with_values = True )

    output = options.output
    outfile = options.stdout
    
    if output == "name":
        outfile.write( "name1\tname2\n" )
        outf = lambda x: x.mFields[0]
    else:
        outf = str
        
    for bed in Bed.iterator( infile1 ):
        try:
            overlaps = idx[bed.contig].find( bed.start, bed.end )
        except (KeyError, IndexError):
            # ignore missing contig and zero length intervals
            continue

        for o in overlaps:
            outfile.write( "\t".join( (outf(bed), outf(o[2]))) + "\n" )

    E.Stop()


if __name__ == "__main__":
    sys.exit( main( sys.argv) )

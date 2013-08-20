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
gtfs2graph.py - compute overlap graph between genes
===================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script compares two :term:`gtf` output files and outputs a bi-partite
graph connecting overlapping genes/exons/...

This script needs some work, only gene-gene output works.

Usage
-----

Example::

   python gtfs2graph.py --help

Type::

   python gtfs2graph.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import optparse
import collections

USAGE="""python %s [OPTIONS] input1 input2


Version: $Id: gtfs2graph.py 2781 2009-09-10 11:33:14Z andreas $
""" % sys.argv[0]


import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import bx.intervals.io
import bx.intervals.intersection
import numpy

class Counter:

    mPercentFormat = "%5.2f"

    def __init__(self, outfile):
        self.mOutfile = outfile
        
    def write( self, *args):
        self.mOutfile.write( "\t".join( args ) + "\n" )

    def getHeader( self ):
        h = []
        for a in ("genes", "exons", "bases" ):
            for b in ("total", "ovl", "unique"):
                for c in ("1","2"):
                    h.append( "n" + a + "_" + b + c )
        for a in ("genes", "exons", "bases" ):
            for b in ("ovl", "unique"):
                for c in ("1","2"):
                    h.append( "p" + a + "_" + b + c )

        return "\t".join(h)

    @E.cachedmethod
    def buildIndex( self, filename ):
        """read and index."""

        idx = {}
        infile = open(filename, "r")
        for e in GTF.readFromFile(infile):
            if e.contig not in idx: idx[e.contig] = bx.intervals.intersection.Intersecter()
            idx[e.contig].add_interval( bx.intervals.Interval(e.start,e.end,value=e) )
        infile.close()
        return idx

    def _count( self, filename, idx ):

        overlapping_genes = set()
        genes = set()
        # iterate over exons
        infile = open( filename, "r" )
        it = GTF.iterator( infile )

        nexons, nexons_overlapping = 0, 0
        nbases, nbases_overlapping = 0, 0
        for this in it:
            nexons += 1
            nbases += this.end - this.start 
            genes.add( this.gene_id )
            
            try:
                intervals = idx[this.contig].find( this.start, this.end )
            except KeyError:
                continue


            if len(intervals) == 0: 
                continue

            overlapping_genes.add( this.gene_id )
            nexons_overlapping += 1
            start, end = this.start, this.end 
            counts = numpy.zeros( end - start, numpy.int )
            for other in intervals:
                for x in range( max(start, other.start ) - start, min(end, other.end ) - start ): 
                    counts[x] += 1
            nbases_overlapping += sum( [1 for x in counts if x > 0 ] )

        infile.close()

        return len(genes), len( overlapping_genes), nexons, nexons_overlapping, nbases, nbases_overlapping

    def count( self, filename1, filename2 ):
        """count overlap between two gtf files."""

        E.info( "counting started for %s versus %s" % (filename1, filename2))

        idx2 = self.buildIndex( filename2 )

        (self.mGenes1, self.mGenesOverlapping1,
         self.mExons1, self.mExonsOverlapping1,
         self.mBases1, self.mBasesOverlapping1 ) = \
         self._count( filename1, idx2 )

        self.mGenesUnique1 = self.mGenes1 - self.mGenesOverlapping1
        self.mExonsUnique1 = self.mExons1 - self.mExonsOverlapping1
        self.mBasesUnique1 = self.mBases1 - self.mBasesOverlapping1

        idx1 = self.buildIndex( filename1 )

        (self.mGenes2, self.mGenesOverlapping2,
         self.mExons2, self.mExonsOverlapping2,
         self.mBases2, self.mBasesOverlapping2 ) = \
         self._count( filename2, idx1 ) 

        self.mGenesUnique2 = self.mGenes2 - self.mGenesOverlapping2
        self.mExonsUnique2 = self.mExons2 - self.mExonsOverlapping2
        self.mBasesUnique2 = self.mBases2 - self.mBasesOverlapping2

    def __str__(self):

        return "\t".join( map( str, (
                        self.mGenes1, self.mGenes2,
                        self.mGenesOverlapping1, self.mGenesOverlapping2,
                        self.mGenesUnique1, self.mGenesUnique2,
                        self.mExons1, self.mExons2,
                        self.mExonsOverlapping1, self.mExonsOverlapping2,
                        self.mExonsUnique1, self.mExonsUnique2,
                        self.mBases1, self.mBases2,
                        self.mBasesOverlapping1, self.mBasesOverlapping2,
                        self.mBasesUnique1, self.mBasesUnique2 ) ) ) + "\t" +\
              "\t".join( map( lambda x: IOTools.prettyPercent( *x), (
                    (self.mGenesOverlapping1, self.mGenes1),
                    (self.mGenesOverlapping2, self.mGenes2),
                    (self.mGenesUnique1, self.mGenes1), 
                    (self.mGenesUnique2, self.mGenes2),
                    (self.mExonsOverlapping1, self.mExons1), 
                    (self.mExonsOverlapping2, self.mExons2),
                    (self.mExonsUnique1, self.mExons1), 
                    (self.mExonsUnique2, self.mExons2),
                    (self.mBasesOverlapping1, self.mBases1), 
                    (self.mBasesOverlapping2, self.mBases2),
                    (self.mBasesUnique1, self.mBases1), 
                    (self.mBasesUnique2, self.mBases2)  ) ) )

class CounterGenes(Counter):
    """output only genes."""
    
    mSeparator = ";"

    def __init__(self, *args, **kwargs ):
        Counter.__init__(self, *args, **kwargs)

    def getHeader( self ):
        h = [ "genes", "gene2" ]
        return "\t".join(h)

    def _run( self, filename, idx ):

        # iterate over exons
        infile = open( filename, "r" )
        it = GTF.iterator( infile )

        keys = set()

        for this in it:
            
            try:
                intervals = idx[this.contig].find( this.start, this.end )
            except KeyError:
                continue

            if len(intervals) == 0: 
                continue

            for i in intervals:
                key = "%s-%s" % (this.gene_id, i.value.gene_id )
                if key not in keys:
                    self.write( this.gene_id, i.value.gene_id )
                    keys.add( key )

        infile.close()


    def run( self, filename1, filename2 ):
        """count overlap between two gtf files."""

        E.info( "counting started for %s versus %s" % (filename1, filename2))

        idx2 = self.buildIndex( filename2 )
        self._run( filename1, idx2 )

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: gtfs2graph.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-s", "--ignore-strand", dest="ignore_strand", action="store_true",
                      help="ignore strand information [default=%default]." )

    parser.add_option("-u", "--update", dest="filename_update", type="string",
                      help="if filename is given, previous results will be read from there and only changed sets will be computed [default=%default]." )

    parser.add_option("-p", "--pattern-id", dest="pattern_id", type="string",
                      help="pattern to convert a filename to an id [default=%default]." )

    parser.add_option("-g", "--genes", dest="genes", action="store_true",
                      help="only output gene stats (includes gene lists) [default=%default]." )

    parser.set_defaults(
        ignore_strand = False,
        filename_update = None,
        pattern_id = "(.*).gtf",
        genes = False,
        )

    (options, args) = E.Start( parser )

    if len(args) != 2:
        print USAGE
        raise ValueError( "two arguments are required" )

    if options.genes:
        counter = CounterGenes( options.stdout )
    else:
        counter = Counter( options.stdout )

    options.stdout.write( counter.getHeader() + "\n" )

    counter.run( args[0], args[1] )
    
    E.Stop()    

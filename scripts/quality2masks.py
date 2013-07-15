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
quality2masks.py - mask low-quality bases in multiple alignment
===============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python quality2masks.py --help

Type::

   python quality2masks.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import optparse
import math
import time
import random
import types
import CGAT.Experiment as E
import CGAT.Blat as Blat
import CGAT.Iterators as Iterators
import alignlib
import CGAT.IndexedFasta as IndexedFasta

def fillAlignment( map_alignment, alignment ):

    i = 0
    for x,c in enumerate( alignment ):
        if c != "-":
            map_alignment.addPair( i, x )
            i += 1


def main():

    parser = E.OptionParser( version = "%prog version: $Id: quality2masks.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("--quality-threshold", dest="quality_threshold", type="int",
                      help="quality threshold for masking positions [default=%default]" )

    parser.add_option("--random", dest="random", action="store_true",
                      help="shuffle quality scores before masking [default=%default]" )

    parser.add_option("--filename-map", dest="filename_map", type="string",
                      help="filename in psl format mapping entries in multiple alignment to the genome [default=%default]" )

    parser.add_option("-q", "--quality-file", dest="quality_file", type="string",
                      help="filename with genomic base quality information [default=%default]."  )


    parser.set_defaults(
        quality_threshold = 40,
        quality_file = "quality",
        filename_map = None,
        frame = 3,
        )

    (options, args) = E.Start( parser )

    ##################################################
    ##################################################
    ##################################################
    ## read map
    ##################################################
    infile = open(options.filename_map) 
    map_genes2genome = {}
    for match in Blat.iterator( infile ):
        assert match.mQueryId not in map_genes2genome, "duplicate entry %s" % match.mQueryId
        map_genes2genome[match.mQueryId] = match
    infile.close()

    ##################################################
    ##################################################
    ##################################################
    ## get quality scores
    ##################################################
    quality = IndexedFasta.IndexedFasta( options.quality_file )
    quality.setTranslator( IndexedFasta.TranslatorBytes() )

    ##################################################
    ##################################################
    ##################################################
    ## main loop
    ##################################################
    ninput, noutput, nmissed = 0, 0, 0

    options.stdout.write( "cluster_id\tstart\tend\n" )

    for line in options.stdin:
        if line.startswith("cluster_id"): continue
        ninput += 1
        cluster_id, gene_id, alignment = line[:-1].split("\t")

        if gene_id not in map_genes2genome:
            nmissed += 1
            E.warn( "gene_id %s not found in map." % gene_id )
            continue
        
        match = map_genes2genome[gene_id]
        map_gene2genome = match.getMapQuery2Target()
        is_negative = match.strand == "-"

        # if strand is negative, the coordinates are 
        # on the negative strand of the gene/query
        # in order to work in the right coordinate system
        # revert the sequence
        if is_negative: 
            alignment = alignment[::-1]

        # get map of gene to alignment
        map_gene2mali = alignlib.makeAlignmentVector()
        fillAlignment( map_gene2mali, alignment )

        # get quality scores
        try:
            quality_scores = quality.getSequence( match.mSbjctId, "+", match.mSbjctFrom, match.mSbjctTo)
        except ValueError, msg:
            nmissed += 1
            E.warn( "could not retrieve quality scores for %s:%i-%i: %s" % (match.mSbjctId, match.mSbjctFrom, match.mSbjctTo, msg) )
            continue

        # print str(alignlib.AlignmentFormatEmissions( map_gene2genome))
        # print str(alignlib.AlignmentFormatEmissions( map_gene2mali))
        # print quality_scores

        map_mali2genome = alignlib.makeAlignmentVector()
        alignlib.combineAlignment( map_mali2genome, map_gene2mali, map_gene2genome, alignlib.RR )
        # print str(alignlib.AlignmentFormatEmissions( map_mali2genome))

        # shuffle quality scores, but only those that are aligned
        if options.random:
            positions = []
            for fp,c in enumerate(alignment):
                if c == "-": continue
                y = map_mali2genome.mapRowToCol( fp ) - match.mSbjctFrom 
                if y < 0: continue
                positions.append( y )
            scores = [ quality_scores[ x ] for x in positions ]
            random.shuffle(scores)
            for p,q in zip( positions,scores): quality_scores[p] = q

        # negative strand
        to_mask = []
        ## reverse position
        rp = len(alignment)
        for fp,c in enumerate(alignment):
            rp -= 1
            if c == "-": continue
            y = map_mali2genome.mapRowToCol( fp ) - match.mSbjctFrom
            if y < 0: continue
            if quality_scores[y] < options.quality_threshold:
                if is_negative: p = rp
                else: p = fp
                E.debug( "low quality base: id=%s, mali=%i, char=%s, contig=%s, strand=%s, pos=%i, quality=%i" % \
                             (cluster_id, p, c, match.mSbjctId, match.strand, map_mali2genome.mapRowToCol( fp ), quality_scores[y] ) )
                if options.frame > 1:
                    start = (p // options.frame) * options.frame
                    to_mask.extend( list( range(start, start + options.frame) ) )
                else:
                    to_mask.append( p ) 

        regions = Iterators.group_by_distance( sorted(to_mask) )
            
        for start,end in regions:
            options.stdout.write( "%s\t%i\t%i\n" % (cluster_id, start, end ) )

        noutput += 1

    E.info( "ninput=%i, noutput=%i, nmissed=%i" % (ninput, noutput, nmissed) )

    E.Stop()
    
if __name__ == "__main__":
    sys.exit(main())

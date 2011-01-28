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
Tophat.py - working with tophat/cufflinks output files
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''

import os, sys, re, collections, itertools

import IOTools

class CuffCompareValues:
    def __init__( self, vals ):
        assert len(vals) == 4
        try: self.sn = float( vals[0] )
        except ValueError: self.sn = None

        try: self.sp = float( vals[1] )
        except ValueError: self.sp = None

        try: self.fsn = float( vals[2] )
        except ValueError: self.fsn = None

        try: self.fsp = float( vals[3] )
        except ValueError: self.fsp = None
    
    def __str__(self):
        return "\t".join( map(IOTools.val2str,
                              ( self.sn,
                                self.sp,
                                self.fsn,
                                self.fsp ) ) )
    @classmethod
    def getHeaders( cls ):
        return ("sn", "sp", "fsn", "fsp" )

class CuffCompareResult:

    def __init__( self ):
        ( self.baselevel,
          self.exonlevel,
          self.intronlevel,
          self.intronchainlevel,
          self.transcriptlevel,
          self.locuslevel,
          self.missedexons_counts,
          self.missedexons_total,
          self.wrongexons_counts,
          self.wrongexons_total,
          self.missedintrons_counts,
          self.missedintrons_total,
          self.wrongintrons_counts,
          self.wrongintrons_total,
          self.missedloci_counts,
          self.missedloci_total,
          self.wrongloci_counts,
          self.wrongloci_total,
          self.query,
          self.query_loci,
          self.query_multi_exon,
          self.reference,
          self.reference_loci,
          self.reference_multi_exon,
          self.loci_multi_exon,
          self.loci_transcripts,
          ) = [None] * 26

    def fromLines( self, lines ):
        '''parse from cuffcompare output.'''
        self.is_empty = True

        for line in lines:
            if line.startswith("#"):
                d = line[1:-1].strip()
                if d.startswith("Query"):
                    self.query, self.query_loci, self.query_multi_exon = map(int, re.match( \
                        "Query mRNAs :\s+(\d+)\s+in\s+(\d+)\s+loci\s+\((\d+)", d ).groups() )
                elif d.startswith("Reference"):
                    self.reference, self.reference_loci, self.reference_multi_exon = map( int, re.match( \
                        "Reference mRNAs :\s+(\d+)\s+in\s+(\d+)\s+loci\s+\((\d+)", d ).groups() )
                elif d.startswith( "("):
                    self.loci_multi_exon, self.loci_transcripts = re.match( \
                        "\((\d+) multi-transcript loci, ~(\S+)", d ).groups()
                    self.loci_multi_exon = int( self.loci_multi_exon )
                    self.loci_transcripts = float( self.loci_transcripts )
                continue

            line = line[:-1].strip()
            if not line: continue

            try:
                tag, data = line.split(":")
            except ValueError:
                raise ValueError("parsing error in line %s" % line )

            tag = re.sub( "\s", "", tag ).lower()
            
            if tag.startswith( "wrong" ) or tag.startswith("missed"):
                counts, total = map(int, re.match( "\s+(\d+)/(\d+)", data).groups() )
                setattr( self, "%s_counts" % tag, counts )
                setattr( self, "%s_total" % tag, total )
                self.is_empty = False
            elif tag.startswith( "total"): 
                # stop at total union across all super-loci
                break
            else:
                values = data.strip().split()
                setattr( self, tag, CuffCompareValues( values ) )
                self.is_empty = False
            
    def __str__(self):
        return "\t".join( map(str, (
                    self.baselevel,
                    self.exonlevel,
                    self.intronlevel,
                    self.intronchainlevel,
                    self.transcriptlevel,
                    self.locuslevel,
                    self.missedexons_counts,
                    self.missedexons_total,
                    self.wrongexons_counts,
                    self.wrongexons_total,
                    self.missedintrons_counts,
                    self.missedintrons_total,
                    self.wrongintrons_counts,
                    self.wrongintrons_total,
                    self.missedloci_counts,
                    self.missedloci_total,
                    self.wrongloci_counts,
                    self.wrongloci_total,
                    self.query,
                    self.query_loci,
                    self.query_multi_exon,
                    self.reference,
                    self.reference_loci,
                    self.reference_multi_exon,
                    self.loci_multi_exon,
                    self.loci_transcripts
                    ) ) )

    @classmethod
    def getHeaders( cls ):
        m = ( "baselevel", "exonlevel", "intronlevel", "intronchainlevel",
              "transcriptlevel", "locuslevel" )
        a = ( "missed", "wrong")
        b = ( "exons", "introns", "loci" )
        c = ( "counts", "total" )
        return [ "%s_%s" % (x,y) for x,y in itertools.product( m, CuffCompareValues.getHeaders() ) ] +\
            [ "%s%s_%s" % (x,y,z) for x,y,z in itertools.product( a, b, c) ] +\
            [ "query", "query_loci", "query_multi_exon",
              "reference", "reference_loci", "reference_multi_exon",
              "loci", "loci_multi_exon" ]
    
def parseTranscriptComparison( infile ):
    '''read cufflinks output in infile stream.

    returns a two-level dictionary mapping with levels track and contig.
    '''

    def __blocker( infile ):
        blocks = {}
        contig, block = None, []
        y = False
        for line in infile:
            if line.startswith("#> Genomic sequence"):
                if block: blocks[contig] = block
                if y:
                    yield dataset, blocks
                    blocks = {}
                    y = False

                contig = re.match("#> Genomic sequence: (\S+)", line).groups()[0]
                block = []
                continue
            elif line.startswith( "#= Summary for dataset:" ):
                dataset = re.match("#= Summary for dataset: (\S+)", line).groups()[0]
                contig = "all"
                block = []
                y = True
                continue
                
            block.append( line )

        blocks[contig] = block
        yield dataset, blocks
    
    result = collections.defaultdict( dict )

    for track, blocks in __blocker( infile ):
        for contig, block in blocks.iteritems():
            r = CuffCompareResult()
            r.fromLines( block )
            result[track][contig] = r

    return result

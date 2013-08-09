#! /bin/env python
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
blast2table.py - output tabular results from BLAST/PSIBLAST runs
================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

A simple BLAST parser. This parser requires BLAST+. BLAST+ should be run with the
following output options::

   -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq"

The option *iterative* describes how results from iterative 
searches (PSIBLAST) should be reported. If set, it will automatically
add an additional field iteration:

all
   Results of all iterations are output.

last
   Results of the last iteration only are output.

first
   For each alignment between query and sbjct, the first 
   alignment is output.

.. note::

   Coordinates are 0-based, open-closed.

The option *alignment_format* selects the alignment format:

emissions

   The alignment is output as a string of positive and negative 
   integers for both query and sbjct sequence. Positive integers
   are emissions, negative integers denotes gaps.

blocks
   Block based alignment format. The alignment is reported as
   three ',' separated lists of block sizes and block positions
   in the query and sbjct sequence. Coordinates are relative to
   the alignment start.

Usage
-----

Type::

   python blast2table.py --help

for command line help.

Blast parsing
+++++++++++++

The following command outputs blast results using the pairsdb alignment format::

   blastp -query <(head -n 10 nrdb.fasta ) -db nrdb -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq" |
   python python blast2table.py --output-format=emissions

The following command outputs blast results using the pairsdb alignment format::

   blastp -query <(head -n 10 nrdb.fasta ) -db nrdb -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq" |
   python python blast2table.py --output-format=blocks

Psiblast parsing
++++++++++++++++

The following command outputs alignments from all iterations::

   psiblast -query <(head -n 2 nrdb.fasta ) -db nrdb -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq" -num_iterations=5
   python blast2table.py --alignment-format=blocks --iterations=all

The following command outputs alignments from the last iteration only::

   psiblast -query <(head -n 2 nrdb.fasta ) -db nrdb -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq" -num_iterations=5
   python blast2table.py --alignment-format=blocks --iterations=all

The following command outputs alignments for the first iteration that they are found::

   psiblast -query <(head -n 2 nrdb.fasta ) -db nrdb -outfmt "6 qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq" -num_iterations=5
   python blast2table.py --alignment-format=blocks --iterations=first

'''

import os
import sys
import re
import optparse
import collections
import CGAT.Experiment as E

BlastResult = collections.namedtuple( "blastresult",
                                      "qseqid qstart qend sseqid sstart send evalue bitscore pident score qseq sseq" )

class Output( object ):

    def __init__(self, line ):

        data = line[:-1].split("\t")
        try:
            self.r = BlastResult._make( data )
        except TypeError:
            raise ValueError( "parsing error in line: '%s'" % line[:-1])
        
        # define some shortcuts for fields used in processing later
        self.query = self.r.qseqid
        self.sbjct = self.r.sseqid
        self.bitscore = float(self.r.bitscore)
        self.iteration = None

class OutputEmissions( Output ):
    '''output blast results in the legacy pairsdb emissions format.'''

    header = ["query_nid", "sbjct_nid", "evalue",
              "query_start", "query_end", "query_ali",
              "sbjct_start", "sbjct_end", "sbjct_ali",
              "bitscore", "pid" ]
    
    def __str__( self ):
        
        def seq2ali( seq ):
            ali = seq.split( "-" )
            result = []
            x = 0
            while x < len(ali)-1:
                result.append( "+%i" % len(ali[x] ) )
                x += 1
                g = x
                while ali[x] == "": x += 1
                result.append( "-%i" % (x-g+1) ) 
            result.append( "+%i" % len(ali[-1] ) )
            return "".join(result)

        return "\t".join( ( self.r.qseqid,
                            self.r.sseqid,
                            self.r.evalue,
                            str(int(self.r.qstart) -1),
                            self.r.qend,
                            seq2ali( self.r.qseq ),
                            str(int(self.r.sstart) -1),
                            self.r.send,
                            seq2ali( self.r.sseq ),
                            self.r.bitscore,
                            self.r.pident ) )

class OutputBlocks( Output ):
    '''output blast alignments as blocks.'''

    header = ["query_nid", "sbjct_nid", "evalue",
              "query_start", "query_end", 
              "sbjct_start", "sbjct_end", 
              "block_sizes", "query_starts", "query_ends",
              "bitscore", "pid" ]
    
    def __str__( self ):

        block_sizes, query_starts, sbjct_starts = [], [], []

        in_block = True
        block_size = 0
        qc, sc = 0, 0
        query_starts.append( qc )
        sbjct_starts.append( sc )
        
        for qa, sa in zip( self.r.qseq, self.r.sseq ):
            qisgap = qa == "-"
            sisgap = sa == "-"
            if qisgap or sisgap:
                if in_block: block_sizes.append( block_size )
                in_block = False
            else:
                if not in_block:
                    # start of new block
                    query_starts.append( qc )
                    sbjct_starts.append( sc )
                    block_size = 0
                block_size += 1
                in_block = True

            if not qisgap: qc += 1
            if not sisgap: sc += 1

        block_sizes.append( block_size)

        return "\t".join( ( self.r.qseqid,
                            self.r.sseqid,
                            self.r.evalue,
                            str(int(self.r.qstart) -1),
                            self.r.qend,
                            str(int(self.r.sstart) -1),
                            self.r.send,
                            ",".join( map(str, block_sizes) ),
                            ",".join( map(str, query_starts) ),
                            ",".join( map(str, sbjct_starts) ),
                            self.r.bitscore,
                            self.r.pident ) )

def addIterations( group ):
    '''for each pair of query and sbjct add the
    iteration.

    There are no markers in psiblast output
    for the end of an iteration. Hence,
    this method uses as a proxy the bitscore.

    returns a dictionary and the number of iterations
    '''

    query2sbjct_pairs = collections.defaultdict( list )
    last_score = group[0].bitscore
    iteration = 1
    for r in group:
        if r.bitscore > last_score:
            iteration += 1
        r.iteration = iteration
        query2sbjct_pairs[ (r.query, r.sbjct) ].append( iteration )
        last_score = r.bitscore

    # normalize
    normed = {}
    for pair, iterations in query2sbjct_pairs.iteritems():
        normed[pair] = sorted(list( set( iterations )))
    return normed, iteration

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-f", "--alignment-format", dest="alignment_format", type="choice",
                      choices = ("emissions", "blocks" ),
                      help="output format options [default=%default]."  )

    parser.add_option( "--no-header", dest="with_header", action="store_false",
                       help="skip output of header [default=%default]."  )

    parser.add_option("-i", "--iterations", dest="iterations", type="choice",
                      choices = ("last", "all", "first", None ),
                      help="output for iterative searches [default=%default]."  )


    parser.set_defaults(
        alignment_format = "emissions",
        with_header = True,
        iterations = None,
        )

    (options, args) = parser.parse_args() 

    if options.alignment_format == "emissions":
        outer = OutputEmissions

    elif options.alignment_format == "blocks":
        outer = OutputBlocks

    if options.with_header:
        header = outer.header[:]
        if options.iterations != None: header.append( "iteration" )
        sys.stdout.write( "\t".join( header) + "\n" )

    def grouper( infile ):
        
        group = []
        
        query = None
        for line in infile:
            if line.startswith("#"): continue
            if not line.strip(): continue
            if line.startswith("Search has CONVERGED!"):
                sys.stdout.write("# %s" % line)
                continue

            r = outer( line )
            if r.query != query:
                if query: yield group
                query, group = r.query, []

            group.append( r )

        if query: yield group

    for group in grouper( sys.stdin ):
        if options.iterations == None:
            for r in group:
                sys.stdout.write( str(r) + "\n" )
        else:
            query2sbjct_pairs, max_iterations = addIterations( group )
            if options.iterations == "all":
                filtered = group
            elif options.iterations == "last":
                filtered = ( r for r in group if r.iteration == max_iterations )
            elif options.iterations == "first":
                output_pairs = {}
                for r in group:
                    key = (r.query,r.sbjct)
                    if key not in output_pairs: output_pairs[key] = r.iteration
                filtered = (r for r in group if r.iteration <= output_pairs[(r.query,r.sbjct)] )
                        
            for r in filtered:
                sys.stdout.write( "%s\t%i\n" % (str(r), r.iteration) )

if __name__ == "__main__":
    sys.exit( main( sys.argv) )
    

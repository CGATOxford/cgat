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
gff2transcripts.py - 
======================================================

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

   python gff2transcripts.py --help

Type::

   python gff2transcripts.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys, string, re, getopt

USAGE="""python %s [OPTIONS] < psl > predictions

Convert GFF exon list to predictions.

Version: $Id: gff2transcripts.py 18 2005-08-09 15:32:24Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-P, --table-predictions=        table name with predictions
-C, --table-contigs=            table name with contig sizes
-m, --file-map=                 file name with mappings of queries to transcripts
""" % sys.argv[0]



param_long_options=["verbose=", "help", "table-predictions=", "table-contigs=", "file-map="]
param_short_options="v:hP:C:m:"

param_tablename_predictions = None
param_tablename_contigs = None
param_filename_map = None

param_connection = "db:andreas"

import pgdb

import Experiment
import PredictionParser

def ProcessChunk( chunk, prediction_id, lcontig  ):
    """process a chunk.
    """
    
    for x in chunk:
        (id, code1, code2, first_res, last_res, score,
         strand, dummy, exon_id, gene_id, transcript_id,
         cdna) = x
        
        if lcontig > 0:
            t = lcontig - first_res
            first_res = lcontig - last_res
            last_res = t
        if first_res < 0 or last_res < 0 or first_res >= last_res:
            print prediction_id
            print lcontig
            print chunk
            raise "error in coordinates"
        
        print string.join( map(str, (
            prediction_id,
            first_res, last_res, cdna )) , "\t")
        
def GetMapTranscript2Prediction( dbhandle, tablename_predictions ):
    """retrieve map of transcripts to predictions."""

    map_transcript2prediction = {}
    
    statement = "SELECT query_token, prediction_id, sbjct_token, sbjct_strand FROM %s" % tablename_predictions
    cc = dbhandle.cursor()

    try:
        cc.execute(statement)
        result = cc.fetchall()
    except pgdb.DatabaseError, msg:
        print "# query failed with message", msg
        return map_transcript2prediction
    cc.close()

    for query_token, prediction_id, sbjct_token, sbjct_strand in result:
        map_transcript2prediction[query_token] = (prediction_id, sbjct_token, sbjct_strand)

    return map_transcript2prediction
    
def GetGenomeLengths( dbhandle, tablename_contigs ):
    """get sizes of contigs."""
    genome_lengths = {}
    
    statement = "SELECT sbjct_token, size FROM %s" % tablename_contigs
    cc = dbhandle.cursor()
    try:
        cc.execute(statement)
        result = cc.fetchall()
    except pgdb.DatabaseError, msg:
        print "# query failed with message", msg
        return genome_lengths

    for x, y in result:
        genome_lengths[x] = y

    cc.close()
    return genome_lengths

if __name__ == "__main__":

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)
        elif o in ("-P", "--table-predictions"):
            param_tablename_predictions = a
        elif o in ("-P", "--table-contigs"):
            param_tablename_contigs = a
        elif o in ("-m", "--file-map"):
            param_filename_map = a

    print Experiment.GetHeader()
    print Experiment.GetParams()

    dbhandle = pgdb.connect( param_connection )

    contig_sizes = GetGenomeLengths( dbhandle, param_tablename_contigs )

    map_transcript2prediction = GetMapTranscript2Prediction( dbhandle, param_tablename_predictions )
    nmissed = 0
    
    map_transcript2query = {}
    if param_filename_map:
        infile = open(param_filename_map, "r")
        for line in infile:
            if line[0] == "#": continue
            transcript, query = string.split(line[:-1], "\t")
            map_transcript2query[transcript] = query
    
    last_transcript_id = None
    for line in sys.stdin:
        if line[0] == "#": continue

        (id, code1, code2, first_res, last_res, score,
         strand, dummy, exon_id, gene_id, transcript_id,
         cdna) = re.split("\s", line[:-1])

        ## coordinates start from 0, but are not in-between
        first_res, last_res = map(int, (first_res, last_res) )
        last_res += 1

        if map_transcript2query:
            transcript_id = map_transcript2query[transcript_id]
        
        if last_transcript_id != transcript_id:

            if last_transcript_id:
                if map_transcript2prediction.has_key(last_transcript_id):
                
                    prediction_id, sbjct_token, sbjct_strand = map_transcript2prediction[last_transcript_id]
                    if sbjct_strand == "+":
                        lgenome = 0
                    else:
                        lgenome = contig_sizes[sbjct_token]
                        
                    ProcessChunk( chunk, prediction_id, lgenome )
                else:
                    print "# missed prediction: %s" % last_transcript_id
                    nmissed += 1
                    
            last_transcript_id = transcript_id
            chunk = []
            
        chunk.append( (id, code1, code2, first_res, last_res, score,
                       strand, dummy, exon_id, gene_id, transcript_id,
                       cdna) )

    if map_transcript2prediction.has_key(last_transcript_id):        
        prediction_id, sbjct_token, sbjct_strand = map_transcript2prediction[last_transcript_id]
        if sbjct_strand == "+":
            lgenome = 0
        else:
            lgenome = contig_sizes[sbjct_token]

        ProcessChunk( chunk, prediction_id, lgenome )
    else:
        print "# missed prediction: %s" % last_transcript_id
        nmissed += 1

    print "# missed=%i" % (nmissed)
    

    

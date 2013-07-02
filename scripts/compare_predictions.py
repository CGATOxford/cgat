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
compare_predictions.py - 
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

   python compare_predictions.py --help

Type::

   python compare_predictions.py --help

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
import getopt
import math

USAGE="""python %s [OPTIONS] < exonerate_output > filtered

Version: $Id: compare_predictions.py 1799 2008-03-28 11:44:19Z andreas $

Evaluate genewise alignments.

Build a file with exon comparisions between genes.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-b, --boundaries=               file with exon boundaries
-e, --exons=                    file with exons (output)
-p, --peptides=                 file with input peptide sequences.
-g, --genome-file=           pattern for filenames with the genomic DNA (FASTA).
-w, --write-notfound           print exons for predictions not found in reference.
-q, --quality-pide=             quality threshold (pide) for exons.
""" % sys.argv[0]

param_long_options=["verbose=", "help", "table-reference=", "table-target="]
param_short_options="v:hR:T:"

import alignlib
import CGAT.Experiment as Experiment
import CGAT.Genomics as Genomics
import CGAT.PredictionParser as PredictionParser
import CGAT.Exons as Exons
import pgdb

param_loglevel = 0

param_connection = "fgu202:andreas"

param_gop = -10.0
param_gep = -2.0

param_tablename_predictions_reference = None
param_tablename_predictions_target = None

##------------------------------------------------------------
if __name__ == '__main__':

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
        elif o in ("-R", "--table-reference"):
            param_tablename_predictions_reference = a
        elif o in ("-R", "--table-target"):
            param_tablename_predictions_target = a

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    print Experiment.GetHeader()
    print Experiment.GetParams()

    dbhandle = pgdb.connect( param_connection )

    statement = """
    SELECT
	prediction_id, 
	query_token, 
	sbjct_token, 
	sbjct_strand, 
	rank, 
	score, 
	query_from, 
	query_to, 
	query_ali, 
	sbjct_from, 
	sbjct_to, 
	sbjct_ali, 
	query_length, 
	query_coverage, 
	ngaps, 
	nframeshifts, 
	nintrons, 
	nsplits, 
	nstopcodons, 
	pidentity, 
	psimilarity, 
	sequence, 
	sbjct_genome_from, 
	sbjct_genome_to, 
	map_query2genome
    FROM %s AS p ORDER BY
    sbjct_token, sbjct_strand, sbjct_genome_from
    """

    cr = dbhandle.cursor()
    cr.execute( statement % param_tablename_predictions_reference )

    statement = """
    SELECT
	prediction_id, 
	query_token, 
	sbjct_token, 
	sbjct_strand, 
	rank, 
	score, 
	query_from, 
	query_to, 
	query_ali, 
	sbjct_from, 
	sbjct_to, 
	sbjct_ali, 
	query_length, 
	query_coverage, 
	ngaps, 
	nframeshifts, 
	nintrons, 
	nsplits, 
	nstopcodons, 
	pidentity, 
	psimilarity, 
	sequence, 
	sbjct_genome_from, 
	sbjct_genome_to, 
	map_query2genome
    FROM %s AS p 
    WHERE p.sbjct_token = '%s' AND
    p.sbjct_strand = '%s' AND 
    OVERLAP( %i, %i, p.sbjct_genome_from, sbjct_genome_to) > 0 
    """

    alignator = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_LOCAL, param_gop, param_gep )
    map_reference2target = alignlib.makeAlignmentVector()    
    assignment_id = 0
    
    for line in cr.fetchall():

        reference = PredictionParser.PredictionParserEntry()
        reference.FillFromTable( line )

        ct = dbhandle.cursor()
        ct.execute( statement % (param_tablename_predictions_target,
                    reference.mSbjctToken, reference.mSbjctStrand,
                    reference.mSbjctGenomeFrom, reference.mSbjctGenomeTo ))

        reference_exons = Exons.Alignment2Exons( reference.mMapPeptide2Genome,
                                                 0, 
                                                 reference.mSbjctFrom)

        for line2 in ct.fetchall():
            target = PredictionParser.PredictionParserEntry()
            target.FillFromTable( line2 )

            target_exons = Exons.Alignment2Exons( target.mMapPeptide2Genome,
                                                  0, 
                                                  target.mSbjctFrom)

            ## check for exon overlap
            rr, tt = 0, 0
            overlap = 0
            while rr < len(reference_exons) and tt < len(target_exons):
                
                r = reference_exons[rr]
                t = target_exons[tt]
                if r.mGenomeTo < t.mGenomeFrom:
                    rr += 1
                    continue
                elif t.mGenomeTo < r.mGenomeFrom:
                    tt += 1
                    continue
                overlap += ( min(r.mGenomeTo, t.mGenomeTo) - max(r.mGenomeFrom, t.mGenomeFrom))
                rr += 1
                tt += 1
                
            if overlap == 0:
                continue
            
            map_reference2target.clear()
            row = alignlib.makeSequence(reference.mTranslation)
            col = alignlib.makeSequence(target.mTranslation)
            alignator.align( map_reference2target, row, col )

            f = alignlib.AlignmentFormatEmissions( map_reference2target )
            row_ali, col_ali = f.mRowAlignment, f.mColAlignment
            pidentity = 100.0 * alignlib.calculatePercentIdentity( map_reference2target, row, col )
            psimilarity = 100.0 * alignlib.calculatePercentSimilarity( map_reference2target )        

            union = max( reference.mSbjctGenomeTo, target.mSbjctGenomeTo) - \
                    min( reference.mSbjctGenomeFrom, target.mSbjctGenomeFrom )
            inter = min( reference.mSbjctGenomeTo, target.mSbjctGenomeTo) - \
                    max( reference.mSbjctGenomeFrom, target.mSbjctGenomeFrom )
            
            assignment_id += 1
            
            print string.join( map(str, (
                assignment_id,
                reference.mPredictionId,
                target.mPredictionId,
                0, 0,
                overlap,
                "%5.2f" % (100.0 * float(overlap) / float(min(len(reference.mTranslation), len(target.mTranslation)) * 3)),
                "%5.2f" % (100.0 * float(overlap) / float(max(len(reference.mTranslation), len(target.mTranslation)) * 3)),
                "%5.2f" % (100.0 * float(inter) / float(union)),
                "%5.2f" % (100.0 * float(inter) / float(reference.mSbjctGenomeTo - reference.mSbjctGenomeFrom)),
                "%5.2f" % (100.0 * float(inter) / float(target.mSbjctGenomeTo - target.mSbjctGenomeFrom)),
                reference.mNIntrons - target.mNIntrons,
                reference.mSbjctGenomeFrom - target.mSbjctGenomeFrom,
                reference.mSbjctGenomeTo - target.mSbjctGenomeTo,                
                map_reference2target.getScore(),
                map_reference2target.getRowFrom(),
                map_reference2target.getRowTo(),
                row_ali,
                map_reference2target.getColFrom(),
                map_reference2target.getColTo(),
                col_ali,
                "%5.2f" % pidentity,
                "%5.2f" % psimilarity,
                map_reference2target.getNumGaps(),
                row.getLength(),
                col.getLength()) ), "\t" )

        ct.close()
        
    cr.close()

    print Experiment.GetFooter()
        
        

        

    

    

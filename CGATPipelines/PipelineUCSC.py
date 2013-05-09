################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
PipelineUCSC.py - utility functions for accessing UCSC data
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This module provides methods for accessing ucsc data.

Usage
-----

This module reads the ``[ucsc]`` section of a configuration file requiring the following variables to be
set:

host
    host to connect to

user
    username to connect with

database
    database to use




Documentation
-------------

Code
----

'''

# for UCSC import
import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections

import CGAT.Experiment as E
import CGAT.Pipeline as P
from ruffus import *

import MySQLdb

try:
    PARAMS = P.getParameters()
except IOError:
    pass

############################################################
############################################################
############################################################
## UCSC tracks
############################################################
def connectToUCSC():
    dbhandle = MySQLdb.Connect( host = PARAMS["ucsc_host"],
                                user = PARAMS["ucsc_user"] )

    cc = dbhandle.cursor()
    cc.execute( "USE %s " %  PARAMS["ucsc_database"] )

    return dbhandle

def getRepeatsFromUCSC( dbhandle, repclasses, outfile ):
    '''select repeats from UCSC and write to *outfile* in gff format.

    If *repclasses* is None, all repeats will be collected.
    '''

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is 
    # queried for tables that end in rmsk.
    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES LIKE '%rmsk'")
    tables = [ x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError( "could not find any `rmsk` tables" )

    # now collect repeats
    tmpfile = P.getTempFile(".")
    
    for table in tables:

        cc = dbhandle.cursor()
        sql = """SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd, '.', strand, '.', 
                      CONCAT('class \\"', repClass, '\\"; family \\"', repFamily, '\\"; repName \\"', repName, '\\";' )
               FROM %(table)s"""
        
        if repclasses:
            sql += " WHERE repClass in ('%(repclasses)s') "
               
        sql = sql % locals()

        E.debug( "executing sql statement: %s" % sql )
        cc.execute( sql )
        for data in cc.fetchall():
            tmpfile.write( "\t".join(map(str,data)) + "\n" )

    tmpfile.close()

    to_cluster = True

    # sort gff and make sure that names are correct
    tmpfilename = tmpfile.name

    statement = [ '''cat %(tmpfilename)s
        | %(scriptsdir)s/gff_sort pos 
        | python %(scriptsdir)s/gff2gff.py 
            --sanitize=genome 
            --skip-missing 
            --genome-file=%(genome_dir)s/%(genome)s
            --log=%(outfile)s.log ''']

    if PARAMS["geneset_remove_contigs"]:
        statement.append( ''' --remove-contigs="%(geneset_remove_contigs)s" ''' )

    statement.append( '''| gzip > %(outfile)s ''' )

    statement = " ".join( statement )
    
    P.run()

    os.unlink( tmpfilename)

def importRefSeqFromUCSC( infile, outfile, remove_duplicates = True ):
    '''import gene set from UCSC database
    based on refseq mappings.

    Outputs a gtf-formatted file a la ENSEMBL.

    Depending on *remove_duplicates*, duplicate mappings are either
    removed or kept.

    Matches to chr_random are ignored (as does ENSEMBL).

    Note that this approach does not work as a gene set, as refseq maps 
    are not real gene builds and unalignable parts cause
    differences that are not reconcilable.
    '''

    import MySQLdb
    dbhandle = MySQLdb.Connect( host = PARAMS["ucsc_host"],
                                user = PARAMS["ucsc_user"] )
        
    cc = dbhandle.cursor()
    cc.execute( "USE %s " %  PARAMS["ucsc_database"] )
        
    duplicates = set()

    if remove_duplicates:
        cc.execute( """SELECT name, COUNT(*) AS c FROM refGene 
                        WHERE chrom NOT LIKE '%_random'
                        GROUP BY name HAVING c > 1""" )
        duplicates = set( [x[0] for x in cc.fetchall() ] )
        E.info( "removing %i duplicates" % len(duplicates ) )

    # these are forward strand coordinates
    statement = '''
        SELECT gene.name, link.geneName, link.name, gene.name2, product, 
               protAcc, chrom, strand, cdsStart, cdsEnd, 
               exonCount, exonStarts, exonEnds, exonFrames 
        FROM refGene as gene, refLink as link 
        WHERE gene.name = link.mrnaAcc 
              AND chrom NOT LIKE '%_random'
        ORDER by chrom, cdsStart 
        '''
    
    outf = gzip.open(outfile, "w")
    
    cc = dbhandle.cursor()
    cc.execute(statement)

    SQLResult = collections.namedtuple('Result', 
        '''transcript_id, gene_id, gene_name, gene_id2, description,
        protein_id, contig, strand, start, end, 
        nexons, starts, ends, frames''')

    counts = E.Counter()
    counts.duplicates = len(duplicates)

    for r in map( SQLResult._make, cc.fetchall() ):

        if r.transcript_id in duplicates: continue

        starts = map( int, r.starts.split(",")[:-1])
        ends = map( int, r.ends.split(",")[:-1])
        frames = map( int, r.frames.split(",")[:-1])

        gtf = GTF.Entry()
        gtf.contig = r.contig
        gtf.source = "protein_coding"
        gtf.strand = r.strand
        gtf.gene_id = r.gene_id
        gtf.transcript_id = r.transcript_id
        gtf.addAttribute( "protein_id", r.protein_id )
        gtf.addAttribute( "transcript_name", r.transcript_id )
        gtf.addAttribute( "gene_name", r.gene_name )

        assert len(starts) == len(ends) == len(frames)

        if gtf.strand == "-":
            starts.reverse()
            ends.reverse()
            frames.reverse()

        counts.transcripts += 1
        i = 0
        for start, end, frame in zip( starts, ends, frames ):
            gtf.feature = "exon"
            counts.exons += 1
            i += 1
            gtf.addAttribute("exon_number", i)
            # frame of utr exons is set to -1 in UCSC
            gtf.start, gtf.end, gtf.frame = start, end, "."
            outf.write( "%s\n" % str(gtf))
            
            cds_start, cds_end = max( r.start, start), min( r.end, end)
            if cds_start >= cds_end: 
                # UTR exons have no CDS
                # do not expect any in UCSC
                continue
            gtf.feature = "CDS"
            # invert the frame
            frame = (3 - frame % 3 ) % 3
            gtf.start, gtf.end, gtf.frame = cds_start, cds_end, frame
            outf.write( "%s\n" % str(gtf))

    outf.close()
    
    E.info("%s" % str(counts))

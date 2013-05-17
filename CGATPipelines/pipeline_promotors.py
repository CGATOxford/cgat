################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
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
==========================================================
pipeline_promotors.py - Defining promotor characteristics
==========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Pipeline to annotate promotors according to certain properties:

1. TATA box
2. G+C Island

Overview
========

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineReporting`). To start with, use the files supplied with the
Example_ data.

Input
-----

Optional inputs
+++++++++++++++

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|                    |                   |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_template.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_template.tgz
   tar -xvzf pipeline_template.tgz
   cd pipeline_template
   python <srcdir>/pipeline_template.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import glob
import gzip
import os
import itertools
import re
import math
import types
import collections
import time
import optparse
import shutil
import sqlite3
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGAT.FastaIterator as FastaIterator
import CGAT.Bed as Bed

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import CGAT.Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineGeneset as PipelineGeneset

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    Use this method to connect to additional databases.

    Returns a database connection.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

###################################################################
###################################################################
###################################################################
## worker tasks
###################################################################

@merge( (os.path.join( PARAMS['annotations_dir'],
                       PARAMS_ANNOTATIONS['interface_tss_bed']),
         os.path.join( PARAMS['annotations_dir'],
                       PARAMS_ANNOTATIONS['interface_contigs'])),
        'tata.bed.gz' )
def findTATABox( infiles, outfile ):
    '''find TATA box in promotors. There are several matrices to choose from:

    M00216 V$TATA_C Retroviral TATA box
    M00252 V$TATA_01 cellular and viral TATA box elements
    M00311 V$ATATA_B Avian C-type TATA box
    M00320 V$MTATA_B Muscle TATA box
    '''
    
    # 1. create fasta file - look for TATA box
    # 
    bedfile, genomefile = infiles

    statement = '''
    slopBed -i %(bedfile)s
            -l %(tata_search_upstream)i
            -r %(tata_search_downstream)i
            -s
            -g %(genomefile)s
    | python %(scriptsdir)s/bed2fasta.py 
       --use-strand
       --genome=%(genome_dir)s/%(genome)s
       --log=%(outfile)s.log
    > %(outfile)s.fasta
    '''
       
    P.run()
    
    match_executable = '/ifs/data/biobase/transfac/match/bin/match_linux64'
    match_matrix = '/ifs/data/biobase/transfac/dat/matrix.dat'
    match_profile = 'minFP_good.prf'
    match_profile = outfile + ".prf"

    prf = '''tata.prf
prf to minimize sum of both errors - derived from minSUM.prf
 MIN_LENGTH 300
0.0
 1.000 0.716 0.780 M00216 V$TATA_C
 1.000 0.738 0.856 M00252 V$TATA_01
 1.000 0.717 0.934 M00311 V$ATATA_B
 1.000 0.711 0.784 M00320 V$MTATA_B
//
'''

    with IOTools.openFile(match_profile, "w") as outf:
        outf.write( prf )

    # -u : uniq - only one best match per sequence
    statement = '''
         %(match_executable)s
         %(match_matrix)s
         %(outfile)s.fasta
         %(outfile)s.match
         %(match_profile)s
         -u
    >> %(outfile)s.log
    '''
    P.run()

    transcript2pos = {}
    for entry in FastaIterator.iterate( IOTools.openFile( outfile + ".fasta" )):
        transcript_id, contig, start, end, strand = re.match( "(\S+)\s+(\S+):(\d+)..(\d+)\s+\((\S)\)", entry.title ).groups()
        transcript2pos[transcript_id] = (contig, int(start), int(end), strand )

    MATCH = collections.namedtuple( "MATCH", "pid transfac_id pos strand core_similarity matrix_similarity sequence" )
    def _grouper( infile ):
        r = []
        keep = False
        for line in infile:
            if line.startswith("Inspecting sequence ID"): 
                keep = True
                if r: yield pid, r
                r = []
                pid = re.match( "Inspecting sequence ID\s+(\S+)", line ).groups()[0]
                continue
            elif line.startswith(" Total"): break

            if not keep: continue
            if line[:-1].strip() == "": continue
            transfac_id, v, core_similarity, matrix_similarity, sequence = [ x.strip() for x in line[:-1].split("|")]
            pos, strand = re.match( "(\d+) \((\S)\)", v ).groups()
            r.append( MATCH._make( (pid, transfac_id, int(pos), strand,
                                    float(core_similarity), float(matrix_similarity), sequence ) ) )

                                    
        yield pid, r
        
    offset = PARAMS["tata_search_upstream"]

    outf = IOTools.openFile( outfile + ".table.gz", "w" )
    outf.write("\t".join( ("transcript_id", "strand",
                           "start", "end",
                           "relative_start", "relative_end",
                           "transfac_id",
                           "core_similarity",
                           "matrix_similarity",
                           "sequence") ) + "\n" )

    bedf = IOTools.openFile( outfile, "w" )

    c = E.Counter()
    found = set()
    for transcript_id, matches in _grouper( IOTools.openFile(outfile + ".match") ):
        contig, seq_start, seq_end, strand = transcript2pos[transcript_id]
        c.promotor_with_matches += 1
        nmatches = 0
        found.add( transcript_id) 
        for match in matches:

            c.matches_total += 1
            lmatch = len( match.sequence )
            if match.strand == "-": 
                c.matches_wrong_strand += 1
                continue

            # get genomic location of match
            if strand == "+":
                genome_start = seq_start + match.pos 
            else:
                genome_start = seq_end - match.pos - lmatch
            
            genome_end = genome_start + lmatch

            # get relative location of match
            if strand == "+":
                tss_start = seq_start + offset
                relative_start = genome_start - tss_start 
            else:
                tss_start = seq_end - offset
                relative_start = tss_start - genome_end
                
            relative_end = relative_start + lmatch
            
            outf.write( "\t".join( map(str, (
                            transcript_id, strand,
                            genome_start, genome_end,
                            relative_start, relative_end,
                            match.transfac_id,
                            match.core_similarity,
                            match.matrix_similarity,
                            match.sequence ) ) ) + "\n" )
            c.matches_output += 1
            nmatches += 1

            bedf.write( "\t".join( map(str, (contig, genome_start, genome_end, transcript_id, strand,
                                             match.matrix_similarity) )) + "\n" )

        if nmatches == 0:
            c.promotor_filtered += 1
        else:
            c.promotor_output += 1

    c.promotor_total = len(transcript2pos)
    c.promotor_without_matches = len(set( transcript2pos.keys() ).difference( found ))

    outf.close()
    bedf.close()

    with IOTools.openFile( outfile + ".summary", "w" ) as outf:
        outf.write ("category\tcounts\n" )
        outf.write( c.asTable() + "\n" )
    
    E.info( c )

@transform( findTATABox, suffix(".bed.gz"), ".load" )
def loadTATABox( infile, outfile ):
    '''load TATA box information.'''
    
    P.load( infile + ".table.gz", outfile, "--index=transcript_id" )

###################################################################
###################################################################
###################################################################
## 
###################################################################

############################################################
############################################################
############################################################
## 
############################################################
@merge( None, "cpg.bed.gz" )
def collectCpGIslands( infile, outfile ):
    '''select repeats from UCSC and write to *outfile* in gff format.
    '''

    dbhandle = PipelineGeneset.connectToUCSC()

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is 
    # queried for tables that end in rmsk.
    cc = dbhandle.cursor()
    table = "cpgIslandExt"
    sql = """SELECT chrom, chromStart, chromEnd, name, obsExp
               FROM %(table)s
    """ % locals() 
    E.debug( "executing sql statement: %s" % sql )
    cc.execute( sql )
    outf = IOTools.openFile( outfile, "w" )
    for data in cc.fetchall():
        outf.write( "\t".join( map(str, data) ) + "\n" )

    outf.close()

############################################################
############################################################
############################################################
## 
############################################################
@merge( (collectCpGIslands,
         os.path.join( PARAMS['annotations_dir'],
                       PARAMS_ANNOTATIONS['interface_tss_bed']) ),
        "cpg.tsv.gz" )
def annotateCpGIslands( infiles, outfile ):
    '''annotate transcript by absence/presence of CpG islands
    '''
    cpgfile, tssfile = infiles
    cpg = Bed.readAndIndex( IOTools.openFile( cpgfile ) )
    
    extension_upstream = PARAMS["cpg_search_upstream"]
    extension_downstream = PARAMS["cpg_search_downstream"]

    c = E.Counter()
    outf = IOTools.openFile( outfile, "w" )
    outf.write("transcript_id\tstrand\tstart\tend\trelative_start\trelative_end\n" )

    for tss in Bed.iterator(IOTools.openFile( tssfile ) ):
        c.tss_total += 1

        if tss.strand == "+":
            start, end = tss.start - extension_upstream, tss.start + extension_downstream
        else:
            start, end = tss.end - extension_downstream, tss.end + extension_upstream

        try:
            matches = list(cpg[tss.contig].find( start, end ))
        except KeyError:
            c.promotor_without_matches += 1
            continue

        if len(matches) == 0:
            c.promotor_without_matches += 1
            continue

        c.promotor_output += 1
        for match in matches:
            c.matches_total += 1
            genome_start, genome_end, x = match

            l = genome_end - genome_start

            # get relative location of match
            if tss.strand == "+":
                relative_start = genome_start - tss.start 
            else:
                relative_start = tss.end - genome_end
            
            relative_end = relative_start + l

            outf.write( "\t".join( map(str, (
                            tss.name, tss.strand,
                            genome_start, genome_end,
                            relative_start, relative_end ))) + "\n" )
            c.matches_output += 1

    outf.close()
            
    with IOTools.openFile( outfile + ".summary", "w" ) as outf:
        outf.write ("category\tcounts\n" )
        outf.write( c.asTable() + "\n" )
    
    E.info( c )

@transform( annotateCpGIslands, suffix(".tsv.gz"), ".load" )
def loadCpGIslands( infile, outfile ):
    '''load CpG Islands information.'''
    
    P.load( infile, outfile, "--index=transcript_id" )

@merge( (loadTATABox, loadCpGIslands), "promotorinfo_transcripts.load" )
def loadTranscriptSummary( infile, outfile ):
    '''summarize binding information per transcript.'''

    dbh = connect()

    table = P.toTable( outfile )

    cc = dbh.cursor()
    # sqlite can not do full outer join
    cc.execute( """DROP TABLE IF EXISTS %(table)s""" % locals() )

    transcripts = [ x[0] for x in cc.execute( "SELECT DISTINCT(transcript_id) FROM annotations.transcript_info" ).fetchall() ]

    tmpf = P.getTempFile()

    tables = ("tata", "cpg" )
    titles = tables

    vals = []
    for table in tables:
        t = set([ x[0] for x in cc.execute( "SELECT DISTINCT(transcript_id) FROM %(table)s" % locals()).fetchall() ])
        vals.append( t )
        
    tmpf.write("transcript_id\t%s\n" % "\t".join( titles ) )
    
    for transcript_id in transcripts:
        tmpf.write( "%s\t%s\n" % (transcript_id,
                                  "\t".join( [ str(int(transcript_id in v)) for v in vals ] ) ) )
        
    tmpf.close()

    P.load( tmpf.name, outfile )
    os.unlink( tmpf.name )

@merge( loadTranscriptSummary, "promotorinfo_genes.load" )
def loadGeneSummary( infile, outfile ):
    '''summarize binding information per gene.'''

    dbh = connect()

    table = P.toTable( outfile )

    cc = dbh.cursor()
    cc.execute( """DROP TABLE IF EXISTS %(table)s """ % locals() )
    cc.execute( """CREATE TABLE %(table)s AS
            SELECT gene_id, SUM( tata ) AS tata, SUM( cpg ) AS cpg 
                   FROM promotorinfo_transcripts AS p,
                        annotations.transcript_info as i
                   WHERE i.transcript_id = p.transcript_id
                   GROUP BY gene_id""" % locals())
    cc.close()

    P.touch( outfile )

@merge( loadGeneSummary, PARAMS["interface_promotor_ontology"] )
def buildGeneOntology( infile, outfile ):
    '''create an output file akin to GO ontology files to be
    used with GO.py
    '''
    
    table = P.toTable( infile )
    columns = ("cpg", "tata" )
    dbh = connect()
    cc = dbh.cursor()

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "go_type\tgene_id\tgo_id\tdescription\tevidence\n" )

    i = 1
    for c in columns:
        cc.execute( "SELECT DISTINCT gene_id FROM %(table)s WHERE %(c)s" % locals() )
        outf.write( "".join( ["promotor\t%s\tGO:%07i\twith_%s\tNA\n" % (x[0],i,c) for x in cc ] ) )
        i += 1
        cc.execute( "SELECT DISTINCT gene_id FROM %(table)s WHERE %(c)s = 0" % locals() )
        outf.write( "".join( ["promotor\t%s\tGO:%07i\twithout_%s\tNA\n" % (x[0],i,c) for x in cc ] ) )
        i += 1

    outf.close()

@follows( loadTranscriptSummary, loadGeneSummary, buildGeneOntology )
def full(): pass

###################################################################
###################################################################
###################################################################
## primary targets
###################################################################
@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting report build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating report" )
    P.run_report( clean = False )

@follows( update_report )
def publish_report():
    '''publish report.'''

    E.info( "publishing report" )
    P.publish_report()

if __name__== "__main__":

    # P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )

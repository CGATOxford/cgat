################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_chipseq.py 2900 2010-04-13 14:38:00Z andreas $
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
:Release: $Id: pipeline_chipseq.py 2900 2010-04-13 14:38:00Z andreas $
:Date: |today|
:Tags: Python

===================
Annotation pipeline
===================

The annotation pipeline imports various annotations that can be used
in other pipelines.

   * a geneset
   * UCSC tracks

The pipeline works on a single genome.

Configuration
=============

Input
=====


Output
======

The results of the computation are all stored in an sqlite relational
database file.

Usage
=====

Type::

   python <script_name>.py --help

for command line help.

Code
----

TODO: currently the bed files and the intervals are inconsistent 
    (due to filtering, there are more intervals in the bed files than
     in the table. The ids do correspond).

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections

import Experiment as E
import Pipeline as P
from ruffus import *

import sqlite3

# for UCSC import
import MySQLdb

import IndexedFasta, IOTools
import PipelineGeneset as PGeneset

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS

if os.path.exists("conf.py"):
    E.info( "reading additional configuration from conf.py" )
    execfile("conf.py")

###################################################################
###################################################################
###################################################################
## genome section
############################################################
############################################################
############################################################
@files( PARAMS["genome"] + ".fasta", PARAMS['interface_contigs'] )
def buildContigSizes( infile, outfile ):
    '''output contig sizes
    '''
    fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )
    outs = IOTools.openFile(outfile, "w" )

    for contig, size in fasta.getContigSizes( with_synonyms = False ).iteritems():
        outs.write( "%s\t%i\n" % ( contig,size) )

    outs.close()

###################################################################
###################################################################
###################################################################
## gene set section
############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], PARAMS['interface_annotation_gff'] )
def buildGeneRegions( infile, outfile ):
    '''annotate genomic regions with reference gene set.

    Only considers protein coding genes. In case of overlapping
    genes, only take the longest (in genomic coordinates).
    Genes not on UCSC contigs are removed.
    '''
    PGeneset.buildGeneRegions( infile, outfile )

############################################################
############################################################
############################################################
@follows( buildGeneRegions )
@files( PARAMS["ensembl_filename_gtf"], PARAMS['interface_genes_gtf'] )
def buildGenes( infile, outfile ):
    '''build a collection of exons from the protein-coding
    section of the ENSEMBL gene set. The exons include both CDS
    and UTR.

    The set is filtered in the same way as in :meth:`buildGeneRegions`.
    '''
    PGeneset.buildProteinCodingGenes( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "gene_info.load" )
def loadGeneInformation( infile, outfile ):
    '''load the transcript set.'''
    PGeneset.loadGeneInformation( infile, outfile )

############################################################
############################################################
############################################################
@files( buildGenes, "gene_stats.load" )
def loadGeneStats( infile, outfile ):
    '''load the transcript set.'''

    PGeneset.loadGeneStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_cds_gtf"] )
def buildCDSTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set.

    Only CDS exons are parts of exons are output - UTR's are removed.
    '''
    PGeneset.buildCDS( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], 
        PARAMS["interface_exons_gtf"] )
def buildExonTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set.

    Only CDS exons are parts of exons are output - UTR's are removed.
    
    '''
    PGeneset.buildExons( infile, outfile )

############################################################
############################################################
############################################################
@transform( buildCDSTranscripts, suffix(".gtf.gz"), "_gtf.load" )
def loadCDSTranscripts( infile, outfile ):
    '''load the transcript set.'''
    PGeneset.loadTranscripts( infile, outfile )

############################################################
############################################################
############################################################
@files( buildCDSTranscripts, "cds_stats.load" )
def loadCDSStats( infile, outfile ):
    '''load the transcript set.'''

    PGeneset.loadTranscriptStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "transcript_info.load" )
def loadTranscriptInformation( infile, outfile ):
    '''load the transcript set.'''
    PGeneset.loadTranscriptInformation( infile, 
                                        outfile,
                                        only_proteincoding = PARAMS["ensembl_only_proteincoding"] )

###################################################################
###################################################################
###################################################################
@files( ( (PARAMS["ensembl_filename_pep"], 
           PARAMS["interface_peptides_fasta"] ), ) )
def buildPeptideFasta( infile, outfile ):
    '''load ENSEMBL peptide file
    
    *infile* is an ENSEMBL .pep.all.fa.gz file.
    '''
    PGeneset.buildPeptideFasta( infile, outfile )

###################################################################
###################################################################
###################################################################
@files( ( (PARAMS["ensembl_filename_cdna"], 
           PARAMS["interface_cdna_fasta"] ), ) )
def buildCDNAFasta( infile, outfile ):
    '''load ENSEMBL peptide file
    
    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''
    PGeneset.buildCDNAFasta( infile, outfile )

###################################################################
###################################################################
###################################################################
@merge( buildCDSTranscripts, PARAMS["interface_cds_fasta"] )
def buildCDSFasta( infile, outfile ):
    '''build cds sequences from peptide and cds file.
    
    *infile* is an ENSEMBL .cdna.all.fa.gz file.
    '''

    PGeneset.buildCDSFasta( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_pep"], "protein_stats.load" )
def loadProteinStats( infile, outfile ):
    '''load the transcript set.'''

    PGeneset.loadProteinStats( infile, outfile )

############################################################
############################################################
############################################################
@merge( (loadProteinStats, loadTranscriptInformation), "seleno.list")
def buildSelenoList( infile, outfile ):
    '''export a list of seleno cysteine transcripts.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    statement = '''
    SELECT DISTINCT transcript_id
    FROM transcript_info as t,
         protein_stats as p
    WHERE p.protein_id = t.protein_id AND
         p.nU > 0
    '''
    outf = open(outfile, "w" )
    outf.write("transcript_id\n")
    outf.write("\n".join( [x[0] for x in cc.execute( statement) ] ) + "\n" )
    outf.close()


############################################################
############################################################
############################################################
@merge( buildCDSTranscripts, PARAMS["interface_promotors_bed"] )
def buildPromotorRegions( infile, outfile ):
    '''annotate promotor regions from reference gene set.'''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing 
                                           --genome-file=%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=%(geneset_promotor_size)s
                              --genome-file=%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | gzip 
        > %(outfile)s
    """

    P.run()

############################################################
############################################################
############################################################
@merge( buildCDSTranscripts, PARAMS["interface_tss_bed"] )
def buildTSSRegions( infile, outfile ):
    '''annotate transcription start sites from reference gene set.

    Similar to promotors, except that the witdth is set to 1.
    '''
    statement = """
        gunzip < %(infile)s 
        | python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome)s --log=%(outfile)s.log 
        | python %(scriptsdir)s/gff2bed.py --is-gtf --name=transcript_id --log=%(outfile)s.log 
        | gzip
        > %(outfile)s
    """

    P.run()

############################################################
############################################################
############################################################
@follows( buildPromotorRegions, buildGeneRegions, buildTSSRegions )
@files( [ ("%s.gtf.gz" % x, "%s.bed.gz" % x ) for x in ("ensembl", "promotors", "tss" ) ] )
def exportRegionAsBed( infile, outfile ):
    '''export a reference gtf file as bed for computing overlap.'''

    bed = Bed.Bed()
    outfile = open( outfile, "w" )
    with open(infile, "r") as inf: 
        for gff in GTF.iterator( inf ):
            bed.contig, bed.start, bed.end = gff.contig, gff.start, gff.end
            bed.mFields = [ gff.gene_id ]
            outfile.write( "%s\n" % str(bed) )
    outfile.close()

############################################################
############################################################
############################################################
## UCSC tracks
############################################################

############################################################
############################################################
############################################################
@files( ((None, PARAMS["interface_repeats_gff"] ), ) )
def importRepeatsFromUCSC( infile, outfile ):
    '''import repeats from a UCSC formatted file.

    The repeats are stored as a :term:`gff` formatted file.
    '''

    repclasses="','".join(PARAMS["ucsc_repeattypes"] )

    # Repeats are either stored in a single ``rmsk`` table (hg19) or in
    # individual ``rmsk`` tables (mm9) like chr1_rmsk, chr2_rmsk, ....
    # In order to do a single statement, the ucsc mysql database is 
    # queried for tables that end in rmsk.
    dbhandle = MySQLdb.Connect( host = PARAMS["ucsc_host"],
                                user = PARAMS["ucsc_user"] )

    cc = dbhandle.cursor()
    cc.execute( "USE %s " %  PARAMS["ucsc_database"] )

    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES LIKE '%rmsk'")
    tables = [ x[0] for x in cc.fetchall()]
    if len(tables) == 0:
        raise ValueError( "could not find any `rmsk` tables" )

    tmpfile = P.getTempFile()
    
    for table in tables:
        E.info( "loading repeats from %s" % table )
        cc = dbhandle.cursor()
        cc.execute("""SELECT genoName, 'repeat', 'exon', genoStart+1, genoEnd, strand, '.', '.', 
                      CONCAT('class \\"', repClass, '\\"; family \\"', repFamily, '\\";')
               FROM %(table)s
               WHERE repClass in ('%(repclasses)s') """ % locals() )
        for data in cc.fetchall():
            tmpfile.write( "\t".join(map(str,data)) + "\n" )

    tmpfile.close()
    tmpfilename = tmpfile.name

    statement = '''cat %(tmpfilename)s
        | %(scriptsdir)s/gff_sort pos 
        | python %(scriptsdir)s/gff2gff.py 
            --sanitize=genome 
            --skip-missing 
            --genome-file=%(genome)s
            --log=%(outfile)s.log 
        | gzip
        > %(outfile)s
    '''
    P.run()
    
    os.unlink( tmpfilename )

if 0:
    ############################################################
    ############################################################
    ############################################################
    ## get UCSC tables
    ############################################################
    def getUCSCTracks(infile = PARAMS["filename_ucsc_encode"]):
        '''return a list of UCSC tracks from infile.'''
        tables = []
        with open(infile) as f:
            for line in f:
                if line.startswith("#"): continue
                tablename = line[:-1].strip()
                if tablename == "": continue
                tables.append( tablename )
        return tables

    ############################################################
    ############################################################
    ############################################################
    ## import UCSC encode tracks
    ############################################################
    @posttask( touch_file("ucsc_encode.import") )
    @files( PARAMS["filename_ucsc_encode"], "ucsc_encode.import") 
    def importUCSCEncodeTracks( infile, outfile ):

        statement = '''
        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -B -e "SELECT * FROM %(tablename)s" %(ucsc_database)s |\
        csv2db.py %(csv2db_options)s \
                  --table=%(tablename)s \
        >> %(outfile)s

        '''

        dbhandle = sqlite3.connect( PARAMS["database"] )

        cc = dbhandle.cursor()
        tables = set( [ x[0] for x in cc.execute( "SELECT name FROM sqlite_master WHERE type='table'") ] )
        cc.close()

        for tablename in getUCSCTracks( infile ):
            if tablename in tables:
                E.info( "skipping %(tablename)s - already exists" % locals())
                continue            

            E.info( "importing %(tablename)s" % locals() )
            P.run()

    ############################################################
    ############################################################
    ############################################################
    ## import UCSC encode tracks
    ############################################################
    @files( PARAMS["filename_mappability"], PARAMS["annotator_mappability"] ) 
    def exportUCSCMappabilityTrackToBed( infile, outfile ):
        '''convert wiggle track with mappability information to a bed-formatted file
        with only the mappable regions in the genome.'''

        infile = gzip.open( infile, "r" )

        outf = open( outfile, "w" )

        for line in infile:
            if line.startswith("fixedStep" ):
                contig, start, step = re.match( "fixedStep chrom=(\S+) start=(\d+) step=(\d+)", line ).groups()
                start, step = int(start)-1, int(step)
                end = start + step
                last_val = None
            else:
                val = int(line)
                if last_val != val:
                    if last_val == 1:
                        outf.write( "\t".join( (contig, str(start), str(end ) ) ) + "\n" )
                    start = end
                    end = start + step
                else:
                    end += step
                last_val = val
        outf.close()
    ############################################################
    ############################################################
    ############################################################
    ## export UCSC encode tracks as bed
    ############################################################
    @transform( importUCSCEncodeTracks, suffix(".import"), ".bed")
    def exportUCSCEncodeTracks( infile, outfile ):

        dbhandle = sqlite3.connect( PARAMS["database"] )

        outs = open(outfile, "w")
        for tablename in getUCSCTracks():
            outs.write( "track name=%s\n" % tablename )

            cc = dbhandle.cursor()
            statement = "SELECT chrom, chrostart, chroend FROM %s ORDER by chrom, chrostart" % (tablename)
            cc.execute( statement )
            for contig, start, end in cc:
                outs.write("%s\t%i\t%i\n" % (contig, start, end) )
        outs.close()


##################################################################
##################################################################
##################################################################
## Primary targets
##################################################################
@follows( buildContigSizes )
def genome():
    '''import information on geneset.'''
    pass

@follows( loadCDSTranscripts,
          loadTranscriptInformation,
          loadGeneStats,
          loadGeneInformation,
          buildExonTranscripts,
          buildSelenoList)
def geneset():
    '''import information on geneset.'''
    pass

@follows( buildPeptideFasta,
          buildCDSFasta,
          buildCDNAFasta )
def fasta():
    '''build fasta files.'''
    pass

@follows( buildPromotorRegions,
          buildTSSRegions )
def promotors():
    '''build promotors.'''
    pass

@follows( geneset, fasta, promotors, genome )
def full():
    '''build all targets.'''
    pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

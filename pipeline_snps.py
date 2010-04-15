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

:Author: Andreas Heger
:Release: $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

snp annotation pipeline.

Input:

Indels in pileup format.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

"""
from ruffus import *
import sys, glob, gzip, os, itertools, CSV, re
import optparse, shutil
import sqlite3
import GFF, GTF
import Experiment as E
import Pipeline as P
import IOTools
import Genomics
import Database
import FastaIterator
import PipelineGeneset as PGeneset
import PipelineGO as PGO

###################################################################
###################################################################
###################################################################
## read global options from configuration file
PARAMS = P.getParameters()

PARAMS["transcripts"] = "transcripts.gtf"
PARAMS["geneset"] = 'ensembl.gtf'
PARAMS["annotation"] = 'ensembl_regions.gff'

if not os.path.exists("conf.py"):
    raise IOError( "could not find configuration file conf.py" )

execfile("conf.py")
SEPARATOR = "|"

###################################################################
###################################################################
###################################################################
## import
###################################################################
if PARAMS["filename_snps"]:
    @split( PARAMS["filename_snps"], "*.pileup.gz" )
    def importSNPs( infile, outfile ):
        '''build samtools pileup formatted files from tabular output
        #CHROM  POS     REF     129P2   129S1   129S5   AKR     A_J     BALB    C3H     C57BL   CAST    CBA     DBA     LP_J    NOD     NZO     PWK     SPRET   WSB
        '''

        outfiles = IOTools.FilePool( "mouse%s.pileup.gz" )

        inf = gzip.open(infile,"r")
        headers = []
        for line in inf:
            data = line[:-1].split("\t")
            if line.startswith("#"):
                if not headers: headers = data[3:]
                continue
            contig, pos, ref = data[:3]
            for h, genotype in zip(headers, data[3:]):
                if genotype == "..": continue
                outfiles.write( h, 
                                "\t".join( map(str, (
                                contig,
                                pos,
                                ref,
                                Genomics.encodeGenotype( genotype ),
                                "0",
                                "0",
                                "0",
                                "0",
                                genotype,
                                "<" * len(genotype) ) ) ) + "\n" )

        outfiles.close()

elif PARAMS["filename_pileup"]:
    pass

###################################################################
###################################################################
###################################################################
@jobs_limit(2)
@transform( "*.pileup.gz", suffix(".pileup.gz"), "_pileup.import")
def importPileup( infile, outfile ):
    '''import pileup information.

    only imports chromosome, pos, genotype (first three columns)
    '''
    
    to_cluster = False
    tablename = outfile[:-len(".import")]
    statement = '''
    gunzip < %(infile)s |
    awk 'BEGIN {printf("contig\\tpos\\treference\\tgenotype\\n")} 
               {printf("%%s\\t%%s\\t%%s\\t%%s\\n",$1,$2,$3,$4)}' |
        csv2db.py %(csv2db_options)s \
        --index=contig,pos \
        --table=%(tablename)s
    > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################
## gene set section
###################################################################
############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], PARAMS['annotation'] )
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
@files( PARAMS["ensembl_filename_gtf"], PARAMS['geneset'] )
def buildGenes( infile, outfile ):
    '''build a collection of exons from the protein-coding
    section of the ENSEMBL gene set. The exons include both CDS
    and UTR.

    The set is filtered in the same way as in :meth:`buildGeneRegions`.
    '''
    PGeneset.buildGenes( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "gene_info.import" )
def importGeneInformation( infile, outfile ):
    '''import the transcript set.'''
    PGeneset.importGeneInformation( infile, outfile )

############################################################
############################################################
############################################################
@files( buildGenes, "gene_stats.import" )
def importGeneStats( infile, outfile ):
    '''import the transcript set.'''

    PGeneset.importGeneStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], PARAMS["transcripts"] )
def buildTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set.

    Only CDS are used.
    '''
    PGeneset.buildTranscripts( infile, outfile )

############################################################
############################################################
############################################################
@transform( buildTranscripts, suffix(".gtf"), "_gtf.import" )
def importTranscripts( infile, outfile ):
    '''import the transcript set.'''
    PGeneset.importTranscripts( infile, outfile )

############################################################
############################################################
############################################################
@files( buildTranscripts, "transcript_stats.import" )
def importTranscriptStats( infile, outfile ):
    '''import the transcript set.'''

    PGeneset.importTranscriptStats( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_gtf"], "transcript_info.import" )
def importTranscriptInformation( infile, outfile ):
    '''import the transcript set.'''
    PGeneset.importTranscriptInformation( infile, 
                                          outfile,
                                          only_proteincoding = PARAMS["ensembl_only_proteincoding"] )

############################################################
############################################################
############################################################
@files( PARAMS["ensembl_filename_pep"], "protein_stats.import" )
def importProteinStats( infile, outfile ):
    '''import the transcript set.'''

    PGeneset.importProteinStats( infile, outfile )

############################################################
############################################################
############################################################
@merge( (importProteinStats, importTranscriptInformation), "seleno.list")
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

###################################################################
###################################################################
###################################################################
@files( PARAMS["ensembl_filename_gtf"], "annotations_bases.fasta" )
def buildBaseAnnotations( infile, outfile ):
    """build base annotations"""       

    to_cluster = True
    job_queue = "server_jobs.q"

    dbname = outfile[:-len(".fasta")]
    statement = """
        gunzip < %(infile)s |\
        python %(scriptsdir)s/gtf2fasta.py \
                --force \
                --genome=%(genome)s \
                --output-filename-pattern=annotations_bases.%%s \
                --log=%(outfile)s.log |\
        python %(toolsdir)s/index_fasta.py \
        --log=%(outfile)s.log \
        %(dbname)s - > %(outfile)s.log
    """

    P.run( **dict( locals().items() + PARAMS.items() ) )

###################################################################
###################################################################
###################################################################
@files( PARAMS["ensembl_filename_gtf"], "annotations_exons.gtf" )
def buildExonAnnotations( infile, outfile ):
    """build exon annotations"""

    to_cluster = True
    job_queue = "server_jobs.q"
    
    statement = """
        gunzip < %(infile)s |\
        awk '$3 == "CDS"' |\
        python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s --log=%(outfile)s.log |\
        python %(scriptsdir)s/gtf2gff.py \
                --method=exons \
                --restrict-source=protein_coding \
                --log=%(outfile)s.log \
        > %(outfile)s
    """

    P.run( **dict( locals().items() + PARAMS.items() ) )

###################################################################
###################################################################
###################################################################
@files( PARAMS["ensembl_filename_gtf"], "annotations_genes.gtf", )
def buildGeneAnnotations( infile, outfile ):
    """build gene annotations.

    Merge exons per gene within the reference set. The
    output includes the UTR and non-coding genes.
    """
    statement = """
        gunzip < %(infile)s |\
        python %(scriptsdir)s/gtf2gtf.py --merge-exons --with-utr --log=%(outfile)s.log |\
        python %(scriptsdir)s/gtf2gtf.py --set-transcript-to-gene --log=%(outfile)s.log |\
        python %(scriptsdir)s/gff2gff.py --skip-missing --sanitize=genome --genome-file=genome --log=%(outfile)s.log |\
        %(scriptsdir)s/gff_sort gene-pos \
        > %(outfile)s
    """
    queue = "server"
    P.run( **dict( locals().items() + PARAMS.items() ) )


###################################################################
###################################################################
###################################################################
@files( buildGeneAnnotations, "annotations_genes.counts" )
def makeGeneCounts( infile, outfile ):
    """coun gene exon statistics.
    """
    
    statement = """
    cat < %(infile)s |\
    python %(scriptsdir)s/gtf2table.py \
        --genome-file=genome \
        --counter=length \
        --log=%(outfile)s.log \
    > %(outfile)s
    """ 
    P.run( **dict( locals().items() + PARAMS.items() ) )

###################################################################
###################################################################
###################################################################
@follows( buildBaseAnnotations, buildExonAnnotations)
@transform(  PARAMS["filename_pileup"],
             suffix(".pileup.gz"), 
             ".annotations.gz" )
def makeAnnotations( infile, outfile ):
    """annotate snps with gene set."""

    to_cluster = True
    
    exons = "annotations_exons.gtf"
    bases = "annotations_bases"

    statement = """
    gunzip < %(infile)s |\
    grep -v "^NT" |\
    python %(scriptsdir)s/snp2table.py \
        --genome-file=genome \
        --filename-exons=%(exons)s \
        --filename-annotations=%(bases)s \
        --filename-junctions=%(bases)s.junctions \
        --log=%(outfile)s.log |\
    gzip > %(outfile)s
    """ 
    P.run( **dict( locals().items() + PARAMS.items() ) )

###################################################################
###################################################################
###################################################################
@follows( buildSelenoList )
@transform(  '*.pileup.gz', 
             suffix(".pileup.gz"), 
             ".effects.gz" )
def makeEffects( infile, outfile ):
    """annotate snps with gene set."""

    to_cluster = True
    
    exons = "annotations_exons.gtf"
    bases = "annotations_bases"
    seleno = "seleno.list"

    statement = """
    gunzip < %(infile)s |\
    grep -v "^NT" |\
    python %(scriptsdir)s/snp2counts.py \
        --genome-file=genome \
        --module=transcript-effects \
        --filename-seleno=%(seleno)s \
        --filename-exons=%(transcripts)s \
        --output-filename-pattern=%(outfile)s.%%s \
        --log=%(outfile)s.log |\
    gzip > %(outfile)s
    """ 
    P.run( **dict( locals().items() + PARAMS.items() ) )

###################################################################
###################################################################
###################################################################
@transform(  makeEffects, 
             suffix(".effects.gz"), 
             "_effects.import" )
def importEffects( infile, outfile ):
    '''import transcript effects into tables.'''

    root = infile[:-len(".effects.gz")]

    statement = '''
    csv2db.py %(csv2db_options)s \
              --from-zipped \
              --index=transcript_id \
              --table=%(root)s_effects \
    < %(infile)s > %(outfile)s
    '''
    P.run()

    for suffix in ("cds", "intron", "splicing", "translation"):
        
        statement = '''
        csv2db.py %(csv2db_options)s \
        --index=transcript_id \
        --table=%(root)s_effects_%(suffix)s \
        < %(infile)s.%(suffix)s >> %(outfile)s
        '''
        P.run()


###################################################################
###################################################################
###################################################################
@transform(  '*_pileup.import', 
             suffix("_pileup.import"), 
             inputs( ("_pileup.import", buildTranscripts, buildSelenoList) ),
             ".alleles" )
def buildAlleles( infiles, outfile ):
    """annotate snps with gene set."""

    to_cluster = True

    infile, transcripts, seleno = infiles
    tablename = infile[:-len(".import")]

    statement = """cat < %(transcripts)s | 
    python %(scriptsdir)s/gtf2alleles.py 
        --genome-file=genome 
        --filename-seleno=%(seleno)s 
        --output-filename-pattern=%(outfile)s.%%s 
        --tablename=%(tablename)s
    > %(outfile)s
    """ 
    P.run( **dict( locals().items() + PARAMS.items() ) )

###################################################################
###################################################################
###################################################################
@transform(  buildAlleles, 
             suffix(".alleles"), 
             "_alleles.import" )
def importAlleles( infile, outfile ):
    '''import allele.'''

    tablename = outfile[:-len(".import")] 

    statement = '''cat
    < %(infile)s.table 
    | perl -p -e "s/False/0/g; s/True/1/g;"
    | csv2db.py %(csv2db_options)s \
              --index=gene_id \
              --index=transcript_id \
              --table=%(tablename)s \
    > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
###################################################################
@transform(importAlleles, 
           suffix("_alleles.import"),
           "_alleles_transcripts.import" )
def summarizeAllelesPerTranscript( infile, outfile ):
    '''summarize effects on a per-gene level.

    The following fields are exclusive:
    is_wildtype
       both alleles are wildtype
    is_knockout
       both alleles knocked out
    is_truncated
       both alleles truncated or truncated and knocked out
    is_affected
       one allele is truncated or knocked out

    The other fields are not necessarily exclusive, for example there
    could be transcripts with one knocked out allele and one wildtype
    allele, such that ``is_nmd_affected``, ``is_affected`` and ``has_wildtype`` 
    are all true.
    '''
    
    tablename = outfile[:-len(".import")]
    track = infile[:-len("_alleles.import")]

    dbhandle = sqlite3.connect( PARAMS["database"] )

    statement = '''
    CREATE TABLE %(tablename)s AS
    SELECT DISTINCT 
           transcript_id,
           COUNT(DISTINCT allele_id) AS nalleles,
           CASE WHEN SUM( is_nmd_knockout) = 2 THEN 1 ELSE 0 END AS is_nmd_knockout,
           CASE WHEN SUM( is_nmd_knockout) >= 1 THEN 1 ELSE 0 END AS is_nmd_affected,
           CASE WHEN SUM( is_splice_truncated) = 2 THEN 1 ELSE 0 END AS is_splice_truncated,
           CASE WHEN SUM( is_splice_truncated) >= 1 THEN 1 ELSE 0 END AS is_splice_affected,
           CASE WHEN SUM( is_stop_truncated) = 2 THEN 1 ELSE 0 END AS is_stop_truncated,
           CASE WHEN SUM( is_stop_truncated) >= 1 THEN 1 ELSE 0 END AS is_stop_affected,
           CASE WHEN SUM( is_wildtype ) = 2 THEN 1 ELSE 0 END AS is_wildtype, 
           CASE WHEN SUM( is_wildtype ) >= 1 THEN 1 ELSE 0 END AS has_wildtype, 
           contig AS contig, 
           strand AS strand, 
           GROUP_CONCAT( reference_first_stop_start ) AS stop_codons_start,
           GROUP_CONCAT( reference_first_stop_end ) AS stop_codons_end,
           0 AS is_knockout,
           0 AS is_truncated,
           0 AS is_affected
    FROM %(track)s_alleles AS a
    GROUP BY transcript_id
    ''' % locals()
    
    Database.executewait( dbhandle, "DROP TABLE IF EXISTS %(tablename)s" % locals() )
    Database.executewait( dbhandle, statement )
    Database.executewait( dbhandle, "CREATE INDEX %(tablename)s_transcript_id ON %(tablename)s (transcript_id)" % locals())
    Database.executewait( dbhandle, "UPDATE %(tablename)s SET is_knockout = is_nmd_knockout" % locals())
    Database.executewait( dbhandle, '''UPDATE %(tablename)s SET is_truncated = 
                                       is_splice_truncated OR is_stop_truncated OR 
                                       (is_splice_affected AND is_stop_affected) OR 
                                       (is_splice_affected AND is_nmd_affected) OR 
                                       (is_stop_affected AND is_nmd_affected)
                                       ''' % locals())
    Database.executewait( dbhandle, 'UPDATE %(tablename)s SET is_affected ='
                          '(is_nmd_affected OR is_splice_affected OR is_stop_affected) AND NOT'
                          '(is_knockout or is_truncated)'% locals())
    dbhandle.commit()

    P.touch(outfile)

###################################################################
###################################################################
###################################################################
@transform(summarizeAllelesPerTranscript, 
           suffix("_alleles_transcripts.import"),
           "_alleles_genes.import" )
def summarizeAllelesPerGene( infile, outfile ):
    '''summarize effects on a per-gene level.'''
    
    tablename = outfile[:-len(".import")]
    track = infile[:-len(".import")]

    dbhandle = sqlite3.connect( PARAMS["database"] )

    statement = '''
    CREATE TABLE %(tablename)s AS
    SELECT DISTINCT 
           i.gene_id AS gene_id,
           COUNT( DISTINCT a.transcript_id ) AS ntranscripts,
           CASE WHEN SUM( is_nmd_knockout ) = COUNT(DISTINCT a.transcript_id) THEN 1 ELSE 0 END AS is_nmd_knockout,
           SUM( is_nmd_knockout ) AS is_nmd_affected,
           CASE WHEN SUM( is_splice_truncated) = COUNT(DISTINCT a.transcript_id) THEN 1 ELSE 0 END AS is_splice_truncated,
           SUM( is_splice_truncated ) AS is_splice_affected,
           CASE WHEN SUM( is_stop_truncated ) = COUNT(DISTINCT a.transcript_id) THEN 1 ELSE 0 END AS is_stop_truncated,
           SUM( is_stop_truncated ) AS is_stop_affected,
           CASE WHEN SUM( is_wildtype ) = COUNT(DISTINCT a.transcript_id) THEN 1 ELSE 0 END AS is_wildtype, 
           SUM( is_wildtype ) AS has_wildtype, 
           contig AS contig, 
           strand AS strand, 
           GROUP_CONCAT( stop_codons_start ) AS stop_codons_start,
           GROUP_CONCAT( stop_codons_end ) AS stop_codons_end,
           0 AS is_knockout,
           0 AS is_truncated,
           0 AS is_affected
    FROM %(track)s AS a, transcript_info AS i
    WHERE i.transcript_id = a.transcript_id
    GROUP BY i.gene_id
    ''' % locals()
    
    Database.executewait( dbhandle, "DROP TABLE IF EXISTS %(tablename)s" % locals() )
    Database.executewait( dbhandle, statement )
    Database.executewait( dbhandle, "CREATE INDEX %(tablename)s_gene_id ON %(tablename)s (gene_id)" % locals())
    Database.executewait( dbhandle, "UPDATE %(tablename)s SET is_knockout = is_nmd_knockout" % locals())
    Database.executewait( dbhandle, '''UPDATE %(tablename)s SET is_truncated = 
                                       is_splice_truncated OR is_stop_truncated OR 
                                       (is_splice_affected + is_stop_affected >= ntranscripts)  
                          ''' % locals())
    Database.executewait( dbhandle, 'UPDATE %(tablename)s SET is_affected ='
                          '(is_nmd_affected OR is_splice_affected OR is_stop_affected) AND NOT'
                          '(is_knockout or is_truncated)'% locals())

    dbhandle.commit()

    P.touch(outfile)

###################################################################
###################################################################
###################################################################
@transform(importEffects, 
           suffix("_effects.import"),
           "_effects_genes.import" )
def summarizeEffectsPerGene( infile, outfile ):
    '''summarize effects on a per-gene level.'''
    
    tablename = outfile[:-len(".import")]
    track = infile[:-len("_effects.import")]

    dbhandle = sqlite3.connect( PARAMS["database"] )

    statement = '''
    CREATE TABLE %(tablename)s AS
    SELECT DISTINCT 
           gene_id, 
           COUNT(*) AS ntranscripts,
           MIN(e.nalleles) AS min_nalleles,
           MAX(e.nalleles) AS max_nalleles,
           MIN(e.stop_min) AS min_stop_min,
           MAX(e.stop_min) AS max_stop_min,
           MIN(e.stop_max) AS min_stop_max,
           MAX(e.stop_max) AS max_stop_max,
           SUM( CASE WHEN stop_min > 0 AND cds_len - stop_min * 3 < last_exon_start THEN 1  
                     ELSE 0 END) AS nmd_knockout,
           SUM( CASE WHEN stop_max > 0 AND cds_len - stop_max * 3 < last_exon_start THEN 1  
                     ELSE 0 END) AS nmd_affected
    FROM transcript_info as i,
         %(track)s_effects AS e
    WHERE i.transcript_id = e.transcript_id
    GROUP BY i.gene_id
    ''' % locals()
    
    Database.executewait( dbhandle, "DROP TABLE IF EXISTS %(tablename)s" % locals() )
    Database.executewait( dbhandle, statement )
    Database.executewait( dbhandle, "CREATE INDEX %(tablename)s_gene_id ON %(tablename)s (gene_id)" % locals())
    dbhandle.commit()

    P.touch(outfile)
    
###################################################################
###################################################################
###################################################################
@transform( makeAnnotations,
            suffix('.annotations.gz'), 
            '.annotations.import' )
def importAnnotations( infile, outfile ):
    '''import annotations'''

    root = infile[:-len(".annotations.gz")]

    statement = '''
    csv2db.py --map gene_id:str \
              -b sqlite \
              --from-zipped \
              --index=gene_id \
              --table=%(root)s_annotations \
              --map=base_qualities:text \
    < %(infile)s > %(outfile)s
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

###################################################################
@follows( buildGeneAnnotations )
@files_re(  glob.glob( '*.pileup.gz'),
            '(.*).pileup.gz', 
            [r'\1.pileup.gz', "annotations_genes.gtf" ],
            r'\1.genecounts.gz' )
def makeSNPCountsPerGene( infiles, outfile ):
    """count snps within genes"""
    
    infile_snps, infile_genes = infiles

    statement = """
    gunzip < %(infile_snps)s |\
    grep -v "^NT" |\
    python %(scriptsdir)s/snp2counts.py \
        --genome-file=genome \
        --filename-exons=%(ensembl_filename_gtf)s \
        --log=%(outfile)s.log |\
    gzip > %(outfile)s
    """ 
    P.run()


###################################################################
@files( [ (None, "assignments.go" ), ] )
def createGO( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    PGO.createGO( infile, outfile )

############################################################
@files_re( createGO, "(.*).go", r"\1.goslim") 
def createGOSlim( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    PGO.createGOSlim( infile, outfile )


###################################################################
@follows(summarizeEffectsPerGene, createGO, createGOSlim)
@files( [ ("%s_effects_genes.import" % x, "%s_%s.%s" % (x,y,z), (x,y,z)) for x,y,z in 
      itertools.product( TRACKS_GO, 
                         ("nmdknockouttranscript",
                          "nmdaffectedtranscript",
                          "nmdknockoutgenes",
                          "nmdaffectedgenes"),
                         ("go", "goslim") ) ] )
def runGOAnalyses( infile, outfile, options ):
    '''run GO analysis on transcripts that have been knocked out
    by premature stop codons.

    ``options`` is a tuple of (``track``, ``analysis``, ``ontology``)

    ``analysis`` can be:

    nmdknockouttranscript
        genes for which one transcript has been knocked out
        due to NMD
    nmdaffectedtranscript
        genes in which one transcript is affected by NMD
    nmdknockoutgenes
        genes in which all transcripts have been knocked out 
        due to NMD

    '''
        
    track, analysis, ontology = options

    # setup foreground set
    if analysis == "nmdknockouttranscript":
        field_where = "e.nmd_knockout > 0"
    elif analysis == "nmdaffectedtranscript":
        field_where = "e.nmd_affected > 0"
    elif analysis == "nmdknockoutgenes":
        field_where = "e.nmd_knockout = e.ntranscripts"

    statement_fg = '''
    SELECT DISTINCT e.gene_id 
        FROM
            %(track)s_effects_genes AS e
        WHERE 
              %(field_where)s
        ORDER BY e.gene_id
    ''' % locals()

    # setup background set
    statement_bg = '''SELECT DISTINCT gene_id FROM gene_info''' % locals()
        
    # choose ontology
    if ontology == "go":
        gofile = "assignments.go"
    elif ontogoly == "goslim":
        gofile = "assignments.goslim"

    # create result directory
    outdir = os.path.abspath( outfile + ".dir" )
    try: os.makedirs( outdir )
    except OSError: pass

    # run
    PGO.runGOFromDatabase( outfile, 
                           outdir, 
                           statement_fg,
                           statement_bg,
                           gofile )

############################################################################
@follows( mkdir( os.path.join( PARAMS["scratchdir"], "malis.dir" ) ) )
@merge( buildAlleles, "malis.map" )
def setupMultipleAlignment( infiles, outfile ):
    '''prepare input files for multiple alignment computations.

    This script does some id-mapping to resolve coordinates.

    Basically, each genome is separated into two alleles. 
    Gene_id's will be suffixed with the allele_id. This ensures
    that exons of a gene with multiple transcipts will be resolved
    correctly with consistent coordinates. 
 
    From an alignment point of view, the two alleles of the genes will be treated
    independently, but transcripts within a gene  will be merged correctly at exon 
    boundaries, again on a per-allele basis.

    Later, when collecting the results, the allele id is moved from the gene to
    the transcript.
    '''

    targetdir = os.path.join( PARAMS["scratchdir"], "malis.dir" )

    filepool_gtf = IOTools.FilePoolMemory( "%(targetdir)s/cluster_%%s.dir/cluster_%%s.gtf" % locals() )
    filepool_pep = IOTools.FilePoolMemory( "%(targetdir)s/cluster_%%s.dir/cluster_%%s_pep.fasta" % locals() )
    filepool_cds = IOTools.FilePoolMemory( "%(targetdir)s/cluster_%%s.dir/cluster_%%s_cds.fasta" % locals() )

    outf = open( outfile, "w")
    outf.write("id\tgroup_id\n")

    map_gene2group = {}
    map_seqid2code = {}
    x = 0
    counts = E.Counter()
    for infile in infiles:
        track = infile[:-len(".alleles")]
        E.info( "adding track %s" % track )

        reader = CSV.DictReader( open(infile+".table","rU"), dialect="excel-tab" )
        for row in reader:
            counts.input += 1
            gene_id, allele_id,transcript_id = row["gene_id"], row["allele_id"], row["transcript_id"]
            if gene_id not in map_gene2group:
                map_gene2group[gene_id] = len(map_gene2group)
            group_id = map_gene2group[gene_id]
            new_gene_id = "-".join( (gene_id, allele_id))
            if row["is_wildtype"] == "1": code = "WT"
            if row["is_nmd_knockout"] == "1": 
                counts.nmd_knockouts += 1
                continue

            else: code = "VA"
            seq_id = SEPARATOR.join( (track, transcript_id, new_gene_id ))
            map_seqid2code[seq_id] = code
            seq_id = SEPARATOR.join( (seq_id, code))
            outf.write( "%s\t%i\n" % (seq_id, group_id))
            filepool_pep.write( str(group_id), ">%s\n%s\n" % (seq_id, row["peptide"] ) )
            filepool_cds.write( str(group_id), ">%s\n%s\n" % (seq_id, row["cds"] ) )
            counts.written += 1

        with inf as open(infile+".gtf"):
            for gtf in GTF.iterator( inf ): 
                group_id = map_gene2group[gtf.gene_id]
                new_gene_id = "-".join( (gtf.gene_id, gtf["allele_id"]))
                seq_id = SEPARATOR.join( (track, gtf.transcript_id, new_gene_id ))
                seq_id = SEPARATOR.join( (seq_id, map_seqid2code[seq_id]))
                gtf.transcript_id = seq_id
                filepool_gtf.write( group_id, str(gtf) + "\n")

        x += 1
        # if x > 2: break

    E.info( "writing data" )
    filepool_gtf.close()
    filepool_pep.close()
    filepool_cds.close()
    outf.close()
    counts.ngroups = len(map_gene2group)
    counts.nsequences = len(map_seqid2code)

    E.info( "%s\n" % (str(counts)) )

@transform( os.path.join( PARAMS["scratchdir"], "malis.dir", "*", "*.gtf"), 
            suffix(".gtf"), 
            ".mali")
def buildMultipleAlignments( infile, outfile ):
    '''build multiple alignments.'''

    track = infile[:-len(".gtf")]
    filename_cds = track + "_cds.fasta"
    filename_pep = track + "_pep.fasta"

    to_cluster = True

    statement = '''
	python %(scriptsdir)s/align_transcripts.py \
		--gtf=%(infile)s \
		--cds=%(filename_cds)s \
		--force-map \
		--verbose=2 \
		--output-filename-pattern=%(track)s_%%s.fasta \
		--output=final_aa \
		--output=final_na \
		--output=aligned_aa \
		--output=aligned_na \
		--output-format="plain-fasta" \
	< %(filename_pep)s > %(outfile)s
      '''

    P.run()

@merge( buildMultipleAlignments, "variants" )
def buildMultipleAlignmentVariantColumns( infile, outfile ):
    '''build multiple alignments.'''

    track = infile[:-len(".gtf")]
    filename_cds = track + "_cds.fasta"
    filename_pep = track + "_pep.fasta"

    to_cluster = True

    statement = '''
	python %(scriptsdir)s/malis2mali.py \
		--gtf=%(infile)s \
		--cds=%(filename_cds)s \
		--force-map \
		--verbose=2 \
		--output-filename-pattern=%(track)s_%%s.fasta \
		--output=final_aa \
		--output=final_na \
		--output=aligned_aa \
		--output=aligned_na \
		--output-format="plain-fasta" \
	< %(filename_pep)s > %(outfile)s
      '''

    P.run()

@merge( buildMultipleAlignments, "malis.result" )
def mergeMultipleAlignments( infiles, outfile ):
    '''collect multiple alignment results into files that
    are compatible with OPTIC.
    '''

    for section in ("final_aa", "final_na", "aligned_aa", "aligned_na"):
        outfilename = outfile + "." + section + ".gz"

        counter = E.Counter()

        E.info("processing %s into %s" % (section, outfilename ))
        outf = gzip.open( outfilename, "w" )
        outf.write("cluster_id\tspecies\ttranscript_id\tgene_id\tcode\tsequence\n")
        for infile in infiles:
            counter.input += 1
            dirname, filename = os.path.split( infile )
            cluster_id = re.match("cluster_(\d+).mali", filename ).groups()[0]
            infilename = os.path.join( dirname, "cluster_%s_%s.fasta" % (cluster_id, section))
            # E.debug( "adding %s - %s from %s" % (filename, cluster_id, infilename) )
            if not os.path.exists(infilename):
                counter.missing += 1
                E.warn("multiple alignment %s missing" % infilename )
                continue
            for entry in FastaIterator.FastaIterator( open( infilename, "r")):
                parts = entry.title.split(SEPARATOR)
                if len(parts) == 4:
                    species, transcript_id, gene_id, code = entry.title.split(SEPARATOR)
                elif len(parts) == 2:
                    species, gene_id = entry.title.split(SEPARATOR)
                    transcipt_id = gene_id
                    code = "CG"
                # transfer the allele_id from the gene to the transcript
                gene_id, allele_id = gene_id.split("-")
                transcript_id += "-" + allele_id

                outf.write( "\t".join(map(str,
                                          (cluster_id,
                                           species,
                                           transcript_id,
                                           gene_id,
                                           code,
                                           entry.sequence))) + "\n")
            counter.output += 1

        outf.close()
        E.info( "%s: %s" % (outfilename, str(counter)))

    P.touch(outfile)

@merge( '*_pileup.import', 
        "genome.maf.gz" )
def buildMAF( infiles, outfile ):

    tracks = " ".join( ["--track=%s" % x[:-len(".import")] for x in infiles] )

    statement = '''
    python %(scriptsdir)s/gtf2gtf.py
           --merge-transcripts --with-utr 
    < transcripts.gtf |
    %(cmd-farm)s --split-at-lines=100 --log=%(outfile)s.log --binary -v 10 
    "python %(scriptsdir)s/snp2maf.py 
          --genome=genome 
          %(tracks)s 
          --reference=mm9 
          --is-gtf 
          --pattern='\(\\\\\\S+\)_pileup'
          --log=%(outfile)s.log" | gzip
    > %(outfile)s
    ''' 

    P.run()

@follows( importTranscripts,
          importTranscriptInformation,
          importGeneStats,
          importGeneInformation )
def prepare():
    pass

@follows( makeEffects, importEffects )
def effects(): pass

@follows( buildAlleles, importAlleles,
          summarizeAllelesPerTranscript,
          summarizeAllelesPerGene )
def alleles(): pass

@follows( prepare, effects )
def full():
    pass

if __name__== "__main__":
    P.checkFiles( ("genome.fasta", "genome.idx" ) )
    P.checkExecutables( ("liftOver",) )
    sys.exit( P.main(sys.argv) )


###############################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Tildon Grant Belgard
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
========================================
RNA-Seq Differential expression pipeline
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The RNA-Seq differential expression pipeline performs differential
expression analysis. It requires three inputs:

   1 one or more genesets in :term:`gtf` formatted files 
   2 mapped reads in :term:`bam` formatted files
   3 design files as :term:`tsv`-separated format

This pipeline works on a single genome.

Overview
========

The pipeline assumes the data derive from multiple tissues/conditions (:term:`experiment`) 
with one or more biological and/or technical replicates (:term:`replicate`). A :term:`replicate`
within each :term:`experiment` is a :term:`track`.

The pipeline performs the following:
   * compute tag counts on an exon, transcript and gene level for each of the genesets
   * estimate FPKM using cufflinks for each of the gene sets
   * perform differential expression analysis. The methods currently implemented are:
      * DESeq
      * EdgeR
      * Cuffdiff

Background
============

The quality of the rnaseq data (read-length, paired-end) determines the
quality of transcript models. For instance, if reads are short (35bp) and/or 
reads are not paired-ended, transcript models will be short and truncated.
In these cases it might be better to concentrate the analysis on only previously
known transcript models. 


Transcripts are the natural choice to measure expression of. However other quantities
might be of interest. Some quantities are biological meaningful, for example differential 
expression from a promotor shared by several trancripts. Other quantities might no biologically
meaningful but are necessary as a technical comprise.
For example, the overlapping transcripts might be hard to resolve and thus might need to be
aggregated per gene. Furthermore, functional annotation is primarily associated with genes
and not individual transcripts. The pipeline attempts to measure transcription and differential
expression for a variety of entities following the classification laid down by :term:`cuffdiff`:

isoform
   Transcript level
gene
   Gene level, aggregates several isoform/transcripts
tss
   Transcription start site. Aggregate all isoforms starting from the same :term:`tss`.
cds
   Coding sequence expression. Ignore reads overlapping non-coding parts of transcripts (UTRs, etc.). Requires
   annotation of the cds and thus only available for :file:`reference.gtf.gz`.

Methods differ in their ability to measure transcription on all levels. 
 
.. todo::
   add promoters and splicing output

Overprediction of differential expression for low-level expressed transcripts with :term:`cuffdiff`
is a `known problem <http://seqanswers.com/forums/showthread.php?t=6283&highlight=fpkm>`_.

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

Reads
+++++

Reads are imported by placing :term:`bam` formatted files are linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`.

Genesets
++++++++

Genesets are imported by placing :term:`gtf` formatted files in the :term:`working directory`.

Design matrices
+++++++++++++++

Design matrices are imported by placing :term:`tsv` formatted files into the :term:`working directory`.
The file has four columns::

      track   include group   pair
      CW-CD14-R1      0       CD14    1
      CW-CD14-R2      0       CD14    1
      CW-CD14-R3      1       CD14    1
      CW-CD4-R1       1       CD4     1
      FM-CD14-R1      1       CD14    2
      FM-CD4-R2       0       CD4     2
      FM-CD4-R3       0       CD4     2
      FM-CD4-R4       0       CD4     2

track
     name of track - should correspond to column header in *infile*
include
     flag to indicate whether or not to include this data  
group
     group indicator - experimental group
pair
     pair that sample belongs to (for paired tests) - set to 0 if the
     design is not paired.

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
|cufflinks_          |>=1.3.0            |transcription levels                            |
+--------------------+-------------------+------------------------------------------------+
|samtools            |>=0.1.16           |bam/sam files                                   |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |                   |working with intervals                          |
+--------------------+-------------------+------------------------------------------------+
|R/DESeq             |                   |differential expression                         |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.42             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+
|bamstats_           |>=1.22             |from CGR, Liverpool                             |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

FPKM values
-----------

The directory :file:`fpkm.dir` contains the output of running cufflinks against each 
set of reads for each geneset. 

The syntax of the output files is ``<geneset>_<track>.cufflinks``. The FPKM values are
uploaded into tables as ``<geneset>_<track>_fpkm`` for per-transcript FPKM values and 
``<geneset>_<track>_genefpkm`` for per-gene FPKM values.

Differential gene expression results
-------------------------------------

Differential expression is estimated for different genesets
with a variety of methods. Results are stored per method in 
subdirectories such as :file:`deseq.dir`, :file:`edger.dir`
or :file:`cuffdiff.dir`.

The syntax of the output files is ``<design>_<geneset>_...``.

:file:`.diff` files:
    The .diff files are :term:`tsv` formatted tables with the result of the differential
    expression analysis. 

    The column names are:

    test_id 
       gene_id/transcript_id/tss_id, etc
    treatment_name  
       name of the treatment
    control_name    
       name of the control
    <>_mean  
       mean expression value of treatment/control
    <>_std   
       standard deviation of treatment/control
    pvalue  
       pvalue of test
    qvalue  
       qvalue of test (FDR if the threshold was set at pvalue)
    l2fold  
       log2 of fold change
    fold    
       fold change, treatment / control
    significant     
       flag: 1 if called significant
    status
       status of test. Values are 
       OK: test ok, FAIL: test failed, NOCALL: no test (insufficient data, etc.).

The results are uploaded into the database as tables called ``<design>_<geneset>_<method>_<level>_<section>``.

Level denotes the biological entity that the differential analysis was performed on.
Level can be one of ``cds``,``isoform``,``tss`` and ``gene``. The later should always be present.

Section is ``diff`` for differential expression results (:file:`.diff`` files) and ``levels`` for expression levels.


Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_rnaseqdiffexpression.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_rnaseqdiffexpression.tgz
   tar -xvzf pipeline_rnaseq.tgz
   cd pipeline_rnaseq
   python <srcdir>/pipeline_rnaseq.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   cufflinks
      cufflinks_ - transcriptome analysis

   tophat
      tophat_ - a read mapper to detect splice-junctions

   deseq
      deseq_ - differential expression analysis

   cuffdiff
      find differentially expressed transcripts. Part of cufflinks_.

   cuffcompare
      compare transcriptomes. Part of cufflinks_.

.. _cufflinks: http://cufflinks.cbcb.umd.edu/index.html
.. _tophat: http://tophat.cbcb.umd.edu/
.. _bamstats: http://www.agf.liv.ac.uk/454/sabkea/samStats_13-01-2011
.. _deseq: http://www-huber.embl.de/users/anders/DESeq/

Code
====

"""

# load modules
from ruffus import *

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
import CGAT.CSV as CSV

import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random

import numpy
import sqlite3
import CGAT.GFF as GFF
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Tophat as Tophat
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError

import CGAT.Expression as Expression

import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineRnaseq as PipelineRnaseq
import CGATPipelines.PipelineMappingQC as PipelineMappingQC
import CGAT.Stats as Stats

# levels of cuffdiff analysis
# (no promotor and splice -> no lfold column)
CUFFDIFF_LEVELS = ("gene", "cds", "isoform", "tss" )

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
import CGATPipelines.PipelineTracks as PipelineTracks

# collect sra nd fastq.gz tracks
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "*.bam" ), "(\S+).bam" )

ALL = PipelineTracks.Sample3()
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

GENESETS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob( "*.gtf.gz" ), "(\S+).gtf.gz" )

###################################################################
###################################################################
###################################################################
def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

#########################################################################
#########################################################################
#########################################################################
## Definition of targets
## Differential expression: design vs geneset
TARGETS_DE = [ ( (x, y, glob.glob("*.bam") ), 
                 "%s_%s.diff" % (P.snip(x, ".tsv"), P.snip(y,".gtf.gz") )) 
               for x,y in itertools.product( glob.glob( "design*.tsv"),
                                             glob.glob( "*.gtf.gz" ) ) ]

## Expression levels: geneset vs bam files
TARGETS_FPKM = [ ( ( "%s.gtf.gz" % x.asFile(), "%s.bam" % y.asFile()),
                   "%s_%s.cufflinks" % (x.asFile(), y.asFile()) )
                 for x,y in itertools.product( GENESETS, TRACKS) ]

#########################################################################
#########################################################################
#########################################################################
## preparation targets

#########################################################################
#########################################################################
#########################################################################
@files(os.path.join(PARAMS["annotations_dir"], "geneset_all.gtf.gz"), "geneset_mask.gtf")
def buildMaskGtf(infile, outfile):
    '''
    This takes ensembl annotations (geneset_all.gtf.gz) and writes out all entries that 
    have a 'source' match to "rRNA" or 'contig' match to "chrM". for use with cufflinks
    '''
    geneset = IOTools.openFile(infile)
    outf = open(outfile, "wb")
    for entry in GTF.iterator(geneset):
        if re.findall("rRNA", entry.source) or re.findall("chrM", entry.contig):
            outf.write("\t".join((map(str,[entry.contig
                              , entry.source
                              , entry.feature
                              , entry.start
                              , entry.end
                              , "."
                              , entry.strand
                              , "."
                              , "transcript_id" + " " + '"' + entry.transcript_id + '"' + ";" + " " + "gene_id" + " " + '"' + entry.gene_id + '"'])))
                               + "\n")

    outf.close()

#########################################################################
#########################################################################
#########################################################################
@transform( "*.gtf.gz", 
            suffix(".gtf.gz"),
            "_geneinfo.load" )
def loadGeneSetGeneInformation( infile, outfile ):
    PipelineGeneset.loadGeneStats( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( "fpkm.dir"  ) )
@files( [ (x, os.path.join( "fpkm.dir",y)) for x,y in TARGETS_FPKM ] )
def runCufflinks(infiles, outfile):
    '''estimate expression levels in each set.
    '''

    gtffile, bamfile = infiles
    to_cluster = True
 
    job_options= "-pe dedicated %i -R y" % PARAMS["cufflinks_threads"]

    track = os.path.basename( P.snip( gtffile, ".gtf.gz" ) )
    
    tmpfilename = P.getTempFilename( "." )
    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )

    gtffile = os.path.abspath( gtffile )
    bamfile = os.path.abspath( bamfile )
    outfile = os.path.abspath( outfile )

    # note: cufflinks adds \0 bytes to gtf file - replace with '.'
    # increase max-bundle-length to 4.5Mb due to Galnt-2 in mm9 with a 4.3Mb intron.
    statement = '''mkdir %(tmpfilename)s; 
    cd %(tmpfilename)s; 
    cufflinks --label %(track)s      
              --GTF <(gunzip < %(gtffile)s)
              --num-threads %(cufflinks_threads)i
              --frag-bias-correct %(bowtie_index_dir)s/%(genome)s.fa
              --library-type %(cufflinks_library_type)s
              %(cufflinks_options)s
              %(bamfile)s 
    >& %(outfile)s;
    perl -p -e "s/\\0/./g" < transcripts.gtf | gzip > %(outfile)s.gtf.gz;
    gzip < isoforms.fpkm_tracking > %(outfile)s.fpkm_tracking.gz;
    gzip < genes.fpkm_tracking > %(outfile)s.genes_tracking.gz;
    '''

    P.run()

    shutil.rmtree( tmpfilename )

#########################################################################
#########################################################################
#########################################################################
@transform( runCufflinks,
            suffix(".cufflinks"), 
            ".load")
def loadCufflinks( infile, outfile ):
    '''load expression level measurements.'''

    track = P.snip( outfile, ".load" )
    P.load( infile + ".genes_tracking.gz",
            outfile = track + "_genefpkm.load",
            options = "--index=gene_id --ignore-column=tracking_id --ignore-column=class_code --ignore-column=nearest_ref_id" )

    track = P.snip( outfile, ".load" )
    P.load( infile + ".fpkm_tracking.gz",
            outfile = track + "_fpkm.load",
            options = "--index=tracking_id --ignore-column=nearest_ref_id --rename-column=tracking_id:transcript_id" )

    P.touch( outfile )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( "cuffdiff.dir" ), buildMaskGtf )
@files( [ (x, os.path.join( "cuffdiff.dir", y)) for x, y in TARGETS_DE ] )
def runCuffdiff( infiles, outfile ):
    '''perform differential expression analysis using deseq.'''
    
    design_file = infiles[0]
    geneset_file = infiles[1]
    bamfiles = infiles[2]

    if PARAMS["cuffdiff_include_mask"]:
        mask_file = os.path.abspath( "geneset_mask.gtf" )
    else:
        mask_file = None

    options = PARAMS["cuffdiff_options"] + " --library-type %s" % PARAMS["cufflinks_library_type"]

    Expression.runCuffdiff( bamfiles, 
                            design_file,
                            geneset_file,
                            outfile,
                            threads = PARAMS.get("cuffdiff_threads",4),
                            cuffdiff_options = options,
                            fdr = PARAMS["cuffdiff_fdr"],
                            mask_file = mask_file )

@transform( runCuffdiff, 
            suffix(".diff"), 
            "_cuffdiff.load" )
def loadCuffdiff( infile, outfile ):
    '''load results from differential expression analysis and produce
    summary plots.

    Note: converts from ln(fold change) to log2 fold change.
   
    The cuffdiff output is parsed. 

    Pairwise comparisons in which one gene is not expressed (fpkm < fpkm_silent)
    are set to status 'NOCALL'. These transcripts might nevertheless be significant.
    '''

    Expression.loadCuffdiff( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
def buildExpressionStats( tables, method, outfile, outdir ):
    '''build expression summary statistics.
    
    Creates some diagnostic plots in

    <exportdir>/<method> directory.
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    def _split( tablename ):
        #parts = tablename.split("_")
        design, geneset = re.match("(.+)_([^_]+)_%s" % method, tablename).groups()
        return design, geneset

        # return re.match("([^_]+)_", tablename ).groups()[0]

    
    keys_status = "OK", "NOTEST", "FAIL", "NOCALL"

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "\t".join( ("design", "geneset", "level", "treatment_name", "control_name", "tested",
                            "\t".join( [ "status_%s" % x for x in keys_status ] ),
                            "significant",
                            "twofold" ) ) + "\n" )

    all_tables = set(Database.getTables( dbhandle ))

    for level in CUFFDIFF_LEVELS:

        for tablename in tables:

            tablename_diff = "%s_%s_diff" % (tablename, level)
            tablename_levels = "%s_%s_diff" % (tablename, level )
            design, geneset = _split( tablename_diff )
            if tablename_diff not in all_tables: continue

            def toDict( vals, l = 2 ):
                return collections.defaultdict( int, [ (tuple( x[:l]), x[l]) for x in vals ] )
            
            tested = toDict( Database.executewait( dbhandle,
                                               """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename_diff)s 
                                    GROUP BY treatment_name,control_name""" % locals() ).fetchall() )
            status = toDict( Database.executewait( dbhandle,
                                                   """SELECT treatment_name, control_name, status, COUNT(*) FROM %(tablename_diff)s 
                                    GROUP BY treatment_name,control_name,status""" % locals() ).fetchall(), 3 )
            signif = toDict( Database.executewait( dbhandle,
                                                   """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename_diff)s 
                                    WHERE significant
                                    GROUP BY treatment_name,control_name""" % locals() ).fetchall() )
            fold2 = toDict( Database.executewait( dbhandle,
                    """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename_diff)s 
                                    WHERE (l2fold >= 1 or l2fold <= -1) AND significant
                                    GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )
            
            for treatment_name, control_name in tested.keys(): 
                outf.write( "\t".join(map(str, (
                                design,
                                geneset,
                                level,
                                treatment_name,
                                control_name,
                                tested[(treatment_name,control_name)],
                                "\t".join( [ str(status[(treatment_name,control_name,x)]) for x in keys_status]),
                                signif[(treatment_name,control_name)],
                                fold2[(treatment_name,control_name)] ) ) ) + "\n" )
                
            ###########################################
            ###########################################
            ###########################################
            # plot length versus P-Value
            data = Database.executewait( dbhandle, 
                                         '''SELECT i.sum, pvalue 
                                 FROM %(tablename_diff)s, 
                                 %(geneset)s_geneinfo as i 
                                 WHERE i.gene_id = test_id AND significant'''% locals() ).fetchall()

            # require at least 10 datapoints - otherwise smooth scatter fails
            if len(data) > 10:
                data = zip(*data)

                pngfile = "%(outdir)s/%(design)s_%(geneset)s_%(level)s_pvalue_vs_length.png" % locals()
                R.png( pngfile )
                R.smoothScatter( R.log10( ro.FloatVector(data[0]) ),
                                 R.log10( ro.FloatVector(data[1]) ),
                                 xlab = 'log10( length )',
                                 ylab = 'log10( pvalue )',
                                 log="x", pch=20, cex=.1 )

                R['dev.off']()

    outf.close()

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "cuffdiff" ) ) )
@transform( loadCuffdiff,
            suffix(".load"), 
            ".plots" )
def buildCuffdiffPlots( infile, outfile ):
    '''create summaries of cufflinks results (including some diagnostic plots)

    Plots are created in the <exportdir>/cuffdiff directory.

    Plots are:

    <geneset>_<method>_<level>_<track1>_vs_<track2>_significance.png
        fold change against expression level
    '''
    ###########################################
    ###########################################
    ## create diagnostic plots
    ###########################################
    outdir = os.path.join( PARAMS["exportdir"], "cuffdiff" )
    
    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    prefix = P.snip( infile, ".load" )

    geneset, method = prefix.split("_")
    
    for level in CUFFDIFF_LEVELS:
        tablename_diff = prefix + "_%s_diff" % level
        tablename_levels = prefix + "_%s_levels" % level
        
        # note that the ordering of EXPERIMENTS and the _diff table needs to be the same
        # as only one triangle is stored of the pairwise results.
        # do not plot "undefined" lfold values (where treatment_mean or control_mean = 0)
        # do not plot lfold values where the confidence bounds contain 0.
        for track1, track2 in itertools.combinations( EXPERIMENTS, 2 ):
            statement = """
                        SELECT CASE WHEN d.treatment_mean < d.control_mean THEN d.treatment_mean 
                                          ELSE d.control_mean END, 
                               d.l2fold, d.significant
                        FROM %(tablename_diff)s AS d
                        WHERE treatment_name = '%(track1)s' AND 
                              control_name = '%(track2)s' AND 
                              status = 'OK' AND
                              treatment_mean > 0 AND 
                              control_mean > 0 
                        """ % locals()
            
            data = zip( *Database.executewait( dbhandle, statement ))
            
            pngfile = "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(track1)s_vs_%(track2)s_significance.png" % locals()
            #ian: Bug fix: moved R.png to after data check so that no plot is started if there is no data
            #     this was leading to R falling over from too many open devices

            if len(data) == 0:
                E.warn( "no plot for %s - %s -%s vs %s" % ( pngfile, level, track1, track2))
                continue

            R.png( pngfile )
            R.plot( ro.FloatVector(data[0]), 
                    ro.FloatVector(data[1]), 
                    xlab = 'min(FPKM)',
                    ylab = 'log2fold',
                    log="x", pch=20, cex=.1,
                    col = R.ifelse( ro.IntVector(data[2]), "red", "black" ) )

            R['dev.off']()

    P.touch( outfile )
        
#########################################################################
#########################################################################
#########################################################################
@follows( loadGeneSetGeneInformation )
@merge( loadCuffdiff,
        "cuffdiff_stats.tsv" )
def buildCuffdiffStats( infiles, outfile ):
    tablenames = [P.toTable(x) for x in infiles ]
    outdir = os.path.dirname( infiles[0] )
    buildExpressionStats( tablenames, "cuffdiff", outfile, outdir )

#########################################################################
#########################################################################
#########################################################################
@transform( buildCuffdiffStats,
            suffix(".tsv"), 
            ".load" )
def loadCuffdiffStats( infile, outfile ):
    '''import cuffdiff results.'''
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@merge( os.path.join( PARAMS["annotations_dir"], 
                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
        "coding_exons.gtf.gz" )
def buildCodingExons( infile, outfile ):
    '''compile set of protein coding exons.

    This set is used for splice-site validation
    '''

    to_cluster = True
    statement = '''
    zcat %(infile)s 
    | awk '$2 == "protein_coding" && $3 == "CDS"'
    | perl -p -e "s/CDS/exon/" 
    | python %(scriptsdir)s/gtf2gtf.py --merge-exons --log=%(outfile)s.log 
    | gzip 
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( "*.gtf.gz",
            suffix(".gtf.gz"),
            ".unionintersection.bed.gz" )
def buildUnionIntersectionExons( infile, outfile ):
    '''build union/intersection genes according to Bullard et al. (2010) BMC Bioinformatics.

    Builds a single-segment bed file.
    '''

    to_cluster = True
    statement = '''
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --intersect-transcripts --with-utr --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2gff.py --is-gtf --crop-unique  --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --log=%(outfile)s.log
    | sort -k1,1 -k2,2n
    | gzip 
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( "*.gtf.gz",
            suffix(".gtf.gz"),
            ".union.bed.gz" )
def buildUnionExons( infile, outfile ):
    '''build union genes.

    Exons across all transcripts of a gene are merged.
    They are then intersected between genes to remove any overlap.

    Builds a single-segment bed file.    
    '''

    to_cluster = True
    statement = '''
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --merge-exons --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2gff.py --is-gtf --crop-unique  --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --log=%(outfile)s.log
    | sort -k1,1 -k2,2n
    | gzip 
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
# note - needs better implementation, currently no dependency checks.
@follows( buildUnionExons, mkdir( "exon_counts.dir" ) )
@files( [ ( ("%s.bam" % x.asFile(), "%s.union.bed.gz" % y.asFile() ),
            ("exon_counts.dir/%s_vs_%s.bed.gz" % (x.asFile(),y.asFile() ) ) )
          for x,y in itertools.product( TRACKS, GENESETS) ] )
def buildExonLevelReadCounts( infiles, outfile ):
    '''compute coverage of exons with reads.
    '''

    infile, exons = infiles

    to_cluster = True

    # note: needs to set flags appropriately for
    # single-end/paired-end data sets
    # set filter options
    # for example, only properly paired reads
    paired = False
    if paired:
        flag_filter = "-f 0x2"
    else:
        flag_filter = ""

    # note: the -split option only concerns the stream in A - multiple
    # segments in B are not split. Hence counting has to proceed via
    # single exons - this can lead to double counting if exon counts
    # are later aggregated.

    statement = '''
    samtools view -b %(flag_filter)s -q %(deseq_min_mapping_quality)s %(infile)s
    | coverageBed -abam stdin -b %(exons)s -split
    | sort -k1,1 -k2,2n
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@collate(buildExonLevelReadCounts,
         regex(r"exon_counts.dir/(.+)_vs_(.+)\.bed.gz"),  
         r"\2.exon_counts.load")
def loadExonLevelReadCounts( infiles, outfile ):
    '''load exon level read counts.
    '''
    
    to_cluster = True

    # aggregate not necessary for bed12 files, but kept in
    # ims: edited so that picks up chromosome, start pos and end pos for downstream use.
    src = " ".join( [ "<( zcat %s | cut -f 1,2,3,4,7 )" % x for x in infiles] )

    tmpfile = P.getTempFilename( "." )
    tmpfile2 = P.getTempFilename( "." )
    
    statement = '''paste %(src)s 
                > %(tmpfile)s'''
    
    P.run()

    tracks = [P.snip(x, ".bed.gz" ) for x in infiles ]
    tracks = [re.match( "exon_counts.dir/(\S+)_vs.*", x).groups()[0] for x in tracks ]

    outf = IOTools.openFile( tmpfile2, "w")
    outf.write( "gene_id\tchromosome\tstart\tend\t%s\n" % "\t".join( tracks ) )
    
    for line in open( tmpfile, "r" ):
        data = line[:-1].split("\t")
        # ims: edit so that now skips five and ens_id is in 3rd index
        genes = list(set([ data[x] for x in range(3,len(data), 5 ) ]))
        # ims: add entries for chromosome, start and ends
        chrom = list(set( [data[x] for x in range(0,len(data), 5 ) ]))
        starts = list(set( [data[x] for x in range(1,len(data), 5 ) ]))
        ends = list(set( [data[x] for x in range(2,len(data), 5 ) ]))
        # ims: edit as value is now in postion 4 and there are 5 columns per line
        values = [ data[x] for x in range(4,len(data), 5 ) ]
        # ims: extra assets for chrom, starts and ends
        assert len(genes) == 1, "paste command failed, wrong number of genes per line"
        assert len(chrom) == 1, "paste command failed, wrong number of chromosomes per line"
        assert len(starts) == 1, "paste command failed, wrong number of starts per line"
        assert len(ends) == 1, "paste command failed, wrong number of ends per line"
        # ims: add extra coloumns into output
        outf.write( "%s\t%s\t%s\t%s\t%s\n" % (genes[0], chrom[0], starts[0], ends[0], "\t".join(map(str, values) ) ) )
    
    outf.close()

    P.load( tmpfile2, outfile )

    os.unlink( tmpfile )
    os.unlink( tmpfile2 )

#########################################################################
#########################################################################
#########################################################################
@follows( buildUnionExons, mkdir( "gene_counts.dir" ) )
@files( [ ( ("%s.bam" % x.asFile(), "%s.gtf.gz" % y.asFile() ),
            ("gene_counts.dir/%s_vs_%s.tsv.gz" % (x.asFile(),y.asFile() ) ) )
          for x,y in itertools.product( TRACKS, GENESETS) ] )
def buildGeneLevelReadCounts( infiles, outfile ):
    '''compute coverage of exons with reads.
    '''

    infile, exons = infiles

    to_cluster = True

    statement = '''
    zcat %(exons)s 
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=genes
          --bam-file=%(infile)s 
          --counter=length
          --prefix="exons_"
          --counter=read-counts 
          --prefix=""
          --counter=read-coverage
          --prefix=coverage_
    | gzip
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildGeneLevelReadCounts,
            suffix(".tsv.gz"),
            "_gene_counts.load" )
def loadGeneLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=gene_id" )
    
#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( "extension_counts.dir" ) )
@transform( "*.bam", 
            regex(r"(\S+).bam"), 
            r"extension_counts.dir/\1.extension_counts.tsv.gz" )
def buildGeneLevelReadExtension( infile, outfile ):
    '''compute extension of cds. 

    Known UTRs are counted as well.
    '''

    to_cluster = True

    cds = os.path.join( PARAMS["annotations_dir"],
                        PARAMS_ANNOTATIONS["interface_geneset_cds_gtf"] )
    
    territories = os.path.join( PARAMS["annotations_dir"],
                                PARAMS_ANNOTATIONS["interface_territories_gff"] )

    utrs = os.path.join( PARAMS["annotations_dir"],
                         PARAMS_ANNOTATIONS["interface_annotation_gff"] )

    if "geneset_remove_contigs" in PARAMS:
        remove_contigs = '''| awk '$1 !~ /%s/' ''' % PARAMS["geneset_remove_contigs"] 
    else:
        remove_contigs = ""

    statement = '''
    zcat %(cds)s 
    %(remove_contigs)s
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=genes
          --bam-file=%(infile)s 
          --counter=position
          --counter=read-extension
          --output-filename-pattern=%(outfile)s.%%s.tsv.gz
          --filename-gff=%(territories)s
          --filename-gff=%(utrs)s
    | gzip
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("transcript_counts.dir") )
@files( [ ( ("%s.bam" % x.asFile(), "%s.gtf.gz" % y.asFile() ),
            ("transcript_counts.dir/%s_vs_%s.tsv.gz" % (x.asFile(),y.asFile() ) ) )
          for x,y in itertools.product( TRACKS, GENESETS) ] )
def buildTranscriptLevelReadCounts( infiles, outfile):
    '''count reads falling into transcripts of protein coding 
       gene models.

    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.
       
    '''
    infile, geneset = infiles
    
    to_cluster = True

    statement = '''
    zcat %(geneset)s 
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=transcripts
          --bam-file=%(infile)s 
          --counter=length
          --prefix="exons_"
          --counter=read-counts 
          --prefix=""
          --counter=read-coverage
          --prefix=coverage_
    | gzip
    > %(outfile)s
    '''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildTranscriptLevelReadCounts,
            suffix(".tsv.gz"),
            ".load" )
def loadTranscriptLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=transcript_id" )

#########################################################################
#########################################################################
#########################################################################
@collate(buildExonLevelReadCounts,
         regex(r"exon_counts.dir/(.+)_vs_(.+)\.bed.gz"),  
         r"exon_counts.dir/\2.exon_counts.tsv.gz")
def aggregateExonLevelReadCounts( infiles, outfile ):
    '''aggregate exon level tag counts for each gene.

    coverageBed adds the following four columns:

    1) The number of features in A that overlapped (by at least one base pair) the B interval.
    2) The number of bases in B that had non-zero coverage from features in A.
    3) The length of the entry in B.
    4) The fraction of bases in B that had non-zero coverage from features in A.

    For bed6: use column 7
    For bed12: use column 13

    This method uses the maximum number of reads found in any exon as the tag count.

    '''
    
    to_cluster = True

    # aggregate not necessary for bed12 files, but kept in
    src = " ".join( [ "<( zcat %s | sort -k4,4 | groupBy -i stdin -g 4 -c 7 -o %s | sort -k1,1)" % (x,PARAMS['counting_aggregate']) for x in infiles ] )

    tmpfile = P.getTempFilename( "." )
    
    statement = '''paste %(src)s 
                > %(tmpfile)s''' % locals()
    
    P.run()

    tracks = [P.snip(x, ".bed.gz" ) for x in infiles ]
    tracks = [re.match( "exon_counts.dir/(\S+)_vs.*", x).groups()[0] for x in tracks ]

    outf = IOTools.openFile( outfile, "w")
    outf.write( "gene_id\t%s\n" % "\t".join( tracks ) )
    
    for line in open( tmpfile, "r" ):
        data = line[:-1].split("\t")
        genes = list(set([ data[x] for x in range(0,len(data), 2 ) ]))
        values = [ data[x] for x in range(1,len(data), 2 ) ]
        assert len(genes) == 1, "paste command failed, wrong number of genes per line"
        outf.write( "%s\t%s\n" % (genes[0], "\t".join(map(str, values) ) ) )
    
    outf.close()

    os.unlink( tmpfile )

#########################################################################
#########################################################################
#########################################################################
@transform( ( aggregateExonLevelReadCounts),
            suffix(".tsv.gz"),
            ".load" )
def loadAggregateExonLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=gene_id" )

TARGETS_DE = [ ( (x, y, glob.glob("*.bam"),
                  "exon_counts.dir/%s.exon_counts.tsv.gz" % P.snip(y,
                                                    ".gtf.gz",path=True)),
                 "%s_%s.diff" % (P.snip(x, ".tsv"), P.snip(y,".gtf.gz") )) 
               for x,y in itertools.product( glob.glob( "design*.tsv"),
                                             glob.glob( "*.gtf.gz" ) ) ]

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("deseq.dir"), loadAggregateExonLevelReadCounts )
@files( [ (x, os.path.join( "deseq.dir", y)) for x, y in TARGETS_DE ] )
def runDESeq( infiles, outfile ):
    '''perform differential expression analysis using deseq.'''

    to_cluster = True 
    design_file, geneset_file, bamfiles, count_file = infiles

    track = P.snip( outfile, ".diff")

    statement = '''python %(scriptsdir)s/Expression.py
              --method=deseq
              --filename-tags=%(count_file)s
              --filename-design=%(design_file)s
              --output-filename-pattern=%(track)s_
              --outfile=%(outfile)s
              --fdr=%(deseq_fdr)f
              --deseq-fit-type=%(deseq_fit_type)s
              --deseq-dispersion-method=%(deseq_dispersion_method)s
              > %(outfile)s.log '''

    P.run()


#########################################################################
#########################################################################
#########################################################################
@transform( runDESeq, suffix(".diff"), "_deseq.load" )
def loadDESeq( infile, outfile ):
    '''load differential expression results.
    '''
    # add gene level follow convention "<level>_diff"
    tablename = P.toTable(outfile) + "_gene_diff" 
    statement = '''cat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=test_id
              --table=%(tablename)s 
            > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( loadGeneSetGeneInformation )
@merge( loadDESeq, "deseq_stats.tsv" )
def buildDESeqStats( infiles, outfile ):
    tablenames = [P.toTable( x ) for x in infiles ] 
    outdir = os.path.dirname( infiles[0] )
    buildExpressionStats( tablenames, "deseq", outfile, outdir )

#########################################################################
#########################################################################
#########################################################################
@transform( buildDESeqStats,
            suffix(".tsv"), 
            ".load" )
def loadDESeqStats( infile, outfile ):
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@follows( aggregateExonLevelReadCounts,mkdir("edger.dir") )
@files( [ (x, os.path.join( "edger.dir", y)) for x, y in TARGETS_DE ] )
def runEdgeR( infiles, outfile ):
    '''perform differential expression analysis using edger.'''

    to_cluster = True 
    design_file, geneset_file, bamfiles, infile = infiles


    track = P.snip( outfile, ".diff")

    statement = '''python %(scriptsdir)s/Expression.py
              --method=edger
              --filename-tags=%(infile)s
              --filename-design=%(design_file)s
              --output-filename-pattern=%(track)s_
              --outfile=%(outfile)s
              --fdr=%(edger_fdr)f
              > %(outfile)s.log '''

    P.run()
    
#########################################################################
#########################################################################
#########################################################################
@transform( runEdgeR, suffix(".diff"), "_edger.load" )
def loadEdgeR( infile, outfile ):
    '''load differential expression results.
    '''
    # add gene level follow convention "<level>_diff"
    tablename = P.toTable(outfile) + "_gene_diff" 
    statement = '''cat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=test_id
              --table=%(tablename)s 
            > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( loadGeneSetGeneInformation )
@merge( loadEdgeR, "edger_stats.tsv" )
def buildEdgeRStats( infiles, outfile ):
    tablenames = [P.toTable( x ) for x in infiles ] 
    outdir = os.path.dirname( infiles[0] )
    buildExpressionStats( tablenames, "edger", outfile, outdir )

#########################################################################
#########################################################################
#########################################################################
@transform( buildEdgeRStats,
            suffix(".tsv"), 
            ".load" )
def loadEdgeRStats( infile, outfile ):
    P.load( infile, outfile )

###################################################################
###################################################################
###################################################################
@follows( loadCufflinks, loadGeneLevelReadCounts )
def expression(): pass

mapToTargets = { 'cuffdiff': loadCuffdiffStats,
                 'deseq': loadDESeqStats,
                 'edger': loadEdgeRStats,
                 }
TARGETS_DIFFEXPRESSION = [ mapToTargets[x] for x in P.asList( PARAMS["methods"] ) ]

@follows( *TARGETS_DIFFEXPRESSION )
def diff_expression(): pass

###################################################################
###################################################################
###################################################################
@jobs_limit(1,"R")
@follows( mkdir("tagplots.dir"), aggregateExonLevelReadCounts )
@files( [ (x, os.path.join( "tagplots.dir", y)) for x, y in TARGETS_DE ] )
def plotRNASEQTagData( infiles, outfile ):
    '''perform differential expression analysis using deseq.'''
    
    design_file, geneset_file = infiles[0], infiles[1]

    infile = os.path.join( "exon_counts.dir", P.snip( geneset_file, ".gtf.gz") + ".exon_counts.tsv.gz" )
    Expression.plotTagStats( infile, design_file, outfile )

    P.touch( outfile )

mapToQCTargets = { 'cuffdiff': runCuffdiff,
                   'deseq': runDESeq,
                   'edger': runEdgeR,
                   }
QCTARGETS = [ mapToQCTargets[x] for x in P.asList( PARAMS["methods"] ) ]

@jobs_limit(1,"R")
@transform( QCTARGETS,
            suffix(".diff"),
            ".plots" )
def plotRNASEQDEData( infile, outfile ):
    '''plot differential expression stats'''
    Expression.plotDETagStats( infile, outfile )
    P.touch(outfile)

###################################################################
###################################################################
###################################################################
@follows( plotRNASEQTagData, plotRNASEQDEData )
def qc(): pass

###################################################################
###################################################################
###################################################################
@follows( expression, diff_expression, qc )
def full(): pass

###################################################################
###################################################################
###################################################################
@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

###################################################################
###################################################################
###################################################################
@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating documentation" )
    P.run_report( clean = False )

###################################################################
###################################################################
###################################################################
@follows( mkdir( "%s/bamfiles" % PARAMS["web_dir"]), 
          mkdir("%s/genesets" % PARAMS["web_dir"]),
          mkdir("%s/classification" % PARAMS["web_dir"]),
          mkdir("%s/differential_expression" % PARAMS["web_dir"]),
          update_report,
          )
def publish():
    '''publish files.'''
    # publish web pages
    P.publish_report()

    # publish additional data
    web_dir = PARAMS["web_dir"]
    project_id = P.getProjectId()

    # directory, files
    exportfiles = {
        "bamfiles" : glob.glob( "*.bam" ) + glob.glob( "*.bam.bai" ),
        "genesets": [ "lincrna.gtf.gz", "abinitio.gtf.gz" ],
        "classification": glob.glob("*.class.tsv.gz") ,
        "differential_expression" : glob.glob( "*.cuffdiff.dir" ),
        }
    
    bams = []

    for targetdir, filenames in exportfiles.iteritems():
        for src in filenames:
            dest = "%s/%s/%s" % (web_dir, targetdir, src)
            if dest.endswith( ".bam"): bams.append( dest )
            dest = os.path.abspath( dest )
            if not os.path.exists( dest ):
                os.symlink( os.path.abspath(src), dest )
    
    # output ucsc links
    for bam in bams: 
        filename = os.path.basename( bam )
        track = P.snip( filename, ".bam" )
        print """track type=bam name="%(track)s" bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/bamfiles/%(filename)s""" % locals()

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

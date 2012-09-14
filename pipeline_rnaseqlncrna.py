################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2012 Nick Ilott
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
==============================
Long non-coding RNA pipeline
==============================


Overview
========

The pipeline_rnaseqLncRNA.py pipeline aims to predict lncRNAs from an ab initio
assembly of transcripts.

It requires that raw reads have been mapped to a reference transcriptome and assembled
into transcript models using cufflinks and therefore makes the assumption that the
transcript building has gone to plan.


Details
========

Prediction of LncRNA are not based on a reference annotation to begin with, although 
downstream comparisons are made to a reference non-coding gene set. The main features
of the pipeline are as follows:

* Build a coding gene set based on an ab initio assembly.
   The ab initio assembly is filtered for protein coding genes. In this step only genes 
   that are compatible with an annotated protein coding transcript are kept - this will
   reduce noise that is associated with a large number of incomplete transfrags. The 
   filtering is based on output from cuffcompare (class code "=").

* build a non-coding gene set.
   A reference non-coding set of transcripts is built by filtering a provided ensembl
   reference set (usually a set that is built from the transcript building pipeline) for 
   transcripts that do not belong to one of the following biotypes

   protein_coding\n
   ambiguous_orf\n
   retained_intron\n
   sense_intronic\n
   antisense\n
   sense_overlapping\n


   This set of non-coding transcripts is required in the filtering of the ab initio geneset 
   for LncRNA prediction. This is because many putative lncRNA have multiple associated
   biotypes. For example MALAT1 is described as both a lncRNA and a processed transcript. To 
   avoid removing known ncRNA we therefore check for existence of putative transcripts in this
   set.

* Build a putative lncRNA gene set.
   The ab initio set of lncRNA are filtered to remove overlapping protein coding exons. This 
   filtering is performed on the level of the transcript - although there may be multiple isoform
   predictions per lncRNA, at this point the sensitivity of lncRNA prediction is increased. 
   Antisense transcripts overlapping protein coding transcripts are retained.
   
* Filter putative lncRNA gene set
   Due to many fragments being produced from RNA-seq data, putative single exon lncRNA are flagged
   with in the lncRNA gtf file so it is easy to filter for the more reliable multi-exonic lncRNA. 
   Although many single exon lncRNA are likely to be artifacts, we assess the overlap of putative
   single exon lncRNA with sets of lncRNA that have been previously identified. If an overlap is
   found with a transcript in the reference set then the reference is added to the lncRNA gene set.
   This means that true single exon lncRNA are still picked up - as long as there is previous evidence
   to support their existence.

* Build final lncRNA gene set
   The putative set of lncRNA are assessed for coding potential using the coding potential calculator
   (CPC). Any lncRNA that are annotated as 'coding' in this analysis are removed from downstream analysis.

* Combine coding and non-coding gene sets
   In order to assess expression levels between genes within samples i.e. protein coding vs. lncRNA, it
   is required that the FPKM estimation be made on a complete geneset. Therefore the lncRNA geneset is 
   concatenated to the protein coding gene set for use in downstream analysis.


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`)


input
-----

The pipeline generally assumes that transcripts have been assembled using the pipeline_rnaseqtranscripts.py
pipeline using cufflinks.

Files are supplied in the working directory. They are specified in the configuration file
and refer to:

* A coding geneset that is the output from a cufflinks transcript assembly.

    abinitio_coding = :file:`<name>.gtf.gz` 

* An abinitio geneset that is the output from a cufflinks transcript assembly. This is to be used for lncRNA prediction. 
  (note that this may be different to the abinitio_coding geneset). 

    abinitio_lncrna = :file:`<name>.gtf.gz`

* A reference geneset containing known protein coding transcripts. This is used for comparisons in the report.

    refcoding = :file:`<name>.gtf.gz`

* A reference geneset from ensembl with all known expressed transcripts

    reference = :file:`<name>.gtf.gz`

* An optional geneset containing previously identified lncRNA. If this is not supplied then the pipeline uses
  a reference non-coding set from the ensembl reference.

    previous = :file:`<name>.gtf.gz`


Pipeline output
----------------

The pipeline produces three main files of interest:

+------------------------------------+--------------------------------------------------+
|           Filename                 |             Description                          |
+------------------------------------+--------------------------------------------------+
|                                    |Ab initio set of lncRNA transcripts filtered for  |
|:file:`lncrna_final.class.gtf.gz`   |single exon status (excl.previously observed) and |
|                                    |classified relative to protein coding transcripts |
+------------------------------------+--------------------------------------------------+
|                                    |Ab inito assembled protein coding transcipts - for|
|:file:`<name>_coding.gtf.gz`        |a comparable set to lncRNA transcripts            |
|                                    |                                                  |
+------------------------------------+--------------------------------------------------+
|                                    |Combined set from the two sets above. to be used  |
|:file:`transcripts.gtf.gz`          |for downstream FPKM estimation and differential   |
|                                    |expression analysis                               |
+------------------------------------+--------------------------------------------------+


code
=====
'''

##########################################################
##########################################################
##########################################################

# load modules                                                                                                                                                                                                                               
from ruffus import *
import Experiment as E
import logging as L
import Database, CSV
import numpy as np
import fnmatch
import sqlite3
import PipelineMappingQC
import GTF
import IOTools
import gzip
import GFF
import PipelineRnaseq
import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random
from rpy2.robjects import r as R
import Expression
import IndexedGenome
import PipelineLncRNA

import Pipeline as P

# get parameters
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")]
     ,"pipeline.ini"],
    defaults = {"annotations_annotations_dir": "",
                "genesets_abinitio_coding": "pruned.gtf.gz",
                "genesets_abinitio_lncrna": "pruned.gtf.gz",
                "genesets_reference": "reference.gtf.gz",
                "genesets_refcoding": "refcoding.gtf.gz",
                "genesets_previous": ""})

PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_annotations_dir"],
                                       "pipeline_annotations.py" )



#########################################################################
#########################################################################
#########################################################################
def connect():
    '''connect to database.

    This method also attaches to helper databases.
    '''

    dbh = sqlite3.connect( PARAMS["database"] )
    statement = '''ATTACH DATABASE '%s' as annotations''' % (PARAMS["annotations_annotations_database"])
    cc = dbh.cursor()
    cc.execute( statement )
    cc.close()

    return dbh

#########################################################################
#########################################################################
#########################################################################
def Rconnect():
    '''
    connect to a database through R
    '''
    R('''library("RSQLite")''')
    R('''library("sciplot")''')
    R('''drv <- dbDriver("SQLite")''')
    R('''con <- dbConnect(drv, dbname = "%s") ''' % PARAMS["database"])
    return R('''con''')

#########################################################################
#########################################################################
#########################################################################
def filenameToTablename(filename):
    '''
    converts filename containing "." to tablename where "." converted to "_"
    '''
    return filename.replace(".", "_" )

#########################################################################
#########################################################################
#########################################################################
def updateFile(filename):
    '''                                                                                                                                                                                                                                      
    create empty file for updating purposes                                                                                                                                                                                                  
    '''

    outf = open(filename, "w")
    outf.write("file created for ruffus update")
    outf.close()

#########################################################################
#########################################################################
#########################################################################
def tabSplit(line):
    '''
    generic split by newline and tab for reading tsv files
    '''
    return line[:-1].split("\t")

#########################################################################
#########################################################################
#########################################################################
@follows(mkdir("gtfs"))
@merge( [PARAMS["genesets_abinitio_coding"], PARAMS["genesets_reference"]]
       , os.path.join("gtfs", P.snip(PARAMS["genesets_abinitio_coding"], ".gtf.gz") + "_coding.gtf.gz") )
def buildCodingGeneSet(infiles, outfile):
    '''
    takes the output from cuffcompare of a transcript
    assembly and filters for annotated protein coding
    genes. 
    
    NB "pruned" refers to nomenclature in the transcript
    building pipeline - transcripts that appear in at least
    two samples.
    
    Because an abinitio assembly will often contain
    fragments of known transcripts and describe them as 
    novel, the default behaviour is to produce a set that
    is composed of 'complete' or 'contained' transcripts
    i.e. nothing novel. This may underestimate the number 
    of transcripts that are actually expressed
    '''
    PipelineLncRNA.buildCodingGeneSet(infiles[0], infiles[1], outfile)

##########################################################
##########################################################
##########################################################
@follows(buildCodingGeneSet)
@transform(PARAMS["genesets_refcoding"], regex(r"(\S+)/(\S+).gtf.gz"), add_inputs(buildCodingGeneSet), r"gtfs/\2.gtf.gz")
def buildRefcodingGeneSet(infiles, outfile):
    '''
    builds a refcoding geneset based on the genes that are present in
    the abinitio assembly
    '''
    PipelineLncRNA.buildRefcodingGeneSet(infiles[1], infiles[0], outfile )

##########################################################
##########################################################
##########################################################
@follows(mkdir("gtfs"))
@files(PARAMS["genesets_reference"], "gtfs/refnoncoding.gtf.gz")
def buildRefnoncodingGeneSet(infile, outfile):
    '''
    filter the refnoncoding geneset for things that are described in ensembl
    as being:
    Ambiguous_orf
    Retained_intron
    Sense_intronic
    antisense
    Sense_overlapping
    Processed transcript
    '''
    PipelineLncRNA.buildRefnoncodingGeneSet(infile, outfile)

##########################################################
##########################################################
##########################################################
@follows(mkdir("gtfs"))
@merge( (PARAMS["genesets_abinitio_lncrna"], PARAMS["genesets_reference"], buildRefnoncodingGeneSet,
         os.path.join(PARAMS["annotations_annotations_dir"], PARAMS_ANNOTATIONS["interface_pseudogenes_gtf"] ),
         os.path.join(PARAMS["annotations_annotations_dir"], PARAMS_ANNOTATIONS["interface_numts_gtf"] ),
         ), "gtfs/lncrna.gtf.gz" )

def buildLncRNAGeneSet(infiles, outfile):
    '''
    build lncRNA gene set. 
    
    This is a set of transcripts in the abinitio set that
    do not overlap at any protein coding or pseudogene transcripts
    or additional biotypes from ensembl that are unwanted
    (exons) in a reference gene set.
    
    Transcripts need to have a length of at least 200 bp.
    '''
    PipelineLncRNA.buildLncRNAGeneSet( infiles[0], infiles[1], infiles[2], infiles[3], infiles[4], outfile, PARAMS["lncrna_min_length"] )        
               
##########################################################################
##########################################################################
##########################################################################
@transform(buildLncRNAGeneSet, suffix(".gtf.gz"), "_flag.gtf.gz")
def flagExonStatus(infile, outfile):
    '''
    Adds an attribute to the gtf entry dependent on whether the lncRNA 
    is multi or single exon
    '''
    PipelineLncRNA.flagExonStatus(infile, outfile)

##########################################################################
##########################################################################
##########################################################################
if PARAMS["genesets_previous"]:
    @transform(flagExonStatus, regex(r"(\S+)_flag.gtf.gz")
               , add_inputs([PARAMS["genesets_previous"], buildRefnoncodingGeneSet]), r"\1_filtered.gtf.gz")
    def buildFilteredLncRNAGeneSet(infiles, outfile):
        '''
        Creates a filtered lncRNA geneset. 
        This geneset will not include any single exon lncRNA 
        unless it has been seen previously i.e. it overlaps
        a previously identified lncRNA
        '''
        PipelineLncRNA.buildFilteredLncRNAGeneSet(infiles[0], infiles[1][1], outfile, geneset_previous = infiles[1][0])
else:
    E.info("no previous lncRNA set provided: Using refnoncoding set")
    @transform(flagExonStatus, regex(r"(\S+)_flag.gtf.gz")
               , add_inputs(["refnoncoding.gtf.gz"]), r"\1_filtered.gtf.gz")
     
    def buildFilteredLncRNAGeneSet(infiles, outfile):
        PipelineLncRNA.buildFilteredLncRNAGeneSet(infiles[0],infiles[1][0], outfile)
    
##########################################################################
##########################################################################
##########################################################################
@transform(buildFilteredLncRNAGeneSet, suffix(".gtf.gz"),add_inputs(PARAMS["genesets_refcoding"]), ".class.gtf.gz")
def classifyFilteredLncRNA(infiles, outfile):
    '''
    classifies all lincRNA before cpc filtering to define any classes that
    are represented in the coding set that are  filtered
    '''
    PipelineLncRNA.classifyLncRNAGenes(infiles[0], infiles[1], outfile, dist = PARAMS["lncrna_dist"])

##########################################################################
##########################################################################
##########################################################################
@follows(mkdir("fasta"))
@transform(buildFilteredLncRNAGeneSet, regex(r"gtfs/(\S+).gtf.gz"), r"fasta/\1.fasta")
def buildLncRNAFasta(infile, outfile):
    '''
    create fasta file from lncRNA geneset for testing coding
    potential of transcripts
    '''
    genome = os.path.join(PARAMS["general_genomedir"], PARAMS["genome"] + ".fasta")
    statement = '''zcat %(infile)s | python %(scriptsdir)s/gff2fasta.py --genome-file=%(genome)s --log=%(outfile)s.log --is-gtf > %(outfile)s'''
    P.run()

##########################################################################
##########################################################################
##########################################################################
@transform(buildLncRNAFasta, regex(r"fasta/(\S+).fasta"), r"cpc/\1.cpc.result")
def runCPC(infile, outfile):
    '''
    run coding potential calculations on lncRNA geneset
    '''
    result_table = P.snip(infile, ".fasta") + ".result"
    result_evidence = P.snip(outfile, ".result") + ".evidence"
    working_dir = "cpc"
    statement = '''%(scriptsdir)s/cpc.sh %(infile)s %(result_table)s %(working_dir)s %(result_evidence)s '''
    P.run()
    statement = '''mv %(result_table)s cpc/'''
    P.run()

##########################################################################
##########################################################################
##########################################################################
@transform(runCPC, suffix(".result"), ".load")
def loadCPCResults(infile, outfile):
    '''
    load the results of the cpc analysis
    '''
    tablename = filenameToTablename(os.path.basename(infile))
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log 
                   --header=transcript_id,feature,C_NC,CP_score --index=transcript_id < %(infile)s > %(outfile)s'''
    P.run()

##########################################################################
##########################################################################
##########################################################################
@follows(loadCPCResults)
@transform(buildFilteredLncRNAGeneSet, regex(r"(\S+)_filtered.gtf.gz"), r"\1_final.gtf.gz")
def buildFinalLncRNAGeneSet(infile, outfile):
    '''
    the final lncRNA gene set consists of transcripts that pass
    the initial filtering stage i.e. are;
    multi-exonic/previously seen single exon transcripts
    display low evidence for coding potential
    '''
    
    # filter based on coding potential
    PipelineLncRNA.buildFinalLncRNAGeneSet(infile, "lncrna_filtered_cpc_result", outfile)

##########################################################################
##########################################################################
##########################################################################
@transform([buildCodingGeneSet, buildFinalLncRNAGeneSet, buildRefcodingGeneSet], regex(r"(\S+)/(\S+).gtf.gz"), r"\2.stats")
def buildGeneSetStats(infile, outfile):
    '''
    counts:
    no. of transcripts
    no. genes
    average number of exons per transcript
    average number of exons per gene
    no. multi-exon transcripts
    no. single exon transcripts
    no. multi-exon genes
    no. single exon genes

    in the coding and lncRNA genesets
    '''
    outf = open(outfile, "w")
    outf.write("\t".join(["no_transcripts"
                , "no_genes"
                , "no_exons_per_transcript"
                , "no_exons_per_gene"
                , "no_single_exon_transcripts"
                , "no_multi_exon_transcripts"
                , "no_single_exon_genes"
                , "no_multi_exon_genes"]) + "\n")
    outf.write("\t".join(map(str, [PipelineLncRNA.CounterTranscripts(infile).count()
                          , PipelineLncRNA.CounterGenes(infile).count()
                          , PipelineLncRNA.CounterExonsPerTranscript(infile).count()
                          , PipelineLncRNA.CounterExonsPerGene(infile).count()
                          , PipelineLncRNA.CounterSingleExonTranscripts(infile).count()
                          , PipelineLncRNA.CounterMultiExonTranscripts(infile).count()
                          , PipelineLncRNA.CounterSingleExonGenes(infile).count()
                          , PipelineLncRNA.CounterMultiExonGenes(infile).count()])))

##########################################################################
##########################################################################
##########################################################################
@transform(buildGeneSetStats, suffix(".stats"), ".load")
def loadGeneSetStats(infile, outfile):
    '''
    load stats on coding and lncRNA gene sets
    '''
    tablename = filenameToTablename(infile)
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log < %(infile)s > %(outfile)s'''
    P.run()

##########################################################################
##########################################################################
##########################################################################
@transform(buildFinalLncRNAGeneSet, regex(r"(\S+).gtf.gz")
           , add_inputs(PARAMS["genesets_refcoding"]), r"\1.class.gtf.gz")
def classifyLncRNA(infiles, outfile):
    '''
    
    Classify lncRNA realtive to protein coding loci
    
    Classify lincRNA in terms of their relationship to 
    protein coding genes - creates indices for intervals on the 
    fly - mayb should be creating additional annotations:

    antisense - transcript overlapping protein coding exons on opposite strand
    antisense_upstream - transcript < 2kb from tss on opposite strand
    antisense_downstream - transcript < 2kb from gene end on opposite strand
    sense_upstream - transcript < 2kb from tss on same strand
    sense_downstream - transcript < 2kb from gene end on same strand
     intergenic - >2kb from anyt protein coding gene
    intronic - overlaps protein coding gene intron on same strand
    antisense_intronic - overlaps protein coding intron on opposite strand
    '''
    
    PipelineLncRNA.classifyLncRNAGenes(infiles[0], infiles[1], outfile, dist = PARAMS["lncrna_dist"])

##########################################################################
##########################################################################
##########################################################################
@transform(classifyLncRNA, suffix(".gtf.gz"), ".load")
def loadLncRNAClass(infile, outfile):
    '''
    load the lncRNA classifications
    '''
    tablename = filenameToTablename(P.snip(infile, ".gtf.gz"))
    
    # just load each transcript with its classification
    temp = P.getTempFile()
    for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(infile))):
        temp.write("%s\t%s\t%s\n" % (transcript[0].transcript_id, transcript[0].gene_id, transcript[0].source))
    temp.close()
        
    inf = temp.name
    statement = '''python %(scriptsdir)s/csv2db.py -t %(tablename)s --log=%(outfile)s.log --header=transcript_id,gene_id,class < %(inf)s > %(outfile)s'''
    P.run()

##########################################################################
##########################################################################
##########################################################################
@merge([buildCodingGeneSet, classifyLncRNA], "gtfs/transcripts.gtf.gz")
def buildFullGeneSet(infiles, outfile):
    '''
    produces a final gene set that can be used for 
    differential expression analysis and comparisons
    between protein coding and lncRNA transcripts
    '''
    # change the source to be in keeping with classification 
    # of trnscripts

    inf = " ".join(infiles)
    statement = '''zcat %(inf)s | sed "s/Cufflinks/protein_coding/g" 
                   | python %(scriptsdir)s/gtf2gtf.py --sort=position --log=%(outfile)s.log | gzip  > %(outfile)s''' 
    P.run()

##########################################################################
##########################################################################
##########################################################################
# targets
@follows(buildCodingGeneSet
         , buildRefnoncodingGeneSet
         , buildFinalLncRNAGeneSet
         , loadGeneSetStats
         , classifyLncRNA
         , buildFullGeneSet)
def GeneSets():
    pass


@follows(loadCPCResults)
def CodingPotential():
    pass
##########################################################################
##########################################################################
##########################################################################

@follows(GeneSets
         ,CodingPotential)
def full():
    pass

##########################################################################
##########################################################################
##########################################################################
@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

@follows( mkdir( "report" ) )
def update_report():
    '''update report.'''

    E.info( "updating documentation" )
    P.run_report( clean = False )

###################################################################
###################################################################
###################################################################
@follows(update_report)
def publish():
    '''publish files.'''
    # publish web pages
    P.publish_report()

    # publish additional data - i.e. the final lncRNA gtf file
    web_dir = PARAMS["web_dir"]

    if not os.path.exists(os.path.join(web_dir), "lncrna_final.class.gtf.gz"):
        os.symlink("lncrna_final.class.gtf.gz", os.path.abspath(os.path.join(os.path.join(web_dir), "lncrna_final.class.gtf.gz")))

##########################################################################
##########################################################################
##########################################################################
if __name__== "__main__":
    sys.exit( P.main(sys.argv) )




"""
=====================
Read mapping pipeline
=====================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The read mapping pipeline imports unmapped reads from one or more
NGS experiments and maps reads against a reference genome.

This pipeline works on a single genome.

Overview
========

The pipeline implements various mappers and QC plots. It can be used for

* Mapping against a genome
* Mapping RNASEQ data against a genome 
* Mapping against a transcriptome

Principal targets
-----------------

mapping
    perform all mappings

qc
    perform all QC steps

full
    compute all mappings and QC

Optional targets
----------------

merge
    merge mapped :term:`bam` formatted files, for example if reads from different lanes were mapped
    separately. After merging, the ``qc`` target can be run again to get qc stats
    for the merged :term:`bam` formatted files.


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

Reads are imported by placing files are linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   filename.<suffix>

The ``suffix`` determines the file type. The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format. The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files. Thus it might be
   difficult to mix different formats.

Optional inputs
+++++++++++++++

Requirements
-------------

The pipeline requires the results from :doc:`pipeline_annotations`. Set the configuration variable 
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|bowtie_             |>=0.12.7           |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|tophat_             |>=1.4.0            |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|gsnap_              |>=2012.07.20       |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|samtools            |>=0.1.16           |bam/sam files                                   |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |                   |working with intervals                          |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.42             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+
|star_               |>=2.2.0c           |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|bamstats_           |>=1.22             |from CGR, Liverpool                             |
+--------------------+-------------------+------------------------------------------------+

Merging bam files
-----------------

The pipeline has the ability to merge data post-mapping. This is useful if data have
been split over several lanes and have been provide as separate fastq files.

To enable merging, set regular expression for the input and output in the [merge] section 
of the configuration file.

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_mapping.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_mapping.tgz
   tar -xvzf pipeline_mapping.tgz
   cd pipeline_mapping
   python <srcdir>/pipeline_mapping.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

Glossary
========

.. glossary::

   tophat
      tophat_ - a read mapper to detect splice-junctions

   bowtie
      bowtie_ - a read mapper

   star
      star_ - a read mapper for RNASEQ data

   
.. _tophat: http://tophat.cbcb.umd.edu/
.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _gsnap: http://research-pub.gene.com/gmap/
.. _bamstats: http://www.agf.liv.ac.uk/454/sabkea/samStats_13-01-2011
.. _star: http://code.google.com/p/rna-star/

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
import CGATPipelines.PipelinePublishing as PipelinePublishing

import CGAT.Stats as Stats

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
     "pipeline.ini" ],
    defaults = {
        'annotations_dir' : "",
        'paired_end' : False } )

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )



###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

# collect sra nd fastq.gz tracks
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob( "*.sra" ), "(\S+).sra" ) +\
    PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
        glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" ) +\
        PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
            glob.glob( "*.fastq.1.gz" ), "(\S+).fastq.1.gz" ) +\
            PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
                glob.glob( "*.csfasta.gz" ), "(\S+).csfasta.gz" )

###################################################################
## Global flags
###################################################################
MAPPERS = P.asList( PARAMS["mappers" ] )
SPLICED_MAPPING = "tophat" in MAPPERS or "gsnap" in MAPPERS or "star" in MAPPERS

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

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("pipeline_conf.py"): 
    L.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

#########################################################################
#########################################################################
#########################################################################
def mergeAndFilterGTF( infile, outfile, logfile ):
    '''sanitize transcripts file for cufflinks analysis.

    Merge exons separated by small introns (< 5bp).

    Transcript will be ignored that
       * have very long introns (max_intron_size) (otherwise, cufflinks complains)
       * are located on contigs to be ignored (usually: chrM, _random, ...)

    This method preserves all features in a gtf file (exon, CDS, ...).

    returns a dictionary of all gene_ids that have been kept.
    '''

    max_intron_size =  PARAMS["max_intron_size"]

    c = E.Counter()

    outf = gzip.open( outfile, "w" )

    E.info( "filtering by contig and removing long introns" )    
    contigs = set(IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"]) ).getContigs())

    rx_contigs = None
    if "geneset_remove_contigs" in PARAMS:
        rx_contigs = re.compile( PARAMS["geneset_remove_contigs"] )
        E.info( "removing contigs %s" % PARAMS["geneset_remove_contigs"] )

    rna_index = None
    if "geneset_remove_repetetive_rna" in PARAMS:
        rna_file = os.path.join( PARAMS["annotations_dir"],
                                 PARAMS_ANNOTATIONS["interface_rna_gff"] )
        if not os.path.exists( rna_file ):
            E.warn( "file '%s' to remove repetetive rna does not exist" % rna_file )
        else:
            rna_index = GTF.readAndIndex( GTF.iterator( IOTools.openFile( rna_file, "r" ) ) )
            E.info( "removing ribosomal RNA in %s" % rna_file )
    
    gene_ids = {}

    logf = IOTools.openFile( logfile, "w" )
    logf.write( "gene_id\ttranscript_id\treason\n" )

    for all_exons in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( infile )) ):

        c.input += 1
        
        e = all_exons[0]
        # filtering 
        if e.contig not in contigs:
            c.missing_contig += 1
            logf.write( "\t".join( (e.gene_id, e.transcript_id, "missing_contig" )) + "\n" )
            continue

        if rx_contigs and rx_contigs.search(e.contig):
            c.remove_contig += 1
            logf.write( "\t".join( (e.gene_id, e.transcript_id, "remove_contig" )) + "\n" )
            continue

        if rna_index and all_exons[0].source != 'protein_coding':
            found = False
            for exon in all_exons:
                if rna_index.contains( e.contig, e.start, e.end ):
                    found = True
                    break
            if found:
                logf.write( "\t".join( (e.gene_id, e.transcript_id, "overlap_rna" )) + "\n" )
                c.overlap_rna += 1
                continue
                
        is_ok = True

        # keep exons and cds separate by grouping by feature
        all_exons.sort( key = lambda x: x.feature )
        new_exons = []

        for feature, exons in itertools.groupby( all_exons, lambda x: x.feature ):

            tmp = sorted( list(exons) , key = lambda x: x.start )

            gene_ids[tmp[0].transcript_id] = tmp[0].gene_id
            
            l, n = tmp[0], []

            for e in tmp[1:]:
                d = e.start - l.end
                if d > max_intron_size:
                    is_ok = False
                    break
                elif d < 5:
                    l.end = max(e.end, l.end)
                    c.merged += 1
                    continue

                n.append( l )
                l = e

            n.append( l )
            new_exons.extend( n )

            if not is_ok: break

        if not is_ok: 
            logf.write( "\t".join( (e.gene_id, e.transcript_id, "bad_transcript" )) + "\n" )
            c.skipped += 1
            continue

        new_exons.sort( key = lambda x: x.start )

        for e in new_exons:
            outf.write( "%s\n" % str(e) )
            c.exons += 1

        c.output += 1


    outf.close()

    L.info( "%s" % str(c) )
    
    return gene_ids

#########################################################################
#########################################################################
#########################################################################
@active_if( SPLICED_MAPPING )
@merge( os.path.join( PARAMS["annotations_dir"], 
                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]),
        "reference.gtf.gz" )
def buildReferenceGeneSet( infile, outfile ):
    '''sanitize ENSEMBL transcripts file for cufflinks analysis.

    Merge exons separated by small introns (< 5bp).

    Removes unwanted contigs according to configuration
    value ``geneset_remove_contigs``.

    Removes transcripts overlapping ribosomal genes if 
    ``geneset_remove_repetitive_rna`` is set. Protein
    coding transcripts are not removed.

    Transcript will be ignored that
       * have very long introns (max_intron_size) (otherwise, cufflinks complains)
       * are located on contigs to be ignored (usually: chrM, _random, ...)
       
    The result is run through cuffdiff in order to add the p_id and tss_id tags
    required by cuffdiff. 

    This will only keep sources of the type 'exon'. It will also remove
    any transcripts not in the reference genome.

    Cuffdiff requires overlapping genes to have different tss_id tags.

    This gene is the source for most other gene sets in the pipeline.
    '''

    tmpfilename = P.getTempFilename( "." )
    tmpfilename2 = P.getTempFilename( "." )
    tmpfilename3 = P.getTempFilename( "." )

    gene_ids = mergeAndFilterGTF( infile, tmpfilename, "%s.removed.gz" % outfile )
    
    #################################################
    E.info( "adding tss_id and p_id" )

    # The p_id attribute is set if the fasta sequence is given.
    # However, there might be some errors in cuffdiff downstream:
    #
    # cuffdiff: bundles.cpp:479: static void HitBundle::combine(const std::vector<HitBundle*, std::allocator<HitBundle*> >&, HitBundle&): Assertion `in_bundles[i]->ref_id() == in_bundles[i-1]->ref_id()' failed.
    #
    # I was not able to resolve this, it was a complex
    # bug dependent on both the read libraries and the input reference gtf files
    statement = '''
    cuffcompare -r <( gunzip < %(tmpfilename)s )
         -T 
         -s %(bowtie_index_dir)s/%(genome)s.fa
         -o %(tmpfilename2)s
         <( gunzip < %(tmpfilename)s )
         <( gunzip < %(tmpfilename)s )
    > %(outfile)s.log
    '''
    P.run()

    #################################################
    E.info( "resetting gene_id and transcript_id" )

    # reset gene_id and transcript_id to ENSEMBL ids
    # cufflinks patch: 
    # make tss_id and p_id unique for each gene id
    outf = IOTools.openFile( tmpfilename3, "w")
    map_tss2gene, map_pid2gene = {}, {}
    inf = IOTools.openFile( tmpfilename2 + ".combined.gtf" )

    def _map( gtf, key, val, m ):
        if val in m:
            while gene_id != m[val]:
                val += "a" 
                if val not in m: break
        m[val] = gene_id

        gtf.setAttribute( key, val )
        
    for gtf in GTF.iterator( inf ):
        transcript_id = gtf.oId
        gene_id = gene_ids[transcript_id] 
        gtf.setAttribute( "transcript_id", transcript_id )
        gtf.setAttribute( "gene_id", gene_id )
        
        # set tss_id
        try: tss_id = gtf.tss_id
        except AttributeError: tss_id = None
        try: p_id = gtf.p_id
        except AttributeError: p_id = None
        
        if tss_id: _map( gtf, "tss_id", tss_id, map_tss2gene)
        if p_id: _map( gtf, "p_id", p_id, map_pid2gene)

        outf.write( str(gtf) + "\n" )
        
    outf.close()

    # sort gtf file
    PipelineGeneset.sortGTF( tmpfilename3, outfile )

    os.unlink( tmpfilename )
    # make sure tmpfilename2 is NEVER empty
    assert tmpfilename2
    for x in glob.glob( tmpfilename2 + "*" ): os.unlink( x )
    os.unlink( tmpfilename3 )

#########################################################################
#########################################################################
#########################################################################
@active_if( SPLICED_MAPPING )
@transform( buildReferenceGeneSet, 
            suffix("reference.gtf.gz"),
            "refcoding.gtf.gz" )
def buildCodingGeneSet( infile, outfile ):
    '''build a gene set with only protein coding 
    transcripts.

    Genes are selected via their gene biotype in the GTF file.
    Note that this set will contain all transcripts of protein
    coding genes, including processed transcripts.

    This set includes UTR and CDS.
    '''
    
    to_cluster = True
    statement = '''
    zcat %(infile)s | awk '$2 == "protein_coding"' | gzip > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@active_if( SPLICED_MAPPING )
@merge( os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_geneset_flat_gtf"]),
        "introns.gtf.gz" )
def buildIntronGeneModels(infile, outfile ):
    '''build protein-coding intron-transcipts.

    Intron-transcripts are the reverse complement of transcripts.

    Only protein coding genes are taken.

    10 bp are truncated on either end of an intron and need
    to have a minimum length of 100.

    Introns from nested genes might overlap, but all exons
    are removed.
    '''

    to_cluster = True
    
    filename_exons = os.path.join( PARAMS["annotations_dir"],
                                   PARAMS_ANNOTATIONS["interface_geneset_exons_gtf"] )

    statement = '''gunzip
        < %(infile)s 
        | awk '$2 == "protein_coding"'
	| python %(scriptsdir)s/gtf2gtf.py --sort=gene 
	| python %(scriptsdir)s/gtf2gtf.py 
               --exons2introns 
               --intron-min-length=100 
               --intron-border=10 
               --log=%(outfile)s.log
        | python %(scriptsdir)s/gff2gff.py
               --crop=%(filename_exons)s
               --log=%(outfile)s.log
	| python %(scriptsdir)s/gtf2gtf.py 
              --set-transcript-to-gene 
              --log=%(outfile)s.log 
	| perl -p -e 's/intron/exon/'
        | gzip
        > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@active_if( SPLICED_MAPPING )
@transform( buildCodingGeneSet,
            suffix(".gtf.gz"),
            "_transcript2gene.load" )
def loadGeneInformation( infile, outfile ):
    PipelineGeneset.loadTranscript2Gene( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@active_if( SPLICED_MAPPING )
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
@active_if( SPLICED_MAPPING )
@transform( buildCodingGeneSet, suffix(".gtf.gz"), ".fa")
def buildReferenceTranscriptome( infile, outfile ):
    '''build reference transcriptome. 

    The reference transcriptome contains all known 
    protein coding transcripts.

    The sequences include both UTR and CDS.

    '''
    to_cluster = True
    gtf_file = P.snip(infile, ".gz") 

    genome_file = os.path.abspath( os.path.join( PARAMS["bowtie_index_dir"], PARAMS["genome"] + ".fa" ) )

    statement = '''
    zcat %(infile)s
    | awk '$3 == "exon"' > %(gtf_file)s;
    gtf_to_fasta %(gtf_file)s %(genome_file)s %(outfile)s;
    checkpoint; 
    samtools faidx %(outfile)s
    ''' 
    P.run()
    
    os.symlink(gtf_file, P.snip(gtf_file,".gtf") + ".gff")
        
    prefix = P.snip( outfile, ".fa" )

    # build raw index
    statement = '''
    bowtie-build -f %(outfile)s %(prefix)s >> %(outfile)s.log 2>&1
    '''

    P.run()

    # build color space index
    statement = '''
    bowtie-build -C -f %(outfile)s %(prefix)s_cs >> %(outfile)s.log 2>&1
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@active_if( SPLICED_MAPPING )
@transform( buildCodingGeneSet, suffix(".gtf.gz"), ".junctions")
def buildJunctions( infile, outfile ):
    '''build file with splice junctions from gtf file.
    
    A junctions file is a better option than supplying a GTF
    file, as parsing the latter often fails. See:

    http://seqanswers.com/forums/showthread.php?t=7563

    '''
    
    outf = IOTools.openFile( outfile, "w" )
    njunctions = 0
    for gffs in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( infile, "r" ) )):
        
        gffs.sort( key = lambda x: x.start )
        end = gffs[0].end
        for gff in gffs[1:]:
            # subtract one: these are not open/closed coordinates but the 0-based coordinates
            # of first and last residue that are to be kept (i.e., within the exon).
            outf.write( "%s\t%i\t%i\t%s\n" % (gff.contig, end - 1, gff.start, gff.strand ) )
            end = gff.end
            njunctions += 1
            
    outf.close()

    if njunctions == 0:
        E.warn( 'no junctions found in gene set' )
        return
    else:
        E.info( 'found %i junctions before removing duplicates' % njunctions )


    # make unique
    statement = '''mv %(outfile)s %(outfile)s.tmp; 
                   cat < %(outfile)s.tmp | sort | uniq > %(outfile)s;
                   rm -f %(outfile)s.tmp; '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@active_if( SPLICED_MAPPING )
@transform( buildCodingGeneSet, suffix(".gtf.gz"), ".splicesites.iit")
def buildGSNAPSpliceSites( infile, outfile ):
    '''build file with know splice sites for GSNAP.
    '''
    
    outfile = P.snip( outfile, ".iit" )
    statement = '''
    zcat %(infile)s | gtf_splicesites | iit_store -o %(outfile)s > %(outfile)s.log'''
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Read mapping
#########################################################################
SEQUENCEFILES=("*.fastq.1.gz", 
               "*.fastq.gz",
               "*.sra",
               "*.export.txt.gz",
               "*.csfasta.gz",
               "*.csfasta.F3.gz",
               )
SEQUENCEFILES_REGEX=regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|csfasta.F3.gz|export.txt.gz)")

###################################################################
###################################################################
###################################################################
## load number of reads
###################################################################
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            r"\1.nreads" )
def countReads( infile, outfile ):
    '''count number of reads in input files.'''
    to_cluster = True
    m = PipelineMapping.Counter()
    statement = m.build( (infile,), outfile ) 
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Map reads with tophat
#########################################################################
@active_if( SPLICED_MAPPING )
@follows( mkdir("tophat.dir" ) )
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            add_inputs( buildJunctions, buildReferenceTranscriptome ), 
            r"tophat.dir/\1.tophat.bam" )
def mapReadsWithTophat( infiles, outfile ):
    '''map reads from .fastq or .sra files.

    A list with known splice junctions is supplied.

    If tophat fails with an error such as::

       Error: segment-based junction search failed with err =-6
       what():  std::bad_alloc
 
    it means that it ran out of memory.

    '''
    job_options= "-pe dedicated %i -R y" % PARAMS["tophat_threads"]

    if "--butterfly-search" in PARAMS["tophat_options"]:
        # for butterfly search - require insane amount of
        # RAM.
        job_options += " -l mem_free=50G"
    else:
        job_options += " -l mem_free=%s" % PARAMS["tophat_memory"]

    to_cluster = True
    m = PipelineMapping.Tophat( executable = P.substituteParameters( **locals() )["tophat_executable"],
                                strip_sequence = PARAMS["strip_sequence"] )
    infile, reffile, transcriptfile = infiles
    tophat_options = PARAMS["tophat_options"] + " --raw-juncs %(reffile)s " % locals()
    
    # Nick - added the option to map to the reference transcriptome first (built within the pipeline)
    if PARAMS["tophat_include_reference_transcriptome"]:
        prefix = os.path.abspath( P.snip( transcriptfile, ".fa" ) )
        tophat_options = tophat_options + " --transcriptome-index=%s -n 2" % prefix

    statement = m.build( (infile,), outfile ) 
    P.run()

@active_if( SPLICED_MAPPING )
@follows( mkdir("tophat2.dir" ) )
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            add_inputs( buildJunctions, buildReferenceTranscriptome ), 
            r"tophat2.dir/\1.tophat2.bam" )
def mapReadsWithTophat2( infiles, outfile ):
    '''map reads from .fastq or .sra files.

    A list with known splice junctions is supplied.

    If tophat fails with an error such as::

       Error: segment-based junction search failed with err =-6
       what():  std::bad_alloc
 
    it means that it ran out of memory.

    '''
    job_options= "-pe dedicated %i -R y" % PARAMS["tophat2_threads"]

    if "--butterfly-search" in PARAMS["tophat2_options"]:
        # for butterfly search - require insane amount of
        # RAM.
        job_options += " -l mem_free=50G"
    else:
        job_options += " -l mem_free=%s" % PARAMS["tophat2_memory"]

    to_cluster = True
    m = PipelineMapping.Tophat2( executable = P.substituteParameters( **locals() )["tophat2_executable"],
                                 strip_sequence = PARAMS["strip_sequence"] )
    infile, reffile, transcriptfile = infiles
    tophat_options = PARAMS["tophat2_options"] + " --raw-juncs %(reffile)s " % locals()
    
    # Nick - added the option to map to the reference transcriptome first (built within the pipeline)
    if PARAMS["tophat2_include_reference_transcriptome"]:
        prefix = os.path.abspath( P.snip( transcriptfile, ".fa" ) )
        tophat_options = tophat_options + " --transcriptome-index=%s -n 2" % prefix

    statement = m.build( (infile,), outfile ) 
    P.run()

############################################################
############################################################
############################################################
@active_if( SPLICED_MAPPING )
@merge( mapReadsWithTophat, "tophat_stats.tsv" )
def buildTophatStats( infiles, outfile ):

    def _select( lines, pattern ):
        x = re.compile(pattern)
        for line in lines:
            r = x.search( line )
            if r: 
                g = r.groups()
                if len(g) > 1: return g
                else: return g[0]

        raise ValueError( "pattern '%s' not found %s" % (pattern, lines ))

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "\t".join( ("track",
                            "reads_in",
                            "reads_removed", 
                            "reads_out", 
                            "junctions_loaded", 
                            "junctions_found", 
                            "possible_splices" ) ) + "\n" )
    
    for infile in infiles:
        
        track = P.snip( infile, ".bam" )
        indir = infile + ".logs" 

        fn = os.path.join( indir, "prep_reads.log" )
        lines = open( fn ).readlines()
        reads_removed, reads_in = map(int, _select( lines, "(\d+) out of (\d+) reads have been filtered out" ) )
        reads_out = reads_in - reads_removed
        prep_reads_version = _select( lines, "prep_reads (.*)$" )
        
        fn = os.path.join( indir, "reports.log" )
        lines = open( fn ).readlines()
        tophat_reports_version = _select( lines, "tophat_reports (.*)$" )
        junctions_loaded = int( _select( lines, "Loaded (\d+) junctions") )
        junctions_found = int( _select( lines, "Found (\d+) junctions from happy spliced reads") )

        fn = os.path.join( indir, "segment_juncs.log" )
        
        
        if os.path.exists(fn):
            lines = open( fn ).readlines()
            if len(lines) > 0:
                segment_juncs_version =  _select( lines, "segment_juncs (.*)$" )
                possible_splices = int( _select( lines, "Reported (\d+) total possible splices") )
            else:
                segment_juncs_version = "na"
                possible_splices = ""
        else:
            segment_juncs_version = "na"
            possible_splices = ""

        # fix for paired end reads - tophat reports pairs, not reads
        if PARAMS["paired_end"]:
            reads_in *= 2
            reads_out *= 2
            reads_removed *= 2

        outf.write( "\t".join( map(str, (track,
                                         reads_in, reads_removed, reads_out, 
                                         junctions_loaded, junctions_found, possible_splices ) ) ) + "\n" )

    outf.close()

############################################################
############################################################
############################################################
@jobs_limit( 1, "db" )
@active_if( SPLICED_MAPPING )
@transform( buildTophatStats, suffix(".tsv"), ".load" )
def loadTophatStats( infile, outfile ):
    P.load( infile, outfile )

############################################################
############################################################
############################################################
@active_if( SPLICED_MAPPING )
@follows( mkdir("gsnap.dir" ) )
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            add_inputs( buildGSNAPSpliceSites ), 
            r"gsnap.dir/\1.gsnap.bam" )
def mapReadsWithGSNAP( infiles, outfile ):
    '''map reads from .fastq or .sra files.

    '''

    infile, infile_splices = infiles
    job_options= "-pe dedicated %i -R y -l mem_free=%s" % (PARAMS["gsnap_threads"],
                                                           PARAMS["gsnap_memory"])
    to_cluster = True
    m = PipelineMapping.GSNAP( executable = P.substituteParameters( **locals() )["gsnap_executable"],
                               strip_sequence = PARAMS["strip_sequence"] )
    
    if PARAMS["gsnap_include_known_splice_sites"]:
        gsnap_options = PARAMS["gsnap_options"] + " --use-splicing=%(infile_splices)s " % locals()

    statement = m.build( (infile,), outfile ) 
    P.run()

############################################################
############################################################
############################################################
@active_if( SPLICED_MAPPING )
@follows( mkdir("star.dir" ) )
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            r"star.dir/\1.star.bam" )
def mapReadsWithSTAR( infile, outfile ):
    '''map reads from .fastq or .sra files.

    '''

    job_options= "-pe dedicated %i -R y -l mem_free=%s" % (PARAMS["star_threads"],
                                                           PARAMS["star_memory"])
    to_cluster = True
    
    star_mapping_genome = PARAMS["star_genome"] or PARAMS["genome"]
    
    m = PipelineMapping.STAR( executable = P.substituteParameters( **locals() )["star_executable"],
                              strip_sequence = PARAMS["strip_sequence"] )
    
    statement = m.build( (infile,), outfile ) 
    P.run()

############################################################
############################################################
############################################################
@active_if( SPLICED_MAPPING )
@merge( mapReadsWithSTAR, "star_stats.tsv" )
def buildSTARStats( infiles, outfile ):
    '''load stats from STAR run.'''

    data = collections.defaultdict( list )
    for infile in infiles:
        fn = infile + ".final.log"
        if not os.path.exists ( fn ):
            raise ValueError( "incomplete run: %s" % infile )
        
        for line in IOTools.openFile( fn ):
            if not "|" in line: continue
            header, value = line.split("|")
            header = re.sub( "%", "percent", header )
            data[header.strip()].append( value.strip() )
    
    keys = data.keys()
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "track\t%s\n" % "\t".join(keys) )
    for x, infile in enumerate(infiles):
        track = P.snip( os.path.basename( infile ), ".bam" )
        outf.write( "%s\t%s\n" % (track, "\t".join( [ data[key][x] for key in keys ] ) ) )
    outf.close()

############################################################
############################################################
############################################################
@jobs_limit( 1, "db" )
@active_if( SPLICED_MAPPING )
@transform( buildSTARStats, suffix(".tsv"), ".load")
def loadSTARStats( infile, outfile ):
    '''load stats from STAR run.'''
    P.load( infile, outfile )
    
############################################################
############################################################
############################################################
@active_if( SPLICED_MAPPING )
@follows( mkdir("transcriptome.dir" ) )
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            add_inputs( buildReferenceTranscriptome ), 
            r"transcriptome.dir/\1.trans.bam" )
def mapReadsWithBowtieAgainstTranscriptome( infiles, outfile ):
    '''map reads using bowtie against transcriptome data.
    '''

    # Mapping will permit up to one mismatches. This is sufficient
    # as the downstream filter in rnaseq_bams2bam requires the
    # number of mismatches less than the genomic number of mismatches.
    # Change this, if the number of permitted mismatches for the genome
    # increases.

    # Output all valid matches in the best stratum. This will 
    # inflate the file sizes due to matches to alternative transcripts
    # but otherwise matches to paralogs will be missed (and such
    # reads would be filtered out).
    job_options= "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
    to_cluster = True
    m = PipelineMapping.BowtieTranscripts( executable = P.substituteParameters( **locals() )["bowtie_executable"],
                                           strip_sequence = PARAMS["strip_sequence"] )
    infile, reffile = infiles
    prefix = P.snip( reffile, ".fa" )
    #IMS: moved reporting options to ini
    #bowtie_options = "%s --best --strata -a" % PARAMS["bowtie_transcriptome_options"] 
    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################
###################################################################
###################################################################
## Map reads with bowtie
###################################################################
@follows( mkdir("bowtie.dir" ) )
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            add_inputs( os.path.join( PARAMS["bowtie_index_dir"], PARAMS["genome"] + ".fa") ),
            r"bowtie.dir/\1.bowtie.bam" )
def mapReadsWithBowtie( infiles, outfile ):
    '''map reads with bowtie'''

    job_options= "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
    to_cluster = True
    m = PipelineMapping.Bowtie( executable = P.substituteParameters( **locals() )["bowtie_executable"],
                                strip_sequence = PARAMS["strip_sequence"] )
    infile, reffile = infiles
    #IMS remove reporting options to the ini
    #bowtie_options = "%s --best --strata -a" % PARAMS["bowtie_options"] 
    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################
###################################################################
###################################################################
## Map reads with bwa
###################################################################
@follows( mkdir("bwa.dir") )
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            r"bwa.dir/\1.bwa.bam" )
def mapReadsWithBWA( infile, outfile ):
    '''map reads with bwa'''


    job_options= "-pe dedicated %i -R y -l mem_free=%s" % (PARAMS["bwa_threads"],
                                                           PARAMS["bwa_memory"] )
    to_cluster = True
    m = PipelineMapping.BWA( remove_non_unique = PARAMS["remove_non_unique"],
                             strip_sequence = PARAMS["strip_sequence"] )

    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################
###################################################################
###################################################################
## Map reads with bwa
###################################################################
@follows( mkdir("stampy.dir") )
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            r"stampy.dir/\1.stampy.bam" )
def mapReadsWithStampy( infile, outfile ):
    '''map reads with stampy'''

    job_options= "-pe dedicated %i -R y -l mem_free=%s" % ( PARAMS["stampy_threads"],
                                                            PARAMS["stampy_memory"] )
    to_cluster = True
    m = PipelineMapping.Stampy( strip_sequence = PARAMS["strip_sequence"] )
    statement = m.build( (infile,), outfile ) 
    P.run()

MAPPINGTARGETS = []
mapToMappingTargets = { 'tophat': (mapReadsWithTophat, loadTophatStats),
                        'tophat2': (mapReadsWithTophat2,),
                        'bowtie': (mapReadsWithBowtie,),
                        'bwa': (mapReadsWithBWA,),
                        'stampy': (mapReadsWithStampy,),
                        'transcriptome': (mapReadsWithBowtieAgainstTranscriptome,),
                        'gsnap' : (mapReadsWithGSNAP,),
                        'star' : (mapReadsWithSTAR,loadSTARStats),
                        }

for x in P.asList( PARAMS["mappers"]):
    MAPPINGTARGETS.extend( mapToMappingTargets[x] )

@follows( *MAPPINGTARGETS )
def mapping(): pass

###################################################################
###################################################################
###################################################################
if "merge_pattern_input" in PARAMS and PARAMS["merge_pattern_input"]:
    if "merge_pattern_output" not in PARAMS or not PARAMS["merge_pattern_output"]:
        raise ValueError("no output pattern 'merge_pattern_output' specificied")
    @collate( MAPPINGTARGETS, 
              regex( "%s\.([^.]+).bam" % PARAMS["merge_pattern_input"].strip() ),
              # the last expression counts number of groups in pattern_input
              r"%s.\%i.bam" % (PARAMS["merge_pattern_output"].strip(), 
                               PARAMS["merge_pattern_input"].count("(")+1),
              )
    def mergeBAMFiles( infiles, outfile ):
        '''merge BAM files from the same experiment.'''
        if len(infiles) == 1:
            E.info("%(outfile)s: only one file for merging - creating softlink" % locals() )
            P.clone( infiles[0], outfile )
            P.clone( infiles[0] + ".bai", outfile + ".bai" )
            return
        
        to_cluster = True

        infiles = " ".join(infiles)
        statement = '''
        samtools merge %(outfile)s %(infiles)s >& %(outfile)s.log;
        checkpoint;
        samtools index %(outfile)s
        '''
        P.run()

    # add to bam files produced
    MAPPINGTARGETS.append( mergeBAMFiles )

    @collate( countReads,
              regex( "%s.nreads" % PARAMS["merge_pattern_input"] ),
              r"%s.nreads" % PARAMS["merge_pattern_output"],
              )
    def mergeReadCounts( infiles, outfile ):
        '''merge BAM files from the same experiment.'''
        
        to_cluster = True

        nreads = 0
        for infile in infiles:
            with IOTools.openFile( infile, "r" ) as inf:
                for line in infiles:
                    if not line.startswith("nreads"): continue
                    nreads += int(line[:-1].split("\t")[1])
                    
        outf = IOTools.openFile( outfile, "w")
        outf.write("nreads\t%i\n" % nreads )
        outf.close()

else:
    @follows( countReads )
    def mergeReadCounts(): pass

###################################################################
###################################################################
###################################################################
## QC targets
###################################################################

############################################################
############################################################
############################################################
@active_if( SPLICED_MAPPING )
@transform( MAPPINGTARGETS,
            suffix(".bam" ),
            ".picard_inserts")
def buildPicardTranscriptomeInsertSize( infiles, outfile ):
    '''build alignment stats using picard.

    Note that picards counts reads but they are in fact alignments.
    '''
    infile, reffile = infiles

    PipelineMappingQC.buildPicardAlignmentStats( infile, 
                                                 outfile,
                                                 reffile )

############################################################
############################################################
############################################################
@transform( MAPPINGTARGETS,
            suffix(".bam" ), 
            add_inputs( os.path.join( PARAMS["bowtie_index_dir"],
                                PARAMS["genome"] + ".fa" ) ), 
            ".picard_stats")
def buildPicardStats( infiles, outfile ):
    '''build alignment stats using picard.

    Note that picards counts reads but they are in fact alignments.
    '''
    infile, reffile = infiles

    # patch for mapping against transcriptome - switch genomic reference
    # to transcriptomic sequences
    if "transcriptome.dir" in infile:
        reffile = "refcoding.fa"

    PipelineMappingQC.buildPicardAlignmentStats( infile, 
                                                 outfile,
                                                 reffile )

############################################################
############################################################
############################################################
@jobs_limit( 1, "db" )
@merge( buildPicardStats, "picard_stats.load" )
def loadPicardStats( infiles, outfile ):
    '''merge alignment stats into single tables.'''

    PipelineMappingQC.loadPicardAlignmentStats( infiles, outfile )

############################################################
############################################################
############################################################
@transform( MAPPINGTARGETS,
            suffix(".bam" ), 
            ".picard_stats.duplication_metrics")
def buildPicardDuplicationStats( infile, outfile ):
    '''Get duplicate stats from picard MarkDuplicates.
    Pair duplication is properly handled, including inter-chromosomal cases. SE data is also handled.
    These stats also contain a histogram that estimates the return from additional sequecing.
    No marked bam files are retained (/dev/null...)
    Note that picards counts reads but they are in fact alignments.
    '''
    PipelineMappingQC.buildPicardDuplicationStats( infile,outfile )

############################################################
############################################################
############################################################
@jobs_limit( 1, "db" )
@merge( buildPicardDuplicationStats, "picard_duplication_stats.load" )
def loadPicardDuplicationStats( infiles, outfile ):
    '''merge alignment stats into single tables.'''
    #separate load function while testing
    PipelineMappingQC.loadPicardDuplicationStats( infiles, outfile )

# ############################################################
# ############################################################
# ############################################################
# @merge( buildBAMs, "mapping_stats.load" )
# def loadMappingStats( infiles, outfile ):

#     header = ",".join( [P.snip( x, ".bam") for x in infiles] )
#     filenames = " ".join( [ "%s.tsv" % x for x in infiles ] )
#     tablename = P.toTable( outfile )

#     statement = """python %(scriptsdir)s/combine_tables.py
#                       --headers=%(header)s
#                       --missing=0
#                       --ignore-empty
#                    %(filenames)s
#                 | perl -p -e "s/bin/track/"
#                 | perl -p -e "s/unique/unique_alignments/"
#                 | python %(scriptsdir)s/table2table.py --transpose
#                 | python %(scriptsdir)s/csv2db.py
#                       --index=track
#                       --table=%(tablename)s 
#                 > %(outfile)s
#             """
#     P.run()


############################################################
############################################################
############################################################
@follows( countReads, mergeReadCounts )
@transform( MAPPINGTARGETS,
            regex("(.*)/(.*)\.(.*).bam"),
            add_inputs( r"\2.nreads" ),
            r"\1/\2.\3.readstats" )
def buildBAMStats( infiles, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = True

    rna_file = os.path.join( PARAMS["annotations_dir"],
                             PARAMS_ANNOTATIONS["interface_rna_gff"] )

    job_options = "-l mem_free=8G"

    bamfile, readsfile = infiles

    nreads = PipelineMappingQC.getNumReadsFromReadsFile( readsfile )
    track = P.snip( readsfile, ".nreads" )

    # if a fastq file exists, submit for counting
    if os.path.exists( track + ".fastq.gz" ):
        fastqfile = track + ".fastq.gz"
    elif os.path.exists( track + ".fastq.1.gz" ):
        fastqfile = track + ".fastq.1.gz"
    else:
        fastqfile = None

    if fastqfile != None:
        fastq_option = "--filename-fastq=%s" % fastqfile
    else:
        fastq_option = ""

    statement = '''python
    %(scriptsdir)s/bam2stats.py
         %(fastq_option)s
         --force
         --filename-rna=%(rna_file)s
         --remove-rna
         --input-reads=%(nreads)i
         --output-filename-pattern=%(outfile)s.%%s
    < %(bamfile)s
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
@jobs_limit( 1, "db" )
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''import bam statisticis.'''
    
    PipelineMappingQC.loadBAMStats( infiles, outfile )

############################################################
############################################################
############################################################
@transform( MAPPINGTARGETS,
            suffix(".bam"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_genomic_context_bed"] ) ),
            ".contextstats" )
def buildContextStats( infiles, outfile ):
    '''build mapping context stats.

    Examines the genomic context to where reads align.

    A read is assigned to the genomic context that it
    overlaps by at least 50%. Thus some reads mapping
    several contexts might be dropped.
    '''

    infile, reffile = infiles

    min_overlap = 0.5

    job_options = "-l mem_free=4G"

    to_cluster = True
    statement = '''
       python %(scriptsdir)s/bam_vs_bed.py
              --min-overlap=%(min_overlap)f
              --log=%(outfile)s.log
              %(infile)s %(reffile)s
       > %(outfile)s
       '''

    P.run()

############################################################
############################################################
############################################################
@jobs_limit( 1, "db" )
@follows( loadBAMStats )
@merge( buildContextStats, "context_stats.load" )
def loadContextStats( infiles, outfile ):
    """
    load context mapping statistics."""

    header = ",".join( [os.path.basename( P.snip( x, ".contextstats") ) for x in infiles] )
    filenames = " ".join( infiles  )
    tablename = P.toTable( outfile )

    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --skip-titles
                   %(filenames)s
                | perl -p -e "s/(bin|category)/track/; s/\?/Q/g"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
                """
    P.run()
    
    dbhandle = sqlite3.connect( PARAMS["database"] )
 
# The following is not necessary any more as context stats now also outputs a "total" column   
#    cc = Database.executewait( dbhandle, '''ALTER TABLE %(tablename)s ADD COLUMN mapped INTEGER''' % locals())
#    statement = '''UPDATE %(tablename)s SET mapped = 
#                                       (SELECT b.alignments_mapped FROM bam_stats AS b 
#                                            WHERE %(tablename)s.track = b.track)''' % locals()#
#
#    cc = Database.executewait( dbhandle, statement )
#    dbhandle.commit()

###################################################################
###################################################################
###################################################################
## QC specific to spliced mapping
###################################################################
###################################################################
###################################################################
@active_if( SPLICED_MAPPING )
@transform( MAPPINGTARGETS,
            suffix(".bam"),
            add_inputs( buildCodingExons ),
            ".exon.validation.tsv.gz" )
def buildExonValidation( infiles, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = True
    infile, exons = infiles
    statement = '''cat %(infile)s
    | python %(scriptsdir)s/bam_vs_gtf.py
         --filename-exons=%(exons)s
         --force
         --log=%(outfile)s.log
         --output-filename-pattern="%(outfile)s.%%s.gz"
    | gzip
    > %(outfile)s
    '''

    P.run()


############################################################
############################################################
############################################################
@jobs_limit( 1, "db" )
@active_if( SPLICED_MAPPING )
@merge( buildExonValidation, "exon_validation.load" )
def loadExonValidation( infiles, outfile ):
    '''merge alignment stats into single tables.'''
    suffix = ".exon.validation.tsv.gz" 
    P.mergeAndLoad( infiles, outfile, suffix = suffix )
    for infile in infiles:
        track = P.snip( infile, suffix )
        o = "%s_overrun.load" % track 
        P.load( infile + ".overrun.gz", o )

#########################################################################
#########################################################################
#########################################################################
@active_if(SPLICED_MAPPING)
@split(buildCodingGeneSet,"refcoding.*.gtf.gz")
def splitCodingGeneSetByChr(infile,outfiles):
    '''split coding geneset by chromosome to allow parallel
    read counting '''

    PipelineMapping.splitGeneSet(infile)

@active_if( SPLICED_MAPPING )
@split( MAPPINGTARGETS, 
        regex("(.+).bam"),
        add_inputs( splitCodingGeneSetByChr ),
        r"\1.*.transcript_counts.tsv.gz" )
def buildTranscriptLevelReadCounts( infiles, outfile):
    '''count reads falling into transcripts of protein coding 
       gene models.
       
    .. note::
       In paired-end data sets each mate will be counted. Thus
       the actual read counts are approximately twice the fragment
       counts.
       
    '''
    infile, genesets = infiles[0],infiles[1:]
    
    to_cluster = True
    statements = []
    
    for geneset in genesets:

        chrom = re.match("refcoding\.(.+)\.gtf.gz",geneset).groups()[0]
        bam = P.snip(infile,".bam")
        outfile = "%s.%s.transcript_counts.tsv.gz" % (bam,chrom)

        statement = '''
               zcat %(geneset)s 
             | python %%(scriptsdir)s/gtf2table.py 
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
         ''' % locals()
        statements.append(statement)

    P.run()

#########################################################################
@active_if( SPLICED_MAPPING )
@collate(buildTranscriptLevelReadCounts, 
         regex("(.+)\..+\.transcript_counts.tsv.gz"),
         r"\1.transcript_counts.tsv.gz")
def collateTranscriptCounts(infiles,outfile):
    ''' pull together the transcript counts over each chromosome '''

    infiles = " ".join(infiles)
    statement = '''python %(scriptsdir)s/combine_tables.py
                                      --cat
                                      --log=%(outfile)s.log
                                      %(infiles)s
                  | cut -f1 --complement
                  | gzip
                  > %(outfile)s '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@jobs_limit( 1, "db" )
@active_if( SPLICED_MAPPING )
@transform( collateTranscriptCounts,
            suffix(".tsv.gz"),
            ".load" )
def loadTranscriptLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=transcript_id --allow-empty" )

#########################################################################
#########################################################################
#########################################################################
@active_if( SPLICED_MAPPING )
@transform( MAPPINGTARGETS,
            suffix(".bam"),
            add_inputs( buildIntronGeneModels ),
            ".intron_counts.tsv.gz" )
def buildIntronLevelReadCounts( infiles, outfile ):
    '''compute coverage of exons with reads.
    '''

    infile, exons = infiles

    if "transcriptome.dir" in infile:
        P.touch(outfile)
        return

    to_cluster = True

    statement = '''
    zcat %(exons)s 
    | python %(scriptsdir)s/gtf2table.py 
          --reporter=genes
          --bam-file=%(infile)s 
          --counter=length
          --prefix="introns_"
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
@jobs_limit( 1, "db" )
@active_if( SPLICED_MAPPING )
@transform(buildIntronLevelReadCounts,
           suffix(".tsv.gz"),
           ".load" )
def loadIntronLevelReadCounts( infile, outfile ):
    P.load( infile, outfile, options="--index=gene_id --allow-empty" )

###################################################################
###################################################################
###################################################################
@merge((countReads, mergeReadCounts), "reads_summary.load" )
def loadReadCounts( infiles, outfile ):
    '''load read counts into database.'''

    outf = P.getTempFile( "." )
    outf.write( "track\ttotal_reads\n")
    for infile in infiles:
        track = P.snip(infile, ".nreads")
        lines = IOTools.openFile( infile ).readlines()
        nreads = int( lines[0][:-1].split("\t")[1])
        outf.write( "%s\t%i\n" % (track,nreads))
    outf.close()
        
    P.load( outf.name, outfile )
    
    os.unlink(outf.name)

###################################################################
###################################################################
###################################################################
## various export functions
###################################################################
@transform( MAPPINGTARGETS,
            regex(".bam"),
            ".bw" )
def buildBigWig( infile, outfile ):
    '''build wiggle files from bam files.'''
    to_cluster = True

    # wigToBigWig observed to use 16G
    job_options = "-l mem_free=16G"

    statement = '''python %(scriptsdir)s/bam2wiggle.py 
                         --output-format=bigwig 
                         %(bigwig_options)s
                         %(infile)s 
                         %(outfile)s
                   > %(outfile)s.log''' 
    P.run()

###################################################################
###################################################################
###################################################################
## 
###################################################################
@merge( buildBigWig,
        "bigwig_stats.load" )
def loadBigWigStats( infiles, outfile ):
    '''load bigwig summary for all wiggle files.'''

    to_cluster = True
    
    data = " ".join( [ '<( bigWigInfo %s | perl -p -e "s/:/\\t/; s/ //g; s/,//g")' % x for x in infiles] )
    headers = ",".join( [P.snip(os.path.basename(x), ".bw") for x in infiles] )

    tablename = P.toTable( outfile)

    statement = '''python %(scriptsdir)s/combine_tables.py 
                         --header=%(headers)s
                         --skip-titles
                         --missing=0
                         --ignore-empty
                         %(data)s 
                | perl -p -e "s/bin/track/" | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
    '''

    P.run()

@transform( MAPPINGTARGETS,
            regex(".bam"),
            ".bed.gz" )
def buildBed( infile, outfile ):
    '''build bed files from bam files.'''
    to_cluster = True

    statement = '''
    cat %(infile)s 
    | python %(scriptsdir)s/bam2bed.py
          %(bed_options)s
          --log=%(outfile)s.log
          -
    | sort -k1,1 -k2,2n
    | bgzip 
    > %(outfile)s
    '''
    P.run()

    statement = '''
    tabix -p bed %(outfile)s
    '''
    P.run()


###################################################################
###################################################################
###################################################################
@follows( loadReadCounts, 
          loadPicardStats, 
          loadBAMStats, 
          loadContextStats )
def general_qc(): pass

@active_if( SPLICED_MAPPING )
@follows( loadExonValidation,
          loadGeneInformation,
          loadTranscriptLevelReadCounts,
          loadIntronLevelReadCounts )
def spliced_qc(): pass

@follows( general_qc, spliced_qc )
def qc(): pass

@follows( loadPicardDuplicationStats )
def duplication(): pass

@follows( buildBigWig, loadBigWigStats )
def wig(): pass

###################################################################
###################################################################
###################################################################
## export targets
###################################################################
@merge( (loadBAMStats, loadPicardStats, loadContextStats), "view_mapping.load" )
def createViewMapping( infile, outfile ):
    '''create view in database for alignment stats.

    This view aggregates all information on a per-track basis.

    The table is built from the following tracks:

       context_stats
       bam_stats
    
    '''

    dbh = connect()

    tablename = P.toTable( outfile )
    view_type = "TABLE"

    tables = (( "bam_stats", "track", ),
              ( "context_stats", "track", ) )

    # do not use: ( "picard_stats_alignment_summary_metrics", "track" ), )
    # as there are multiple rows per track for paired-ended data.

    P.createView( dbh, tables, tablename, outfile, view_type )

###################################################################
###################################################################
###################################################################
@follows( createViewMapping )
def views():
    pass

###################################################################
###################################################################
###################################################################
@follows( mapping, qc, views)
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
          mkdir( "%s/bigwigfiles" % PARAMS["web_dir"]), 
          update_report,
          )
def publish():
    '''publish files.'''

    # directory, files
    export_files = {
        "bamfiles" : glob.glob( "*/*.bam" ) + glob.glob( "*/*.bam.bai" ),
        "bigwigfiles" : glob.glob( "*/*.bw" ),
        }

    # publish web pages
    E.info( "publishing report")
    P.publish_report( export_files = export_files)

    E.info( "publishing UCSC data hub" )
    PipelinePublishing.publish_tracks( export_files )

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

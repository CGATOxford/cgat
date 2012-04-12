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

The pipeline assumes the data derive from multiple tissues/conditions (:term:`experiment`) 
with one or more biological and/or technical replicates (:term:`replicate`). A :term:`replicate`
within each :term:`experiment` is a :term:`track`.

Mapping strategy
----------------

The best strategy for mapping and transcriptome assembly depends on the length of your reads.
With short reads, detecting novel splice-junctions is a difficult task. In this case it will be
best to rely on a set of known splice-junctions. Longer reads map more easily across splice-junctions.

From the tophat manual::

   TopHat finds splice junctions without a reference annotation. TopHat version 1.4 maps RNA-seq reads
   first to a reference transcriptome. Only those reads that don't map in this initial process are 
   mapped against the genome.

   Through the second stage of genome mapping, TopHat identifies novel splice junctions and then confirms
   these through mapping to known junctions.
   
   Short read sequencing machines can currently produce reads 100bp or longer, but many exons are 
   shorter than this, and so would be missed in the initial mapping. TopHat solves this problem 
   by splitting all input reads into smaller segments, and then mapping them independently. The segment 
   alignments are "glued" back together in a final step of the program to produce the end-to-end read alignments.

   TopHat generates its database of possible splice junctions from three sources of evidence. The first 
   source is pairings of "coverage islands", which are distinct regions of piled up reads in the 
   initial mapping. Neighboring islands are often spliced together in the transcriptome, so 
   TopHat looks for ways to join these with an intron. The second source is only used when 
   TopHat is run with paired end reads. When reads in a pair come from different exons of a 
   transcript, they will generally be mapped far apart in the genome coordinate space. When 
   this happens, TopHat tries to "close" the gap between them by looking for subsequences of 
   the genomic interval between mates with a total length about equal to the expected distance 
   between mates. The "introns" in this subsequence are added to the database. The third, and 
   strongest, source of evidence for a splice junction is when two segments from the same read #
   are mapped far apart, or when an internal segment fails to map. With long (>=75bp) reads, 
   "GT-AG", "GC-AG" and "AT-AC" introns be found ab initio. With shorter reads, TopHat only 
   reports alignments across "GT-AG" introns

Thus, in order to increase the sensitivity of splice-site detection, it might be best to derive a set of 
splice-junctions using all reads. This is not done automatically, but can be done manually by 
adding a file with junctions to the ``tophat_options`` entry in the configuration file.

The pipeline supplies tophat with a list of all coding exons to facilitate mapping across known
splice-junctions. If they are prioritized, I do not know.

Transcripts are built individually for each :term:`track`. This seems to be the most rigorous way
as there might be conflicting transcripts between replicates and merging the sets might confuse transcript
reconstruction. Also, conflicting transcripts between replicates give an idea of the variability of the data. 
However, if there are only few reads,  there might be a case for building transcript models using reads 
from all replicates of an experiment. However, there is no reason to merge reads between experiments.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

Reads
+++++

Reads are imported by placing files are linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. The ``suffix`` determines the file type.
The following suffixes/file types are possible:

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
------------

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
|samtools            |>=0.1.16           |bam/sam files                                   |
+--------------------+-------------------+------------------------------------------------+
|bedtools            |                   |working with intervals                          |
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


   
.. _tophat: http://tophat.cbcb.umd.edu/
.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml
.. _bamstats: http://www.agf.liv.ac.uk/454/sabkea/samStats_13-01-2011

Code
====

"""

# load modules
from ruffus import *

import Experiment as E
import logging as L
import Database, CSV

import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random

import numpy, sqlite3
import GFF, GTF, IOTools, IndexedFasta
import Tophat
from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.vectors as rovectors
from rpy2.rinterface import RRuntimeError

import Expression

import PipelineGeneset
import PipelineMapping
import PipelineRnaseq
import PipelineMappingQC
import Stats

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % __file__[:-len(".py")],
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
import PipelineTracks

# collect sra nd fastq.gz tracks
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
    glob.glob( "*.sra" ), "(\S+).sra" ) +\
    PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
        glob.glob( "*.fastq.gz" ), "(\S+).fastq.gz" ) +\
        PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
            glob.glob( "*.fastq.1.gz" ), "(\S+).fastq.1.gz" ) +\
            PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
                glob.glob( "*.csfasta.gz" ), "(\S+).csfasta.gz" )

ALL = PipelineTracks.Sample3()
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

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
            rna_index = GFF.readAndIndex( GFF.iterator( IOTools.openFile( rna_file, "r" ) ) )
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
## Read mapping
#########################################################################
SEQUENCEFILES=("*.fastq.1.gz", 
               "*.fastq.gz",
               "*.sra",
               "*.csfasta.gz",
               "*.csfasta.F3.gz",
               )
SEQUENCEFILES_REGEX=regex( r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz|csfasta.F3.gz)")

#########################################################################
#########################################################################
#########################################################################
## Map reads with tophat
#########################################################################
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            add_inputs( buildJunctions, buildReferenceTranscriptome ), 
            r"\1.tophat.bam" )
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
        job_options += " -l mem_free=8G"
    else:
        job_options += " -l mem_free=2G"

    to_cluster = True
    m = PipelineMapping.Tophat( executable = P.substituteParameters( **locals() )["tophat_executable"] )
    infile, reffile, transcriptfile = infiles
    tophat_options = PARAMS["tophat_options"] + " --raw-juncs %(reffile)s " % locals()
    
    # Nick - added the option to map to the reference transcriptome first (built within the pipeline)
    if PARAMS["tophat_include_reference_transcriptome"]:
        prefix = os.path.abspath( P.snip( transcriptfile, ".fa" ) )
        tophat_options = tophat_options + " --transcriptome-index=%s -n 2" % prefix

    statement = m.build( (infile,), outfile ) 
    P.run()

############################################################
############################################################
############################################################
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
@transform( buildTophatStats, suffix(".tsv"), ".load" )
def loadTophatStats( infile, outfile ):
    P.load( infile, outfile )

###################################################################
###################################################################
###################################################################
## Map reads with bowtie
###################################################################
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            add_inputs( os.path.join( PARAMS["bowtie_index_dir"], PARAMS["genome"] + ".fa") ),
            r"\1.bowtie.bam" )
def mapReadsWithBowtie( infiles, outfile ):
    '''map reads with bowtie'''

    job_options= "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
    to_cluster = True
    m = PipelineMapping.Bowtie()
    infile, reffile = infiles
    bowtie_options = "%s --best --strata -a" % PARAMS["bowtie_options"] 
    statement = m.build( (infile,), outfile ) 
    P.run()

###################################################################
###################################################################
###################################################################
## Map reads with bwa
###################################################################
@transform( SEQUENCEFILES,
            SEQUENCEFILES_REGEX,
            r"\1.bwa.bam" )
def mapReadsWithBWA( infile, outfile ):
    '''map reads with shrimp'''

    job_options= "-pe dedicated %i -R y -l mem_free=16G" % PARAMS["bwa_threads"] 
    to_cluster = True
    m = PipelineMapping.BWA()
    statement = m.build( (infile,), outfile ) 
    P.run()

MAPPINGTARGETS = []
mapToMappingTargets = { 'tophat': loadTophatStats,
                        'bowtie': mapReadsWithBowtie,
                        'bwa': mapReadsWithBWA,
                        }
for x in P.asList( PARAMS["mappers"]):
    MAPPINGTARGETS.append( mapToMappingTargets[x] )

@follows( *MAPPINGTARGETS )
def mapping(): pass

###################################################################
###################################################################
###################################################################
## QC targets
###################################################################

# The following is specific to spliced mapping:
# ###################################################################
# ###################################################################
# ###################################################################
# @transform( MAPPINGTARGETS,
#             suffix(".bam"),
#             add_inputs( buildCodingExons ),
#             ".exon.validation.tsv.gz" )
# def buildExonValidation( infiles, outfile ):
#     '''count number of reads mapped, duplicates, etc.
#     '''

#     to_cluster = True
#     infile, exons = infiles
#     statement = '''cat %(infile)s
#     | python %(scriptsdir)s/rnaseq_bam_vs_exons.py
#          --filename-exons=%(exons)s
#          --force
#          --log=%(outfile)s.log
#          --output-filename-pattern="%(outfile)s.%%s.gz"
#     | gzip
#     > %(outfile)s
#     '''

#     P.run()

# ############################################################
# ############################################################
# ############################################################
# @merge( buildExonValidation, "exon_validation.load" )
# def loadExonValidation( infiles, outfile ):
#     '''merge alignment stats into single tables.'''
#     suffix = suffix = ".exon.validation.tsv.gz" 
#     P.mergeAndLoad( infiles, outfile, suffix = suffix )
#     for infile in infiles:
#         track = P.snip( infile, suffix )
#         o = "%s_overrun.load" % track 
#         P.load( infile + ".overrun.gz", o )


############################################################
############################################################
############################################################
@transform( MAPPINGTARGETS,
            suffix(".bam" ),
            add_inputs( buildReferenceTranscriptome ), 
            ".picard_inserts")
def buildPicardInsertSize( infiles, outfile ):
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
            suffix(".bam" ), ".picard_stats")
def buildPicardStats( infile, outfile ):
    '''build alignment stats using picard.

    Note that picards counts reads but they are in fact alignments.
    '''
    PipelineMappingQC.buildPicardAlignmentStats( infile, outfile,
                                                 os.path.join( PARAMS["bowtie_index_dir"],
                                                               PARAMS["genome"] + ".fa" ) )

############################################################
############################################################
############################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "bamstats" ) ) )
@transform( MAPPINGTARGETS,
            suffix(".bam" ), 
             ".bam.report")
def buildBAMReports( infile, outfile ):
    '''build alignment stats using bamstats

    '''
    to_cluster = True

    # requires a large amount of memory to run.
    # only use high-mem machines
    job_options = "-l mem_free=32G"

    # xvfb-run  -f ~/.Xauthority -a 
    track = P.snip( infile, ".bam" )

    # use a fake X display in order to avoid problems with
    # no X connection on the cluster
    xvfb_command = P.which("xvfb-run" )

    # permit multiple servers using -a option
    if xvfb_command: xvfb_command+= " -a "

    # bamstats can not accept a directory as output, hence cd to exportdir
    statement = '''
    cd %(exportdir)s/bamstats;
    %(xvfb_command)s
    %(cmd-run)s bamstats -i ../../%(infile)s -v html -o %(track)s.html 
                         --qualities --mapped --lengths --distances --starts
    >& ../../%(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
@merge( buildPicardStats, "picard_stats.load" )
def loadPicardStats( infiles, outfile ):
    '''merge alignment stats into single tables.'''

    PipelineMappingQC.loadPicardAlignmentStats( infiles, outfile )

        
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
@transform( MAPPINGTARGETS,
            suffix(".bam"),
            ".readstats" )
def buildBAMStats( infile, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = True

    rna_file = os.path.join( PARAMS["annotations_dir"],
                             PARAMS_ANNOTATIONS["interface_rna_gff"] )

    statement = '''python
    %(scriptsdir)s/bam2stats.py
         --force
         --filename-rna=%(rna_file)s
         --remove-rna
         --output-filename-pattern=%(outfile)s.%%s
    < %(infile)s
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''import bam statisticis.'''

    header = ",".join( [P.snip( x, ".readstats") for x in infiles] )
    # filenames = " ".join( [ "<( cut -f 1,2 < %s)" % x for x in infiles ] )
    filenames = " ".join( infiles )
    tablename = P.toTable( outfile )
    E.info( "loading bam stats - summary" )
    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                      --take=2
                   %(filenames)s
                | perl -p -e "s/bin/track/"
                | perl -p -e "s/unique/unique_alignments/"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
            """
    P.run()

    for suffix in ("nm", "nh"):
        E.info( "loading bam stats - %s" % suffix )
        filenames = " ".join( [ "%s.%s" % (x, suffix) for x in infiles ] )
        tname = "%s_%s" % (tablename, suffix)
        
        statement = """python %(scriptsdir)s/combine_tables.py
                      --header=%(header)s
                      --skip-titles
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/%(suffix)s/"
                | python %(scriptsdir)s/csv2db.py
                      --allow-empty
                      --table=%(tname)s 
                >> %(outfile)s
                """
    
        P.run()

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

    to_cluster = True
    statement = '''
       python %(scriptsdir)s/rnaseq_bam_vs_bed.py
              --min-overlap=%(min_overlap)f
              --log=%(outfile)s.log
              %(infile)s %(reffile)s
       > %(outfile)s
       '''

    P.run()

############################################################
############################################################
############################################################
@follows( loadBAMStats )
@merge( buildContextStats, "context_stats.load" )
def loadContextStats( infiles, outfile ):
    """load context mapping statistics."""

    header = ",".join( [P.snip( x, ".contextstats") for x in infiles] )
    filenames = " ".join( infiles  )
    tablename = P.toTable( outfile )

    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --skip-titles
                   %(filenames)s
                | perl -p -e "s/bin/track/; s/\?/Q/g"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
                """
    P.run()
    
    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    cc = Database.executewait( dbhandle, '''ALTER TABLE %(tablename)s ADD COLUMN mapped INTEGER''' % locals())
    statement = '''UPDATE %(tablename)s SET mapped = 
                                       (SELECT b.mapped FROM bam_stats AS b 
                                            WHERE %(tablename)s.track = b.track)''' % locals()

    cc = Database.executewait( dbhandle, statement )
    dbhandle.commit()

###################################################################
###################################################################
###################################################################
## export targets
###################################################################
@merge( mapping,  "view_mapping.load" )
def createViewMapping( infile, outfile ):
    '''create view in database for alignment stats.

    This view aggregates all information on a per-track basis.

    The table is built from the following tracks:

    always present:
       context_stats
       bam_stats
       picard_stats_alignment_summary_metrics
    
    '''

    tablename = P.toTable( outfile )
    # can not create views across multiple database, so use table
    view_type = "TABLE"
    
    dbhandle = connect()
    Database.executewait( dbhandle, "DROP %(view_type)s IF EXISTS %(tablename)s" % locals() )

    statement = '''
    CREATE %(view_type)s %(tablename)s AS
    SELECT b.track, *
    FROM bam_stats AS b,
          context_stats AS c,          
          picard_stats_alignment_summary_metrics AS a
    WHERE 
      b.track = c.track
      AND b.track = a.track
    '''

    Database.executewait( dbhandle, statement % locals() )

    nrows = Database.executewait( dbhandle, "SELECT COUNT(*) FROM view_mapping" ).fetchone()[0]
    
    if nrows == 0:
        raise ValueError( "empty view mapping, check statement = %s" % (statement % locals()) )

    E.info( "created view_mapping with %i rows" % nrows )

    P.touch( outfile )

###################################################################
###################################################################
###################################################################
@follows( createViewMapping )
def views():
    pass

###################################################################
###################################################################
###################################################################
@follows( loadPicardStats, loadBAMStats, loadContextStats )
def qc(): pass

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
        "bamfiles" : glob.glob( "*.accepted.bam" ) + glob.glob( "*.accepted.bam.bai" ),
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

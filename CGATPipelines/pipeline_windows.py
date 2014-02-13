"""
================
Windows pipeline
================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

This pipeline takes mapped reads from ChIP-Seq experiments
such has chromatin marks, MeDIP and performs a window based
enrichment analysis.

   1. Identify differentially occupied regions
   4. Filter DMRs
   5. Calculate DMR statistics
   6. Produce report (SphinxReport)

Methods
=======

Read processing
---------------

For medip-seq analysis, the following filtering steps are typically applied to the mapped data:

   1. removing duplicate reads
   2. removing reads with a mapping quality of less than 10

Medip analysis
--------------

The medip analysis makes use of the MEDIPS R package by `Chavez et al. <http://medips.molgen.mpg.de/>`_ 
(see :pmid:`20802089`).

Briefly, the data is processed in the following way:

1. Quality control
   1. Saturation analysis
   2. Computing CpG coverage
   3. Computing CpG enrichment
   
2. Normalization
   1. Output data normalized by total read depth (rpm - reads per million)
   2. Output normalized relative methylation scores (rms)
   3. Output normalized absolute methylation scores (ams)

Tiling strategies
-----------------

The pipeline implements different tiling strategies.

variable width
   variable width tiles. Tiles are defined based on regions that contain
   short reads and are present in a minimum number of samples.
   
fixwidth_nooverlap
   tiles of size ``tiling_window_size`` with adjacent tiles not overlapping.

fixwidth_overlap
   tiles of size ``tiling_window_size`` with adjacent tiles overlapping by 
   by 50%.

cpg
   windows of size ``tiling_window_size`` are defined around CpG sites.
   Overlapping windows are merged and only windows with a minimum number
   (``tiling_min_cpg``) of CpG sites are kept.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

Input
-----

Reads are imported by placing files or linking to files in the :term:`working directory`.

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

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|Stampy              |>=0.9.0            |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|BWA                 |                   |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|SAMtools            |                   |filtering, SNV / indel calling                  |
+--------------------+-------------------+------------------------------------------------+
|BEDTools            |                   |filtering, SNV / indel calling                  |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.38             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

?

Example
=======

ToDo: make exome sequencing example


Code
====

"""

# load modules
from ruffus import *

import logging as L
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
import csv
import numpy
import sqlite3
import pandas

import CGAT.Experiment as E
import CGAT.Database as Database
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.Stats as Stats
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Pipeline as P
import CGAT.Expression as Expression
import CGAT.BamTools as BamTools
import CGATPipelines.PipelineGeneset as PipelineGeneset
import CGATPipelines.PipelineWindows as PipelineWindows
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelineTracks as PipelineTracks
import CGATPipelines.PipelineMappingQC as PipelineMappingQC

from rpy2.robjects import r as R
import rpy2.robjects as ro

#########################################################################
#########################################################################
#########################################################################
# load options from the config file
P.getParameters( ["%s/pipeline.ini" % os.path.splitext(__file__)[0], 
                  "../pipeline.ini", 
                  "pipeline.ini" ] )

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

###################################################################
###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
# load all tracks - exclude input/control tracks
Sample = PipelineTracks.AutoSample

METHODS = P.asList( PARAMS["methods" ] )

###################################################################
###################################################################
###################################################################
# if conf.py exists: execute to change the above assignmentsn
if os.path.exists("pipeline_conf.py"):
    L.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

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
#@transform( "*CD4*/bam/*.genome.bam",
@follows( mkdir("tags.dir") )
@transform( '*.bam',
            regex("(.*).bam"),
            r"tags.dir/\1.bed.gz" )
def prepareTags( infile, outfile ):
    '''prepare tag files from bam files for medip-seq analysis.

    Optional steps include:

    * deduplication - remove duplicate reads
    * quality score filtering - remove reads below a certain quality score.
    * paired ended data - merge pairs
    * paired ended data - filter by insert size

    '''
    to_cluster = True

    track = P.snip( outfile, ".bed.gz" )

    tmpdir = P.getTempFilename()

    is_paired = BamTools.isPaired( infile )

    current_file = infile

    nfiles = 0
    statement = [ "mkdir %(tmpdir)s" ]

    if "filtering_quality" in PARAMS and PARAMS["filtering_quality"] > 0:
        next_file = "%(tmpdir)s/bam_%(nfiles)i.bam" % locals()
        statement.append( '''samtools view -q %%(filtering_quality)i -b 
                             %(current_file)s 
                             2>> %%(outfile)s.log 
                             > %(next_file)s ''' % locals())
        nfiles += 1
        current_file = next_file

    if "filtering_dedup" in PARAMS and PARAMS["filtering_dedup"]:
        # Picard's MarkDuplicates requries an explicit bam file.
        next_file = "%(tmpdir)s/bam_%(nfiles)i.bam" % locals()

        dedup_method = PARAMS["filtering_dedup_method"]

        if dedup_method == 'samtools':
            statement.append( '''samtools rmdup - - ''' )

        elif dedup_method == 'picard':
            statement.append('''MarkDuplicates INPUT=%(current_file)s
                                               OUTPUT=%(next_file)s
                                               ASSUME_SORTED=true 
                                               METRICS_FILE=%(outfile)s.duplicate_metrics
                                               REMOVE_DUPLICATES=TRUE 
                                               VALIDATION_STRINGENCY=SILENT
                                               2>> %%(outfile)s.log ''' % locals() )
        nfiles += 1
        current_file = next_file

    if is_paired:
        statement.append( '''cat %(current_file)s 
            | python %(scriptsdir)s/bam2bed.py
              --merge-pairs
              --min-insert-size=%(filtering_min_insert_size)i
              --max-insert-size=%(filtering_max_insert_size)i
              --log=%(outfile)s.log
              -
            | python %(scriptsdir)s/bed2bed.py
              --method=sanitize-genome
              --genome-file=%(genome_dir)s/%(genome)s
              --log=%(outfile)s.log
            | cut -f 1,2,3,4
            | sort -k1,1 -k2,2n
            | bgzip > %(outfile)s''')
    else:
        statement.append( '''cat %(current_file)s 
            | python %(scriptsdir)s/bam2bed.py
              --log=%(outfile)s.log
              -
            | python %(scriptsdir)s/bed2bed.py
              --method=sanitize-genome
              --genome-file=%(genome_dir)s/%(genome)s
              --log=%(outfile)s.log
            | cut -f 1,2,3,4
            | sort -k1,1 -k2,2n
            | bgzip > %(outfile)s''')
        

    statement.append( "tabix -p bed %(outfile)s" )
    statement.append( "rm -rf %(tmpdir)s" )

    statement = " ; ".join( statement )

    P.run()

    os.unlink( tmpdir )

#########################################################################
#########################################################################
#########################################################################
@merge( prepareTags, "picard_duplicates.load" )
def loadPicardDuplicateStats( infiles, outfile ):
    '''Merge Picard duplicate stats into single table and load into SQLite.
    '''
    PipelineMappingQC.loadPicardDuplicateStats( infiles, outfile, pipeline_suffix = ".bed.gz" )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( "background.dir" ))
@transform( "*[Ii]nput*.bw",
            regex("(.*).bw"),
            r"background.dir/\1.bed.gz" )
def buildBackgroundWindows( infile, outfile ):
    '''compute regions with high background count in input 
    '''

    job_options = "-l mem_free=16G"

    statement = '''
    python %(scriptsdir)s/wig2bed.py 
             --bigwig-file=%(infile)s 
             --genome-file=%(genome_dir)s/%(genome)s
             --threshold=%(filtering_background_density)f
             --method=threshold
             --log=%(outfile)s.log
    | bgzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( buildBackgroundWindows, "background.dir/background.bed.gz" )
def mergeBackgroundWindows( infiles, outfile ):
    '''build a single bed file of regions with elevated background.'''

    if len(infiles) == 0: 
        # write a dummy file with a dummy chromosome
        # an empty background file would otherwise cause
        # errors downstream in bedtools intersect
        outf = IOTools.openFile( outfile, "w" )
        outf.write("chrXXXX\t1\t2\n" )
        outf.close()
        return

    infiles = " ".join(infiles)
    genomefile = os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS['interface_contigs'])
    statement = '''
    zcat %(infiles)s 
    | bedtools slop -i stdin -b %(filtering_background_extension)i -g %(genomefile)s
    | mergeBed 
    | bgzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( os.path.join( PARAMS["annotations_dir"], 
                          PARAMS_ANNOTATIONS["interface_cpg_bed"] ),
            regex(".*/([^/]*).bed.gz"), 
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_genomic_context_bed"]) ),
            "cpg_context.tsv.gz" )
def buildCpGAnnotation( infiles, outfile ):
    '''annotate the location of CpGs within the genome.'''

    cpg_bed, context_bed = infiles

    statement = '''
    python %(scriptsdir)s/bam_vs_bed.py --min-overlap=0.5 %(cpg_bed)s %(context_bed)s
    | gzip 
    > %(outfile)s'''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildCpGAnnotation, suffix( ".tsv.gz" ), ".load" )
def loadCpGAnnotation( infile, outfile ):
    '''load CpG annotations.'''
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( prepareTags, suffix(".bed.gz"), ".covered.bed.gz" )
def buildCoverageBed( infile, outfile ):
    '''build bed file with regions covered by reads.

    Intervals containing only few reads (tiling_min_reads) are removed.
    '''

    to_cluster = True

    statement = '''
    zcat %(infile)s 
    | cut -f 1,2,3
    | python %(scriptsdir)s/bed2bed.py
          --method=merge
          --merge-distance=%(medips_extension)i
          --log=%(outfile)s.log
          --merge-min-intervals=%(tiling_min_reads)i
    | gzip
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildCoverageBed, suffix(".bed.gz"), ".tsv.gz" )
def buildCpGComposition( infile, outfile ):
    '''compute CpG density across regions covered by reads.
    '''

    statement = '''
    zcat %(infile)s 
    | python %(scriptsdir)s/bed2table.py
          --counter=composition-cpg
          --genome-file=%(genome_dir)s/%(genome)s
    | gzip
    > %(outfile)s
    '''
    P.run()

@merge( buildCoverageBed, "tags.dir/genomic.covered.tsv.gz" )
def buildReferenceCpGComposition( infiles, outfile ):
    '''compute CpG densities across reference windows across
    the genome. 

    This will take the first file of the input and
    shuffle the intervals, and then compute.

    Using fixed size windows across the genome results in
    a very discretized distribution compared to the other 
    read coverage tracks which have intervals of different size.
    '''

    infile = infiles[0]
    contig_sizes = os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_contigs"] )
    gaps_bed = os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_gaps_bed"] )

    # remove windows which are more than 50% N - column 17
    statement = '''bedtools shuffle 
                      -i %(infile)s
                      -g %(contig_sizes)s
                      -excl %(gaps_bed)s
                      -chromFirst
    | python %(scriptsdir)s/bed2table.py
          --counter=composition-cpg
          --genome-file=%(genome_dir)s/%(genome)s
    | gzip
    > %(outfile)s
    '''
    P.run()

    # | awk '$1 !~ /%(tiling_remove_contigs)s/'
    # | awk '$1 == "contig" || $17 < 0.5'

@transform( (buildCpGComposition,buildReferenceCpGComposition), 
            suffix(".tsv.gz"), 
            ".composition.load" )
def loadCpGComposition( infile, outfile ):
    '''load CpG Composition data.'''
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( prepareTags, 
            suffix(".bed.gz"), 
            add_inputs( os.path.join( PARAMS["annotations_dir"], 
                                      PARAMS_ANNOTATIONS["interface_cpg_bed"] ) ),
            ".cpg_coverage.gz" )
def buildCpGCoverage( infiles, outfile ):
    '''count number times certain CpG are covered by reads.

    Reads are processed in the same way as by buildCoverageBed.
    '''

    # coverageBed is inefficient. If bedfile and cpgfile
    # were sorted correspondingly the overlap analysis
    # could be done in very little memory.

    infile, cpg_file = infiles
    to_cluster = True

    job_options = "-l mem_free=16G"

    statement = '''
    zcat %(infile)s 
    | coverageBed -a stdin -b %(cpg_file)s -counts
    | cut -f 6
    | python %(scriptsdir)s/data2histogram.py
    | gzip
    > %(outfile)s
     '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( buildCpGCoverage, "cpg_coverage.load" )
def loadCpGCoverage( infiles, outfile ):
    '''load cpg coverag data.'''
    P.mergeAndLoad( infiles, outfile, 
                    regex="/(.*).cpg_coverage.gz",
                    row_wise = False )

#########################################################################
#########################################################################
#########################################################################
@merge( (buildCoverageBed, mergeBackgroundWindows), "windows.bed.gz")
def buildWindows( infiles, outfile ):
    '''build tiling windows according to parameter tiling_method.

    Remove windows in background.
    '''

    tiling_method = PARAMS["tiling_method"]

    coverage_bed, background_bed = infiles[:-1], infiles[-1]

    to_cluster = True

    coverage_bed = " ".join( coverage_bed )

    if tiling_method == "varwidth":

        infiles = " ".join( infiles )

        statement = '''
        zcat %(coverage_bed)s
        | sort -k1,1 -k2,2n
        | python %(scriptsdir)s/bed2bed.py
              --method=merge
              --merge-distance=0
              --log=%(outfile)s.log
        '''

    elif tiling_method == "fixwidth_nooverlap":

        statement = '''python %(scriptsdir)s/genome_bed.py
                      -g %(genome_dir)s/%(genome)s
                      --window=%(tiling_nonoverlapping_window)i
                      --shift=%(tiling_window_size)i
                      --log=%(outfile)s.log'''

    elif tiling_method == "fixwidth_overlap":

        assert PARAMS["tiling_window_size"] % 2 == 0
        shift = PARAMS["tiling_window_size"] // 2

        statement = '''python %(scriptsdir)s/genome_bed.py
                      -g %(genome_dir)s/%(genome)s
                      --window=%(tiling_window_size)i
                      --shift=%(shift)i
                      --log=%(outfile)s.log'''

    elif tiling_method == "cpg":

        statement = '''cat %(genome_dir)s/%(genome)s.fasta
                       | python %(scriptsdir)s/fasta2bed.py
                      --method=windows-cpg
                      --window-size=%(tiling_window_size)i
                      --min-cpg=%(tiling_min_cpg)i
                      --log=%(outfile)s.log'''
        
    elif os.path.exists( tiling_method ):
        # existing file
        statement = '''mergeBed -i %(tiling_method)s'''

    else:
        raise ValueError("unknow tiling method '%s'" % tiling_method)

    statement += '''
        | awk '$1 !~ /%(tiling_remove_contigs)s/'
        | bedtools intersect -v -wa -a stdin -b %(background_bed)s
        | gzip
        > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildWindows,
            suffix(".bed.gz"),
            ".stats")
def buildWindowStats( infile, outfile ):
    '''compute tiling window size statistics from bed file.'''

    use_cluster = True

    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/gff2histogram.py 
                   --force
                   --format=bed 
                   --data=size
                   --method=hist
                   --method=stats
                   --output-filename-pattern=%(outfile)s.%%s.tsv
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildWindowStats,
            suffix(".stats"),
            "_stats.load")
def loadWindowStats( infile, outfile ):
    '''load window statistics.'''
    P.load( infile + ".hist.tsv", P.snip(infile,".stats") + "_hist" + ".load")
    P.load( infile + ".stats.tsv", outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( buildWindows,
            suffix(".bed.gz"), 
            ".bigbed")
def buildBigBed( infile, outfile ):
    '''bed file with intervals that are covered by reads in any of the experiments.
    '''

    to_cluster = True
    to_cluster = False

    tmpfile = P.getTempFilename()

    contig_sizes = os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_contigs"] )

    statement = '''
    zcat %(infile)s > %(tmpfile)s;
    bedToBigBed %(tmpfile)s %(contig_sizes)s %(outfile)s;
    rm -f %(tmpfile)s
    '''
    P.run()

    try: os.unlink( tmpfile )
    except OSError: pass

#########################################################################
#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( "counts.dir" ) )
@transform( prepareTags,
            regex(".*/(.*).bed.gz"), 
            add_inputs( buildWindows ),
            r"counts.dir/\1.counts.bed.gz" )
def countReadsWithinWindows(infiles, outfile ):
    '''build read counds for windows.'''
    bedfile, windowfile = infiles
    PipelineWindows.countReadsWithinWindows( bedfile,
                                             windowfile,
                                             outfile,
                                             counting_method = PARAMS['tiling_counting_method'] )

#########################################################################
#########################################################################
#########################################################################
@merge( countReadsWithinWindows,
        r"counts.dir/counts.tsv.gz")
def aggregateWindowsReadCounts( infiles, outfile ):
    '''aggregate tag counts for each window.
    '''
    PipelineWindows.aggregateWindowsReadCounts( infiles, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( aggregateWindowsReadCounts, suffix(".tsv.gz"), ".load")
def loadWindowsReadCounts( infile, outfile ):
    '''load a sample of window composition data for QC purposes.'''
    P.load( infile, outfile, limit=10000, shuffle = True)

###################################################################
###################################################################
###################################################################
# 
###################################################################
def getInput( track ):
    '''return a list of input tracks associated with track.

    Associations can be defined in the .ini file in the section
    [input]. For example, the following snippet associates track
    track1 with the bamfiles :file:`track1.bam` and :file:`track2.bam`::
    
       [input]
       track1.bam=input1.bam,input2.bam

    Glob expressions are permitted.

    Default tracks can be specified using a placeholder ``%``. The
    following will associate all tracks with the same bam file::

        [bams]
        %=all.bam


    '''

    input_files = []
 
    # configparser by default converts option names to lower case
    fn = track.asFile()
    fn = fn.lower()

    if "input_%s" % fn in PARAMS:
        input_files.extend( P.asList( PARAMS["input_%s" % fn ] ) )
    elif P.CONFIG.has_section("input"):
        for pattern, value in P.CONFIG.items( "input" ):
            if "%" in pattern:
                pattern = re.sub( "%", "\S+", pattern )
            if re.search( pattern, fn ):
                input_files.extend( P.asList( value ) )

    return input_files

def mapTrack2Input( tracks ):
    '''given a list of tracks, return a dictionary mapping a track to its input
    '''
    
    # select columns in foreground and background
    map_track2input = {}
    for idx, track in enumerate(tracks):

        if track == "interval_id": continue

        try: 
            t = Sample( tablename = track )
        except ValueError, msg:
            print msg
            continue

        input_files = getInput( t )

        ## currently only implement one input file per track
        assert len(input_files) <= 1, "%s more than input: %s" % (track, input_files)

        if len(input_files) == 0:
            map_track2input[track] = None
        else:
            map_track2input[track] = Sample( filename = input_files[0] ).asTable()

    return map_track2input

#########################################################################
#########################################################################
#########################################################################
@transform( loadWindowsReadCounts, suffix(".load"), "_l2foldchange_input.tsv")
def buildWindowsFoldChangesPerInput( infile, outfile ):
    '''Compute fold changes for each sample compared to appropriate input.

    If no input is present, simply divide by average.

    '''
    
    # get all data
    dbh = connect()
    cc = dbh.cursor()
    cc.execute("SELECT * FROM counts" )
    data = cc.fetchall()

    # transpose, remove interval_id column
    data = zip( *data )
    columns = [x[0] for x in cc.description]

    map_track2input = mapTrack2Input( columns )
    take_tracks = [ x for x,y in enumerate( columns ) if y in map_track2input ]
    take_input = [ x for x,y in enumerate( columns ) if y in map_track2input.values() and y != None]

    # build data frame
    dataframe = pandas.DataFrame( dict( [ (columns[x], data[x]) for x in take_tracks ] ) )
    dataframe = dataframe.astype('float64')
    dataframe_input = pandas.DataFrame( dict( [ (columns[x], data[x]) for x in take_input ] ) )


    # add pseudocounts
    pseudocount = 1
    for column in dataframe.columns: dataframe[column] += pseudocount
    for column in dataframe_input.columns: dataframe_input[column] += pseudocount

    # compute normalization ratios
    # total_input / total_column
    ratios = {}
    for column in dataframe.columns:
        i = map_track2input[column]
        if i != None:
            # ratios[column] = float(sum(dataframe_input[i])) / sum( dataframe[column] )
            ratios[column] = dataframe_input[i].median() / dataframe[column].median()
        else:
            ratios[column] = None

    for column in dataframe.columns:
        if ratios[column] != None:
            # normalize by input
            dataframe[column] *= ratios[column] / dataframe_input[map_track2input[column]] 
        else:
            # normalize by median
            dataframe[column] /= dataframe[column].median()

    dataframe = numpy.log2( dataframe )

    dataframe.to_csv( outfile, sep="\t", index=False)

#########################################################################
#########################################################################
#########################################################################
@transform( loadWindowsReadCounts, suffix(".load"), "_l2foldchange_median.tsv")
def buildWindowsFoldChangesPerMedian( infile, outfile ):
    '''Compute l2fold changes for each sample compared to the median count
    in sample.

    '''
    
    # get all data
    dbh = connect()
    cc = dbh.cursor()
    cc.execute("SELECT * FROM counts" )
    data = cc.fetchall()

    # transpose, remove interval_id column
    data = zip( *data )
    columns = [x[0] for x in cc.description]

    take_tracks = [ x for x,y in enumerate( columns ) if y != "interval_id" ]
    # build data frame
    dataframe = pandas.DataFrame( dict( [ (columns[x], data[x]) for x in take_tracks ] ) )
    dataframe = dataframe.astype('float64')

    # add pseudocounts
    pseudocount = 1
    for column in dataframe.columns: dataframe[column] += pseudocount

    for column in dataframe.columns:
        dataframe[column] /= dataframe[column].median()

    dataframe = numpy.log2( dataframe )

    dataframe.to_csv( outfile, sep="\t", index=False)

#########################################################################
#########################################################################
#########################################################################
@transform((buildWindowsFoldChangesPerMedian, buildWindowsFoldChangesPerInput),
           suffix(".tsv"), ".load")
def loadWindowsFoldChanges( infile,outfile ):
    '''load fold change stats'''
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( aggregateWindowsReadCounts,
            suffix(".tsv.gz"),
            "_stats.tsv" )
def summarizeAllWindowsReadCounts( infile, outfile ):
    '''perform summarization of read counts'''

    prefix = P.snip(outfile, ".tsv")
    job_options = "-l mem_free=32G"
    statement = '''python %(scriptsdir)s/runExpression.py
              --method=summary
              --filename-tags=%(infile)s
              --output-filename-pattern=%(prefix)s_
              --log=%(outfile)s.log
              > %(outfile)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( "design*.tsv",
            regex( "(.*).tsv" ),
            add_inputs( aggregateWindowsReadCounts ),
            r"counts.dir/\1_stats.tsv" )
def summarizeWindowsReadCounts( infiles, outfile ):
    '''perform summarization of read counts within experiments.
    '''

    design_file, counts_file = infiles    
    prefix = P.snip(outfile, ".tsv")
    statement = '''python %(scriptsdir)s/runExpression.py
              --method=summary
              --filename-design=%(design_file)s
              --filename-tags=%(counts_file)s
              --output-filename-pattern=%(prefix)s_
              --log=%(outfile)s.log
              > %(outfile)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("dump.dir") )
@transform( "design*.tsv",
            regex( "(.*).tsv" ),
            add_inputs( aggregateWindowsReadCounts ),
            r"dump.dir/\1.tsv.gz" )
def dumpWindowsReadCounts( infiles, outfile ):
    '''output tag tables used for analysis.

    This is for debugging purposes. The tables
    can be loaded into R for manual analysis.
    '''
    design_file, counts_file = infiles

    statement = '''python %(scriptsdir)s/runExpression.py
              --method=dump
              --filename-design=%(design_file)s
              --filename-tags=%(counts_file)s
              --log=%(outfile)s.log
              > %(outfile)s'''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( (summarizeWindowsReadCounts, summarizeAllWindowsReadCounts), 
            suffix("_stats.tsv"), "_stats.load" )
def loadTagCountSummary( infile, outfile ):
    '''load windows summary.'''
    P.load(infile, outfile )
    P.load( P.snip(infile, ".tsv")+ "_correlation.tsv",
            P.snip( outfile, "_stats.load") + "_correlation.load",
            options = "--first-column=track")

#########################################################################
def loadMethylationData( infile, design_file ):
    '''load methylation data for deseq/edger analysis.

    This method creates various R objects:

    countsTable : data frame with counts. 
    groups : vector with groups

    '''

    E.info( "reading data")
    R( '''counts_table = read.delim( '%(infile)s', header = TRUE, 
                                                   row.names = 1, 
                                                   stringsAsFactors = TRUE )''' % locals() )

    E.info( "read data: %i observations for %i samples" % tuple(R('''dim(counts_table)''')))

    # Load comparisons from file
    R('''pheno = read.delim( '%(design_file)s', header = TRUE, stringsAsFactors = TRUE )''' % locals() )

    # Make sample names R-like - substitute - for . and add the .prep suffix
    R('''pheno[,1] = gsub('-', '.', pheno[,1]) ''')

    # Ensure pheno rows match count columns 
    R('''pheno2 = pheno[match(colnames(counts_table),pheno[,1]),,drop=FALSE]''' ) 

    # Subset data & set conditions
    R('''includedSamples <- pheno2$include == '1' ''')
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''conds <- pheno2$group[ includedSamples ]''')

    # Subset data & set conditions
    R('''includedSamples <- pheno2$include == '1' ''')
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''groups <- factor(pheno2$group[ includedSamples ])''')
    R('''pairs = factor(pheno2$pair[ includedSamples ])''')

    groups = R('''levels(groups)''')
    pairs = R('''levels(pairs)''')

    E.info( "filtered data: %i observations for %i samples" % tuple( R('''dim(countsTable)''') ) )

    return groups, pairs

@follows( mkdir("deseq.dir"), mkdir("deseq.dir/plots") )
@transform( "design*.tsv",
            regex( "(.*).tsv" ),
            add_inputs( aggregateWindowsReadCounts ),
            r"deseq.dir/\1.tsv.gz" )
def runDESeq( infiles, outfile ):
    '''estimate differential expression using DESeq.

    The final output is a table. It is slightly edited such that
    it contains a similar output and similar fdr compared to cuffdiff.
    '''

    spike_file = os.path.join("spike.dir", infiles[0]) + ".gz"
    if os.path.exists( spike_file):
        outfile_spike = P.snip(outfile, '.tsv.gz') + '.spike.gz'

        PipelineWindows.runDE( infiles, 
                               outfile_spike, 
                               "deseq.dir", 
                               method = "deseq",
                               spike_file = spike_file )

    PipelineWindows.runDE( infiles, 
                           outfile, 
                           "deseq.dir", 
                           method = "deseq" )

#########################################################################
#########################################################################
#########################################################################
@transform( runDESeq, suffix(".tsv.gz"), ".load" )
def loadDESeq( infile, outfile ):
    '''load DESeq per-chunk summary stats.'''

    prefix = P.snip( outfile, ".load" )

    if os.path.exists( infile + "_size_factors.tsv" ):
        P.load( infile + "_size_factors.tsv", 
                prefix + "_deseq_size_factors.load", 
                collapse = True,
                transpose = "sample")

    for fn in glob.glob( infile + "*_summary.tsv" ):
        prefix = P.snip(fn[len(infile)+1:], "_summary.tsv")

        P.load( fn, 
                prefix + ".deseq_summary.load", 
                collapse = 0,
                transpose = "sample")

    P.touch( outfile )


#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("spike.dir"))
@transform( "design*.tsv",
            regex( "(.*).tsv" ),
            add_inputs( aggregateWindowsReadCounts ),
            r"spike.dir/\1.tsv.gz" )
def buildSpikeIns( infiles, outfile ):
    '''build a table with counts to spike into the original count
    data sets.
    '''

    design_file, counts_file = infiles
    design = P.snip( design_file, ".tsv")
    statement = '''
    zcat %(counts_file)s
    | python %(scriptsdir)s/runExpression.py
            --log=%(outfile)s.log          
            --filename-design=%(design_file)s
            --filename-tags=-
            --method=spike
            --output-filename-pattern=%(outfile)s_
    | gzip 
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("edger.dir") )
@transform( "design*.tsv",
            regex( "(.*).tsv" ),
            add_inputs( aggregateWindowsReadCounts ),
            r"edger.dir/\1.tsv.gz" )
def runEdgeR( infiles, outfile ):
    '''estimate differential methylation using EdgeR
    
    This method applies a paired test. The analysis follows
    the example in chapter 11 of the EdgeR manual.
    '''

    spike_file = os.path.join("spike.dir", infiles[0]) + ".gz"
    if os.path.exists( spike_file):
        outfile_spike = P.snip(outfile, '.tsv.gz') + '.spike.gz'
    
        PipelineWindows.runDE( infiles, 
                               outfile_spike,
                               "edger.dir", 
                               method = "edger",
                               spike_file = spike_file )

    PipelineWindows.runDE( infiles, 
                           outfile, 
                           "edger.dir", 
                           method = "edger" )

#########################################################################
#########################################################################
#########################################################################
@transform( runEdgeR, suffix(".tsv.gz"), ".load" )
def loadEdgeR( infile, outfile ):
    '''load EdgeR per-chunk summary stats.'''

    prefix = P.snip( outfile, ".load" )

    for fn in glob.glob( infile + "*_summary.tsv" ):
        prefix = P.snip(fn[len(infile)+1:], "_summary.tsv")

        P.load( fn, 
                prefix + ".deseq_summary.load", 
                collapse = 0,
                transpose = "sample")

    P.touch( outfile )


#########################################################################
#########################################################################
#########################################################################
@follows( mkdir("roi.dir"))
@transform( "design*.tsv",
            regex( "(.*).tsv" ),
            add_inputs( aggregateWindowsReadCounts ),
            r"roi.dir/\1.tsv.gz" )
def runFilterAnalysis( infiles, outfile ):
    '''output windows applying a filtering criterion.
    
    Does not apply a threshold.
    '''
    PipelineWindows.outputRegionsOfInterest( infiles, outfile )

DIFFTARGETS = []
mapToTargets = { 'deseq': (loadDESeq,runDESeq,),
                 'edger' : (runEdgeR,),
                 'filter' : (runFilterAnalysis,)
                 }
for x in METHODS:
    DIFFTARGETS.extend( mapToTargets[x] )
        
@follows( loadTagCountSummary, 
          loadWindowStats,
          loadWindowsReadCounts,
          loadWindowsFoldChanges,
          *DIFFTARGETS )
def diff_windows(): pass

#########################################################################
@transform( DIFFTARGETS, suffix(".gz"), ".cpg.tsv.gz" )
def computeWindowComposition( infile, outfile ):
    '''for the windows returned from differential analysis, compute CpG content
    for QC purposes.
    '''

    to_cluster = True

    statement = '''
    zcat %(infile)s
    | grep -v "^spike"
    | perl -p -e "s/:/\\t/; s/-/\\t/; s/test_id/contig\\tstart\\tend/"
    | python %(scriptsdir)s/bed2table.py
          --counter=composition-cpg
          --genome-file=%(genome_dir)s/%(genome)s
          --has-header
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( computeWindowComposition, suffix(".tsv.gz"), ".load")
def loadWindowComposition( infile, outfile ):
    '''load a sample of window composition data for QC purposes.'''
    P.load( infile, outfile, limit=10000)

#########################################################################
#########################################################################
#########################################################################
@transform( DIFFTARGETS, suffix(".gz"), ".merged.gz" )
def mergeDMRWindows( infile, outfile ):
    '''merge overlapping windows.

    Sample/control labels are by default inverted to reflect
    that unmethylated windows are of principal interest.
    '''

    to_cluster = True

    statement = '''
    zcat %(infile)s
    | grep -v "^spike"
    | python %(scriptsdir)s/medip_merge_intervals.py
          --log=%(outfile)s.log
          --invert
          --output-filename-pattern=%(outfile)s.%%s.bed.gz
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
def outputSpikeCounts( outfile, infile_name, max_expression, max_fold ):

    df = pandas.read_csv( infile_name, 
                          sep="\t", 
                          index_col = 0 )
    
    E.debug( "read %i rows and %i columns of data" % df.shape)

    if "edger" in outfile.lower():
        # edger: treatment_mean and control_mean do not exist
        # use supplied values directly.
        l10average = numpy.log( df['treatment_mean'] )
        l2fold = numpy.log2( df['fold'] )
    else:
        # use pseudocounts to compute fold changes
        treatment_mean = df['treatment_mean'] + 1
        control_mean = df['control_mean'] + 1
        # build log2 average values
        l10average = numpy.log( (treatment_mean + control_mean) / 2 )
        l2fold = numpy.log2( treatment_mean / control_mean)

    # compute expression bins
    expression_bins = numpy.arange( 0, max_expression, 0.5 )
    fold_change_bins = numpy.arange( -max_fold, max_fold, 0.5 )

    d2hist_counts, xedges, yedges = numpy.histogram2d( l10average, l2fold, 
                                                       bins = ( expression_bins, fold_change_bins ))
    
    dd = pandas.DataFrame( d2hist_counts )
    dd.index = list(xedges[:-1]) 
    dd.columns = list(yedges[:-1])
    dd.to_csv( IOTools.openFile( outfile, "w"),
               sep="\t")

    return df, d2hist_counts, xedges, yedges, l10average, l2fold

@transform( DIFFTARGETS, suffix(".gz"), ".power.gz" )
def buildSpikeResults( infile, outfile ):
    '''build matrices with results from spike-in and upload
    into database.
    '''

    max_expression = 5.0
    max_fold = 4.0
    
    ########################################
    # output and load spiked results
    tmpfile_name = P.getTempFilename( "." )

    spikefile = P.snip( infile, '.tsv.gz') + '.spike.gz'

    if not os.path.exists( spikefile ):
        E.warn('no spike data: %s' % spikefile)
        P.touch( outfile )
        return

    statement = '''zcat %(spikefile)s
    | grep -e "^spike" -e "^test_id"
    > %(tmpfile_name)s
    '''
    P.run()

    spiked, spiked_d2hist_counts, xedges, yedges, spiked_l10average, spiked_l2fold = \
        outputSpikeCounts( P.snip(outfile, ".power.gz") + ".spiked.gz", 
                                  tmpfile_name,
                                  max_expression,
                                  max_fold )

    ########################################
    # output and load unspiked results
    tmpfile_name = P.getTempFilename( "." )

    statement = '''zcat %(infile)s
    | grep -v -e "^spike" 
    > %(tmpfile_name)s
    '''
    P.run()
    
    unspiked, unspiked_d2hist_counts, unspiked_xedges, unspiked_yedges, unspiked_l10average, unspiked_l2fold = \
        outputSpikeCounts( P.snip(outfile, ".power.gz") + ".unspiked.gz", 
                           tmpfile_name,
                           max_expression,
                           max_fold )

    assert xedges.all() == unspiked_xedges.all()


    tmpfile = IOTools.openFile( tmpfile_name, "w")
    tmpfile.write( "expression\tfold\tfdr\tcounts\tpercent\n" )

    fdr_thresholds = [0.01, 0.05 ] + list(numpy.arange(0.1, 1.0, 0.1))
    power_thresholds = numpy.arange(0.1, 1.1, 0.1)

    spiked_total = float(spiked_d2hist_counts.sum().sum())
    unspiked_total = float(unspiked_d2hist_counts.sum().sum())

    outf = IOTools.openFile( outfile, "w")
    outf.write( "fdr\tpower\tintervals\tintervals_percent\n" )

    # significant results
    for fdr in fdr_thresholds:
        take = spiked['qvalue'] < fdr
        
        spiked_d2hist_fdr, xedges, yedges = numpy.histogram2d( spiked_l10average[take], 
                                                        spiked_l2fold[take], 
                                                        bins = ( xedges, yedges ))
        
        spiked_d2hist_fdr_normed = spiked_d2hist_fdr / spiked_d2hist_counts
        spiked_d2hist_fdr_normed = numpy.nan_to_num( spiked_d2hist_fdr_normed )

        # set values without data to 0
        spiked_d2hist_fdr_normed[ spiked_d2hist_counts == 0] = -1.0

        for x,y in itertools.product(range( len(xedges)-1), range(len(yedges)-1)):
            tmpfile.write( "\t".join( map(str, (xedges[x], yedges[y], 
                                                fdr,
                                                spiked_d2hist_fdr[x,y], 
                                                100.0 * spiked_d2hist_fdr_normed[x,y])) ) + "\n")

        
        # take elements in spiked_hist_fdr above a certain threshold
        for power in power_thresholds:
            power_take = spiked_d2hist_fdr_normed >= power
            power_counts = unspiked_d2hist_counts[power_take]

            outf.write( "\t".join( map( str, (fdr, power,
                                              power_counts.sum().sum(),
                                              100.0 * power_counts.sum().sum() / unspiked_total) ) ) + "\n" )
                                              

    tmpfile.close()
    outf.close()

    # upload into table
    method = P.snip(os.path.dirname( outfile ), ".dir")
    tablename = P.toTable(P.snip( outfile, "power.gz" ) + method + ".spike.load" )

    statement = '''cat %(tmpfile_name)s
    | python %(scriptsdir)s/csv2db.py 
           --table=%(tablename)s
           --index=fdr
    > %(outfile)s.log'''

    P.run()
    os.unlink(tmpfile_name)


@transform( buildSpikeResults, suffix('.tsv.power.gz'), '.power.load')
def loadSpikeResults( infile, outfile ):
    '''load work results.'''
    method = P.snip(os.path.dirname( outfile ), '.dir')
    tablename = P.toTable( outfile )
    tablename = '_'.join( (tablename, method) )

    P.load( infile, outfile, options = '--index=fdr,power --allow-empty',
            tablename=tablename)

#########################################################################
@transform( mergeDMRWindows, suffix(".merged.gz"), ".stats" )
def buildDMRStats( infile, outfile ):
    '''compute differential methylation stats.'''
    method = os.path.dirname(infile)
    method = P.snip( method, ".dir")
    PipelineWindows.buildDMRStats( infile, outfile, method=method )

#########################################################################
@transform( mergeDMRWindows, suffix(".merged.gz"), ".fdr" )
def buildFDRStats( infile, outfile ):
    '''compute differential methylation stats.'''
    method = os.path.dirname(infile)
    method = P.snip( method, ".dir")
    PipelineWindows.buildFDRStats( infile, outfile, method=method )

#########################################################################
@transform( mergeDMRWindows, suffix(".merged.gz"), ".merged.gz.all.bed.gz" )
def outputAllWindows( infile, outfile ):
    '''output all bed windows.'''
    PipelineWindows.outputAllWindows( infile, outfile )

#########################################################################
@transform( outputAllWindows, suffix(".all.bed.gz"), 
            (".top.bed.gz", ".bottom.bed.gz" ) )
def outputTopWindows( infile, outfiles ):
    '''output bed with largest/smallest l2fold changes.
    '''
    outfile = outfiles[0]

    ignore_pipe_errors = True

    statement = '''zcat %(infile)s
    | awk '$4 !~ /inf/'
    | sort -k4,4n 
    | tail -n %(bed_export)i 
    | gzip > %(outfile)s || true
    '''
    P.run()

    outfile = outfiles[1]

    statement = '''zcat %(infile)s
    | awk '$4 !~ /inf/'
    | sort -k4,4n 
    | head -n %(bed_export)i 
    | gzip 
    > %(outfile)s || true
    '''
    P.run()

#########################################################################
@merge( buildDMRStats, "dmr_stats.load" )
def loadDMRStats( infiles, outfile ):
    '''load DMR stats into table.'''
    P.concatenateAndLoad( infiles, outfile,
                          missing_value = 0,
                          regex_filename = ".*\/(.*).tsv.stats")

#########################################################################
#########################################################################
#########################################################################
@transform( mergeDMRWindows,
            suffix(".merged.gz"),
            ".stats")
def buildDMRWindowStats( infile, outfile ):
    '''compute window size statistics of DMR from bed file.'''

    to_cluster = True

    statement = '''
    zcat %(infile)s
    | grep -v 'contig'
    | python %(scriptsdir)s/gff2histogram.py 
                   --force
                   --format=bed 
                   --data=size
                   --method=hist
                   --method=stats
                   --output-filename-pattern=%(outfile)s.%%s.tsv
    > %(outfile)s
    '''
    P.run()


###################################################################
###################################################################
###################################################################
@follows( mkdir( "transcriptprofiles.dir" ) )
@transform( prepareTags,
            regex(r".*/([^/].*)\.bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_geneset_coding_exons_gtf"] )),
            r"transcriptprofiles.dir/\1.transcriptprofile.tsv.gz" )
def buildIntervalProfileOfTranscripts( infiles, outfile ):
    '''build a table with the overlap profile
    with protein coding exons.
    '''
    
    to_cluster = True

    bedfile, gtffile = infiles
    
    track = P.snip( os.path.basename( outfile ), '.transcriptprofile.tsv.gz')
    try: 
        t = Sample( filename = track )
    except ValueError, msg:
        print msg
        return

    # no input normalization, this is done later.
    options = ''
    # input_files = getInput( t )

    # ## currently only implement one input file per track
    # assert len(input_files) <= 1, "%s more than input: %s" % (track, input_files)
    

    # if len(input_files) == 1:
    #     options = '--controlfile=%s' % \
    #         (os.path.join( os.path.dirname( bedfile ),
    #                        input_files[0] + '.bed.gz') )

    statement = '''zcat %(gtffile)s
                   | python %(scriptsdir)s/gtf2gtf.py 
                     --filter=representative-transcript 
                     --log=%(outfile)s.log 
                   | python %(scriptsdir)s/bam2geneprofile.py
                      --output-filename-pattern="%(outfile)s.%%s"
                      --force
                      --reporter=transcript
                      --method=geneprofile 
                      --method=tssprofile 
                      --normalize-profile=all
                      --output-all-profiles
                      --resolution-upstream=1000
                      --resolution-downstream=1000
                      --resolution-cds=1000
                      --extension-upstream=5000
                      --extension-downstream=5000
                      %(options)s
                      %(bedfile)s -
                   > %(outfile)s
                '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
# @merge( buildDMRBed, "dmr_overlap.tsv.gz" )
# def computeDMROverlap( infiles, outfile ):
#     '''compute overlap between bed sets.'''
    
#     to_cluster = True

#     if os.path.exists(outfile): 
#         # note: update does not work due to quoting
#         os.rename( outfile, "orig." + outfile )
#         options = "--update=orig.%s" % outfile
#     else:
#         options = ""
    
#     infiles = " ".join( infiles )

#     # note: need to quote track names
#     statement = '''
#         python %(scriptsdir)s/diff_bed.py 
#               --pattern-id=".*/(.*).dmr.bed.gz"
#               --log=%(outfile)s.log 
#               %(options)s %(infiles)s 
#         | awk -v OFS="\\t" '!/^#/ { gsub( /-/,"_", $1); gsub(/-/,"_",$2); } {print}'
#         | gzip
#         > %(outfile)s
#         '''

#     P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( mergeDMRWindows, regex(  "(.*)\.(.*).merged.gz"), r"\1_\2.bed.gz" )
def buildMRBed( infile, outfile ):
    '''output bed6 file with methylated regions.

    All regions are output, even the insignificant ones.

    The score is the log fold change.
    '''
    
    outf = IOTools.openFile( outfile, "w" )
    c = E.Counter()
    for row in csv.DictReader( IOTools.openFile( infile ),
                               dialect = "excel-tab" ):
        c.input += 1

        contig, start, end = re.match("(.*):(\d+)-(\d+)", row["interval_id"] ).groups()
        c.output += 1
        outf.write( "\t".join( (contig, start, end, str(c.input), row["lfold"] ) ) + "\n" )
        
    outf.close()
    
    E.info( "%s" % str(c) )

@follows( loadDMRStats, loadSpikeResults,
          outputAllWindows, outputTopWindows )
def dmr(): pass

@follows( diff_windows, dmr) 
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
          mkdir("%s/medips" % PARAMS["web_dir"]),
          )
def publish():
    '''publish files.'''

    # directory : files
    export_files = { #"bamfiles": glob.glob("*.bam") + glob.glob("*.bam.bai"),
                     #"bigwigfiles": glob.glob("*.bw"),
                     "bedfiles": glob.glob("deseq.dir/*.bed.gz") + glob.glob("edger.dir/*.bed.gz"),
                     }

    # publish web pages
    P.publish_report( export_files = export_files )

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )



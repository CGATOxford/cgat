################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_capseq.py 2900 2011-05-24 14:38:00Z david $
#
#   Copyright (C) 2012 David Sims
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
===================
CAPseq pipeline
===================

:Author: David Sims 
:Release: $Id: pipeline_capseq.py 2900 2011-05-24 14:38:00Z david $
:Date: |today|
:Tags: Python

The CAPseq pipeline imports reads from one or more CAPseq experiments and
performs the following tasks:

   * Align reads to the genome using Bowtie
   * Call peaks using MACS
   * Compare intervals across tracks
   * Compare intervals with external bed files 
   * Annotate intervals with respect to a reference gene set (Ensembl or RNAseq-derived)

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline_capseq.ini` file. The pipeline looks for a configuration file in several places:

   1. The default configuration in the :term:`code directory`.
   2. A shared configuration file :file:`../pipeline.ini`.
   3. A local configuration :file:`pipeline.ini`.

The order is as above. Thus, a local configuration setting will
override a shared configuration setting and a default configuration
setting.

Configuration files follow the ini format (see the python
`ConfigParser <http://docs.python.org/library/configparser.html>` documentation).
The configuration file is organized by section and the variables are documented within 
the file. In order to get a local configuration file in the current directory, type::

    python <codedir>/pipeline_cpg.py config

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.


Input
-----

Reads
++++++

Input are :file:`.fastq.gz`-formatted files. The files should be
labeled in the following way::

   sample-condition-replicate.fastq.gz

Note that neither ``sample``, ``condition`` or ``replicate`` should contain 
``_`` (underscore) and ``.`` (dot) characters as these are used by the pipeline
to delineate tasks.

Requirements
------------

The pipeline requires the information from the following pipelines:

:doc:`pipeline_annotations`

set the configuration variables:
   :py:data:`annotations_database` 
   :py:data:`annotations_dir`

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|bowtie              |>=0.12.7           |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|MACS                |14                 |peak finding                                    |
+--------------------+-------------------+------------------------------------------------+
|Picard              |>=1.4              |Dupliate removal, mapping stats                 |
+--------------------+-------------------+------------------------------------------------+
|BEDTools            |                   |interval comparison                             |
+--------------------+-------------------+------------------------------------------------+


Pipline Output
==============

The results of the computation are all stored in an sqlite relational
database :file:`csvdb`.

Glossary
========

.. glossary::

   bowtie
      bowtie_ - a read mapper

.. _bowtie: http://bowtie-bio.sourceforge.net/index.shtml

Code
====

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections, gzip
import sqlite3
import pysam
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import IOTools
import MAST, GTF, GFF, Bed
import cStringIO
import numpy
import Masker
import fileinput
import gff2annotator
import Experiment as E
import logging as L
import PipelineChipseq as PIntervals
import PipelineTracks
import PipelineMapping
import PipelineGO
from ruffus import *
from rpy2.robjects import r as R
import rpy2.robjects as ro

USECLUSTER = True

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
import Pipeline as P
P.getParameters(  ["%s.ini" % __file__[:-len(".py")],  "../pipeline.ini", "pipeline.ini" ] )
PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["ensembl_annotation_dir"],"pipeline_annotations.py" )

###################################################################
###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
# load all tracks - exclude input/control tracks
Sample = PipelineTracks.Sample3

TRACKS = PipelineTracks.Tracks( Sample ).loadFromDirectory( 
    [ x for x in glob.glob( "*.export.txt.gz" ) if PARAMS["tracks_control"] not in x ],
      "(\S+).export.txt.gz" ) +\
      PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
          [ x for x in glob.glob( "*.sra" ) if PARAMS["tracks_control"] not in x ], 
          "(\S+).sra" ) +\
          PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
              [x for x in glob.glob( "*.fastq.gz" ) if PARAMS["tracks_control"] not in x], 
              "(\S+).fastq.gz" ) +\
              PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
                  [x for x in glob.glob( "*.fastq.1.gz" ) if PARAMS["tracks_control"] not in x], 
                  "(\S+).fastq.1.gz" ) +\
                  PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
                      [ x for x in glob.glob( "*.csfasta.gz" ) if PARAMS["track_control"] not in x], 
                        "(\S+).csfasta.gz" )
for X in TRACKS:
    print "TRACK=", X, "\n"

TRACKS_CONTROL = PipelineTracks.Tracks( Sample ).loadFromDirectory( 
    [ x for x in glob.glob( "*.export.txt.gz" ) if PARAMS["tracks_control"] in x ],
      "(\S+).export.txt.gz" ) +\
      PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
          [ x for x in glob.glob( "*.sra" ) if PARAMS["tracks_control"] in x ], 
          "(\S+).sra" ) +\
          PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
              [x for x in glob.glob( "*.fastq.gz" ) if PARAMS["tracks_control"] in x], 
              "(\S+).fastq.gz" ) +\
              PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
                  [x for x in glob.glob( "*.fastq.1.gz" ) if PARAMS["tracks_control"] in x], 
                  "(\S+).fastq.1.gz" ) +\
                  PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
                      [ x for x in glob.glob( "*.csfasta.gz" ) if PARAMS["track_control"] in x], 
                        "(\S+).csfasta.gz" )
for X in TRACKS_CONTROL:
    print "TRACK_CONTROL=", X, "\n"

def getControl( track ):
    '''return appropriate control for a track'''
    n = track.clone()
    n.condition = PARAMS["tracks_control"]
    return n

###################################################################
###################################################################
###################################################################
# aggregate per experiment
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue") )
# aggregate per condition
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition",) )
# aggregate per tissue
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue",) )
# compound targets : all experiments
TRACKS_MASTER = EXPERIMENTS.keys() + CONDITIONS.keys()
# compound targets : correlation between tracks
TRACKS_CORRELATION = TRACKS_MASTER + list(TRACKS)

print "Expts=", EXPERIMENTS, "\n"

###################################################################
###################################################################
###################################################################
## Section 1: MAP READS
@follows(mkdir("bam"))
@transform( ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz" ),
            regex( r"(\S+).(export.txt.gz|fastq.1.gz|fastq.gz|sra|csfasta.gz)"), 
            r"bam/\1.bam" )
def buildBAM( infile, outfile ):
    '''map reads with bowtie'''
    to_cluster = True
    track = P.snip( os.path.basename(outfile), ".bam" )
    job_options= "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
    m = PipelineMapping.Bowtie()
    reffile = PARAMS["samtools_genome"]
    statement = m.build( (infile,), outfile ) 
    P.run()

#########################################################################
@transform( buildBAM, suffix(".bam"), ".alignstats" )
def buildPicardAlignStats( infile, outfile ):
    '''Gather BAM file alignment statistics using Picard '''
    to_cluster = True
    track = P.snip( os.path.basename(infile), ".bam" )
    statement = '''CollectAlignmentSummaryMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(samtools_genome)s ASSUME_SORTED=true OUTPUT=%(outfile)s VALIDATION_STRINGENCY=SILENT ''' % locals()
    P.run()

############################################################
@merge( buildPicardAlignStats, "picard_align_stats.load" )
def loadPicardAlignStats( infiles, outfile ):
    '''Merge Picard alignment stats into single table and load into SQLite.'''
    # Join data for all tracks into single file
    outf = P.getTempFile()
    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".alignstats" )
        if not os.path.exists( f ): 
            E.warn( "File %s missing" % f )
            continue
        lines = [ x for x in open( f, "r").readlines() if not x.startswith("#") and x.strip() ]
        if first: outf.write( "%s\t%s" % ("track", lines[0] ) )
        first = False
        for i in range(1, len(lines)):
            outf.write( "%s\t%s" % (track,lines[i] ))
    outf.close()
    tmpfilename = outf.name

    # Load into database
    tablename = P.toTable( outfile )
    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s'''
    P.run()
    os.unlink( tmpfilename )

#########################################################################
@transform( buildBAM, suffix(".bam"), ".gcstats" )
def buildPicardGCStats( infile, outfile ):
    '''Gather BAM file GC bias stats using Picard '''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".bam" )
    statement = '''CollectGcBiasMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(samtools_genome)s 
                   OUTPUT=%(outfile)s CHART_OUTPUT=%(outfile)s.pdf SUMMARY_OUTPUT=%(outfile)s.summary VALIDATION_STRINGENCY=SILENT ''' % locals()
    P.run()

############################################################
@merge( buildPicardGCStats, "picard_gcbias_stats.load" )
def loadPicardGCStats( infiles, outfile ):
    '''Merge Picard insert size stats into single table and load into SQLite.'''
    # Join data for all tracks into single file
    outf = P.getTempFile()
    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".dedup.gcstats" )
        if not os.path.exists( f ): 
            E.warn( "File %s missing" % f )
            continue
        lines = [ x for x in open( f, "r").readlines() if not x.startswith("#") and x.strip() ]
        if first: outf.write( "%s\t%s" % ("track", lines[0] ) )
        first = False
        outf.write( "%s\t%s" % (track,lines[1] ))
    outf.close()
    tmpfilename = outf.name

    # Load into database
    tablename = P.toTable( outfile )
    statement = '''cat %(tmpfilename)s
                   | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                   > %(outfile)s '''
    P.run()
    os.unlink( tmpfilename )

#########################################################################
@transform( buildBAM, suffix(".bam"), ".readstats" )
def buildBAMStats( infile, outfile ):
    '''Count number of reads mapped, duplicates, etc. using bam2stats.py'''
    to_cluster = USECLUSTER
    scriptsdir = PARAMS["general_scriptsdir"]
    statement = '''python %(scriptsdir)s/bam2stats.py --force 
                   --output-filename-pattern=%(outfile)s.%%s < %(infile)s > %(outfile)s'''
    P.run()

#########################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''Import bam statistics into SQLite'''

    scriptsdir = PARAMS["general_scriptsdir"]
    header = ",".join( [P.snip( os.path.basename(x), ".readstats") for x in infiles] )
    filenames = " ".join( [ "<( cut -f 1,2 < %s)" % x for x in infiles ] )
    tablename = P.toTable( outfile )
    E.info( "loading bam stats - summary" )
    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --ignore-empty
                   %(filenames)s
                | perl -p -e "s/bin/track/"
                | perl -p -e "s/unique/unique_alignments/"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --allow-empty
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s"""
    P.run()

    for suffix in ("nm",):
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
                      --table=%(tname)s 
                      --allow-empty
                >> %(outfile)s """
        P.run()

#########################################################################
@transform( buildBAM, suffix( ".bam"), ".dedup.bam")
def dedup(infiles, outfile):
        '''Remove duplicate alignments from BAM files.'''
        to_cluster = USECLUSTER
        track = P.snip( outfile, ".bam" )
        dedup_method = PARAMS["dedup_method"]
        if dedup_method == 'samtools':
            statement = '''samtools rmdup %(infiles)s %(outfile)s; ''' % locals()    
        elif dedup_method == 'picard':
            statement = '''MarkDuplicates INPUT=%(infiles)s  ASSUME_SORTED=true OUTPUT=%(outfile)s 
                           METRICS_FILE=%(track)s.dupstats REMOVE_DUPLICATES=true 
                           VALIDATION_STRINGENCY=SILENT; ''' % locals()
        statement += '''samtools index %(outfile)s; ''' % locals()
        P.run()

#########################################################################
@merge( dedup, "picard_duplicate_stats.load" )
def loadPicardDuplicateStats( infiles, outfile ):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''
    # Join data for all tracks into single file
    outf = open('dupstats.txt','w')
    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".dedup.bam" )
        statfile = P.snip(f, ".bam" )  + ".dupstats"
        if not os.path.exists( statfile ): 
            E.warn( "File %s missing" % statfile )
            continue
        lines = [ x for x in open( statfile, "r").readlines() if not x.startswith("#") and x.strip() ]
        if first: outf.write( "%s\t%s" % ("track", lines[0] ) )
        first = False
        outf.write( "%s\t%s" % (track,lines[1] ))

    outf.close()
    tmpfilename = outf.name

    # Load into database
    tablename = P.toTable( outfile )
    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s '''
    P.run()

############################################################
@transform( dedup, suffix( ".dedup.bam"), ".dedup.mapped.bam")
def removeUnmapped(infile, outfile):
    ''' Remove unmapped reads from BAM file'''
    to_cluster = True
    statement = '''samtools view -F 4 -bh %(infile)s > %(outfile)s; 
                   samtools index %(outfile)s; ''' % locals()
    P.run()

############################################################
@follows( removeUnmapped )
@files( [ (("bam/%s.dedup.mapped.bam" % x, 
            "bam/%s.dedup.mapped.bam" % getControl(x)), 
            "bam/%s.norm.bam" % x ) for x in TRACKS ] )
def normaliseBAMs( infiles, outfile ):
    '''Sample reads from larger of two BAM files so CAP and control samples have aprox same aligned read number.'''
    infile, controlfile = infiles
    controlfilenorm = controlfile.replace(".dedup.mapped",".norm")
    to_cluster = True

    # Count reads in chip file
    countfile1 = infile.replace(".bam",".count")
    statement= ''' samtools idxstats %s | awk '{s+=$3} END {print s}' > %s ''' % ( infile,countfile1 )
    P.run()
    fh = open(countfile1,"r")
    chip_reads = int(fh.read())
    fh.close()

    # Count reads in input file
    countfile2 = controlfile.replace(".bam",".count")
    statement= ''' samtools idxstats %s | awk '{s+=$3} END {print s}' > %s ''' % ( controlfile,countfile2 )
    P.run()
    fh = open(countfile2,"r")
    input_reads = int(fh.read())
    fh.close()

    # If chip > input then sample chip reads
    if chip_reads > input_reads:
        PIntervals.buildSimpleNormalizedBAM( (infile,countfile1), outfile, input_reads)
        statement = '''cp %(controlfile)s %(controlfilenorm)s; cp %(controlfile)s.bai %(controlfilenorm)s.bai; '''
        P.run()
    else:
        PIntervals.buildSimpleNormalizedBAM( (controlfile,countfile2), controlfilenorm, chip_reads)
        statement = '''cp %(infile)s %(outfile)s; cp %(infile)s.bai %(outfile)s.bai; '''
        P.run()


############################################################
@follows(normaliseBAMs, mkdir("merged_bams"))
@files( [( [ "bam/%s.norm.bam" % y for y in EXPERIMENTS[x]], 
           "merged_bams/%s.merge.bam" % str(x).replace("-agg","")) for x in EXPERIMENTS ] )
def mergeReplicateBAMs( infiles, outfile ):
    '''Merge normalised BAM files for all replicates, then sort and index. '''
    track = P.snip( outfile, ".merge.bam" )
    in_list = " ".join(infiles)
    statement = '''samtools merge %(track)s.rep.bam %(in_list)s;
                   samtools sort  %(track)s.rep.bam %(track)s.merge; 
                   samtools index %(outfile)s;
                   rm %(track)s.rep.bam; ''' % locals()
    P.run()

############################################################   
@transform( mergeBams, suffix(".bam"), ".bw" )
def getMergedBigWig( infile, outfile ):
    '''Merge multiple BAM files per replicate to produce a single non peak-shifted bigwig file'''

    statement = '''python %(scriptsdir)s/bam2wiggle.py --output-format=bigwig %(infile)s %(outfile)s '''
    P.run()
    
############################################################   
@follows( normaliseBAMs, mkdir("merged_bigwigs") )
@files( [( [ "bam/%s.norm.bam" % y for y in EXPERIMENTS[x]], 
           "merged_bigwigs/%s.merged.bw" % str(x).replace("-agg","")) for x in EXPERIMENTS ] )
def getMergedBigWigPeakShift( infiles, outfile ):
    '''Merge multiple BAM files per replicate to produce a single peak-shifted bigwig file'''
    expt = P.snip( os.path.basename( outfile ), ".merged.bw").replace("-agg","")
    in_list = " --bamfile=".join(infiles)
    
    offsets = []
    for t in infiles:
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    shifts = " --shift=".join(offsets)
    statement = '''python %(scriptsdir)s/bam2wiggle.py 
                      --output-format=bigwig
                      %(in_list)s
                      %(shifts)s '''
    P.run()
    
############################################################
############################################################
############################################################
## BUILD INTERVALS USING CONTROL SAMPLE
@follows( normaliseBAMs, mkdir("macs"), mkdir("macs/with_input") )
@files( [ (("bam/%s.norm.bam" % x, 
            "bam/%s.norm.bam" % getControl(x)), 
            "macs/with_input/%s.macs" % x ) for x in TRACKS ] )
def runMACS( infiles, outfile ):
    '''Run MACS for peak detection using control sample.'''
    infile, controlfile = infiles
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".norm.bam" )
    dir = os.path.dirname(infile)

    # change to macs directory and run MACS from there so that wig files end up in right place
    statement = '''cd macs/with_input; 
                   macs14 -t ../../%(infile)s 
                          -c ../../%(controlfile)s
                          --name=%(track)s
                          --format=BAM
                          --diag
                          --wig -S
                          %(macs_options)s 
                   >& ../../%(outfile)s;
                   cd ../..;''' 
    P.run() 
    
############################################################
@transform( runMACS, regex(r"macs/with_input/(\S+).macs"),
            inputs( (r"macs/with_input/\1.macs", r"bam/\1.norm.bam")), 
            r"macs/with_input/\1_macs.load" )
def loadMACS( infiles, outfile ):
    '''Load MACS intervals into database and filter on fold change and p-value'''
    infile, bamfile = infiles
    PIntervals.loadMACS( infile, outfile, bamfile )
    
############################################################
@merge( runMACS, "macs.summary" )
def summarizeMACS( infiles, outfile ):
    '''Parse MACS summary statistics from log file'''
    PIntervals.summarizeMACS( infiles, outfile )

############################################################
@transform( summarizeMACS, suffix(".summary"), "_summary.load" )
def loadMACSSummary( infile, outfile ):
    '''load macs summary into database'''
    P.load( infile, outfile, "--index=track" )

############################################################
@follows( mkdir("intervals") )
@transform( loadMACS, regex(r"macs/with_input/(\S+)_macs.load"), r"intervals/\1.bed" )
def exportIntervalsAsBed( infile, outfile ):
    '''Export MACS intervals from database as BED file and filter on fold change'''
    fc = PARAMS["intervals_min_fc"]
    PIntervals.exportMacsIntervalsAsBed( infile, outfile, fc )

############################################################
############################################################
############################################################
## BUILD INTERVALS WITHOUT CONTROL SAMPLE
@follows( normaliseBAMs, mkdir("macs"), mkdir("macs/no_input") )
@files( [ ("bam/%s.norm.bam" % x, 
           "macs/no_input/%s.solo.macs" % x ) for x in TRACKS ] )
def runMACSsolo( infile, outfile ):
    '''Run MACS for peak detection.'''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".norm.bam" )
    statement = '''cd macs/no_input/; 
                   macs14 -t ../../%(infile)s 
                          --name=%(track)s.solo
                          --format=BAM
                          --diag
                          --wig -S
                          %(macs_options)s 
                   >& ../../%(outfile)s;
                   cd ../..;''' 
    P.run() 
    
############################################################
@transform( runMACSsolo, regex(r"macs/no_input/(\S+).solo.macs"),
            inputs( (r"macs/no_input/\1.solo.macs", r"bam/\1.bam")), 
            r"macs/no_input/\1.solo_macs.load" )
def loadMACSsolo( infiles, outfile ):
    infile, bamfile = infiles
    PIntervals.loadMACS( infile, outfile, bamfile )
    
############################################################
@merge( runMACSsolo, "macs_solo.summary" )
def summarizeMACSsolo( infiles, outfile ):
    '''run MACS for peak detection.'''
    PIntervals.summarizeMACSsolo( infiles, outfile )

############################################################
@transform( summarizeMACSsolo, suffix(".summary"), "_summary.load" )
def loadMACSsoloSummary( infile, outfile ):
    '''load macs summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
@transform( loadMACSsolo, regex(r"macs/no_input/(\S+).macs.load"), r"intervals/\1.bed" )
def exportIntervalsAsBedsolo( infile, outfile ):
    '''Export list of intervals passing fold change threshold to file '''
    fc = PARAMS["intervals_min_fc"]
    PIntervals.exportMacsIntervalsAsBed( infile, outfile, fc )

############################################################
############################################################
############################################################
## Merge nearby intervals
@transform( (exportIntervalsAsBed), suffix(".bed"), ".merged.bed")
def mergeIntervals(infile, outfile):
    '''Merge intervals less than n bases apart in each dataset and update foldchange scores'''
    d = PARAMS["intervals_merge_dist"]
    method = PARAMS["intervals_foldchange_merge_method"]
    PIntervals.mergeIntervalsWithScores( infile, outfile, d, method )

############################################################
@transform( mergeIntervals, suffix(".merged.bed"), ".merged.cleaned.bed")
def sanitiseIntervals(infile, outfile):
    '''sanatise so that intervals do not exceed contig length'''
    statement = '''cat %(infile)s | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s -L %(outfile)s.log > %(outfile)s;'''
    P.run()

############################################################
@transform( (sanitiseIntervals), suffix(".merged.cleaned.bed"), ".merged.cleaned.bed.load" )
def loadMergedIntervals( infile, outfile ):
    '''load combined intervals.

    Also, re-evaluate the intervals by counting reads within
    the interval. In contrast to the initial pipeline, the
    genome is not binned. In particular, the meaning of the
    columns in the table changes to:

    nProbes: number of reads in interval
    PeakCenter: position with maximum number of reads in interval
    AvgVal: average coverage within interval

    If *replicates* is true, only replicates will be considered
    for the counting. Otherwise the counts aggregate both replicates
    and conditions.
    '''

    # Write header to output file
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    headers = ( "contig","start","end","interval_id","nPeaks","PeakCenter","Length","AvgVal","PeakVal","nProbes", "Fold" )
    tmpfile.write( "\t".join(headers) + "\n" )
    contig,start,end,interval_id,npeaks,peakcenter,length,avgval,peakval,nprobes = "",0,0,0,0,0,0,0,0,0

    # Get SAM file and Macs offset
    samfiles, offsets = [], []
    track = P.snip( os.path.basename(infile), ".merged.cleaned.bed")
    base_track = track.replace(".solo","")

    fn = "bam/%s.norm.bam" % track
    assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, track)
    samfiles.append( pysam.Samfile( fn,  "rb" ) )
    if track.find("solo") > -1:
        fn = "macs/no_input/%s.macs" % track
    else:
        fn = "macs/with_input/%s.macs" % track
    if os.path.exists( fn ):
        offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    # Loop over input Bed file and calculate stats for merged intervals
    c = E.Counter()
    for line in open(infile, "r"):
        c.input += 1
        contig, start, end, int_id, fc = line[:-1].split()[:5]
        start, end = int(start), int(end)
        interval_id = c.input

        npeaks, peakcenter, length, avgval, peakval, nprobes = PIntervals.countPeaks( contig, start, end, samfiles, offsets )

        # nreads can be 0 if the intervals overlap only slightly
        # and due to the binning, no reads are actually in the overlap region.
        # However, most of these intervals should be small and have already be deleted via 
        # the merge_min_interval_length cutoff.
        # do not output intervals without reads.
        if nprobes == 0:
            c.skipped_reads += 1
            
        c.output += 1
        tmpfile.write( "\t".join( map( str, (contig,start,end,int_id,npeaks,peakcenter,length,avgval,peakval,nprobes,fc) )) + "\n" )
 
    tmpfile.close()

    tmpfilename = tmpfile.name
    tablename = "%s_macs_merged_intervals" % track
    
    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --index=interval_id
                       --index=contig,start 
                       --table=%(tablename)s
                   < %(tmpfilename)s > %(outfile)s '''
    P.run()
    os.unlink( tmpfile.name )
    L.info( "%s\n" % str(c) )

############################################################
############################################################
############################################################
## Assess background (non-peak) binding
@follows( sanitiseIntervals, mkdir("background") )
@files( [ (("bam/%s.norm.bam" % x, 
            "intervals/%s.merged.cleaned.bed" % x), 
            "background/%s.bg" % x ) for x in TRACKS ] )
def getBackground(infiles, outfile):
    '''Count the number of reads in the bamfile used for MACS that do not overlap an interval'''
    bam, bed = infiles
    to_cluster = True
    statement = '''intersectBed -abam %(bam)s -b %(bed)s -v -bed | wc -l > %(outfile)s; '''
    statement += '''intersectBed -abam %(bam)s -b %(bed)s -u -bed | wc -l >> %(outfile)s;'''
    statement += '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()

############################################################
@transform( getBackground, suffix(".bg"), ".bg.load")
def loadBackground(infile, outfile):
    '''load background into database'''
    track = P.snip( os.path.basename( infile ), ".bg" )
    header = "out_peaks,in_peaks"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                        --table=%(track)s_background
                        --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## Assess effect of altering fold change threshold
@follows( exportIntervalsAsBed, mkdir("foldchange") )
@files( [ ("macs/with_input/%s.macs" % x, 
           "foldchange/%s.foldchange" % x ) for x in TRACKS ] )
def thresholdFoldChange(infile, outfile):
    '''Assess interval overlap between conditions at different fold change thresholds '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Make bed files for different fold change thresholds
    track = P.snip( os.path.basename( infile ), ".macs" ).replace("-","_")
    macsdir = os.path.dirname(infile)
    foldchange = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    for fc in foldchange:
        cc = dbhandle.cursor()
        query = '''SELECT contig, start, end, interval_id, fold FROM %(track)s_macs_intervals WHERE fold > %(fc)i ORDER by contig, start;''' % locals()
        cc.execute( query )
        outbed = "foldchange/" + track + ".fc" + str(fc) + ".bed"
        outs = open( outbed, "w")

        for result in cc:
            contig, start, end, interval_id,fold = result
            outs.write( "%s\t%i\t%i\t%s\t%i\n" % (contig, start, end, str(interval_id), fold) )
        cc.close()
        outs.close()

        statement = '''echo %(fc)i >> %(outfile)s; cat %(outbed)s | wc -l >> %(outfile)s; '''
        P.run()
    
    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

############################################################
@transform( thresholdFoldChange, suffix(".foldchange"), ".foldchange.load")
def loadFoldChangeThreshold(infile, outfile):
    '''Load intervals overlapping chipseq into database '''
    track = P.snip( os.path.basename( infile ), ".foldchange" ).replace(".","_").replace("-","_")
    header = "threshold,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_foldchange
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
@transform( thresholdFoldChange, suffix(".foldchange"), ".shared.foldchange")
def sharedIntervalsFoldChangeThreshold(infile, outfile):
    '''identify shared intervals between datasets at different foldchange thresholds'''

    # Open foldchange file and read
    fc_file = open(infile, "r")
    foldchange = []
    for line in fc_file:
        threshold, interval_count = line.split()
        foldchange.append(threshold)

    in_track = P.snip( os.path.basename( infile ), ".foldchange" )
    in_dir = os.path.dirname(infile)
    out_dir = in_dir.replace("macs","fc")
    try: os.mkdir( out_dir )
    except OSError: pass

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    # for each foldchange 
    for fc in foldchange:

        in_bed = "foldchange/" + in_track.replace("-","_") + ".fc" + str(fc) + ".bed"

        # For each track
        for track in TRACKS:
            if (str(track) != in_track):
                compare_bed = "foldchange/" + str(track).replace("-","_") + ".fc" + str(fc) + ".bed"
                statement = '''echo %(track)s %(fc)s >> %(outfile)s; intersectBed -a %(in_bed)s -b %(compare_bed)s -u | wc -l >> %(outfile)s; ''' 
                P.run()

    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; sed -i '{s/ /\\t/g}' %(outfile)s; '''
    P.run()

############################################################
@transform( sharedIntervalsFoldChangeThreshold, suffix(".shared.foldchange"), ".shared.foldchange.load")
def loadSharedIntervalsFoldChangeThreshold(infile, outfile):
    '''Load intervals overlapping other tracks into database '''
    track = P.snip( os.path.basename( infile ), ".shared.foldchange" ).replace(".","_").replace("-","_")
    header = "track,threshold,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_foldchange_shared
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()


############################################################
############################################################
############################################################
## Calculate replicated intervals
@follows( sanitiseIntervals, mkdir("replicated_intervals") )
@files( [( [ "intervals/%s.merged.cleaned.bed" % y for y in EXPERIMENTS[x]], 
           "replicated_intervals/%s.rep.bed" % str(x).replace("-agg","")) for x in EXPERIMENTS ] )
def replicatedIntervals( infiles, outfile ):
    '''Combine replicates between experiments.
       First all intervals are merged across replicates 
       and then merged intervals that do not overlap an interval in all replicates are removed. '''

    expt = P.snip( os.path.basename( outfile ), ".rep.bed").replace("-agg","")
    outdir = os.path.dirname( outfile )
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name

    # merge all files to get total set of intervals
    in_list = " ".join(infiles)
    statement = '''cat %(in_list)s | mergeBed -i stdin -scores max > %(outdir)s/%(expt)s.merged.bed; 
                   cat %(outdir)s/%(expt)s.merged.bed > %(outfile)s;''' % locals()
    P.run()

    # Remove intervals that do not overlap intervals in all reps
    for i in infiles:
        statement = '''intersectBed -a %(outfile)s -b %(i)s -u > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s; ''' 
        P.run()
   
    
############################################################
@transform( replicatedIntervals, suffix(".rep.bed"), ".rep.bed.load")
def loadReplicatedIntervals(infile, outfile):
    '''load replicated intervals.

    Also, re-evaluate the intervals by counting reads within
    the interval. In contrast to the initial pipeline, the
    genome is not binned. In particular, the meaning of the
    columns in the table changes to:

    nProbes: number of reads in interval
    PeakCenter: position with maximum number of reads in interval
    AvgVal: average coverage within interval    '''

    track = P.snip( os.path.basename(infile), ".rep.bed" )
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]

    # Write header to output file
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    headers = ( "contig","start","end","interval_id","nPeaks","PeakCenter","Length","AvgVal","PeakVal","nProbes", "Fold" )
    tmpfile.write( "\t".join(headers) + "\n" )
    contig,start,end,interval_id,npeaks,peakcenter,length,avgval,peakval,nprobes = "",0,0,0,0,0,0,0,0,0

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( pysam.Samfile( fn,  "rb" ) )
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    # Loop over input Bed file and calculate stats for replicated intervals
    c = E.Counter()
    for line in open(infile, "r"):
        c.input += 1
        contig, start, end, fc = line[:-1].split()[:4]
        start, end = int(start), int(end)
        interval_id = c.input

        npeaks, peakcenter, length, avgval, peakval, nprobes = PIntervals.countPeaks( contig, start, end, samfiles, offsets )

        # nreads can be 0 if the intervals overlap only slightly
        # and due to the binning, no reads are actually in the overlap region.
        # However, most of these intervals should be small and have already be deleted via 
        # the merge_min_interval_length cutoff.
        # do not output intervals without reads.
        if nprobes == 0:
            c.skipped_reads += 1
            
        c.output += 1
        tmpfile.write( "\t".join( map( str, (contig,start,end,interval_id,npeaks,peakcenter,length,avgval,peakval,nprobes,fc) )) + "\n" )
 
    tmpfile.close()
    tmpfilename = tmpfile.name

    # Load into database
    tablename = "%s_replicated_intervals" % track
    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --index=interval_id
                       --index=contig,start 
                       --table=%(tablename)s
                   < %(tmpfilename)s > %(outfile)s '''
    P.run()
    os.unlink( tmpfile.name )
    L.info( "%s\n" % str(c) )

############################################################
@transform( loadReplicatedIntervals, suffix(".rep.bed.load"), ".replicated.bed" )
def exportReplicatedIntervalsAsBed( infile, outfile ):
    '''export locations for all intervals.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    track = P.snip( os.path.basename(infile), ".rep.bed.load" ).replace("-","_")

    cc = dbhandle.cursor()
    statement = "SELECT contig, start, end, interval_id FROM %s_replicated_intervals ORDER by contig, start" % track
    cc.execute( statement )

    # Write to bed file
    outs = open( outfile, "w")
    for result in cc:
        contig, start, stop, interval_id = result
        outs.write( "%s\t%i\t%i\t%i\n" % (contig, start, stop, interval_id) )
    cc.close()
    outs.close()
    

############################################################
############################################################
############################################################
## Compare intervals
@transform( sanitiseIntervals, suffix(".merged.cleaned.bed"), ".merged.overlap")
def pairwiseIntervals(infile, outfile):
    '''identify overlapping intervals for each pair of datasets'''

    to_cluster = True
    in_track = P.snip( os.path.basename( infile ), ".merged.cleaned.bed")
    statement = '''echo "track" > %(outfile)s; echo "overlap" >> %(outfile)s;'''
    P.run()

    for track in TRACKS:
       if (str(track) <> in_track):
           statement = '''echo %(track)s >> %(outfile)s; intersectBed -a %(infile)s -b intervals/%(track)s.merged.cleaned.bed -u | wc -l >> %(outfile)s; '''
           P.run()
    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

############################################################
@transform( pairwiseIntervals, suffix(".merged.overlap"), ".merged.overlap.load")
def loadPairwiseIntervals(infile, outfile):
    '''Load overlapping intervals into database '''
    track = P.snip( os.path.basename( infile ), ".overlap" ).replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_overlap
                   > %(outfile)s '''
    P.run()

############################################################
@transform( sanitiseIntervals, suffix(".merged.cleaned.bed"), ".merged.unique.bed")
def uniqueIntervals(infile, outfile):
    '''identify unique intervals for each dataset'''

    to_cluster = True
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name
    statement = '''cat %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip( os.path.basename( infile ),".merged.cleaned.bed")
    for track in TRACKS:
       if str(track) <> in_track: 
           statement = '''intersectBed -a %(outfile)s -b intervals/%(track)s.merged.cleaned.bed -v > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s ''' 
           P.run()

############################################################
@transform( uniqueIntervals, suffix(".merged.unique.bed"), ".merged.unique.bed.load")
def loadUniqueIntervals(infile, outfile):
    '''Load unique intervals into database '''
    track = P.snip( os.path.basename( infile ), ".unique.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop,interval_id,fold"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_unique_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()

############################################################
@transform( sanitiseIntervals, suffix(".merged.cleaned.bed"), ".merged.shared.bed")
def sharedIntervals(infile, outfile):
    '''identify shared intervals between datasets'''

    to_cluster = True
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name
    statement = '''cat %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip( os.path.basename( infile ), ".merged.cleaned.bed")
    for track in TRACKS:
       if str(track) != in_track:
           statement = '''intersectBed -a %(outfile)s -b intervals/%(track)s.merged.cleaned.bed -u > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s; ''' 
           P.run()

############################################################
@transform( sharedIntervals, suffix(".merged.shared.bed"), ".merged.shared.bed.load")
def loadSharedIntervals(infile, outfile):
    '''Load shared intervals into database '''
    track = P.snip( os.path.basename( infile ), ".shared.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop,interval_id,fold"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_shared_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()

############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.unique.bed" )
def uniqueReplicatedIntervals(infile, outfile):
    '''identify unique intervals for each dataset'''

    to_cluster = True
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name
    statement = '''cat %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip( os.path.basename( infile ),".replicated.bed")
    for track in EXPERIMENTS:
       track = str(track).replace("-agg","")
       if track <> in_track: 
           statement = '''intersectBed -a %(outfile)s -b replicated_intervals/%(track)s.replicated.bed -v > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s ''' 
           P.run()


###########################################################
@transform( uniqueReplicatedIntervals, suffix(".replicated.unique.bed"), ".replicated.unique.bed.load")
def loadUniqueReplicatedIntervals(infile, outfile):
    '''Load unique intervals into database '''
    track = P.snip( os.path.basename( infile ), ".unique.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop,interval_id,fold"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_unique_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()

############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.shared.bed" )
def sharedReplicatedIntervals(infile, outfile):
    '''identify shared intervals between datasets'''

    to_cluster = True
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name
    statement = '''cat %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip( os.path.basename( infile ), ".replicated.bed")
    for track in EXPERIMENTS:
       track = str(track).replace("-agg","")
       if track != in_track:
           statement = '''intersectBed -a %(outfile)s -b replicated_intervals/%(track)s.replicated.bed -u > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s; ''' 
           P.run()

############################################################
@transform( sharedReplicatedIntervals, suffix(".replicated.shared.bed"), ".replicated.shared.bed.load" )
def loadSharedReplicatedIntervals(infile, outfile):
    '''Load shared intervals into database '''
    track = P.snip( os.path.basename( infile ), ".shared.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop,interval_id"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_shared_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()
    
############################################################
@follows(exportReplicatedIntervalsAsBed)
@files( ("replicated_intervals/*liver*.replicated.bed", "replicated_intervals/*testes*.replicated.bed"), "replicated_intervals/liver.testes.venn" )
def liverTestesVenn(infiles, outfile):
    '''identify interval overlap between liver and testes. Merge intervals first.'''
    liver, testes = infiles
    liver_name = P.snip( os.path.basename(liver), ".replicated.bed" )
    testes_name = P.snip( os.path.basename(testes), ".replicated.bed" )
    to_cluster = True
    
    statement = '''cat %(liver)s %(testes)s | mergeBed -i stdin | awk 'OFS="\\t" {print $1,$2,$3,"CAPseq"NR}' > replicated_intervals/liver.testes.merge.bed;
                   echo "Total merged intervals" > %(outfile)s; cat replicated_intervals/liver.testes.merge.bed | wc -l >> %(outfile)s; 
                   echo "Liver & testes" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.merge.bed -b %(liver)s -u | intersectBed -a stdin -b %(testes)s -u > replicated_intervals/liver.testes.shared.bed; cat replicated_intervals/liver.testes.shared.bed | wc -l >> %(outfile)s; 
                   echo "Testes only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.merge.bed -b %(liver)s -v > replicated_intervals/%(testes_name)s.liver.testes.unique.bed; cat replicated_intervals/%(testes_name)s.liver.testes.unique.bed | wc -l >> %(outfile)s; 
                   echo "Liver only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.merge.bed -b %(testes)s -v > replicated_intervals/%(liver_name)s.liver.testes.unique.bed; cat replicated_intervals/%(liver_name)s.liver.testes.unique.bed | wc -l >> %(outfile)s;                   
                   sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()
    
############################################################    
@transform( liverTestesVenn, suffix(".venn"), ".venn.load" )
def loadLiverTestesVenn(infile, outfile):
    '''Load liver testes venn overlap into database '''
    header = "category,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=liver_testes_venn
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()
    
############################################################   
@follows(liverTestesVenn) 
@files( "replicated_intervals/liver.testes.shared.bed", "replicated_intervals/liver.testes.shared.bed.load" )
def loadLiverTestesShared(infile, outfile):
    '''Load liver testes shared intervals into database '''
    header = "contig,start,end,interval_id"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=liver_testes_shared_intervals
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################       
@follows(liverTestesVenn) 
@transform( "replicated_intervals/*.liver.testes.unique.bed", suffix(".liver.testes.unique.bed"),  ".liver.testes.unique.bed.load" )
def loadLiverTestesUnique(infile, outfile):
    '''Load liver testes unique intervals into database '''
    header = "contig,start,end,interval_id"
    track = P.snip(os.path.basename(infile), ".liver.testes.unique.bed")
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_liver_testes_unique_intervals
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()


############################################################       
@follows(liverTestesVenn) 
@files( "replicated_intervals/liver.testes.merge.bed", "replicated_intervals/liver.testes.merge.bed.load" )
def loadLiverTestesMerge(infile, outfile):
    '''Load liver testes shared intervals into database '''
    header = "contig,start,end,interval_id"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=liver_testes_merged_intervals
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
@follows(loadLiverTestesShared, loadLiverTestesUnique, loadLiverTestesMerge)
@merge( "replicated_intervals/*.liver.testes.unique.bed", "replicated_intervals/liver.testes.merge.sort.bed")
def exportLiverTestesMergeWithSort( infiles, outfile):
    ''' query database to produce a bed file which can be sorted by liver testes unique category and then length'''
    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    tracks = []
    for infile in infiles:
        t = P.snip(os.path.basename(infile), ".liver.testes.unique.bed").replace("-","_")
        tracks.append(t)


    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT m.contig, m.start, m.end, m.interval_id, 
               "sh_" || substr('000000' || (m.end-m.start), -6, 6)  as sort
               FROM liver_testes_merged_intervals m, liver_testes_shared_intervals s 
               WHERE m.interval_id=s.interval_id ''' % locals()
    for i, t in enumerate(tracks):
        query += '''UNION 
                    SELECT m.contig, m.start, m.end, m.interval_id, 
                    "a%(i)s_" || substr('000000' || (m.end-m.start), -6, 6)  as sort
                    FROM liver_testes_merged_intervals m, %(t)s_liver_testes_unique_intervals u%(i)s 
                    WHERE m.interval_id=u%(i)s.interval_id ''' % locals()
    print query
    cc.execute( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
            outs.write("%s%s" % (pre, str(r)) )
            pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()
    
############################################################
@follows( liverTestesVenn )
@files( "replicated_intervals/*.liver.testes.unique.bed", "replicated_intervals/liver.testes.chromatin.log" )    
def liverTestesUniqueChromatinProfile(infiles, outfile):
    '''plot chromatin mark profiles over tissue-specific CAPseq intervals'''
    chromatin = P.asList(PARAMS["bigwig_chromatin"])
    
    for infile in infiles:
        track = P.snip( os.path.basename(infile), ".liver.testes.unique.bed" )
        
        outtemp = P.getTempFile()
        tmpfilename = outtemp.name
    
        for bw in chromatin:
            chromatin_track = P.snip( os.path.basename(bw), ".bw" )
            ofp = "replicated_intervals/" + track + "." + chromatin_track + ".profile"    
            statement = '''cat %(infile)s | python %(scriptsdir)s/bed2gff.py --as-gtf | gzip > %(tmpfilename)s.gtf.gz;
                           python %(scriptsdir)s/bam2geneprofile.py 
                           --bamfile=%(bw)s 
                           --gtffile=%(tmpfilename)s.gtf.gz
                           --output-filename-pattern=%(ofp)s
                           --reporter=interval
                           --method=intervalprofile
                           --log=%(outfile)s'''
            P.run()      

############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".bed"), ".long.bed" )
def getLongIntervalBed( infile, outfile ):
    '''Generate bed file of top 500 longest intervals'''
    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename( infile ), ".replicated.bed" ).replace("-","_").replace(".","_")

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT i.contig, i.start, i.end, i.interval_id
               FROM %(track)s_replicated_intervals i
               ORDER BY length desc LIMIT 500;''' % locals()
    cc.execute( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
            outs.write("%s%s" % (pre, str(r)) )
            pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()
    
############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".bed"), ".long.genelist" )
def getLongIntervalGenes( infile, outfile ):
    '''Generate bed file of top 500 longest intervals'''
    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename( infile ), ".replicated.bed" ).replace("-","_").replace(".","_")
    cc = dbhandle.cursor()
    statement = "ATTACH DATABASE '%s' AS annotations; "  % (PARAMS["ensembl_annotation_database"])
    cc.execute(statement)

    # Extract data from db
    query = '''SELECT distinct t.gene_id
               FROM %(track)s_replicated_intervals i, %(track)s_replicated_tss s, %(track)s_replicated_ensembl_gene_overlap o, annotations.transcript_info t
               WHERE (substr(s.closest_id,1,18)=t.transcript_id
               or substr(s.closest_id,20,18)=t.transcript_id
               or substr(s.closest_id,39,18)=t.transcript_id
               or substr(s.closest_id,58,18)=t.transcript_id
               or substr(s.closest_id,77,18)=t.transcript_id)
               AND i.interval_id=s.gene_id
               AND o.gene_id=i.interval_id
               AND s.is_overlap > 0
               AND t.gene_biotype='protein_coding'
               AND i.length > 3000
               AND o.genes_pover2 > 20
               ORDER BY i.length desc
               LIMIT 500''' % locals()
    cc.execute( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
            outs.write("%s%s" % (pre, str(r)) )
            pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()
    
############################################################
@transform( getLongIntervalGenes, suffix(".genelist"), ".gtf.gz" )
def getLongIntervalGeneGTF( infile, outfile ):
    '''Filter Ensembl GTF file using list of Ensembl gene ids associated with long CAPseq intervals '''
    gene_file = os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_tss_profile"])
    statement = '''zcat %(gene_file)s | python %(scriptsdir)s/gtf2gtf.py --filter=gene --apply=%(infile)s | gzip > %(outfile)s; '''
    P.run()
    
    
############################################################
@transform( getLongIntervalGeneGTF, suffix(".gtf.gz"), ".chromatin_profile.log" )    
def longIntervalGeneChromatinProfile(infile, outfile):
    '''plot chromatin mark profiles over tissue-specific CAPseq intervals'''
    chromatin = P.asList(PARAMS["bigwig_chromatin"])
    
    track = P.snip( os.path.basename(infile), ".gtf.gz" )
        
    outtemp = P.getTempFile()
    tmpfilename = outtemp.name
    
    for bw in chromatin:
        chromatin_track = P.snip( os.path.basename(bw), ".bw" )
        ofp = "replicated_intervals/" + track + ".genes." + chromatin_track + ".profile"    
        statement = '''python %(scriptsdir)s/bam2geneprofile.py 
                       --bamfile=%(bw)s 
                       --gtffile=%(infile)s
                       --output-filename-pattern=%(ofp)s
                       --reporter=interval
                       --method=intervalprofile
                       --log=%(outfile)s'''
        P.run() 
            
    
############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".bed"), ".short.bed" )
def getShortIntervalBed( infile, outfile ):
    '''Generate bed file intervals < 2kb in length'''
    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename( infile ), ".replicated.bed" ).replace("-","_").replace(".","_")

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT i.contig, i.start, i.end, i.interval_id
               FROM %(track)s_replicated_intervals i
               WHERE length < 2000;''' % locals()
    cc.execute( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
            outs.write("%s%s" % (pre, str(r)) )
            pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()
    
############################################################
@merge( (getLongIntervalBed, getShortIntervalBed), "long_interval_chromatin_profile.log" )    
def longIntervalChromatinProfile(infiles, outfile):
    '''plot chromatin mark profiles over tissue-specific CAPseq intervals'''
    chromatin = P.asList(PARAMS["bigwig_chromatin"])
    
    for infile in infiles:
        track = P.snip( os.path.basename(infile), ".bed" )
        
        outtemp = P.getTempFile()
        tmpfilename = outtemp.name
    
        for bw in chromatin:
            chromatin_track = P.snip( os.path.basename(bw), ".bw" )
            ofp = "replicated_intervals/" + track + "." + chromatin_track + ".profile"    
            statement = '''cat %(infile)s | python %(scriptsdir)s/bed2gff.py --as-gtf | gzip > %(tmpfilename)s.gtf.gz;
                           python %(scriptsdir)s/bam2geneprofile.py 
                           --bamfile=%(bw)s 
                           --gtffile=%(tmpfilename)s.gtf.gz
                           --output-filename-pattern=%(ofp)s
                           --reporter=interval
                           --method=intervalprofile
                           --log=%(outfile)s'''
            P.run() 
            
############################################################
@transform( getMergedBigWig, suffix(".bw"), ".profile_long_interval_genes" )    
def longIntervalGeneCAPseqProfile(infiles, outfile):
    '''plot CAPseq profiles over long intervals'''
    
    track = P.snip( os.path.basename(infile), ".bed" )
    capseq = os.path.join ("replicated_intervals", track, ".long.gtf.gz")
        
    ofp = "replicated_intervals/" + track + "." + chromatin_track + ".profile"    
    statement = '''python %(scriptsdir)s/bam2geneprofile.py 
                           --bamfile=%(infile)s 
                           --gtffile=%(capseq)s
                           --output-filename-pattern=%(outfile)s
                           --reporter=interval
                           --method=intervalprofile
                           --log=%(outfile)s.log'''
    P.run() 


############################################################
############################################################
############################################################
## Compare peak shape across intervals
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.peakshape" )
def getPeakShape(infile, outfile):
    '''Cluster intervals based on peak shape '''
    track = P.snip( os.path.basename( infile ), ".replicated.bed" )
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]
    ofp = "replicated_intervals/" + track + ".replicated.peakshape"

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( fn )
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    bamfiles = " ".join( ("--bamfile=%s" % x) for x in samfiles )
    shifts =  " ".join( ("--shift=%s" % y)  for y in offsets )
    
    statement = '''python %(scriptsdir)s/bam2peakshape.py 
                       %(bamfiles)s 
                       --bedfile=%(infile)s
                       --output-filename-pattern=%(ofp)s
                       %(shifts)s
                       --sort=peak-width
                       --window-size=5000
                       --bin-size=10
                       --force
                       --log=%(outfile)s.log
                   -S %(outfile)s '''
    P.run()

############################################################    
@transform( uniqueReplicatedIntervals, suffix(".replicated.unique.bed"), ".replicated.unique.peakshape" )
def getPeakShapeTissueSpecificIntervals(infile, outfile):
    '''Cluster intervals based on peak shape '''
    track = P.snip( os.path.basename( infile ), ".replicated.unique.bed" )

    
    for expt in EXPERIMENTS:

        expt_name = str(expt).replace("-agg", "")
        replicates = EXPERIMENTS[expt]
        ofp = "replicated_intervals/" + expt_name + "." + track +".replicated.unique.peakshape"

        # setup files
        samfiles, offsets = [], []
        for t in replicates:
            fn = "bam/%s.norm.bam" % t.asFile()
            assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
            samfiles.append( fn )
            fn = "macs/with_input/%s.macs" % t.asFile()
            if os.path.exists( fn ):
                offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

        bamfiles = " ".join( ("--bamfile=%s" % x) for x in samfiles )
        shifts =  " ".join( ("--shift=%s" % y)  for y in offsets )
        
        statement = '''python %(scriptsdir)s/bam2peakshape.py 
                           %(bamfiles)s 
                           --bedfile=%(infile)s
                           --output-filename-pattern=%(ofp)s
                           %(shifts)s
                           --sort=peak-width
                           --window-size=5000
                           --bin-size=10
                           --force
                           --log=%(outfile)s.log
                       -S %(outfile)s '''
        P.run()
    
############################################################    
@follows( liverTestesVenn )
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".liver.testes.merge.peakshape.gz" )
def getPeakShapeLiverTestes(infile, outfile):
    '''Cluster intervals based on peak shape '''
    track = P.snip( os.path.basename( infile ), ".replicated.bed" )
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]
    #ofp = "replicated_intervals/" + track + ".liver.testes.merge.peakshape"
    bedfile = "replicated_intervals/liver.testes.merge.sort.bed"
    
    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( fn )
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    bamfiles = " ".join( ("--bamfile=%s" % x) for x in samfiles )
    shifts =  " ".join( ("--shift=%s" % y)  for y in offsets )
    
    statement = '''python %(scriptsdir)s/bam2peakshape.py 
                       %(bamfiles)s 
                       --bedfile=%(bedfile)s
                       --output-filename-pattern=%(outfile)s.%%s
                       %(shifts)s
                       --sort=peak-width
                       --sort=peak-height
                       --sort=interval-width
                       --sort=interval-score
                       --window-size=5000
                       --bin-size=10
                       --normalization=sum
                       --force
                       --log=%(outfile)s.log
                   | gzip
                   > %(outfile)s '''
    P.run()    
    
############################################################
############################################################
############################################################
## Compare replicated intervals with external bed files
@files( PARAMS["bed_ucsc_cgi"], "ucsc.bed.load")
def loadUCSCPredictedCGIIntervals(infile, outfile):
    '''load CGI intervals'''

    header = "contig,start,stop,id"
    statement = '''cat %(infile)s | awk 'OFS="\\t" {print $1,$2,$3,$4NR}' | python %(scriptsdir)s/csv2db.py
                      --table=cgi_intervals
                      --index=contig,start
                      --index=id
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".cgi_cap.bed")
def getCGIIntervals(infile, outfile):
    '''identify intervals overlapping CGI for each datasets'''

    CGI = PARAMS["bed_ucsc_cgi"]
    dataset_name =  P.snip( os.path.basename( CGI ), ".bed")
    statement = '''intersectBed -a %(infile)s -b %(CGI)s -u > %(outfile)s; '''
    P.run()

############################################################
@transform( getCGIIntervals, suffix(".cgi_cap.bed"), ".cgi_cap.bed.load")
def loadCGIIntervals(infile, outfile):
    '''Load intervals overlapping CGI into database '''
    track = P.snip( os.path.basename( infile ), ".cgi_cap.bed" ).replace(".cleaned","").replace(".","_").replace("-","_")
    header = "contig,start,stop,interval_id"
    statement = '''cat %(infile)s | awk 'OFS="\\t" {print $1,$2,$3,$4}' | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_predicted_cgi_and_cap
                      --index=contig,start
                      --index=interval_id
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".cap_only.bed")
def getNonCGIIntervals(infile, outfile):
    '''identify CApseq intervals not overlapping predicted CGI for each dataset'''

    CGI = PARAMS["bed_ucsc_cgi"]
    dataset_name =  P.snip( os.path.basename( CGI ), ".bed")
    statement = '''intersectBed -a %(infile)s -b %(CGI)s -v > %(outfile)s; '''
    P.run()

############################################################
@transform( getNonCGIIntervals, suffix(".cap_only.bed"), ".cap_only.bed.load")
def loadNonCGIIntervals(infile, outfile):
    '''Load intervals not overlapping CGI into database '''
    track = P.snip( os.path.basename( infile ), ".cap_only.bed" ).replace(".cleaned","").replace(".","_").replace("-","_")
    header = "contig,start,stop,interval_id"
    statement = '''cat %(infile)s | awk 'OFS="\\t" {print $1,$2,$3,$4}' | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_cap_not_predicted_cgi
                      --index=contig,start
                      --index=interval_id
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".cgi_only.bed")
def getPredictedCGIIntervals(infile, outfile):
    '''identify predicted CGI intervals not overlapping CAPseq intervals for each dataset'''
    CGI = PARAMS["bed_ucsc_cgi"]
    #dataset_name =  P.snip( os.path.basename( CGI ), ".bed")
    statement = '''cat %(CGI)s | awk 'OFS="\\t" {print $1,$2,$3,$4NR}' | intersectBed -a stdin -b %(infile)s -v > %(outfile)s; '''
    P.run()

############################################################
@transform( getPredictedCGIIntervals, suffix(".cgi_only.bed"), ".cgi_only.bed.load")
def loadPredictedCGIIntervals(infile, outfile):
    '''Load predicted CGI intervals not overlapping CAP-seq intervals into database '''
    track = P.snip( os.path.basename( infile ), ".cgi_only.bed" ).replace(".cleaned","").replace(".replicated","").replace(".merged","")
    replicated = False
    if infile.find("replicated") > -1:
        replicated = True
    table = P.snip( os.path.basename( infile ), ".cgi_only.bed" ).replace(".cleaned","").replace(".","_").replace("-","_")

    if replicated:
        expt_track = track + "-agg"
        replicates = EXPERIMENTS[expt_track]
    else:
        replicates = [track, ]

    # Write header to output file
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    headers = ( "contig","start","stop","interval_id","nPeaks","PeakCenter","Length","AvgVal","PeakVal","nProbes" )
    tmpfile.write( "\t".join(headers) + "\n" )
    contig,start,end,interval_id,npeaks,peakcenter,length,avgval,peakval,nprobes = "",0,0,0,0,0,0,0,0,0

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( pysam.Samfile( fn,  "rb" ) )
        fn = "macs/with_input/%s.macs" % t
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    # Loop over input Bed file and calculate stats for merged intervals
    c = E.Counter()
    for line in open(infile, "r"):
        c.input += 1
        contig, start, end, interval_id = line[:-1].split()[:4]
        start, end = int(start), int(end)
        #interval_id = c.input
        npeaks, peakcenter, length, avgval, peakval, nprobes = PIntervals.countPeaks( contig, start, end, samfiles, offsets )
        if nprobes == 0:
            c.skipped_reads += 1
        c.output += 1
        tmpfile.write( "\t".join( map( str, (contig,start,end,interval_id,npeaks,peakcenter,length,avgval,peakval,nprobes) )) + "\n" )
    tmpfile.close()
    tmpfilename = tmpfile.name
    tablename = "%s_predicted_cgi_not_cap" % table
    
    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                       --index=contig,start 
                       --table=%(tablename)s
                       --allow-empty
                   < %(tmpfilename)s > %(outfile)s '''
    P.run()
    os.unlink( tmpfile.name )
    L.info( "%s\n" % str(c) )

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".cgi_overlap")
def getCGIOverlapCount(infile, outfile):
    '''identify intervals overlapping CGI for each datasets'''

    CGI = P.asList(PARAMS["bed_cgi"])

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    for island in CGI:
       dataset_name =  P.snip( os.path.basename( island ), ".bed")
       statement = '''echo %(dataset_name)s >> %(outfile)s; intersectBed -a %(infile)s -b %(island)s -u | wc -l >> %(outfile)s; '''
       P.run()
    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

############################################################
@transform( getCGIOverlapCount, suffix(".cgi_overlap"), ".cgi_overlap.load")
def loadCGIOverlapCount(infile, outfile):
    '''Load intervals overlapping CGI into database '''
    track = P.snip( os.path.basename( infile ), ".cgi_overlap" ).replace(".cleaned","").replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_cgi
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform((sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".chipseq")
def getChipseqOverlap(infile, outfile):
    '''identify intervals overlapping chipseq intervals for each datasets'''

    chipseq = P.asList(PARAMS["bed_chipseq"])

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    if len(chipseq[0]) > 0:
        for tf in chipseq:
           dataset_name =  P.snip( os.path.basename( tf ), ".bed")
           statement = '''echo %(dataset_name)s >> %(outfile)s; intersectBed -a %(infile)s -b %(tf)s -u | wc -l >> %(outfile)s; '''
           P.run()
        statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
        P.run()
    else:
        statement = '''touch %(outfile)s '''
        P.run()

############################################################
@transform( getChipseqOverlap, suffix(".chipseq"), ".chipseq.load")
def loadChipseqIntervals(infile, outfile):
    '''Load intervals overlapping chipseq into database '''
    track = P.snip( os.path.basename( infile ), ".chipseq" ).replace(".cleaned","").replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_chipseq
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".capseq")
def getCapseqOverlap(infile, outfile):
    '''identify intervals overlapping capseq intervals for each datasets'''

    capseq = P.asList(PARAMS["bed_capseq"])

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    if len(capseq[0]) > 0:
        for x in capseq:
            dataset_name =  P.snip( os.path.basename( x ), ".bed")
            statement = '''echo %(dataset_name)s >> %(outfile)s; intersectBed -a %(infile)s -b %(x)s -u | wc -l >> %(outfile)s; '''
            P.run()
        statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
        P.run()
    else:
        statement = '''touch %(outfile)s '''
        P.run()

############################################################
@transform( getCapseqOverlap, suffix(".capseq"), ".capseq.load")
def loadCapseqIntervals(infile, outfile):
    '''Load intervals overlapping capseq into database '''
    track = P.snip( os.path.basename( infile ), ".capseq" ).replace(".cleaned","").replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_capseq
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".chromatin")
def getChromatinMarkOverlap(infile, outfile):
    '''identify intervals overlapping chromatin mark intervals for each datasets'''

    chromatin = P.asList(PARAMS["bed_chromatin"])

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    if len(chromatin[0]) > 0:
        for mark in chromatin:
           dataset_name =  P.snip( os.path.basename( mark ), ".bed")
           statement = '''echo %(dataset_name)s >> %(outfile)s; intersectBed -a %(infile)s -b %(mark)s -u | wc -l >> %(outfile)s; '''
           P.run()
        statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
        P.run()
    else:
        statement = '''touch %(outfile)s '''
        P.run()

############################################################
@transform( getChromatinMarkOverlap, suffix(".chromatin"), ".chromatin.load")
def loadChromatinMarkIntervals(infile, outfile):
    '''Load intervals overlapping chromatin marks into database '''
    track = P.snip( os.path.basename( infile ), ".chromatin" ).replace(".cleaned","").replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_chromatin
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@files( ( os.path.join( PARAMS["ensembl_annotation_dir"],PARAMS["ensembl_annotation_transcript_tss"] ), PARAMS["bed_ucsc_cgi"]), "overlap.tss.cgi" )
def getCGIEnsemblTranscriptTSSOverlap( infiles, outfile ):
    '''Establish overlap between cgi and tss intervals'''
    tss, cgi = infiles
    tss_extend = PARAMS["ensembl_annotation_tss_extend"]
    to_cluster = True

    outtemp1 = P.getTempFile()
    tmpfilename1 = outtemp1.name

    if os.path.exists( outfile):
        statement = '''rm %s''' % outfile
        P.run()

    statement = """zcat %(tss)s | slopBed -i stdin -g %(samtools_genome)s.fai -b %(tss_extend)s > tss_extended.bed; """
    statement += """echo "Predicted CGIs overlapping 1 or more TSS" >> %(outfile)s; intersectBed -a %(cgi)s -b tss_extended.bed -u | wc -l >> %(outfile)s; """
    statement += """echo "Predicted CGIs not overlapping any TSS" >> %(outfile)s; intersectBed -a %(cgi)s -b tss_extended.bed -v | wc -l >> %(outfile)s; """
    statement += """echo "TSS overlapped by 1 or more CGI" >> %(outfile)s; intersectBed -a tss_extended.bed -b %(cgi)s -u | wc -l >> %(outfile)s; """
    statement += """echo "TSS not overlapped by any predicted CGI" >> %(outfile)s; intersectBed -a tss_extended.bed -b %(cgi)s -v | wc -l >> %(outfile)s; """

    statement += '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()

############################################################
@transform( getCGIEnsemblTranscriptTSSOverlap, regex(r"overlap.tss.cgi"), r"overlap.tss.cgi.load")
def loadCGIEnsemblTranscriptTSSOverlap(infile, outfile):
    '''load TSS CGI overlap into database'''

    header = "track,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                        --table=cgi_ensembl_transcript_tss_venn
                        --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
@files( ( os.path.join( PARAMS["ensembl_annotation_dir"],PARAMS["ensembl_annotation_gene_tss_interval"] ), PARAMS["bed_ucsc_cgi"]), "cgi.ensembl_gene_tss.overlap" )
def getCGIEnsemblGeneTSSOverlap( infiles, outfile ):
    '''Establish overlap between cgi and tss intervals'''
    tss, cgi = infiles
    to_cluster = True

    outtemp1 = P.getTempFile()
    tmpfilename1 = outtemp1.name

    if os.path.exists( outfile):
        statement = '''rm %s''' % outfile
        P.run()

    statement = """echo "Predicted CGIs overlapping 1 or more TSS" >> %(outfile)s; intersectBed -a %(cgi)s -b %(tss)s -u | wc -l >> %(outfile)s; """
    statement += """echo "Predicted CGIs not overlapping any TSS" >> %(outfile)s; intersectBed -a %(cgi)s -b %(tss)s -v | wc -l >> %(outfile)s; """
    statement += """echo "TSS overlapped by 1 or more CGI" >> %(outfile)s; intersectBed -a %(tss)s -b %(cgi)s -u | wc -l >> %(outfile)s; """
    statement += """echo "TSS not overlapped by any predicted CGI" >> %(outfile)s; intersectBed -a %(tss)s -b %(cgi)s -v | wc -l >> %(outfile)s; """
    statement += '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()

############################################################
@transform( getCGIEnsemblGeneTSSOverlap, regex(r"cgi.ensembl_gene_tss.overlap"), r"cgi.ensembl_gene_tss.overlap.load")
def loadCGIEnsemblGeneTSSOverlap(infile, outfile):
    '''load TSS CGI overlap into database'''

    header = "track,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                        --table=cgi_ensembl_gene_tss_venn
                        --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
@merge( "external_bed/*.bed", "external_interval_sets.stats" )
def getExternalBedStats(infiles, outfile):
    '''Calculate statistics for external bed files '''
    chromatin = P.asList(PARAMS["bed_chromatin"])
    capseq = P.asList(PARAMS["bed_capseq"])
    chipseq = P.asList(PARAMS["bed_chipseq"])
    CGI = P.asList(PARAMS["bed_cgi"])
    extBed = chromatin + capseq + chipseq + CGI

    if os.path.exists(outfile):
        statement = '''rm %(outfile)s'''
        P.run()

    for f in extBed:
        if len(f) > 0:
            track = P.snip( os.path.basename(f),".bed" )
            statement = """echo '%(track)s' >> %(outfile)s; cat %(f)s | wc -l >> %(outfile)s; """
            P.run()
    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

############################################################
@transform( getExternalBedStats, regex(r"(\S+).stats"), r"\1.load" )
def loadExternalBedStats(infile, outfile):
    '''Load statistics for external bed files into database '''
    statement = """cat %(infile)s | python %(scriptsdir)s/csv2db.py 
                         --header=bed,intervals
                         --table=external_interval_sets 
                    > %(outfile)s"""
    P.run()

############################################################
############################################################
############################################################
## Compare intervals to external bed files using GAT
@follows( mkdir("gat") )
@files(PARAMS["samtools_genome"]+".fai", "gat/"+PARAMS["genome"]+".bed.gz")
def buildGATWorkspace(infile, outfile ):
    '''Build genomic workspace file for GAT '''
    statement = '''cat %(infile)s | awk 'OFS="\\t" {print $1,0,$2,"workspace"}' | gzip > %(outfile)s '''
    P.run()

############################################################
@follows( buildGATWorkspace )
@merge( (sanitiseIntervals,exportReplicatedIntervalsAsBed), "gat/external_dataset_gat.tsv" )
def runExternalDatasetGAT(infiles, outfile):
    '''Run genome association tester on bed files '''

    to_cluster = True
    segfiles = ""
    for x in infiles:
        track = P.snip(os.path.basename(x), ".bed")
        statement = """cat %(x)s | awk 'OFS="\\t" {print $1,$2,$3,"%(track)s"}' > gat/%(track)s.bed; """
        P.run()
        segfiles += " --segment-file=gat/%s.bed " % track 

    # External datasets
    chromatin = P.asList(PARAMS["bed_chromatin"])
    capseq = P.asList(PARAMS["bed_capseq"])
    chipseq = P.asList(PARAMS["bed_chipseq"])
    CGI = P.asList(PARAMS["bed_cgi"])
    extBed = chromatin + capseq + chipseq + CGI
    annofiles = " ".join( [ "--annotation-file=%s" % x for x in extBed ] )
    statement = """gatrun.py %(segfiles)s %(annofiles)s --workspace=gat/%(genome)s.bed.gz --num-samples=1000 --nbuckets=120000 --force > %(outfile)s"""
    P.run()

############################################################
@transform( runExternalDatasetGAT, suffix(".tsv"), ".load" )
def loadExternalDatasetGAT(infile, outfile):
    '''Load genome association tester results into database '''
    statement = """cat %(infile)s | grep -v "^#" | python %(scriptsdir)s/csv2db.py 
                         --table=external_dataset_gat_results gat/external_dataset_gat.tsv
                    > %(outfile)s"""
    P.run()


############################################################
############################################################
############################################################
## ANNOTATE INTERVAL GENOMIC LOCATION (ENSEMBL)
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".ensembl_transcript_overlap" )
def annotateEnsemblTranscriptOverlap( infile, outfile ):
    '''classify intervals according to their base pair overlap with respect to different genomic features (genes, TSS, upstream/downstream flanks) '''
    to_cluster = True

    feature_list = P.asList( PARAMS["ensembl_annotation_features_transcript"] )

    outfiles = ""
    first = True
    for feature in feature_list:
        feature_name = P.snip( os.path.basename( feature ), ".gtf" ).replace(".","_")
        outfiles += " %(outfile)s.%(feature_name)s " % locals()
        if first:
            cut_command = "cut -f1,4,5,6,8 "
            first = False
        else:
            cut_command = "cut -f4,5,6 "
        statement = """
                cat < %(infile)s 
                | python %(scriptsdir)s/bed2gff.py --as-gtf 
                | python %(scriptsdir)s/gtf2table.py 
		                --counter=overlap  
		                --counter=length  
		                --log=%(outfile)s.log 
		                --filename-gff=%(feature)s 
		                --genome-file=%(genome_dir)s/%(genome)s
                | %(cut_command)s 
                | sed s/nover/%(feature_name)s_nover/g 
                | sed s/pover/%(feature_name)s_pover/g 
                | sed s/min/length/
                > %(outfile)s.%(feature_name)s"""
        P.run()

    # Paste output together
    statement = '''paste  %(outfiles)s > %(outfile)s'''
    P.run()

############################################################
@transform( annotateEnsemblTranscriptOverlap, suffix(".ensembl_transcript_overlap"), ".ensembl_transcript_overlap.load" )
def loadEnsemblTranscriptOverlap( infile, outfile ):
    '''load interval annotations: genome architecture '''
    track= P.snip( os.path.basename(infile), ".ensembl_transcript_overlap").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_ensembl_transcript_overlap
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".ensembl_gene_overlap" )
def annotateEnsemblGeneOverlap( infile, outfile ):
    '''classify intervals according to their base pair overlap with respect to different genomic features (genes, TSS, upstream/downstream flanks) '''
    to_cluster = True

    feature_list = P.asList( PARAMS["ensembl_annotation_features_gene"] )

    outfiles = ""
    first = True
    for feature in feature_list:
        feature_name = P.snip( os.path.basename( feature ), ".gtf" ).replace(".","_")
        outfiles += " %(outfile)s.%(feature_name)s " % locals()
        if first:
            cut_command = "cut -f1,4,5,6,8 "
            first = False
        else:
            cut_command = "cut -f4,5,6 "
        statement = """
                cat < %(infile)s 
                | python %(scriptsdir)s/bed2gff.py --as-gtf 
                | python %(scriptsdir)s/gtf2table.py 
		                --counter=overlap  
		                --counter=length  
		                --log=%(outfile)s.log 
		                --filename-gff=%(feature)s 
		                --genome-file=%(genome_dir)s/%(genome)s
                | %(cut_command)s 
                | sed s/nover/%(feature_name)s_nover/g 
                | sed s/pover/%(feature_name)s_pover/g 
                | sed s/min/length/
                > %(outfile)s.%(feature_name)s"""
        P.run()

    # Paste output together
    statement = '''paste  %(outfiles)s > %(outfile)s'''
    P.run()

############################################################
@transform( annotateEnsemblGeneOverlap, suffix(".ensembl_gene_overlap"), ".ensembl_gene_overlap.load" )
def loadEnsemblGeneOverlap( infile, outfile ):
    '''load interval annotations: genome architecture '''
    track= P.snip( os.path.basename(infile), ".ensembl_gene_overlap").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_ensembl_gene_overlap
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@files( PARAMS["bed_ucsc_cgi"], "cgi.ensembl_transcript" )
def annotateCGIEnsemblTranscriptOverlap( infile, outfile ):
    '''classify predicted CGI intervals according to their base pair overlap 
       with respect to different genomic features (genes, TSS, upstream/downstream flanks) '''
    to_cluster = True

    feature_list = P.asList( PARAMS["ensembl_annotation_features_transcript"] )

    outfiles = ""
    first = True
    for feature in feature_list:
        feature_name = P.snip( os.path.basename( feature ), ".gtf" ).replace(".","_")
        outfiles += " %(outfile)s.%(feature_name)s " % locals()
        if first:
            cut_command = "cut -f1,4,5,6,8 "
            first = False
        else:
            cut_command = "cut -f4,5,6 "
        statement = """
                cat %(infile)s 
                | awk '{print $1,$2,$3,$4-NR}'
                | python %(scriptsdir)s/bed2gff.py --as-gtf 
                | python %(scriptsdir)s/gtf2table.py 
		                --counter=overlap  
		                --counter=length  
		                --log=%(outfile)s.log 
		                --filename-gff=%(feature)s 
		                --genome-file=%(genome_dir)s/%(genome)s
                | %(cut_command)s 
                | sed s/nover/%(feature_name)s_nover/g 
                | sed s/pover/%(feature_name)s_pover/g 
                | sed s/min/length/
                > %(outfile)s.%(feature_name)s"""
        P.run()

    # Paste output together
    statement = '''paste  %(outfiles)s > %(outfile)s'''
    P.run()

############################################################
@transform( annotateCGIEnsemblTranscriptOverlap, suffix(".ensembl_transcript"), ".ensembl_transcript.load" )
def loadCGIEnsemblTranscriptOverlap( infile, outfile ):
    '''load interval annotations: genome architecture '''
    track= P.snip( os.path.basename(infile), ".ensembl_transcript").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_ensembl_transcript_overlap 
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@files( PARAMS["bed_ucsc_cgi"], "cgi.ensembl_gene" )
def annotateCGIEnsemblGeneOverlap( infile, outfile ):
    '''classify predicted CGI intervals according to their base pair overlap 
       with respect to different genomic features (genes, TSS, upstream/downstream flanks) '''
    to_cluster = True

    feature_list = P.asList( PARAMS["ensembl_annotation_features_gene"] )

    outfiles = ""
    first = True
    for feature in feature_list:
        feature_name = P.snip( os.path.basename( feature ), ".gtf" ).replace(".","_")
        outfiles += " %(outfile)s.%(feature_name)s " % locals()
        if first:
            cut_command = "cut -f1,4,5,6,8 "
            first = False
        else:
            cut_command = "cut -f4,5,6 "
        statement = """
                cat %(infile)s 
                | awk '{print $1,$2,$3,$4-NR}'
                | python %(scriptsdir)s/bed2gff.py --as-gtf 
                | python %(scriptsdir)s/gtf2table.py 
		                --counter=overlap  
		                --counter=length  
		                --log=%(outfile)s.log 
		                --filename-gff=%(feature)s 
		                --genome-file=%(genome_dir)s/%(genome)s
                | %(cut_command)s 
                | sed s/nover/%(feature_name)s_nover/g 
                | sed s/pover/%(feature_name)s_pover/g 
                | sed s/min/length/
                > %(outfile)s.%(feature_name)s"""
        P.run()

    # Paste output together
    statement = '''paste  %(outfiles)s > %(outfile)s'''
    P.run()

############################################################
@transform( annotateCGIEnsemblGeneOverlap, suffix(".ensembl_gene"), ".ensembl_gene.load" )
def loadCGIEnsemblGeneOverlap( infile, outfile ):
    '''load interval annotations: genome architecture '''
    track= P.snip( os.path.basename(infile), ".ensembl_gene").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_ensembl_gene_overlap 
                         --index=gene_id
                 > %(outfile)s; """
    P.run()


############################################################
############################################################
############################################################
## Non-coding transcript overlap
@transform( loadEnsemblTranscriptOverlap, suffix(".replicated.ensembl_transcript_overlap.load"), ".replicated.intergenic.bed" )
def getIntergenicCapseqBed( infile, outfile ):
    '''Make bed file of all CAPseq peaks that do not overlap genes/TSS'''
    to_cluster = True
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track= P.snip( os.path.basename(infile), ".replicated.ensembl_transcript_overlap.load")
    table= track.replace(".","_").replace("-","_")

    cc = dbhandle.cursor()
    statement =  '''SELECT i.contig, i.start, i.end, i.interval_id
                    FROM %(table)s_replicated_ensembl_transcript_overlap o, %(table)s_replicated_intervals i
                    WHERE tss_extended_pover1=0
                    AND genes_pover1=0
                    AND downstream_flank_pover1=0
                    AND upstream_flank_pover1=0
                    AND i.interval_id=o.gene_id
                    order by contig, start''' % locals()
    cc.execute( statement )

    # Write to bed file
    outs = open( outfile, "w")
    for result in cc:
        contig, start, stop, interval_id = result
        outs.write( "%s\t%i\t%i\t%i\n" % (contig, start, stop, interval_id) )
    cc.close()
    outs.close()

############################################################
@transform( (getIntergenicCapseqBed, exportReplicatedIntervalsAsBed), suffix(".bed"), ".noncoding.tss" )
def getEnsemblNoncodingTSSDistance( infile, outfile ):
    '''Get distance of intergenic CAPseq peaks to nearest non-coding transcript'''
    to_cluster = True

    annotation_file = os.path.join( PARAMS["ensembl_annotation_dir"],
                                    PARAMS["ensembl_annotation_noncoding_tss"] )

    statement = """cat < %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
	                 | python %(scriptsdir)s/gtf2table.py 
		                   --counter=distance-tss 
		                   --log=%(outfile)s.log 
                       --filename-gff=%(annotation_file)s 
                       --filename-format="bed" 
                   > %(outfile)s"""
    P.run()

############################################################
@transform( getEnsemblNoncodingTSSDistance, suffix( ".tss"), ".tss.load" )
def loadEnsemblNoncodingTSSDistance( infile, outfile ):
    '''load interval annotations: distance to non-coding transcription start sites '''
    track= P.snip( os.path.basename(infile), ".noncoding.tss").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_noncoding 
                         --index=gene_id
                         --index=closest_id
                         --index=id5
                         --index=id3
                 > %(outfile)s; """
    P.run()

############################################################
############################################################
############################################################
## Ensembl TTS/TTS overlap/distance
############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.ensembl_gene_tss.overlap" )
def getCapseqEnsemblGeneTSSOverlap( infile, outfile ):
    '''Establish overlap between capseq and tss intervals'''
    tss = os.path.join( PARAMS["ensembl_annotation_dir"],PARAMS["ensembl_annotation_gene_tss_interval"] )
    to_cluster = True

    outtemp1 = P.getTempFile()
    tmpfilename1 = outtemp1.name

    if os.path.exists( outfile):
        statement = '''rm %s''' % outfile
        P.run()

    statement = """echo "CAPseq intervals overlapping 1 or more TSS" >> %(outfile)s; intersectBed -a %(infile)s -b %(tss)s -u | wc -l >> %(outfile)s; """
    statement += """echo "CAPseq intervals not overlapping any TSS" >> %(outfile)s; intersectBed -a %(infile)s -b %(tss)s -v | wc -l >> %(outfile)s; """
    statement += """echo "TSS overlapped by 1 or more CAPseq interval" >> %(outfile)s; intersectBed -a %(tss)s -b %(infile)s -u | wc -l >> %(outfile)s; """
    statement += """echo "TSS not overlapped by any CAPseq intervals" >> %(outfile)s; intersectBed -a %(tss)s -b %(infile)s -v | wc -l >> %(outfile)s; """
    statement += '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()

############################################################
@transform( getCapseqEnsemblGeneTSSOverlap, suffix("replicated.ensembl_gene_tss.overlap"), "replicated.ensembl_gene_tss.overlap.load")
def loadCapseqEnsemblGeneTSSOverlap(infile, outfile):
    '''load TSS Capseq overlap into database'''

    header = "track,intervals"
    track = P.snip( os.path.basename( infile), ".ensembl_gene_tss.overlap" )
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                        --table=%(track)s_ensembl_gene_tss_venn
                        --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".tss" )
def annotateEnsemblTranscriptTSS( infile, outfile ):
    '''compute distance to TSS'''

    to_cluster = True
    annotation_file = os.path.join( PARAMS["ensembl_annotation_dir"],
                                    PARAMS["ensembl_annotation_transcript_tss"] )

    statement = """cat < %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
	                 | python %(scriptsdir)s/gtf2table.py 
		                   --counter=distance-tss 
		                   --log=%(outfile)s.log 
                       --filename-gff=%(annotation_file)s 
                       --filename-format="bed" 
                   > %(outfile)s"""
    P.run()

############################################################
@transform( annotateEnsemblTranscriptTSS, suffix( ".tss"), ".tss.load" )
def loadEnsemblTranscriptTSS( infile, outfile ):
    '''load interval annotations: distance to transcription start sites '''
    track= P.snip( os.path.basename(infile), ".tss").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_tss 
                         --index=gene_id
                         --index=closest_id
                         --index=id5
                         --index=id3
                 > %(outfile)s; """
    P.run()
    
############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".gene.tss" )
def annotateEnsemblGeneTSS( infile, outfile ):
    '''compute distance to TSS for single TSS per gene (longest transcript)'''

    to_cluster = True

    annotation_file = os.path.join( PARAMS["ensembl_annotation_dir"],
                                    PARAMS["ensembl_annotation_gene_tss"] )

    statement = """cat < %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
	                 | python %(scriptsdir)s/gtf2table.py 
		                   --counter=distance-tss 
		                   --log=%(outfile)s.log 
                       --filename-gff=%(annotation_file)s 
                       --filename-format="bed" 
                   > %(outfile)s"""
    P.run()

############################################################
@transform( annotateEnsemblGeneTSS, suffix( ".gene.tss"), ".gene.tss.load" )
def loadEnsemblGeneTSS( infile, outfile ):
    '''load interval annotations: distance to transcription start sites '''
    track= P.snip( os.path.basename(infile), ".gene.tss").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_gene_tss 
                         --index=gene_id
                         --index=closest_id
                         --index=id5
                         --index=id3
                 > %(outfile)s; """
    P.run()
        
############################################################
@transform( loadEnsemblTranscriptTSS, suffix( ".replicated.tss.load"), ".replicated.tss.bed" )
def exportCapseqEnsemblTSSBed( infile, outfile ):
    '''export bed file of all CAPseq intervals within 1kb of annotated Ensembl transcript TSS '''
    track= P.snip( os.path.basename(infile), ".replicated.tss.load").replace(".cleaned","").replace(".","_").replace("-","_")

    dbhandle = sqlite3.connect( PARAMS["database"] )

    cc = dbhandle.cursor()
    statement = '''SELECT i.contig, i.start, i.end, i.interval_id 
                   FROM %s_replicated_intervals i, %s_replicated_tss t
                   WHERE i.interval_id=t.gene_id 
                   AND t. closest_dist < 1000
                   ORDER by contig, start''' % (track, track)
    cc.execute( statement )

    # Write to bed file
    outs = open( outfile, "w")
    for result in cc:
        contig, start, stop, interval_id = result
        outs.write( "%s\t%i\t%i\t%i\n" % (contig, start, stop, interval_id) )
    cc.close()
    outs.close()
    
############################################################
@follows(exportReplicatedIntervalsAsBed, exportCapseqEnsemblTSSBed)
@files( ("replicated_intervals/*liver*.replicated.tss.bed", "replicated_intervals/*testes*.replicated.tss.bed"), "replicated_intervals/liver.testes.tss.venn" )
def liverTestesTSSVenn(infiles, outfile):
    '''identify interval overlap between liver and testes for TSS associated intervals. Merge intervals first.'''
    liver, testes = infiles
    to_cluster = True
    
    statement = '''cat %(liver)s %(testes)s | mergeBed -i stdin > replicated_intervals/liver.testes.tss.merge.bed;
                   echo "Total merged intervals" > %(outfile)s; cat replicated_intervals/liver.testes.tss.merge.bed | wc -l >> %(outfile)s; 
                   echo "Liver & testes" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.tss.merge.bed -b %(liver)s -u | intersectBed -a stdin -b %(testes)s -u | wc -l >> %(outfile)s; 
                   echo "Testes only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.tss.merge.bed -b %(liver)s -v | wc -l >> %(outfile)s; 
                   echo "Liver only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.tss.merge.bed -b %(testes)s -v | wc -l >> %(outfile)s;                   
                   sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()
    
############################################################    
@transform( liverTestesTSSVenn, suffix(".venn"), ".venn.load" )
def loadLiverTestesTSSVenn(infile, outfile):
    '''Load liver testes venn overlap into database '''
    header = "category,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=liver_testes_tss_venn
                      --header=%(header)s
                   > %(outfile)s '''
    P.run() 
    
############################################################
@transform( loadEnsemblTranscriptTSS, suffix( ".replicated.tss.load"), ".replicated.no.tss.bed" )
def exportCapseqNoEnsemblTSSBed( infile, outfile ):
    '''export bed file of all CAPseq intervals not within 1kb of annotated Ensembl transcript TSS '''
    track= P.snip( os.path.basename(infile), ".replicated.tss.load").replace(".cleaned","").replace(".","_").replace("-","_")

    dbhandle = sqlite3.connect( PARAMS["database"] )

    cc = dbhandle.cursor()
    statement = '''SELECT i.contig, i.start, i.end, i.interval_id 
                   FROM %s_replicated_intervals i, %s_replicated_tss t
                   WHERE i.interval_id=t.gene_id 
                   AND t. closest_dist >= 1000
                   ORDER by contig, start''' % (track, track)
    cc.execute( statement )

    # Write to bed file
    outs = open( outfile, "w")
    for result in cc:
        contig, start, stop, interval_id = result
        outs.write( "%s\t%i\t%i\t%i\n" % (contig, start, stop, interval_id) )
    cc.close()
    outs.close()
    
############################################################
@follows(exportReplicatedIntervalsAsBed, exportCapseqNoEnsemblTSSBed)
@files( ("replicated_intervals/*liver*.replicated.no.tss.bed", "replicated_intervals/*testes*.replicated.no.tss.bed"), "replicated_intervals/liver.testes.intergenic.venn" )
def liverTestesIntergenicVenn(infiles, outfile):
    '''identify interval overlap between liver and testes for non-TSS associated intervals. Merge intervals first.'''
    liver, testes = infiles
    to_cluster = True
    
    statement = '''cat %(liver)s %(testes)s | mergeBed -i stdin > replicated_intervals/liver.testes.intergenic.merge.bed;
                   echo "Total merged intervals" > %(outfile)s; cat replicated_intervals/liver.testes.intergenic.merge.bed | wc -l >> %(outfile)s; 
                   echo "Liver & testes" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.intergenic.merge.bed -b %(liver)s -u | intersectBed -a stdin -b %(testes)s -u | wc -l >> %(outfile)s; 
                   echo "Testes only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.intergenic.merge.bed -b %(liver)s -v | wc -l >> %(outfile)s; 
                   echo "Liver only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.intergenic.merge.bed -b %(testes)s -v | wc -l >> %(outfile)s;                   
                   sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()
    
############################################################    
@transform( liverTestesIntergenicVenn, suffix(".venn"), ".venn.load" )
def loadLiverTestesIntergenicVenn(infile, outfile):
    '''Load liver testes venn overlap into database '''
    header = "category,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=liver_testes_intergenic_venn
                      --header=%(header)s
                   > %(outfile)s '''
    P.run() 

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".tts" )
def annotateEnsemblTranscriptTTS( infile, outfile ):
    '''Compute distance to TTS'''
    to_cluster = USECLUSTER
    annotation_file = os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_transcript_tts"] )
    statement = """cat < %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
	                 | python %(scriptsdir)s/gtf2table.py 
		                   --counter=distance-tss 
		                   --log=%(outfile)s.log 
                       --filename-gff=%(annotation_file)s 
                       --filename-format="bed" 
                   > %(outfile)s"""
    P.run()

############################################################
@transform( annotateEnsemblTranscriptTTS, suffix( ".tts"), ".tts.load" )
def loadEnsemblTranscriptTTS( infile, outfile ):
    '''load interval annotations: distance to transcription termination sites '''
    track= P.snip( os.path.basename(infile), ".tts").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_tts 
                         --index=gene_id
                         --index=closest_id
                         --index=id5
                         --index=id3
                 > %(outfile)s; """
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".gene.tts" )
def annotateEnsemblGeneTTS( infile, outfile ):
    '''compute distance to TTS'''
    to_cluster = True
    annotation_file = os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_gene_tts"] )
    statement = """cat < %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
	                 | python %(scriptsdir)s/gtf2table.py 
		                   --counter=distance-tss 
		                   --log=%(outfile)s.log 
                       --filename-gff=%(annotation_file)s 
                       --filename-format="bed" 
                   > %(outfile)s"""
    P.run()

############################################################
@transform( annotateEnsemblGeneTTS, suffix( ".gene.tts"), ".gene.tts.load" )
def loadEnsemblGeneTTS( infile, outfile ):
    '''load interval annotations: distance to transcription termination sites '''
    track= P.snip( os.path.basename(infile), ".tts").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_gene_tts 
                         --index=gene_id
                         --index=closest_id
                         --index=id5
                         --index=id3
                 > %(outfile)s; """
    P.run()

############################################################
############################################################
############################################################
## TSS profile
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.transcript.tss-profile.png" )
def getReplicatedEnsemblTranscriptTSSProfile(infile, outfile):
    '''Build TSS profile from BAM files'''
    to_cluster = USECLUSTER

    track = P.snip( os.path.basename(infile), ".replicated.bed" )
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]
    tss_file = os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_tss_profile"])
    ofp = "replicated_intervals/" + track + ".replicated.transcript.tss-profile"

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( fn )
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    bamfiles = " ".join( ("--bamfile=%s" % x) for x in samfiles )
    shifts =  " ".join( ("--shift=%s" % y)  for y in offsets )

    statement = '''python %(scriptsdir)s/bam2geneprofile.py 
                       %(bamfiles)s 
                       --gtffile=%(tss_file)s
                       --output-filename-pattern=%(ofp)s
                       %(shifts)s
                       --reporter=transcript
                       --method=tssprofile
                       --normalization=total-sum
                       --normalize-profile=area
                       --normalize-profile=counts
                       --normalize-profile=none'''
    P.run()

############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.transcript.tss-profile.capseq.png" )
def getReplicatedEnsemblTranscriptTSSProfileCapseq(infile,outfile):
    '''Build TSS profile from BAM files'''
    to_cluster = USECLUSTER

    track = P.snip( os.path.basename(infile), ".replicated.bed" )
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]
    ofp = "replicated_intervals/" + track + ".replicated.transcript.tss-profile.capseq"

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( fn )
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    bamfiles = " ".join( ("--bamfile=%s" % x) for x in samfiles )
    shifts =  " ".join( ("--shift=%s" % y)  for y in offsets )

    gene_file = os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_tss_profile"])
    tss_file = os.path.join( PARAMS["ensembl_annotation_dir"],  PARAMS["ensembl_annotation_transcript_tss"])

    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name

    statement = '''intersectBed -a %(tss_file)s -b %(infile)s -u | cut -f4 > %(tmpfilename)s; 
                   zcat %(gene_file)s | python %(scriptsdir)s/gtf2gtf.py --filter=transcript --apply=%(tmpfilename)s | gzip > %(tmpfilename)s.gtf.gz; 
                   python %(scriptsdir)s/bam2geneprofile.py 
                       %(bamfiles)s 
                       --gtffile=%(tmpfilename)s.gtf.gz
                       --output-filename-pattern=%(ofp)s
                       %(shifts)s
                       --reporter=transcript
                       --method=tssprofile
                       --perInterval
                       --normalization=interval'''
    P.run()

############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.transcript.tss-profile.nocapseq.png" )
def getReplicatedEnsemblTranscriptTSSProfileNoCapseq(infile,outfile):
    '''Build TSS profile from BAM files'''
    to_cluster = USECLUSTER

    track = P.snip( os.path.basename(infile), ".replicated.bed" )
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]
    ofp = "replicated_intervals/" + track + ".replicated.transcript.tss-profile.nocapseq"

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( fn )
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    bamfiles = " ".join( ("--bamfile=%s" % x) for x in samfiles )
    shifts =  " ".join( ("--shift=%s" % y)  for y in offsets )

    gene_file = os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_tss_profile"])
    tss_file = os.path.join( PARAMS["ensembl_annotation_dir"],  PARAMS["ensembl_annotation_transcript_tss"])

    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name

    statement = '''intersectBed -a %(tss_file)s -b %(infile)s -v | cut -f4 > %(tmpfilename)s; 
                   zcat %(gene_file)s | python %(scriptsdir)s/gtf2gtf.py --filter=transcript --apply=%(tmpfilename)s | gzip > %(tmpfilename)s.gtf.gz; 
                   python %(scriptsdir)s/bam2geneprofile.py 
                       %(bamfiles)s 
                       --gtffile=%(tmpfilename)s.gtf.gz
                       --output-filename-pattern=%(ofp)s
                       %(shifts)s
                       --reporter=transcript
                       --method=tssprofile
                       --perInterval
                       --normalization=interval'''
    P.run()

############################################################
############################################################
## Per gene
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.gene.tss-profile.png" )
def getReplicatedEnsemblGeneTSSProfile(infile, outfile):
    '''Build TSS profile from BAM files'''
    to_cluster = USECLUSTER

    track = P.snip( os.path.basename(infile), ".replicated.bed" )
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]
    tss_file = os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_tss_profile"])
    ofp = "replicated_intervals/" + track + ".replicated.gene.tss-profile"

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( fn )
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    bamfiles = " ".join( ("--bamfile=%s" % x) for x in samfiles )
    shifts =  " ".join( ("--shift=%s" % y)  for y in offsets )

    statement = '''python %(scriptsdir)s/bam2geneprofile.py 
                       %(bamfiles)s 
                       --gtffile=%(tss_file)s
                       --output-filename-pattern=%(ofp)s
                       %(shifts)s
                       --perInterval
                       --reporter=gene
                       --method=tssprofile
                       --normalization=interval'''
    P.run()

############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.gene.tss-profile.capseq.png" )
def getReplicatedEnsemblGeneTSSProfileCapseq(infile,outfile):
    '''Build TSS profile from BAM files'''
    to_cluster = USECLUSTER

    track = P.snip( os.path.basename(infile), ".replicated.bed" )
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]
    ofp = "replicated_intervals/" + track + ".replicated.gene.tss-profile.capseq"

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( fn )
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    bamfiles = " ".join( ("--bamfile=%s" % x) for x in samfiles )
    shifts =  " ".join( ("--shift=%s" % y)  for y in offsets )

    gene_file = os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_tss_profile"])
    tss_file = os.path.join( PARAMS["ensembl_annotation_dir"],  PARAMS["ensembl_annotation_transcript_tss"])

    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name

    statement = '''intersectBed -a %(tss_file)s -b %(infile)s -u | cut -f4 > %(tmpfilename)s; 
                   zcat %(gene_file)s | python %(scriptsdir)s/gtf2gtf.py --filter=transcript --apply=%(tmpfilename)s | gzip > %(tmpfilename)s.gtf.gz; 
                   python %(scriptsdir)s/bam2geneprofile.py 
                       %(bamfiles)s 
                       --gtffile=%(tmpfilename)s.gtf.gz
                       --output-filename-pattern=%(ofp)s
                       %(shifts)s
                       --reporter=gene
                       --method=tssprofile
                       --perInterval
                       --normalization=interval'''
    P.run()

############################################################
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.gene.tss-profile.nocapseq.png" )
def getReplicatedEnsemblGeneTSSProfileNoCapseq(infile,outfile):
    '''Build TSS profile from BAM files'''
    to_cluster = USECLUSTER

    track = P.snip( os.path.basename(infile), ".replicated.bed" )
    expt_track = track + "-agg"
    replicates = EXPERIMENTS[expt_track]
    ofp = "replicated_intervals/" + track + ".replicated.gene.tss-profile.nocapseq"

    # setup files
    samfiles, offsets = [], []
    for t in replicates:
        fn = "bam/%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( fn )
        fn = "macs/with_input/%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShiftFromMacs( fn ) )

    bamfiles = " ".join( ("--bamfile=%s" % x) for x in samfiles )
    shifts =  " ".join( ("--shift=%s" % y)  for y in offsets )

    gene_file = os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_tss_profile"])
    tss_file = os.path.join( PARAMS["ensembl_annotation_dir"],  PARAMS["ensembl_annotation_transcript_tss"])

    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name

    statement = '''intersectBed -a %(tss_file)s -b %(infile)s -v | cut -f4 > %(tmpfilename)s; 
                   zcat %(gene_file)s | python %(scriptsdir)s/gtf2gtf.py --filter=transcript --apply=%(tmpfilename)s | gzip > %(tmpfilename)s.gtf.gz; 
                   python %(scriptsdir)s/bam2geneprofile.py 
                       %(bamfiles)s 
                       --gtffile=%(tmpfilename)s.gtf.gz
                       --output-filename-pattern=%(ofp)s
                       %(shifts)s
                       --reporter=gene
                       --method=tssprofile
                       --perInterval
                       --normalization=interval'''
    P.run()

############################################################
############################################################
############################################################
## Genomic repeats
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".repeats" )
def annotateRepeats( infile, outfile ):
    '''count the overlap between intervals and repeats.'''
    to_cluster = True
    annotation_file = os.path.join( PARAMS["ensembl_annotation_dir"],PARAMS_ANNOTATIONS["interface_repeats_gff"] )
    statement = """cat < %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=overlap
                      	 --log=%(outfile)s.log
                         --filename-gff=%(annotation_file)s
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s"""
    P.run()

############################################################
@transform( annotateRepeats, suffix(".repeats"), ".repeats.load" )
def loadRepeats( infile, outfile ):
    '''load interval annotations: repeats'''

    track= P.snip( os.path.basename(infile), ".repeats").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_repeats 
                         --index=gene_id
                         --allow-empty
                 > %(outfile)s; """
    P.run()

############################################################
############################################################
############################################################
## RNAseq/CAGE based annotation
@transform( exportReplicatedIntervalsAsBed, suffix(".replicated.bed"), ".replicated.rnaseq_tss.overlap" )
def getCapseqRNAseqTSSOverlap( infile, outfile ):
    '''Establish overlap between capseq and tss intervals'''
    tss = os.path.join( PARAMS["rnaseq_annotation_dir"],PARAMS["rnaseq_annotation_tss_extended"] )
    to_cluster = True

    outtemp1 = P.getTempFile()
    tmpfilename1 = outtemp1.name

    if os.path.exists( outfile):
        statement = '''rm %s''' % outfile
        P.run()

    statement = """echo "CAPseq intervals overlapping 1 or more RNAseq TSS" >> %(outfile)s; intersectBed -a %(infile)s -b %(tss)s -u | wc -l >> %(outfile)s; """
    statement += """echo "CAPseq intervals not overlapping any RNAseq TSS" >> %(outfile)s; intersectBed -a %(infile)s -b %(tss)s -v | wc -l >> %(outfile)s; """
    statement += """echo "RNAseq TSS overlapped by 1 or more CAPseq interval" >> %(outfile)s; intersectBed -a %(tss)s -b %(infile)s -u | wc -l >> %(outfile)s; """
    statement += """echo "RNAseq TSS not overlapped by any CAPseq intervals" >> %(outfile)s; intersectBed -a %(tss)s -b %(infile)s -v | wc -l >> %(outfile)s; """
    statement += '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()
    
############################################################
@transform( getCapseqRNAseqTSSOverlap, suffix("replicated.rnaseq_tss.overlap"), "replicated.rnaseq_tss.overlap.load")
def loadCapseqRNAseqTSSOverlap(infile, outfile):
    '''load TSS Capseq overlap into database'''

    header = "track,intervals"
    track = P.snip( os.path.basename( infile), ".rnaseq_tss.overlap" )
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                        --table=%(track)s_rnaseq_tss_venn
                        --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## ANNOTATE INTERVAL NUCLEOTIDE COMPOSITION
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".composition" )
def annotateCapseqComposition( infile, outfile ):
    '''Establish the nucleotide composition of intervals'''

    to_cluster = True
    statement = """cat %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( annotateCapseqComposition, suffix( ".composition"), ".composition.load" )
def loadCapseqComposition( infile, outfile ):
    '''Load the nucleotide composition of intervals'''

    track= P.snip( os.path.basename(infile), ".composition").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_composition 
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".control.composition" )
def annotateControlComposition( infile, outfile ):
    '''Establish the nucleotide composition of control intervals'''

    to_cluster = True
    track= P.snip( os.path.basename(infile), ".bed")
    dirname= os.path.dirname(infile)

    statement = """cat %(infile)s | python %(scriptsdir)s/bed2bed.py -m shift -g %(genome_dir)s/%(genome)s --offset=-10000 -S %(dirname)s/%(track)s.control.bed;
                   cat %(dirname)s/%(track)s.control.bed
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s.log
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( annotateControlComposition, suffix( ".control.composition"), ".control.composition.load" )
def loadControlComposition( infile, outfile ):
    '''Load the nucleotide composition of intervals'''

    track= P.snip( os.path.basename(infile), ".control.composition").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_composition_control
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".flanking5.composition" )
def annotateFlankingCompositionLeft( infile, outfile ):
    '''Establish the nucleotide composition of intervals immediately upstream'''

    to_cluster = True
    track= P.snip( os.path.basename(infile), ".bed")
    dirname= os.path.dirname(infile)
    flank_size = PARAMS["intervals_flank_size"]

    # Exclude intervals with length < 100bp
    statement = """flankBed -i %(infile)s -l %(flank_size)s -r 0 -g %(samtools_genome)s.fai 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s -L %(dirname)s/%(track)s.flanking5.log 
                   | awk 'OFS="\\t" {if ($3-$2>100) print $1,$2,$3,$4}' > %(dirname)s/%(track)s.flanking5.bed;
                   cat %(dirname)s/%(track)s.flanking5.bed
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s.log
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( annotateFlankingCompositionLeft, suffix( ".flanking5.composition"), ".flanking5.composition.load" )
def loadFlankingCompositionLeft( infile, outfile ):
    '''Load the nucleotide composition of regions flanking intervals'''

    track= P.snip( os.path.basename(infile), ".flanking5.composition").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_composition_flanking5
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@transform( (sanitiseIntervals,exportReplicatedIntervalsAsBed), suffix(".bed"), ".flanking3.composition" )
def annotateFlankingCompositionRight( infile, outfile ):
    '''Establish the nucleotide composition of intervals immediately downstream'''

    to_cluster = True
    track= P.snip( os.path.basename(infile), ".bed")
    dirname= os.path.dirname(infile)
    flank_size = PARAMS["intervals_flank_size"]

    # Exclude intervals with length < 100bp
    statement = """flankBed -i %(infile)s -l 0 -r 1000 -g %(samtools_genome)s.fai 
                   | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s -L %(dirname)s/%(track)s.flanking3.log 
                   | awk 'OFS="\\t" {if ($3-$2>100) print $1,$2,$3,$4}' > %(dirname)s/%(track)s.flanking3.bed;
                   cat %(dirname)s/%(track)s.flanking3.bed
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s.log
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( annotateFlankingCompositionRight, suffix( ".flanking3.composition"), ".flanking3.composition.load" )
def loadFlankingCompositionRight( infile, outfile ):
    '''Load the nucleotide composition of regions flanking intervals'''
    track= P.snip( os.path.basename(infile), ".flanking3.composition").replace(".cleaned","").replace(".","_").replace("-","_")
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_composition_flanking3
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@files( os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_transcript_tss"] ), "tss.transcript.composition" )
def annotateEnsemblTranscriptTSSComposition( infile, outfile ):
    '''Establish the nucleotide composition of tss intervals'''
    to_cluster = True
    tss_extend = PARAMS["ensembl_annotation_tss_extend"]
    statement = """zcat %(infile)s 
                   | slopBed -i stdin -g %(samtools_genome)s.fai -b %(tss_extend)s
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s.log
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( annotateEnsemblTranscriptTSSComposition, suffix( ".composition"), ".composition.load" )
def loadEnsemblTranscriptTSSComposition( infile, outfile ):
    '''Load the nucleotide composition of tss intervals'''
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=tss_transcript_comp
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@files( os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_gene_tss"] ), "tss.gene.composition" )
def annotateEnsemblGeneTSSComposition( infile, outfile ):
    '''Establish the nucleotide composition of tss intervals'''
    to_cluster = True
    tss_extend = PARAMS["ensembl_annotation_tss_extend"]
    statement = """zcat %(infile)s 
                   | slopBed -i stdin -g %(samtools_genome)s.fai -b %(tss_extend)s
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s.log
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( annotateEnsemblGeneTSSComposition, suffix( ".composition"), ".composition.load" )
def loadEnsemblGeneTSSComposition( infile, outfile ):
    '''Load the nucleotide composition of tss intervals'''
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=tss_gene_comp
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@files( os.path.join( PARAMS["ensembl_annotation_dir"], PARAMS["ensembl_annotation_gene_tss_interval"] ), "tss.gene.interval.composition" )
def annotateEnsemblGeneTSSIntervalComposition( infile, outfile ):
    '''Establish the nucleotide composition of tss intervals'''
    to_cluster = True
    statement = """zcat %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s.log
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( annotateEnsemblGeneTSSIntervalComposition, suffix( ".composition"), ".composition.load" )
def loadEnsemblGeneTSSIntervalComposition( infile, outfile ):
    '''Load the nucleotide composition of tss intervals'''
    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=tss_gene_interval_comp
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@files( PARAMS["bed_ucsc_cgi"], "cgi.composition" )
def annotateCGIComposition( infile, outfile ):
    '''Establish the nucleotide composition of CGI intervals'''
    to_cluster = True

    # Give each row a unique identifier
    statement = """cat %(infile)s 
                   | awk '{print $1,$2,$3,$4NR}'
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                         --counter=composition-cpg
                      	 --log=%(outfile)s.log
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( annotateCGIComposition, suffix( ".composition"), ".composition.load" )
def loadCGIComposition( infile, outfile ):
    '''Load the nucleotide composition of CGI intervals'''

    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=cgi_comp
                         --index=gene_id
                 > %(outfile)s; """
    P.run()


############################################################
############################################################
############################################################
## Compare intervals genomic features using GAT
@follows(buildGATWorkspace)
@merge( exportReplicatedIntervalsAsBed, "gat/genomic_features_gat.tsv" )
def runEnsemblGenomicFeaturesGAT(infiles, outfile):
    '''Run genome association tester on bed files '''

    to_cluster = True
    segfiles = ""
    for x in infiles:
        track = P.snip(os.path.basename(x), ".replicated.bed")
        statement = """cat %(x)s | awk 'OFS="\\t" {print $1,$2,$3,"%(track)s"}' > gat/%(track)s.bed; """
        P.run()
        segfiles += " --segment-file=gat/%s.bed " % track 
    annofile = PARAMS["gat_genomic_features_file"]
    annopath = os.path.dirname(annofile)
    annotrack = P.snip(os.path.basename(annofile), ".gff.gz")
    # Convert annotation file to bed
    statement = """zcat %(annofile)s | python %(scriptsdir)s/gff2bed.py --name='feature' --is-gtf -S gat/%(annotrack)s.bed; """
    P.run()
    statement = """gatrun.py %(segfiles)s --annotation-file=gat/%(annotrack)s.bed --workspace=gat/%(genome)s.bed.gz --num-samples=1000 --force --nbuckets=120000 > %(outfile)s"""
    P.run()

############################################################
@transform( runEnsemblGenomicFeaturesGAT, regex(r"gat/(\S+).tsv"), r"\1.load" )
def loadEnsemblGenomicFeaturesGAT(infile, outfile):
    '''Load genome association tester results into database '''
    statement = """cat %(infile)s | grep -v "^#" | python %(scriptsdir)s/csv2db.py 
                       --table=gat_genomic_features_results 
                   > %(outfile)s"""
    P.run()

############################################################
############################################################
############################################################
## EXPORT
@transform( normaliseBAMs, regex("(\S+)/bam/(\S+).norm.bam"), r"%s/\1.bigwig" % PARAMS["exportdir"])
def exportBigwig( infile, outfile ):
    '''convert BAM to bigwig file.'''

    # no bedToBigBed on the 32 bit cluster
    to_cluster = True
    statement = '''python %(scriptsdir)s/bam2wiggle.py
                       --output-format=bigwig 
                       --output-filename=%(outfile)s 
                   %(infile)s > %(outfile)s.log '''
    P.run()


############################################################
@transform( "*/macs/*.macs_model.pdf", regex( r"(\S+)/macs/(\S+)(.macs_model.pdf)"), r"%s/MACS/\2.macs_model.pdf" % PARAMS["exportdir"])
def exportMacsModel( infile, outfile ):
    '''copy MACS model files to export directory.'''
    try: os.mkdir( PARAMS["exportdir"] )
    except OSError: pass
    try: os.mkdir( '''%s/MACS''' % PARAMS["exportdir"] )
    except OSError: pass
    to_cluster = False
    statement = '''cp -p %(infile)s %(outfile)s'''
    P.run()

############################################################
@transform( annotateEnsemblTranscriptTSS, suffix(".replicated.tss"), r".capseq_tss_overlap")
def exportCapseqTSSOverlapGeneList( infile, outfile):
    '''Export list of genes where one or more transcript TSS is overlapped by a replicated CAPseq interval'''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename( infile ), ".replicated.tss" ).replace("-","_").replace(".","_")

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT closest_id FROM %(track)s_replicated_tss where closest_dist=0;''' % locals()
    cc.execute( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        ids = str(result[0])
        genes = ids.split(",")
        for g in genes:
            outs.write( "%s\n" % g )
    cc.close()
    outs.close()

############################################################
@transform( annotateEnsemblTranscriptTSS, suffix(".replicated.tss"), ".capseq_tss_1kb")
def exportCapseqTSS1kbGeneList( infile, outfile):
    '''Export list of genes where one or more transcript TSS is within 1kb of a replicated CAPseq interval'''

    max_gene_dist = 1000

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename( infile ), ".replicated.tss" ).replace("-","_").replace(".","_")

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT closest_id FROM %(track)s_replicated_tss where closest_dist < %(max_gene_dist)s ORDER BY closest_dist;''' % locals()
    cc.execute( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        ids = str(result[0])
        genes = ids.split(",")
        for g in genes:
            outs.write( "%s\n" % g )
    cc.close()
    outs.close()

############################################################
############################################################
@transform( loadCapseqComposition, suffix(".replicated.composition.load"), ".replicated.gc.export" )
def exportCapseqGCProfiles( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename( infile ), ".replicated.composition.load" ).replace("-","_").replace(".","_")

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.pGC, cc.pGC, c3.pGC, c5.pGC 
               FROM %(track)s_replicated_composition c
               left join %(track)s_replicated_composition_control cc on c.gene_id=cc.gene_id
               left join %(track)s_replicated_composition_flanking3 c3 on c.gene_id=c3.gene_id
               left join %(track)s_replicated_composition_flanking5 c5 on c.gene_id=c5.gene_id;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()

############################################################
@transform( loadCGIComposition, suffix("cgi.composition.load"), "cgi.gc.export" )
def exportCGIGCProfiles( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.pGC
               FROM cgi_comp c;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()
    
############################################################
@transform( loadUniqueReplicatedIntervals, suffix(".replicated.unique.bed.load"), ".replicated.unique.genes.export" )
def exportTissueSpecificCAPseqGenes( infile, outfile ):
    '''Export tissue specific CAPseq genes '''

    track = P.snip( os.path.basename( infile ), ".replicated.unique.bed.load" ).replace("-","_").replace(".","_")
    
    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    
    statement = "ATTACH DATABASE '%s' AS annotations; "  % (PARAMS["ensembl_annotation_database"])
    cc.execute(statement)

    # Extract data from db
    query = '''SELECT a.gene_id
               FROM %(track)s_replicated_unique_intervals u, %(track)s_replicated_tss t, annotations.transcript_info a
               WHERE u.interval_id=t.gene_id
               AND t.closest_dist < 1000
               AND t.closest_id=a.transcript_id
               AND a.gene_biotype='protein_coding';''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    outs.write("gene_id\n")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()
    
############################################################
@transform( exportTissueSpecificCAPseqGenes, suffix(".replicated.unique.genes.export"), ".replicated.unique.genes.go" )
def runGOOnGeneLists( infile, outfile ):
    PipelineGO.runGOFromFiles( outfile = outfile,
                               outdir = "go",
                               fg_file = infile,
                               bg_file = None,
                               go_file = os.path.join(PARAMS["ensembl_annotation_dir"], PARAMS_ANNOTATIONS["interface_go"] ),
                               ontology_file = os.path.join(PARAMS["ensembl_annotation_dir"], PARAMS_ANNOTATIONS["interface_go_obo"] ),
                               minimum_counts = PARAMS["go_minimum_counts"] )

############################################################
@transform( exportTissueSpecificCAPseqGenes, suffix(".replicated.unique.genes.export"), ".replicated.unique.genes.goslim" )
def runGOSlimOnGeneLists( infile, outfile ):
    PipelineGO.runGOFromFiles( outfile = outfile,
                               outdir = "go",
                               fg_file = infile,
                               bg_file = None,
                               go_file = os.path.join(PARAMS["ensembl_annotation_dir"], PARAMS_ANNOTATIONS["interface_goslim"] ),
                               ontology_file = os.path.join(PARAMS["ensembl_annotation_dir"], PARAMS_ANNOTATIONS["interface_goslim_obo"]),
                               minimum_counts = PARAMS["go_minimum_counts"] )

############################################################
@transform( (runGOOnGeneLists, runGOSlimOnGeneLists), regex( "(.*)"), r"\1.load" )
def loadGOResults( infile, outfile ):
    '''load GO results.'''
    P.load( os.path.join( "go", "all.biol_process.fold"), outfile )

############################################################
@transform( loadEnsemblTranscriptTSSComposition, suffix("tss.transcript.composition.load"), "tss.transcript.gc.export" )
def exportEnsemblTranscriptTSSGCProfiles( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.pGC
               FROM tss_transcript_comp c;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()

############################################################
@transform( loadEnsemblGeneTSSComposition, suffix("tss.gene.composition.load"), "tss.gene.gc.export" )
def exportEnsemblGeneTSSGCProfiles( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.pGC
               FROM tss_gene_comp c;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()

############################################################
@transform( loadCapseqComposition, suffix(".replicated.composition.load"), ".replicated.cpg.export" )
def exportCapseqCpGObsExp( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename( infile ), ".replicated.composition.load" ).replace("-","_").replace(".","_")

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.CpG_ObsExp, cc.CpG_ObsExp, c3.CpG_ObsExp, c5.CpG_ObsExp 
               FROM %(track)s_replicated_composition c
               left join %(track)s_replicated_composition_control cc on c.gene_id=cc.gene_id
               left join %(track)s_replicated_composition_flanking3 c3 on c.gene_id=c3.gene_id
               left join %(track)s_replicated_composition_flanking5 c5 on c.gene_id=c5.gene_id;''' % locals()
    cc.execute( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
            outs.write("%s%s" % (pre, str(r)) )
            pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()

############################################################
@transform( loadCGIComposition, suffix("cgi.composition.load"), "cgi.cpg.export" )
def exportCGICpGObsExp( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.CpG_ObsExp
               FROM cgi_comp c;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()

############################################################
@transform( loadEnsemblTranscriptTSSComposition, suffix("tss.transcript.composition.load"), "tss.transcript.cpg.export" )
def exportEnsemblTranscriptTSSCpGObsExp( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.CpG_ObsExp
               FROM tss_transcript_comp c;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()

############################################################
@transform( loadEnsemblGeneTSSComposition, suffix("tss.gene.composition.load"), "tss.gene.cpg.export" )
def exportEnsemblGeneTSSCpGObsExp( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.CpG_ObsExp
               FROM tss_gene_comp c;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()



############################################################
@transform( loadCapseqComposition, suffix(".replicated.composition.load"), ".replicated.cpg_density.export" )
def exportCapseqCpGDensity( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename( infile ), ".replicated.composition.load" ).replace("-","_").replace(".","_")

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.pCpG, cc.pCpG, c3.pCpG, c5.pCpG 
               FROM %(track)s_replicated_composition c
               left join %(track)s_replicated_composition_control cc on c.gene_id=cc.gene_id
               left join %(track)s_replicated_composition_flanking3 c3 on c.gene_id=c3.gene_id
               left join %(track)s_replicated_composition_flanking5 c5 on c.gene_id=c5.gene_id;''' % locals()
    cc.execute( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
            outs.write("%s%s" % (pre, str(r)) )
            pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()

############################################################
@transform( loadCGIComposition, suffix("cgi.composition.load"), "cgi.cpg_density.export" )
def exportCGICpGDensity( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.pCpG
               FROM cgi_comp c;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()

############################################################
@transform( loadEnsemblTranscriptTSSComposition, suffix("tss.transcript.composition.load"), "tss.transcript.cpg_density.export" )
def exportEnsemblTranscriptTSSCpGDensity( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.pCpG
               FROM tss_transcript_comp c;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()

############################################################
@transform( loadEnsemblGeneTSSComposition, suffix("tss.gene.composition.load"), "tss.gene.cpg_density.export" )
def exportEnsemblGeneTSSCpGDensity( infile, outfile ):
    '''Export file of GC content '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT c.gene_id, c.pCpG
               FROM tss_gene_comp c;''' % locals()
    cc.execute( query )
    E.info( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        pre = ""
        for r in result:
          outs.write("%s%s" % (pre, str(r)) )
          pre = "\t"
        outs.write("\n")

    cc.close()
    outs.close()
    
############################################################    
@transform( (exportReplicatedIntervalsAsBed, sharedReplicatedIntervals, uniqueReplicatedIntervals), suffix(".bed"), ".length" )
def exportCAPseqReplicatedLength( infile, outfile ):
    '''Export length of CAPseq intervals'''
    statement = '''cat %(infile)s | awk '{print $3-$2}' > %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################
## REPORTS
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

@files( "report.log", "publish.log")
def publish_report(infile, outfile):
    '''Link bed, bam, wig and report files to web '''
    publish_dir = PARAMS["publish_dir"]
    species = PARAMS["genome"]
    report_dir = os.path.join(publish_dir, species)
    bam_dir = os.path.join(publish_dir, "bam")
    bed_dir = os.path.join(publish_dir, "bed")
    wig_dir = os.path.join(publish_dir, "wig")  
    tss_dir = os.path.join(publish_dir, "tss")
    tss_dist_dir = os.path.join(publish_dir, "tss_distance")
    gc_dir = os.path.join(publish_dir, "gc")
    cgi_dir = os.path.join(publish_dir, "cpg")
    cpg_density_dir = os.path.join(publish_dir, "cpg_density")
    length_dir = os.path.join(publish_dir, "length")
    chromatin_dir = os.path.join(publish_dir, "chromatin")
    #report_dir = PARAMS["publish_report"]   
    print(infile)
    working_dir = os.getcwd()
    statement =  '''cp -rf report/html/* %(report_dir)s > %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/bam/*.norm.bam %(bam_dir)s >> %(outfile)s;'''
    statement += '''cp -sn %(working_dir)s/macs/with_input/*/*/*.wig.gz %(wig_dir)s >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.replicated.bed %(bed_dir)s >> %(outfile)s;'''    
    statement += '''cp -sn %(working_dir)s/intervals/*solo*.bed %(bed_dir)s/no_input >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/intervals/*.merged.cleaned.bed %(bed_dir)s/replicates >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.replicated.unique.bed %(bed_dir)s/tissue_specific >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.tss-profile*.tsv.gz %(tss_dir)s >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.replicated.gc.export %(gc_dir)s >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/tss.gene.gc.export %(gc_dir)s/%(species)s.tss.gene.gc.export >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/tss.transcript.gc.export %(gc_dir)s/%(species)s.tss.transcript.gc.export >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/cgi.gc.export %(gc_dir)s/%(species)s.cgi.gc.export >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.replicated.cpg.export %(cgi_dir)s >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/tss.gene.cpg.export %(cgi_dir)s/%(species)s.tss.gene.cpg.export >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/tss.transcript.cpg.export %(cgi_dir)s/%(species)s.tss.transcript.cpg.export >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/cgi.cpg.export %(cgi_dir)s/%(species)s.cgi.cpg.export >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.replicated.cpg_density.export %(cpg_density_dir)s >> %(outfile)s; '''   
    statement += '''cp -sn %(working_dir)s/tss.gene.cpg_density.export %(cpg_density_dir)s/%(species)s.tss.gene.cpg_density.export >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/tss.transcript.cpg_density.export %(cpg_density_dir)s/%(species)s.tss.transcript.cpg_density.export >> %(outfile)s; '''
    statement += '''cp -sn %(working_dir)s/cgi.cpg_density.export %(cpg_density_dir)s/%(species)s.cgi.cpg_density.export >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.length %(length_dir)s >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.gene.tss %(tss_dist_dir)s >> %(outfile)s; '''     
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.profile.tsv.gz %(chromatin_dir)s >> %(outfile)s; ''' 
    statement += '''cp -sn %(working_dir)s/replicated_intervals/*.liver.testes.unique.bed %(bed_dir)s/liver_testes >> %(outfile)s; ''' 
        
    P.run()

############################################################
############################################################
############################################################
## Replication
@files("pipeline.ini", "conf.py")
def replicateMappingAndCalling(infile, outfile):
    '''Link BAM files and MACS results from previous run'''
    # must copy and modify pipeline.ini first
    prev_run = PARAMS["previous_run"]
    statement = '''cp -ps %(prev_run)s/*.fastq.gz .;'''
    statement += '''cp %(prev_run)s/conf.py %(prev_run)s/sphinxreport.ini .; '''
    statement += '''mkdir bam; cp -psr %(prev_run)s/bam/* bam/; '''
    statement += '''mkdir macs; cp -psr %(prev_run)s/macs/* macs/; '''
    statement += '''mkdir external_bed; cp %(prev_run)s/external_bed/*.bed external_bed/;'''
    P.run()

############################################################
############################################################
############################################################
## Pipeline organisation
@follows( buildBAM, 
          buildPicardAlignStats, loadPicardAlignStats,
          buildBAMStats, loadBAMStats,
          dedup, loadPicardDuplicateStats)
def mapReads():
    '''Align reads to target genome.'''
    pass

@follows( removeUnmapped, normaliseBAMs )
def NormaliseBAMFiles():
    '''Normalise BAMS'''
    pass

@follows( runMACS, loadMACS, 
          summarizeMACS, loadMACSSummary, 
          exportIntervalsAsBed )
def buildIntervalsMacs():
    '''Find peaks using MACS using input as control'''
    pass

@follows( runMACSsolo, loadMACSsolo, 
          summarizeMACSsolo, loadMACSsoloSummary, 
          exportIntervalsAsBedsolo )
def buildIntervalsMacsNoControl():
    '''Find peaks using MACS without control sample'''
    pass

@follows( mergeIntervals, loadMergedIntervals )
def mergePeaks():
    '''Merge nearby intervals within a single track'''
    pass

@follows( getBackground, loadBackground )
def background():
    '''Assess level of background (non-peak) binding'''
    pass

@follows( thresholdFoldChange, loadFoldChangeThreshold,
          sharedIntervalsFoldChangeThreshold,
          loadSharedIntervalsFoldChangeThreshold)
def FoldChangeThreshold():
    '''Assess number of intervals above different fold change thresholds'''
    pass

@follows( replicatedIntervals, loadReplicatedIntervals )
def findReplicatedIntervals():
    '''Intersect replicates to identify replicated intervals '''
    pass

@follows( pairwiseIntervals, loadPairwiseIntervals, 
          uniqueIntervals, loadUniqueIntervals, 
          sharedIntervals, loadSharedIntervals,
          uniqueReplicatedIntervals, loadUniqueReplicatedIntervals, 
          sharedReplicatedIntervals, loadSharedReplicatedIntervals )
def comparePeaks():
    '''Compare intervals across tracks'''
    pass

@follows( getPeakShape )
def comparePeakShape():
    '''Compare interval shape within tracks'''
    pass

@follows( getCGIIntervals, loadCGIIntervals,
          getNonCGIIntervals, loadNonCGIIntervals,
          getPredictedCGIIntervals, loadPredictedCGIIntervals,
          getCGIOverlapCount, loadCGIOverlapCount,
          getChipseqOverlap, loadChipseqIntervals,
          getCapseqOverlap, loadCapseqIntervals,
          getChromatinMarkOverlap, loadChromatinMarkIntervals,
          getExternalBedStats, loadExternalBedStats,
          getCGIEnsemblTranscriptTSSOverlap, loadCGIEnsemblTranscriptTSSOverlap,
          getCGIEnsemblGeneTSSOverlap, loadCGIEnsemblGeneTSSOverlap,
          buildGATWorkspace, runExternalDatasetGAT, loadExternalDatasetGAT )
def compareExternal():
    '''Compare intervals external bed files'''
    pass

@follows( annotateEnsemblTranscriptOverlap, loadEnsemblTranscriptOverlap,
          annotateEnsemblGeneOverlap, loadEnsemblGeneOverlap,
          annotateEnsemblTranscriptTSS, loadEnsemblTranscriptTSS, 
          annotateEnsemblGeneTSS, loadEnsemblGeneTSS, 
          annotateEnsemblTranscriptTTS, loadEnsemblTranscriptTTS,
          annotateEnsemblGeneTTS, loadEnsemblGeneTTS,
          annotateCGIEnsemblTranscriptOverlap, loadCGIEnsemblTranscriptOverlap,
          annotateCGIEnsemblGeneOverlap, loadCGIEnsemblGeneOverlap,
          getCapseqEnsemblGeneTSSOverlap, loadCapseqEnsemblGeneTSSOverlap,
#          runEnsemblGenomicFeaturesGAT, loadEnsemblGenomicFeaturesGAT,
          annotateRepeats, loadRepeats )
def annotateEnsemblGenomicLocation():
    '''Annotate interval genomic location using data from Ensembl'''
    pass

@follows( getIntergenicCapseqBed,
          getEnsemblNoncodingTSSDistance,
          loadEnsemblNoncodingTSSDistance )
def noncoding():
    '''annotate distance of intergeneic CAPseq intervals from non-coding transcripts TSSs'''
    pass

@follows( getReplicatedEnsemblTranscriptTSSProfile, 
          getReplicatedEnsemblTranscriptTSSProfileCapseq,
          getReplicatedEnsemblTranscriptTSSProfileNoCapseq,
          getReplicatedEnsemblGeneTSSProfile, 
          getReplicatedEnsemblGeneTSSProfileCapseq,
          getReplicatedEnsemblGeneTSSProfileNoCapseq )
def annotateEnsemblTSSProfile():
    ''' Annotate TSS profile for all Ensembl transcripts'''
    pass
    
@follows( exportCapseqEnsemblTSSBed, exportCapseqNoEnsemblTSSBed )
def exportCapseqTSSIntervals():
    ''' Export CAPseq intervals that overlap with Ensembl TSS as bed file'''
    pass
    
@follows( liverTestesVenn, loadLiverTestesVenn,
          liverTestesTSSVenn, loadLiverTestesTSSVenn,
          liverTestesIntergenicVenn, loadLiverTestesIntergenicVenn )
def liverVsTestes():
    ''' Compare liver vs testes CAPseq intervals to make venn diagram'''
    pass
    
#@follows( annotateRnaseqGenomicFeatureOverlap, loadRnaseqGenomicFeatureOverlap,
#          annotateRnaseqTranscriptTSS, loadRnaseqTranscriptTSS, 
#          annotateRnaseqGeneTSS, loadRnaseqGeneTSS, 
#          annotateRnaseqTranscriptTTS, loadRnaseqTranscriptTTS,
#          annotateRnaseqGeneTTS, loadRnaseqGeneTTS,
#          annotateCGIRnaseqGenomicFeatureOverlap, loadCGIRnaseqGenomicFeatureOverlap,
#          annotateCGIRnaseqGenomicFeatures, loadCGIRnaseqGenomicFeatures,
#          runRnaseqGenomicFeaturesGAT, loadRnaseqGenomicFeaturesGAT )
#def annotateRnaseqGenomicLocation():
#    '''Annotate interval genomic location using data from RNAseq or CAGE datasets'''
#    pass

@follows( annotateCapseqComposition, loadCapseqComposition,
          annotateControlComposition, loadControlComposition,
          annotateFlankingCompositionLeft, loadFlankingCompositionLeft,
          annotateFlankingCompositionRight, loadFlankingCompositionRight,
          annotateCGIComposition, loadCGIComposition )
def annotateCapseqNucleotideComposition():
    '''Annotate interval nucleotide composition'''
    pass

@follows( annotateEnsemblTranscriptTSSComposition, loadEnsemblTranscriptTSSComposition,
          annotateEnsemblGeneTSSComposition, loadEnsemblGeneTSSComposition,
          annotateEnsemblGeneTSSIntervalComposition, loadEnsemblGeneTSSIntervalComposition )
def annotateEnsemblTSSNucleotideComposition():
    '''Annotate Ensembl TSS nucleotide composition'''
    pass

#@follows( annotateRnaseqTranscriptTSSComposition, loadRnaseqTranscriptTSSComposition,
#          annotateRnaseqGeneTSSComposition, loadRnaseqGeneTSSComposition,
#          annotateRnaseqGeneIntervalTSSComposition, loadRnaseqGeneIntervalTSSComposition )
#def annotateRnaseqTSSNucleotideComposition():
#    '''Annotate Rnaseq TSS nucleotide composition'''
#    pass

@follows( exportBigwig, exportMacsModel )
def export():
    '''export files'''
    pass

@follows( exportCapseqTSSOverlapGeneList, exportCapseqTSS1kbGeneList )
def exportGeneLists():
    '''export files'''
    pass

@follows( exportCapseqGCProfiles, exportCGIGCProfiles,
          exportEnsemblTranscriptTSSGCProfiles, exportEnsemblGeneTSSGCProfiles,
          exportCapseqCpGObsExp, exportCGICpGObsExp,
          exportEnsemblTranscriptTSSCpGObsExp, exportEnsemblGeneTSSCpGObsExp,
          exportCapseqCpGDensity, exportCGICpGDensity,
          exportEnsemblTranscriptTSSCpGDensity, exportEnsemblGeneTSSCpGDensity)
def exportCpGData():
    '''export files'''
    pass

@follows( mapReads,
          NormaliseBAMFiles,
          buildIntervalsMacs,
          buildIntervalsMacsNoControl,
          mergePeaks,
          background,
          FoldChangeThreshold,
          findReplicatedIntervals,
          comparePeaks,
          comparePeakShape,
          compareExternal,
          annotateEnsemblGenomicLocation,
          noncoding,
          annotateEnsemblTSSProfile,
          annotateCapseqNucleotideComposition,
          annotateEnsemblTSSNucleotideComposition,
          exportGeneLists )
def full():
    '''run the full pipeline.'''
    pass

@follows( mapReads,
          NormaliseBAMFiles,
          buildIntervalsMacs,
          buildIntervalsMacsNoControl,
          mergePeaks )
def mapAndCallPeaks():
    '''Map reads to genome and call peaks with MACS'''
    pass

@follows( background,
          FoldChangeThreshold,
          findReplicatedIntervals,
          comparePeaks,
          comparePeakShape,
          compareExternal )
def compareIntervals():
    '''Compare MACS intervals from different tracks with each other and with external bed files'''
    pass

@follows( annotateEnsemblGenomicLocation,
          noncoding,
          annotateEnsemblTSSProfile,
          annotateCapseqNucleotideComposition,
          annotateEnsemblTSSNucleotideComposition, )
def annotateEnsembl():
    '''Annotate CAPseq intervals using the Ensembl gene set'''
    pass

#@follows( annotateRnaseqGenomicLocation,
#          annotateCapseqNucleotideComposition,
#          annotateRnaseqTSSNucleotideComposition, )
#def annotateRnaseq():
#    '''Annotate CAPseq intervals using an RNAseq/CAGE based gene set'''
#    pass


if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
    

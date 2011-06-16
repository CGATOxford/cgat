################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_cpg.py 2900 2011-05-24 14:38:00Z david $
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
===================
CpG-Island pipeline
===================

:Author: David Sims 
:Release: $Id: pipeline_cpg.py 2900 2011-05-24 14:38:00Z david $
:Date: |today|
:Tags: Python

The CpG Island pipeline imports reads from one or more CpG island pulldown experiments and
performs the following tasks:

   * align reads to the genome
   * call peaks 
   * annotate intervals with respect to a reference gene set

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline_cpg.ini` file. The pipeline looks for a configuration file in several places:

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
|bowtie_             |>=0.12.7           |read mapping                                    |
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
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections
import Experiment as E
import logging as L
from ruffus import *
import csv
import sqlite3
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import IOTools
import MAST, GTF, GFF, Bed
import cStringIO
import pysam
import numpy
import gzip
import Masker
import fileinput
import gff2annotator

import pipeline_chipseq_intervals as PIntervals
import PipelineGeneset as PGeneset
import PipelineTracks
import PipelineMapping

USECLUSTER = True

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
import Pipeline as P
P.getParameters(  ["%s.ini" % __file__[:-len(".py")],  "../pipeline.ini", "pipeline.ini" ] )
PARAMS = P.PARAMS
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],"pipeline_annotations.py" )

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

def getControl( track ):
    '''return appropriate control for a track'''
    n = track.clone()
    n.condition = PARAMS["tracks_control"]
    return n

###################################################################
###################################################################
###################################################################
# define aggregates
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

# print "EXP=", EXPERIMENTS
# print "COND=", CONDITIONS
# print "TISSUES=", TISSUES
# print "TOSUBTRACT=", TOSUBTRACT
# print "MASTER=", TRACKS_MASTER
# print "CORRELATION=", TRACKS_CORRELATION

###################################################################
###################################################################
###################################################################
## MAP READS
@transform( ("*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra",
             "*.csfasta.gz" ),
            regex( r"(\S+).(export.txt.gz|fastq.1.gz|fastq.gz|sra|csfasta.gz)"), 
            r"\1/bam/\1.bam" )
def buildBAM( infile, outfile ):
    '''map reads with bowtie'''
    to_cluster = True
    track = P.snip( os.path.basename(outfile), ".bam" )
    try: os.mkdir( track )
    except OSError: pass
    try: os.mkdir( '''%(track)s/bam''' % locals() )
    except OSError: pass
    job_options= "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
    m = PipelineMapping.Bowtie()
    reffile = PARAMS["samtools_genome"]
    statement = m.build( (infile,), outfile ) 
    P.run()

#########################################################################
@transform( buildBAM,
            regex( r"(\S+)/bam/(\S+).bam"),
            r"\1/bam/\2.dedup.bam")
def dedup(infiles, outfile):
        '''Remove duplicate alignments from BAM files.'''
        to_cluster = USECLUSTER
        track = P.snip( outfile, ".bam" )
        dedup_method = PARAMS["dedup_method"]
        if dedup_method == 'samtools':
            statement = '''samtools rmdup %(infiles)s %(outfile)s; ''' % locals()    
        elif dedup_method == 'picard':
            statement = '''MarkDuplicates INPUT=%(infiles)s  ASSUME_SORTED=true OUTPUT=%(outfile)s METRICS_FILE=%(track)s.dupstats VALIDATION_STRINGENCY=SILENT; ''' % locals()
        statement += '''samtools index %(outfile)s; ''' % locals()
        #print statement
        P.run()

#########################################################################
@merge( dedup, "picard_duplicate_stats.load" )
def loadPicardDuplicateStats( infiles, outfile ):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''

    tablename = P.toTable( outfile )

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

    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
               '''
    P.run()

#########################################################################
@transform( dedup, 
            regex( r"(\S+)/bam/(\S+).bam"),
            r"\1/bam/\2.alignstats" )
def buildPicardAlignStats( infile, outfile ):
    '''Gather BAM file alignment statistics using Picard '''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".bam" )
    statement = '''CollectAlignmentSummaryMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(samtools_genome)s ASSUME_SORTED=true OUTPUT=%(outfile)s VALIDATION_STRINGENCY=SILENT ''' % locals()
    P.run()

############################################################
@merge( buildPicardAlignStats, "picard_align_stats.load" )
def loadPicardAlignStats( infiles, outfile ):
    '''Merge Picard alignment stats into single table and load into SQLite.'''

    tablename = P.toTable( outfile )

    outf = P.getTempFile()

    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".dedup.alignstats" )
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

    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
               '''
    P.run()

    os.unlink( tmpfilename )

#########################################################################
@transform( dedup, 
            regex(r"(\S+)/bam/(\S+).bam"),
            r"\1/bam/\2.readstats" )
def buildBAMStats( infile, outfile ):
    '''Count number of reads mapped, duplicates, etc. '''
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
    header = ",".join( [P.snip( os.path.basename(x), ".dedup.readstats") for x in infiles] )
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
                      --table=%(tname)s 
                      --allow-empty
                >> %(outfile)s """
        P.run()


############################################################
############################################################
############################################################
## BUILD INTERVALS USING CONTROL SAMPLE
@follows( dedup )
@files( [ (("%s/bam/%s.dedup.bam" % (x, x.asFile()), "%s/bam/%s.dedup.bam" % (getControl(x), getControl(x).asFile())), 
           "%s/macs/%s.macs" % (x, x.asFile()) ) for x in TRACKS ] )
def runMACS( infiles, outfile ):
    '''Run MACS for peak detection.'''
    infile, controlfile = infiles

    to_cluster = True

    track = P.snip( os.path.basename(infile), ".dedup.bam" )
    try: os.mkdir( track )
    except OSError: pass
    try: os.mkdir( '''%(track)s/macs''' % locals() )
    except OSError: pass

    statement = '''macs14 -t %(infile)s 
                          -c %(controlfile)s
                          --name=%(outfile)s
                          --format=BAM
                          --diag
                          %(macs_options)s 
                   >& %(outfile)s''' 
    P.run() 
    
############################################################
@transform( runMACS,
            regex(r"(\S+)/macs/(\S+).macs"),
            inputs( (r"\1/macs/\2.macs", r"\1/bam/\2.bam")), 
            r"\1/macs/\2.macs.load" )
def loadMACS( infiles, outfile ):
    infile, bamfile = infiles
    PIntervals.loadMACS( infile, outfile, bamfile )
    
############################################################
@merge( runMACS, "macs.summary" )
def summarizeMACS( infiles, outfile ):
    '''run MACS for peak detection.'''
    PIntervals.summarizeMACS( infiles, outfile )

############################################################
@transform( summarizeMACS, suffix(".summary"), "_summary.load" )
def loadMACSSummary( infile, outfile ):
    '''load macs summary.'''
    PIntervals.loadMACSSummary( infile, outfile )

############################################################
@transform( loadMACS, regex(r"(\S+)/macs/(\S+).macs.load"), r"\1/macs/\2.bed" )
def exportIntervalsAsBed( infile, outfile ):
    PIntervals.exportMacsAsBed( infile, outfile )

############################################################
############################################################
############################################################
## BUILD INTERVALS WITHOUT CONTROL SAMPLE
@follows( dedup )
@files( [ ("%s/bam/%s.dedup.bam" % (x, x.asFile()), 
           "%s/macs/%s.solo.macs" % (x, x.asFile()) ) for x in TRACKS ] )
def runMACSsolo( infile, outfile ):
    '''Run MACS for peak detection.'''

    to_cluster = True

    track = P.snip( os.path.basename(infile), ".dedup.bam" )
    try: os.mkdir( track )
    except OSError: pass
    try: os.mkdir( '''%(track)s/macs''' % locals() )
    except OSError: pass

    statement = '''macs14 -t %(infile)s 
                          --name=%(outfile)s
                          --format=BAM
                          --diag
                          %(macs_options)s 
                   >& %(outfile)s''' 
    P.run() 
    
############################################################
@transform( runMACSsolo,
            regex(r"(\S+)/macs/(\S+).solo.macs"),
            inputs( (r"\1/macs/\2.solo.macs", r"\1/bam/\2.bam")), 
            r"\1/macs/\2.solo.macs.load" )
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
    PIntervals.loadMACSSummary( infile, outfile )

############################################################
@transform( loadMACSsolo, regex(r"(\S+)/macs/(\S+).macs.load"), r"\1/macs/\2.bed" )
def exportIntervalsAsBedsolo( infile, outfile ):
    PIntervals.exportMacsAsBed( infile, outfile )

############################################################
############################################################
############################################################
## COMBINE INTERVALS
@follows( exportIntervalsAsBed )
@files( [( [ "%s.bed" % y.asFile() for y in EXPERIMENTS[x]], 
           "%s.bed" % x.asFile()) 
         for x in EXPERIMENTS ] )
def combineExperiment( infiles, outfile ):
    '''combine replicates between experiments.
    The replicates are combined using intersection.'''
    PIntervals.intersectBedFiles( infiles, outfile )

############################################################
@follows( exportIntervalsAsBed )    
@files( [( [ "%s.bed" % y.asFile() for y in CONDITIONS[x]], 
           "%s.bed" % x.asFile()) 
         for x in CONDITIONS ] )
def combineCondition( infiles, outfile ):
    '''combine conditions between cell lines. 
    The conditions are merged via intersection.'''
    PIntervals.intersectBedFiles( infiles, outfile )

############################################################
@follows( exportIntervalsAsBed )    
@files( [( [ "%s.bed" % y.asFile() for y in TISSUES[x]], 
           "%s.bed" % x.asFile()) 
         for x in TISSUES ] )
def combineTissue( infiles, outfile ):
    '''combine conditions between cell lines. 
    The conditions are merged via intersection. '''
    PIntervals.intersectBedFiles( infiles, outfile )

############################################################
@transform( (combineExperiment,
             combineCondition),
            suffix(".bed"),
            "_bed.load" )
def loadCombinedIntervals( infile, outfile ):
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

    tmpfile = tempfile.NamedTemporaryFile(delete=False)

    headers = ("AvgVal","DisttoStart","GeneList","Length","PeakCenter","PeakVal","Position","interval_id","nCpGs","nGenes","nPeaks","nProbes","nPromoters", "contig","start","end" )

    tmpfile.write( "\t".join(headers) + "\n" )

    avgval,contig,disttostart,end,genelist,length,peakcenter,peakval,position,start,interval_id,ncpgs,ngenes,npeaks,nprobes,npromoters = \
        0,"",0,0,"",0,0,0,0,0,0,0,0,0,0,0,

    samfiles, offsets = [], []

    track = Sample( filename = P.snip( infile, ".bed") )

    # get replicates / aggregated tracks associated with track
    # remove subtraction as not relevant for tag counting
    unsubtracted_track = getUnsubtracted ( track )
    replicates = PipelineTracks.getSamplesInTrack( unsubtracted_track, TRACKS )

    assert len(replicates) > 0
    
    # setup files
    for t in replicates:
        fn = "%s.norm.bam" % t.asFile()
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( pysam.Samfile( fn,  "rb" ) )
        fn = "%s.macs" % t.asFile()
        if os.path.exists( fn ):
            offsets.append( PIntervals.getPeakShift( fn ) )

    mlength = int(PARAMS["calling_merge_min_interval_length"])

    c = E.Counter()
    # count tags
    for line in open(infile, "r"):
        c.input += 1
        contig, start, end, interval_id = line[:-1].split()[:4]

        start, end = int(start), int(end)

        # remove very short intervals
        if end-start < mlength: 
            c.skipped_length += 1
            continue

        npeaks, peakcenter, length, avgval, peakval, nprobes = \
            PIntervals.countPeaks( contig, start, end, samfiles, offsets )

        # nreads can be 0 if the intervals overlap only slightly
        # and due to the binning, no reads are actually in the overlap region.
        # However, most of these intervals should be small and have already be deleted via 
        # the merge_min_interval_length cutoff.
        # do not output intervals without reads.
        if nprobes == 0:
            c.skipped_reads += 1
            
        c.output += 1
        tmpfile.write( "\t".join( map( str, (avgval,disttostart,genelist,length,
                                             peakcenter,peakval,position,interval_id,
                                             ncpgs,ngenes,npeaks,nprobes,npromoters, 
                                             contig,start,end) )) + "\n" )
 
    tmpfile.close()

    tmpfilename = tmpfile.name
    tablename = "%s_intervals" % track.asTable()
    
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --index=interval_id 
              --table=%(tablename)s
    < %(tmpfilename)s 
    > %(outfile)s
    '''

    P.run()
    os.unlink( tmpfile.name )

    L.info( "%s\n" % str(c) )

############################################################
@merge( exportIntervalsAsBed, "merged.bed" )
def buildMergedIntervals( infiles, outfile ):
    '''combine all experiments.
    The replicates are combined using a merge.'''
    PIntervals.mergeBedFiles( infiles, outfile )

############################################################
@merge( (exportIntervalsAsBed, 
         combineExperiment, 
         combineCondition,),
        "intervals.overlap" )
def buildOverlap( infiles, outfile ):
    '''compute overlap between intervals.'''

    if os.path.exists(outfile): 
        # note: update does not work due to quoting
        os.rename( outfile, outfile + ".orig" )
        options = "--update=%s.orig" % outfile
    else:
        options = ""

    infiles = " ".join( infiles )

    # note: need to quote track names
    statement = '''
        python %(scriptsdir)s/diff_bed.py %(options)s %(infiles)s 
        | awk -v OFS="\\t" '!/^#/ { gsub( /-/,"_", $1); gsub(/-/,"_",$2); } {print}'
        > %(outfile)s
        '''
    P.run()

############################################################
@transform( buildOverlap, suffix(".overlap"), "_overlap.load" )
def loadOverlap( infile, outfile ):
    '''load overlap results. '''

    tablename = "overlap"

    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --index=set1 
                       --index=set2 
                       --table=%(tablename)s 
                   < %(infile)s > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## PEAK CORRELATION
@follows( exportIntervalsAsBed )
@files_re( ["%s.bed" % x.asFile() for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "peakval.correlation" )
def makePeakvalCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.  '''
    PIntervals.makeIntervalCorrelation( infiles, outfile, "peakval", "merged.bed" )

@follows( exportIntervalsAsBed )
@files_re( ["%s.bed" % x.asFile() for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "avgval.correlation" )
def makeAvgvalCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field Avgval. '''
    PIntervals.makeIntervalCorrelation( infiles, outfile, "avgval", "merged.bed" )

@follows( exportIntervalsAsBed )
@files_re( ["%s.bed" % x.asFile() for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "length.correlation" )
def makeLengthCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field length.  '''
    PIntervals.makeIntervalCorrelation( infiles, outfile, "length", "merged.bed" )

@transform( ( makePeakvalCorrelation, makeAvgvalCorrelation, makeLengthCorrelation ),
            suffix(".correlation"),
            "_correlation.load")
def loadCorrelation( infile, outfile ):
    '''load correlation data.'''
    P.load( infile, outfile, "--index=id --map=default:float --map=id:int" )


############################################################
############################################################
############################################################
## ANNOTATE INTERVALS
@transform( (exportIntervalsAsBed,exportIntervalsAsBedsolo), suffix(".bed"), ".annotations" )
def annotateIntervals( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set.
    '''
    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_annotation_gff"] )

    statement = """
    cat < %(infile)s 
    | python %(scriptsdir)s/bed2gff.py --as-gtf 
    | python %(scriptsdir)s/gtf2table.py 
		--counter=position 
		--counter=classifier-chipseq 
		--section=exons 
		--counter=length 
		--log=%(outfile)s.log 
		--filename-gff=%(annotation_file)s 
		--genome-file=%(genome_dir)s/%(genome)s
    > %(outfile)s"""
    
    P.run()

############################################################
@transform( annotateIntervals, suffix(".annotations"), "_annotations.load" )
def loadAnnotations( infile, outfile ):
    '''load interval annotations: genome architecture '''
    P.load( infile, outfile, "--index=gene_id" )

############################################################
@transform( (exportIntervalsAsBed,exportIntervalsAsBedsolo), suffix(".bed"), ".tss" )
def annotateTSS( infile, outfile ):
    '''compute distance to TSS'''

    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_tss_bed"] )

    statement = """cat < %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
	                 | python %(scriptsdir)s/gtf2table.py 
		                   --counter=distance-tss 
		                   --log=%(outfile)s.log 
                       --filename-gff=%(annotation_file)s 
                       --filename-format="bed" 
                       --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s"""
    P.run()

############################################################
@transform( annotateTSS, suffix( ".tss"), "_tss.load" )
def loadTSS( infile, outfile ):
    '''load interval annotations: distance to transcription start sites '''
    P.load( infile, outfile, "--index=gene_id --index=closest_id --index=id5 --index=id3" )

############################################################
@transform( (exportIntervalsAsBed,exportIntervalsAsBedsolo), suffix(".bed"), ".repeats" )
def annotateRepeats( infile, outfile ):
    '''count the overlap between intervals and repeats.'''

    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],PARAMS_ANNOTATIONS["interface_repeats_gff"] )

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
@transform( annotateRepeats, suffix(".repeats"), "_repeats.load" )
def loadRepeats( infile, outfile ):
    '''load interval annotations: repeats'''
    P.load( infile, outfile, "--index=gene_id --allow-empty" )

############################################################
@transform( (exportIntervalsAsBed,exportIntervalsAsBedsolo), suffix(".bed"), ".composition" )
def annotateComposition( infile, outfile ):
    '''Establish the nucleotide composition of intervals'''

    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],PARAMS_ANNOTATIONS["interface_repeats_gff"] )

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
############################################################
############################################################
## EXPORT
@files_re( (buildBAM),
           "(.*).bam",
           "%s/\1.bigwig" % PARAMS["exportdir"])
def exportBigwig( infile, outfile ):
    '''convert BAM to bigwig file.'''

    # no bedToBigBed on the 32 bit cluster
    to_cluster = True

    statement = '''python %(scriptsdir)s/bam2wiggle.py \
                --genome-file=%(genome_dir)s/%(genome)s \
                --output-format=bigwig \
                --output-filename=%(outfile)s \
                %(infile)s \
                > %(outfile)s.log'''
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

@files_re(exportBigwig,combine("(.*).bigwig"), "bigwig.view")
def viewBigwig( infiles, outfile ):

    outs = open( outfile, "w" )
    outs.write( "# paste the following into the UCSC browser: \n" )

    try:
        os.makedirs( PARAMS["ucsc_dir"] )
    except OSError:
        pass
    
    for src in infiles:
        dest = os.path.join( PARAMS["ucsc_dir"], src ) 
        if not os.path.exists( dest ) or \
                os.path.getmtime(src) > os.path.getmtime(dest):
            shutil.copyfile( src, dest )
        track = src[:-len(".bigwig")]
        url = PARAMS["ucsc_url"] % src 
        outs.write( '''track type=bigWig name="%(track)s" description="%(track)s" bigDataUrl=%(url)s\n''' \
            % locals() )
    outs.close()

############################################################
@follows ( mkdir( 'export' ))
@merge( "run*.bed", ("export/intervals_%s.bed" % PARAMS["version"], "intervals.view") )
def viewIntervals( infiles, outfiles ):

    outfile_bed, outfile_code = outfiles
    outs = open( outfile_bed, "w" )
    version = PARAMS["version"]
    for infile in infiles:

        track = infile[:-len(".bed")]
        
        outs.write( '''track name="interval_%(track)s_%(version)s" description="Intervals in %(track)s - version %(version)s" visibility=2\n''' % locals() )

        with open(infile,"r") as f:
            for bed in Bed.iterator( f ):
                # MACS intervals might be less than 0
                if bed.start <= 0: 
                    bed.start = 0
                outs.write(str(bed) + "\n")

    outs.close()
    basename, filename = os.path.split( outfile_bed )

    dest = os.path.join( PARAMS["ucsc_dir"], filename )
    try:
        os.makedirs( PARAMS["ucsc_dir"] )
    except OSError:
        pass
    
    shutil.copyfile( outfile_bed, dest )

    filename = re.sub( "^.*/ucsc_tracks/", "", dest )

    outs = open( outfile_code, "w" )
    outs.write( "#paste the following into the UCSC browser:\n" )
    outs.write( "http://wwwfgu.anat.ox.ac.uk/~andreas/ucsc_tracks/%(filename)s\n" % locals())
    outs.close()

############################################################
############################################################
############################################################
## Pipeline organisation
@follows( buildBAM, dedup, loadPicardDuplicateStats, buildPicardAlignStats, loadPicardAlignStats, 
          buildBAMStats, loadBAMStats)
def mapReads():
    '''Align reads to target genome.'''
    pass

@follows( runMACS, loadMACS, summarizeMACS, loadMACSSummary, exportIntervalsAsBed )
def buildIntervals():
    pass

@follows( runMACSsolo, loadMACSsolo, summarizeMACSsolo, loadMACSsoloSummary, exportIntervalsAsBedsolo )
def buildIntervalsNoControl():
    pass

@follows( loadCombinedIntervals, 
          buildMergedIntervals, )
def combineIntervals():
    pass

@follows( loadCorrelation, 
          loadOverlap,)
def correlation():
    '''run the correlation targets.'''
    pass

@follows( annotateIntervals, loadAnnotations, 
          annotateTSS, loadTSS, 
          annotateRepeats, loadRepeats )
def annotation():
    '''run the annotation targets.'''
    pass

@follows( exportBigwig, viewBigwig, viewIntervals )
def export():
    '''export files'''
    pass

@follows( mapReads,
          buildIntervals, 
          annotation )
def full():
    '''run the full pipeline.'''
    pass

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

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )


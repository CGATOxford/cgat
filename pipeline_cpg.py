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
## Merge nearby intervals
@transform( exportIntervalsAsBed, regex(r"(\S+)/macs/(\S+).bed"), r"\1/macs/\2.merged.bed")
def mergeIntervals(infile, outfile):
    '''Merge intervals less than n bases apart in each dataset'''

    d = PARAMS["intervals_merge_dist"]
    statement = '''mergeBed -d %(d)s -i %(infile)s > %(outfile)s;'''
    P.run()

############################################################
@transform( mergeIntervals, regex(r"(\S+)/macs/(\S+).merged.bed"), r"\1/macs/\2.merged.load")
def loadMergedIntervals(infile, outfile):
    '''Load merged intervals into database '''
    track = P.snip( os.path.basename( infile ), ".merged.bed" )
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_merged_intervals
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## Compare intervals
@transform( exportIntervalsAsBed, regex(r"(\S+)/macs/(\S+).bed"), r"\1/macs/\2.unique.bed")
def uniqueIntervals(infile, outfile):
    '''identify unique intervals for each dataset'''

    tmpfile = P.getTempFilename(".")
    statement = '''mergeBed -i %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip( os.path.basename( infile ), ".bed")
    for track in TRACKS:
       if (str(track) <> in_track):
           statement = '''mergeBed -i %(track)s/macs/%(track)s.bed | intersectBed -a %(outfile)s -b stdin -v > %(tmpfile)s; mv %(tmpfile)s %(outfile)s ''' 
           P.run()
    os.unlink( tmpfile )

############################################################
@transform( uniqueIntervals, regex(r"(\S+)/macs/(\S+).unique.bed"), r"\1/macs/\2.unique.load")
def loadUniqueIntervals(infile, outfile):
    '''Load unique intervals into database '''
    track = P.snip( os.path.basename( infile ), ".unique.bed" )
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_unique_intervals
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
@transform( uniqueIntervals, regex(r"(\S+)/macs/(\S+).unique.bed"), r"\1/macs/\2.unique.coverage")
def analyseUniqueIntervals(infile, outfile):
    '''Analyse coverage of unique intervals in other datasets'''
    track = P.snip( os.path.basename( infile ), ".unique.bed" )
    header = "contig,start,stop"

    tmpfile = P.getTempFilename(".")
    statement = '''cat %(infile)s | sort -k1 -k2 > %(outfile)s;'''
    P.run()
    in_track = P.snip( os.path.basename( infile ), ".bed")
    for track in TRACKS:
       if (str(track) != in_track):
           header = header + "," + str(track)
           statement = '''coverageBed -abam %(track)s/bam/%(track)s.dedup.bam -b %(infile)s | sort -k1 -k2 | cut -f 4 | paste %(outfile)s - > %(tmpfile)s; mv %(tmpfile)s %(outfile)s'''
           P.run()

    # Load into database
    statement = '''cat %(outfile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_unique_coverage
                      --header=%(header)s
                   > %(outfile)s.load '''
    P.run()
    os.unlink( tmpfile )

############################################################
@transform( exportIntervalsAsBed, regex(r"(\S+)/macs/(\S+).bed"), r"\1/macs/\2.shared.bed")
def sharedIntervals(infile, outfile):
    '''identify shared intervals between datasets'''

    tmpfile = P.getTempFilename(".")
    statement = '''mergeBed -i %(infile)s > %(outfile)s;'''
    P.run()
    in_track = P.snip( os.path.basename( infile ), ".bed")
    for track in TRACKS:
       if (str(track) != in_track):
           statement = '''mergeBed -i %(track)s/macs/%(track)s.bed | intersectBed -a %(outfile)s -b stdin -u > %(tmpfile)s; mv %(tmpfile)s %(outfile)s; ''' 
           P.run()
    os.unlink( tmpfile )

############################################################
@transform( sharedIntervals, regex(r"(\S+)/macs/(\S+).shared.bed"), r"\1/macs/\2.shared.load")
def loadSharedIntervals(infile, outfile):
    '''Load shared intervals into database '''
    track = P.snip( os.path.basename( infile ), ".shared.bed" )
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_shared_intervals
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## Assess background (non-peak) binding
@follows( exportIntervalsAsBed )
@files( [ (("%s/bam/%s.dedup.bam" % (x, x.asFile()), "%s/macs/%s.bed" % (x, x.asFile())), "%s/macs/%s.bg" % (x, x.asFile()) ) for x in TRACKS ] )
def getBackground(infiles, outfile):
    '''Count the number of reads in the bamfile used for MACS that do not overlap an interval'''
    bam, bed = infiles
    statement = '''intersectBed -abam %(bam)s -b %(bed)s -v -bed | wc -l > %(outfile)s; '''
    statement += '''intersectBed -abam %(bam)s -b %(bed)s -u -bed | wc -l >> %(outfile)s;'''
    statement += '''sed -i '{N;s/\n/\t/g}' %(outfile)s; '''
    P.run()

############################################################
@transform( getBackground, regex(r"(\S+)/macs/(\S+).bg"), r"\1/macs/\2.bg.load")
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
def addBackground( infile, outfile):
    '''Add random reads to Chipseq input file to increase background '''

    # Identify sample to add reads to (lowest background)

    # Calculate number of reads to add (increase to highest background level)

    # Generate random reads and add to bam
    statement = '''java -jar -Xmx2048m SimSeq.jar 
                       -1 50 -2 50 \
                       --error  /ifs/apps/bio/simseq-72ce499/examples/hiseq_mito_default_bwa_mapping_mq10_1.txt \
                       --error2 /ifs/apps/bio/simseq-72ce499/examples/hiseq_mito_default_bwa_mapping_mq10_2.txt \
                       --insert_size 250 \
                       --insert_stdev 30 \
                       --read_number 10000 \
                       --read_prefix simseq_ \
                       --reference %(samtools_genome)s \
                       --duplicate_probability 0.0 \
                       --out simseq.sam; '''
    statement += '''samtools view -bS -t %(samtools_genome)s.fai -o simseq.bam simseq.sam; 
                    samtools sort out.bam out.srt`
                    samtools index out.sorted.bam`'''
    P.run()


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

    outtemp1 = P.getTempFile()
    tmpfilename1 = outtemp1.name
    outtemp2 = P.getTempFile()
    tmpfilename2 = outtemp2.name

    track= P.snip( os.path.basename(infile), ".bed")

    statement = """cat %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-na
                      	 --log=%(outfile)s
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(tmpfilename1)s; """
    statement += """cat %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(tmpfilename2)s; """
    statement += """python ~/src/combine_tables.py %(tmpfilename1)s %(tmpfilename2)s 
                    | python ~/src/csv2db.py 
                         --table=%(track)s_composition 
                         --index=gene_id
                   >> %(outfile)s"""
    P.run()

############################################################
@transform( (exportIntervalsAsBed), suffix(".bed"), ".control.composition" )
def controlComposition( infile, outfile ):
    '''Establish the nucleotide composition of control intervals'''

    to_cluster = True

    outtemp1 = P.getTempFile()
    tmpfilename1 = outtemp1.name
    outtemp2 = P.getTempFile()
    tmpfilename2 = outtemp2.name

    track= P.snip( os.path.basename(infile), ".bed")
    dirname= os.path.dirname(infile)

    statement = """slopBed -i %(infile)s -g %(samtools_genome)s.fai -l 10000 -r -10000 > %(dirname)s/%(track)s.control.bed; """

    statement += """cat %(dirname)s/%(track)s.control.bed
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-na
                      	 --log=%(outfile)s
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(tmpfilename1)s; """
    statement += """cat %(dirname)s/%(track)s.control.bed
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(tmpfilename2)s; """
    statement += """python ~/src/combine_tables.py %(tmpfilename1)s %(tmpfilename2)s 
                    | python ~/src/csv2db.py 
                         --table=%(track)s_control_composition 
                         --index=gene_id
                   >> %(outfile)s"""
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
## Compare bed files using GAT
@merge( "gat/*.bed", "gat/gat.tsv" )
def runGAT(infiles, outfile):
    '''Run genome association tester on bed files '''
    segfiles = " ".join( [ "--segment-file=%s" % x for x in infiles ] )
    annofiles = " ".join( [ "--annotation-file=%s" % x for x in infiles ] )
    statement = """gatrun.py %(segfiles)s %(annofiles)s --workspace=gat/%(genome)s.bed.gz --num-samples=1000 > %(outfile)s"""
    P.run()

############################################################
@transform( runGAT, regex(r"gat/*.tsv"), "*.load" )
def loadGAT(infile, outfile):
    '''Load genome association tester results into database '''
    statement = """cat %(infile)s | grep -v "^#" | python %(scriptsdir)s/csv2db.py 
                         --table=gat_results 
                    > %(outfile)s"""
    P.run()

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

@follows( annotateIntervals, loadAnnotations, 
          annotateTSS, loadTSS, 
          annotateRepeats, loadRepeats,
          annotateComposition, controlComposition )
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


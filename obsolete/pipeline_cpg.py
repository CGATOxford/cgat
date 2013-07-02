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
|bowtie              |>=0.12.7           |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|MACS                |14                 |peak finding                                    |
+--------------------+-------------------+------------------------------------------------+
|Picard              |>=1.4              |Dupliate removal, mapping stats                 |
+--------------------+-------------------+------------------------------------------------+
|BEDTools            |                   |interval comparison                             |
+--------------------+-------------------+------------------------------------------------+
|Zinba               |                   |peak finding                                    |
+--------------------+-------------------+------------------------------------------------+
|SICER               |                   |peak finding                                    |
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
import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections
import gzip
import sqlite3
import pysam
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.MAST as MAST
import CGAT.GTF as GTF
import CGAT.GFF as GFF
import CGAT.Bed as Bed
import cStringIO
import numpy
import CGAT.Masker as Masker
import fileinput
import gff2annotator

import CGAT.Experiment as E
import logging as L
from ruffus import *
import PipelineChipseq as PipelineChipseq
import CGATPipelines.PipelineTracks as PipelineTracks
import CGATPipelines.PipelineMapping as PipelineMapping

from rpy2.robjects import r as R
import rpy2.robjects as ro

USECLUSTER = True

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
import CGAT.Pipeline as P
P.getParameters(  ["%s/pipeline.ini" % os.path.splitext(__file__)[0],  "../pipeline.ini", "pipeline.ini" ] )
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

print "Expts=", EXPERIMENTS, "\n"

###################################################################
###################################################################
###################################################################
## Section 1: MAP READS
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
            r"\1/bam/\2.alignstats" )
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

    tablename = P.toTable( outfile )
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

    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
               '''
    P.run()
    os.unlink( tmpfilename )

#########################################################################
@transform( buildBAM, 
            regex( r"(\S+)/bam/(\S+).bam"),
            r"\1/bam/\2.gcstats" )
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

    tablename = P.toTable( outfile )
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

    statement = '''cat %(tmpfilename)s
                   | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                   > %(outfile)s '''
    P.run()
    os.unlink( tmpfilename )

#########################################################################
@transform( buildBAM, 
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
            statement = '''MarkDuplicates INPUT=%(infiles)s  ASSUME_SORTED=true OUTPUT=%(outfile)s 
                           METRICS_FILE=%(track)s.dupstats REMOVE_DUPLICATES=true 
                           VALIDATION_STRINGENCY=SILENT; ''' % locals()
        statement += '''samtools index %(outfile)s; ''' % locals()
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
                > %(outfile)s '''
    P.run()

############################################################
@transform( dedup,
            regex( r"(\S+)/bam/(\S+).dedup.bam"),
            r"\1/bam/\2.dedup.mapped.bam")
def removeUnmapped(infile, outfile):
    ''' Remove unmapped reads from BAM file'''
    to_cluster = True
    statement = '''samtools view -F 4 -bh %(infile)s > %(outfile)s; 
                   samtools index %(outfile)s; ''' % locals()
    P.run()

############################################################
@follows( removeUnmapped )
@files( [ (("%s/bam/%s.dedup.mapped.bam" % (x, x.asFile()), "%s/bam/%s.dedup.mapped.bam" % (getControl(x), getControl(x).asFile())), 
           "%s/bam/%s.norm.bam" % (x, x.asFile()) ) for x in TRACKS ] )
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
############################################################
############################################################
## BUILD INTERVALS USING CONTROL SAMPLE
@follows( normaliseBAMs )
@files( [ (("%s/bam/%s.norm.bam" % (x, x.asFile()), "%s/bam/%s.norm.bam" % (getControl(x), getControl(x).asFile())), 
           "%s/macs/%s.macs" % (x, x.asFile()) ) for x in TRACKS ] )
def runMACS( infiles, outfile ):
    '''Run MACS for peak detection using control sample.'''
    infile, controlfile = infiles
    to_cluster = True
    track = P.snip( os.path.basename(infile), ".norm.bam" )
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
            inputs( (r"\1/macs/\2.macs", r"\1/bam/\2.norm.bam")), 
            r"\1/macs/\2.macs.load" )
def loadMACS( infiles, outfile ):
    '''Load MACS intervals into database and filter on fold change and p-value'''
    infile, bamfile = infiles
    PipelineChipseq.loadMACS( infile, outfile, bamfile )
    
############################################################
@merge( runMACS, "macs.summary" )
def summarizeMACS( infiles, outfile ):
<<<<<<< local
    '''Parse MACS summary statistics from log file'''
    PIntervals.summarizeMACS( infiles, outfile )
=======
    '''run MACS for peak detection.'''
    PipelineChipseq.summarizeMACS( infiles, outfile )
>>>>>>> other

############################################################
@transform( summarizeMACS, suffix(".summary"), "_summary.load" )
def loadMACSSummary( infile, outfile ):
<<<<<<< local
    '''load macs summary into database'''
    PIntervals.loadMACSSummary( infile, outfile )
=======
    '''load macs summary.'''
    PipelineChipseq.loadMACSSummary( infile, outfile )
>>>>>>> other

############################################################
@transform( loadMACS, regex(r"(\S+)/macs/(\S+).macs.load"), r"\1/macs/\2.bed" )
def exportIntervalsAsBed( infile, outfile ):
<<<<<<< local
    '''Export MACS intervals from database as BED file and filter on fold change'''
    fc = PARAMS["intervals_min_fc"]
    PIntervals.exportMacsIntervalsAsBed( infile, outfile, fc )
=======
    PipelineChipseq.exportMacsAsBed( infile, outfile )
>>>>>>> other

############################################################
############################################################
############################################################
## BUILD INTERVALS WITHOUT CONTROL SAMPLE
@follows( normaliseBAMs )
@files( [ ("%s/bam/%s.norm.bam" % (x, x.asFile()), 
           "%s/macs/%s.solo.macs" % (x, x.asFile()) ) for x in TRACKS ] )
def runMACSsolo( infile, outfile ):
    '''Run MACS for peak detection.'''

    to_cluster = True
    track = P.snip( os.path.basename(infile), ".norm.bam" )
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
    PipelineChipseq.loadMACS( infile, outfile, bamfile )
    
############################################################
@merge( runMACSsolo, "macs_solo.summary" )
def summarizeMACSsolo( infiles, outfile ):
    '''run MACS for peak detection.'''
    PipelineChipseq.summarizeMACSsolo( infiles, outfile )

############################################################
@transform( summarizeMACSsolo, suffix(".summary"), "_summary.load" )
def loadMACSsoloSummary( infile, outfile ):
    '''load macs summary.'''
    PipelineChipseq.loadMACSSummary( infile, outfile )

############################################################
@transform( loadMACSsolo, regex(r"(\S+)/macs/(\S+).macs.load"), r"\1/macs/\2.bed" )
def exportIntervalsAsBedsolo( infile, outfile ):
<<<<<<< local
    fc = PARAMS["intervals_min_fc"]
    PIntervals.exportMacsIntervalsAsBed( infile, outfile, fc )
=======
    PipelineChipseq.exportMacsAsBed( infile, outfile )
>>>>>>> other

############################################################
############################################################
############################################################
## Find intervals using SICER with input
@follows( normaliseBAMs )
@files( [ (("%s/bam/%s.norm.bam" % (x, x.asFile()), "%s/bam/%s.norm.bam" % (getControl(x), getControl(x).asFile())), 
           "%s/sicer/%s.sicer" % (x, x.asFile()) ) for x in TRACKS ] )
def runSICER( infiles, outfile ):
    '''Run SICER for peak detection.'''
    infile, controlfile = infiles
    to_cluster = False

    track = P.snip( os.path.basename(infile), ".norm.bam" )
    control = P.snip( os.path.basename(controlfile), ".norm.bam" )
    inputdir = os.path.dirname(outfile)
    try: os.mkdir( track )
    except OSError: pass
    try: os.mkdir( '''%(track)s/sicer''' % locals() )
    except OSError: pass

    # convert bam to bed
    statement = '''bamToBed -i %(infile)s > %(track)s/sicer/%(track)s.bed; 
                   bamToBed -i %(controlfile)s > %(track)s/sicer/%(control)s.bed; '''

    # Run SICER
    statement += '''cd %(inputdir)s; SICER.sh . %(track)s.bed %(control)s.bed . %(genome)s %(sicer_params)s >& %(track)s.sicer''' 
    P.run() 
    
############################################################
@transform(runSICER, regex(r"(\S+)/sicer/(\S+).sicer"), r"\1/sicer/\2.sicer.load" )
def loadSICER( infile, outfile ):
    track = P.snip( os.path.basename(infile), ".sicer" )
    sicerdir = os.path.dirname(infile)
    window = PARAMS["sicer_window"]
    gap = PARAMS["sicer_gap"]
    FDR = PARAMS["sicer_fdr"]
    bedfile = sicerdir + "/" + track + "-W" + str(window) + "-G" + str(gap) + "-islands-summary-FDR" + str(FDR)
    tablename = "%s_sicer_intervals" % track
    headers="contig,start,stop,chip_reads,control_reads,pvalue,fold,fdr"

    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --header=%(headers)s
                       --index=contig,start
                       --table=%(tablename)s
                       --allow-empty 
                   < %(bedfile)s > %(outfile)s'''
    P.run()
    
    
############################################################
@merge( runSICER, "sicer.summary" )
def summarizeSICER( infiles, outfile ):
    '''run SICER for peak detection.'''
    def __get( line, stmt ):
        x = line.search(stmt )
        if x: return x.groups() 

    map_targets = [
        ("Window average: (\d+)", "window_mean",()),
        ("Minimum num of tags in a qualified window:  (\d+)", "window_min",()),
        ("The score threshold is:  (\d+)", "score_threshold",()),
        ("Total number of islands:  (\d+)","total_islands", ()),
        ("chip library size   (\d+)", "chip_library_size", ()),
        ("control library size   (\d+)", "control_library_size", ()),
        ("Total number of chip reads on islands is:  (\d+)", "chip_island_reads", ()),
        ("Total number of control reads on islands is:  (\d+)", "control_island_reads", ()),
        ("Given significance 0.01 ,  there are (\d+) significant islands",  "significant_islands", ()) ]

    mapper, mapper_header = {}, {}
    for x,y,z in map_targets: 
        mapper[y] = re.compile( x )
        mapper_header[y] = z

    keys = [ x[1] for x in map_targets ]

    outs = open(outfile,"w")

    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend( ["%s_%s" % (k,x) for x in mapper_header[k] ])
        else:
            headers.append( k )
    outs.write("track\t%s" % "\t".join(headers) + "\n" )

    for infile in infiles:
        results = collections.defaultdict(list)
        with open( infile ) as f:
            for line in f:
                if "diag:" in line: break
                for x,y in mapper.items():
                    s = y.search( line )
                    if s: 
                        results[x].append( s.groups()[0] )
                        break
                
        row = [ P.snip( os.path.basename(infile), ".sicer" ) ]
        for key in keys:
            val = results[key]
            if len(val) == 0: v = "na"
            else: 
                c = len(mapper_header[key])
                if c >= 1: assert len(val) == c, "key=%s, expected=%i, got=%i, val=%s, c=%s" %\
                   (key,
                    len(val),
                    c,
                    str(val), mapper_header[key])
                v = "\t".join( val )
            row.append(v)
        outs.write("\t".join(row) + "\n" )

    outs.close()


############################################################
@transform( summarizeSICER, suffix(".summary"), "_summary.load" )
def loadSICERSummary( infile, outfile ):
    '''load sicer summary.'''
    
    table = P.snip( os.path.basename(outfile), ".load" )
    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                      --index=track 
                      --table=%(table)s
                   < %(infile)s > %(outfile)s'''
    P.run()

############################################################
@transform( loadSICER, regex(r"(\S+)/sicer/(\S+).sicer.load"), r"\1/sicer/\2.sicer.bed" )
def exportSicerAsBed( infile, outfile ):
    '''export locations for all intervals.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename(infile), ".sicer.load" ).replace("-","_")

    cc = dbhandle.cursor()
    statement = "SELECT contig, start, stop FROM %s_sicer_intervals ORDER by contig, start" % track
    cc.execute( statement )

    outs = open( outfile, "w")
    for result in cc:
        contig, start, stop = result
        outs.write( "%s\t%i\t%i\n" % (contig, start, stop) )
    cc.close()
    outs.close()

############################################################
############################################################
############################################################
## Run Zinba to call peaks
@follows( dedup )
@files( [ (("%s/bam/%s.norm.bam" % (x, x.asFile()), "%s/bam/%s.norm.bam" % (getControl(x), getControl(x).asFile())), 
           "%s/zinba/%s.peaks" % (x, x.asFile()) ) for x in TRACKS ] )
def runZinba( infiles, outfile ):
    '''Run Zinba for peak detection.'''
    infile, controlfile = infiles
    to_cluster = False

    track = P.snip( os.path.basename(infile), ".norm.bam" )
    control = P.snip( os.path.basename(controlfile), ".norm.bam" )
    inputdir = os.path.dirname(outfile)
    frag_len = PARAMS['zinba_fragment_size']
    mappability = PARAMS['zinba_mappability']
    genome = PARAMS['zinba_genome']

    try: os.mkdir( track )
    except OSError: pass
    try: os.mkdir( '''%(track)s/zinba''' % locals() )
    except OSError: pass
    try: os.mkdir( '''%(track)s/zinba/map_ext%(frag_len)s''' % locals() )
    except OSError: pass

    # convert bam to bed
    statement = '''bamToBed -i %(infile)s > %(track)s/zinba/%(track)s.bed; 
                   bamToBed -i %(controlfile)s > %(track)s/zinba/%(control)s.bed; '''
    P.run()

    # Run Zinba
    R.library( 'zinba' )
    R( '''generateAlignability( mapdir='%(mappability)s', outdir='%(track)s/zinba/map_ext%(frag_len)s', athresh=1, extension=%(frag_len)s, twoBitFile='%(genome)s' )''' % locals() )
    R( '''basealigncount( inputfile='%(track)s/zinba/%(track)s.bed', outputfile='%(track)s/zinba/%(track)s.basecount', extension=%(frag_len)s, filetype='bed', twoBitFile='%(genome)s' )'''  % locals() )
    R( '''zinba( refinepeaks=1, seq='%(track)s/zinba/%(track)s.bed', input='%(track)s/zinba/%(control)s.bed', filetype='bed',  align='%(track)s/zinba/map_ext%(frag_len)s', twoBit='%(genome)s', outfile='%(track)s/zinba/%(track)s', extension=%(frag_len)s, basecountfile='%(track)s/zinba/%(track)s.basecount', numProc=4, threshold=0.01, broad=FALSE, printFullOut=1, interaction=FALSE, mode='peaks', FDR=TRUE) '''  % locals() )

############################################################
@transform(runZinba, regex(r"(\S+)/zinba/(\S+).peaks"), r"\1/zinba/\2.zinba.load" )
def loadZinba( infile, outfile ):
    track = P.snip( os.path.basename(infile), ".peaks" )
    #zinbadir = os.path.dirname(infile)

    tablename = "%s_zinba_intervals" % track
    headers="peakid,contig,start,stop,strand,sig,maxloc,maxval,pstart,pstop,median,qvalue"

    statement = '''cat | sed -v 'PEAKID'
                   | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --header=%(headers)s
                       --index=contig,start
                       --table=%(tablename)s
                       --allow-empty 
                   > %(outfile)s'''
    P.run()

############################################################
@transform( loadZinba, regex(r"(\S+)/zinba/(\S+).zinba.load"), r"\1/zinba/\2.zinba.bed" )
def exportZinbaAsBed( infile, outfile ):
    '''export locations for all intervals.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    track = P.snip( os.path.basename(infile), ".zinba.load" ).replace("-","_")

    cc = dbhandle.cursor()
    statement = "SELECT contig, pstart, pstop FROM %s_zinba_intervals ORDER by contig, start" % track
    cc.execute( statement )

    outs = open( outfile, "w")

    for result in cc:
        contig, start, stop = result
        outs.write( "%s\t%i\t%i\n" % (contig, start, stop) )
    cc.close()
    outs.close()

############################################################
############################################################
############################################################
## Merge nearby intervals
@transform( (exportIntervalsAsBed), regex(r"(\S+)/macs/(\S+).bed"), r"\1/macs/\2.merged.bed")
def mergeIntervals(infile, outfile):
    '''Merge intervals less than n bases apart in each dataset and update foldchange scores'''

    d = PARAMS["intervals_merge_dist"]
    method = PARAMS["intervals_merge_method"]
    PIntervals.mergeIntervalsWithScores( infile, outfile, d, method )

############################################################
@transform( mergeIntervals, regex(r"(\S+)/macs/(\S+).merged.bed"), r"\1/macs/\2.merged.cleaned.bed")
def sanitiseIntervals(infile, outfile):
    '''sanatise so that intervals do not exceed contig length'''
    statement = '''cat %(infile)s | python %(scriptsdir)s/bed2bed.py --method=filter-genome --genome-file=%(genome_dir)s/%(genome)s -L %(outfile)s.log > %(outfile)s;'''
    P.run()

############################################################
@transform( (sanitiseIntervals), regex(r"(\S+)/macs/(\S+).merged.cleaned.bed"), r"\1/macs/\2.merged.cleaned.bed.load" )
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

    fn = "%s/bam/%s.norm.bam" % (base_track,track)
    assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, track)
    samfiles.append( pysam.Samfile( fn,  "rb" ) )
    fn = "%s/macs/%s.macs" % (base_track,track)
    if os.path.exists( fn ):
        offsets.append( PIntervals.getPeakShift( fn ) )

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
@follows( sanitiseIntervals )
@files( [ (("%s/bam/%s.norm.bam" % (x, x.asFile()), "%s/macs/%s.merged.cleaned.bed" % (x, x.asFile())), "%s/macs/%s.bg" % (x, x.asFile()) ) for x in TRACKS ] )
def getBackground(infiles, outfile):
    '''Count the number of reads in the bamfile used for MACS that do not overlap an interval'''
    bam, bed = infiles
    to_cluster = True
    statement = '''intersectBed -abam %(bam)s -b %(bed)s -v -bed | wc -l > %(outfile)s; '''
    statement += '''intersectBed -abam %(bam)s -b %(bed)s -u -bed | wc -l >> %(outfile)s;'''
    statement += '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
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
############################################################
############################################################
## Assess effect of altering fold change threshold
@follows( exportIntervalsAsBed )
@files( [ ("%s/macs/%s.bed" % (x, x.asFile()), "%s/macs/%s.foldchange" % (x, x.asFile()) ) for x in TRACKS ] )
def thresholdFoldChange(infile, outfile):
    '''Assess interval overlap between conditions at different fold change thresholds '''

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )

    # Make bed files for different fold change thresholds
    track = P.snip( os.path.basename( infile ), ".bed" ).replace("-","_")
    macsdir = os.path.dirname(infile)
    try: os.mkdir( macsdir+"/foldchange" )
    except OSError: pass
    foldchange = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

    for fc in foldchange:
        cc = dbhandle.cursor()
        query = '''SELECT contig, start, end, interval_id, fold FROM %(track)s_macs_intervals WHERE fold > %(fc)i ORDER by contig, start;''' % locals()
        cc.execute( query )

        outbed = macsdir + "/foldchange/" + track + ".fc" + str(fc) + ".bed"
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
@transform( thresholdFoldChange, regex(r"(\S+)/macs/(\S+).foldchange"), r"\1/macs/\2.foldchange.load")
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
@transform( thresholdFoldChange, regex(r"(\S+)/macs/(\S+).foldchange"), r"\1/fc/\2.foldchange")
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

    # for each foldchange 
    for fc in foldchange:

        in_bed = in_dir + "/foldchange/" + in_track.replace("-","_") + ".fc" + str(fc) + ".bed"

        # For each track
        for track in TRACKS:
            if (str(track) != in_track):
                compare_bed = str(track) + "/macs/foldchange/" + str(track).replace("-","_") + ".fc" + str(fc) + ".bed"
                statement = '''echo %(track)s %(fc)s >> %(outfile)s; intersectBed -a %(in_bed)s -b %(compare_bed)s -u | wc -l >> %(outfile)s; ''' 
                P.run()

    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; sed -i '{s/ /\\t/g}' %(outfile)s; '''
    P.run()

############################################################
@transform( sharedIntervalsFoldChangeThreshold, regex(r"(\S+)/fc/(\S+).foldchange"), r"\1/fc/\2.foldchange.load")
def loadSharedIntervalsFoldChangeThreshold(infile, outfile):
    '''Load intervals overlapping other tracks into database '''
    track = P.snip( os.path.basename( infile ), ".foldchange" ).replace(".","_").replace("-","_")
    header = "track,threshold,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_foldchange_shared
                      --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## Compare intervals
@transform( sanitiseIntervals, regex(r"(\S+)/macs/(\S+).merged.cleaned.bed"), r"\1/macs/\2.overlap")
def pairwiseIntervals(infile, outfile):
    '''identify overlapping intervals for each pair of datasets'''

    to_cluster = True
    in_track = P.snip( os.path.basename( infile ), ".merged.cleaned.bed")
    statement = '''echo "track" > %(outfile)s; echo "overlap" >> %(outfile)s;'''
    P.run()

    for track in TRACKS:
       if (str(track) <> in_track):
           statement = '''echo %(track)s >> %(outfile)s; intersectBed -a %(infile)s -b %(track)s/macs/%(track)s.merged.cleaned.bed -u | wc -l >> %(outfile)s; '''
           P.run()
    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

############################################################
@transform( pairwiseIntervals, regex(r"(\S+)/macs/(\S+).overlap"), r"\1/macs/\2.overlap.load")
def loadPairwiseIntervals(infile, outfile):
    '''Load overlapping intervals into database '''
    track = P.snip( os.path.basename( infile ), ".overlap" ).replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_overlap
                   > %(outfile)s '''
    P.run()

############################################################
@transform( sanitiseIntervals, regex(r"(\S+)/macs/(\S+).merged.cleaned.bed"), r"\1/macs/\2.unique.bed")
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
           statement = '''intersectBed -a %(outfile)s -b %(track)s/macs/%(track)s.merged.cleaned.bed -v > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s ''' 
           P.run()

############################################################
@transform( uniqueIntervals, regex(r"(\S+)/macs/(\S+).unique.bed"), r"\1/macs/\2.unique.load")
def loadUniqueIntervals(infile, outfile):
    '''Load unique intervals into database '''
    track = P.snip( os.path.basename( infile ), ".unique.bed" )
    header = "contig,start,stop,interval_id,fold"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_unique_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()

############################################################
@transform( uniqueIntervals, regex(r"(\S+)/macs/(\S+).unique.bed"), r"\1/macs/\2.unique.coverage")
def analyseUniqueIntervals(infile, outfile):
    '''Analyse coverage of unique intervals in other datasets'''

    header = "contig,start,stop"
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name
    in_track = P.snip( os.path.basename( infile ), ".unique.bed")
    # Sort infile
    statement = '''cat %(infile)s | sort -k1 -k2 > %(outfile)s;'''
    P.run()

    for track in TRACKS:
       if (str(track) != in_track):
           header = header + "," + str(track)
           statement = '''coverageBed -abam %(track)s/bam/%(track)s.dedup.bam -b %(infile)s | sort -k1 -k2 | cut -f 4 | paste %(outfile)s - > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s'''
           P.run()

    # Load into database
    statement = '''cat %(outfile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(in_track)s_unique_coverage
                      --header=%(header)s
                   > %(outfile)s.load '''
    P.run()

############################################################
@transform( sanitiseIntervals, regex(r"(\S+)/macs/(\S+).merged.cleaned.bed"), r"\1/macs/\2.shared.bed")
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
           statement = '''intersectBed -a %(outfile)s -b %(track)s/macs/%(track)s.merged.cleaned.bed -u > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s; ''' 
           P.run()

############################################################
@transform( sharedIntervals, regex(r"(\S+)/macs/(\S+).shared.bed"), r"\1/macs/\2.shared.load")
def loadSharedIntervals(infile, outfile):
    '''Load shared intervals into database '''
    track = P.snip( os.path.basename( infile ), ".shared.bed" )
    header = "contig,start,stop,interval_id,fold"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_shared_intervals
                      --header=%(header)s
                      --index=contig,start
                      --index=interval_id
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## Calculate replicated intervals
@follows( sanitiseIntervals, mkdir("replicated_intervals") )
@files( [( [ "%s/macs/%s.merged.cleaned.bed" % (y,y.asFile()) for y in EXPERIMENTS[x]], 
           "replicated_intervals/%s.replicated.bed" % str(x).replace("-agg","")) for x in EXPERIMENTS ] )
def replicatedIntervals( infiles, outfile ):
    '''Combine replicates between experiments.
       First all intervals are merged across replicates 
       and then merged intervals that do not overlap an interval in all replicates are removed. '''

    expt = P.snip( os.path.basename( outfile ), ".replicated.bed").replace("-agg","")
    outdir = os.path.dirname( outfile )
    tmpfile = P.getTempFile()
    tmpfilename = tmpfile.name

    # merge all files to get total set of intervals
    in_list = " ".join(infiles)
    statement = '''cat %(in_list)s | mergeBed -i stdin > %(outdir)s/%(expt)s.merged.bed; 
                   cat %(outdir)s/%(expt)s.merged.bed > %(outfile)s;''' % locals()
    P.run()

    # Remove intervals that do not overlap intervals in all reps
    for i in infiles:
        statement = '''intersectBed -a %(outfile)s -b %(i)s -u > %(tmpfilename)s; mv %(tmpfilename)s %(outfile)s; ''' 
        P.run()

############################################################
@transform( replicatedIntervals, regex(r"replicated_intervals/(\S+).replicated.bed"), r"replicated_intervals/\1.replicated.bed.load")
def loadReplicatedIntervals(infile, outfile):
    '''Load replicated intervals into database '''
    expt = P.snip( os.path.basename( infile ), ".replicated.bed" ).replace("-","_").replace(".","_")
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(expt)s_replicated_intervals
                      --header=%(header)s
                      --index=contig,start
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## Compare replicated intervals with external bed files
@transform( replicatedIntervals, regex(r"replicated_intervals/(\S+).replicated.bed"), r"replicated_intervals/\1.cgi_cap.bed")
def getCGIIntervals(infile, outfile):
    '''identify intervals overlapping CGI for each datasets'''

    CGI = PARAMS["bed_ucsc_cgi"]
    dataset_name =  P.snip( os.path.basename( CGI ), ".bed")
    statement = '''intersectBed -a %(infile)s -b %(CGI)s -u > %(outfile)s; '''
    P.run()

############################################################
@transform( getCGIIntervals, regex(r"replicated_intervals/(\S+).cgi_cap.bed"), r"replicated_intervals/\1.cgi_cap.bed.load")
def loadCGIIntervals(infile, outfile):
    '''Load intervals overlapping CGI into database '''
    track = P.snip( os.path.basename( infile ), ".cgi_cap.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_cgi_cap_bed
                      --index=contig,start
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( replicatedIntervals, regex(r"replicated_intervals/(\S+).replicated.bed"), r"replicated_intervals/\1.cap_only.bed")
def getNonCGIIntervals(infile, outfile):
    '''identify intervals not overlapping CGI for each datasets'''

    CGI = PARAMS["bed_ucsc_cgi"]
    dataset_name =  P.snip( os.path.basename( CGI ), ".bed")
    statement = '''intersectBed -a %(infile)s -b %(CGI)s -v > %(outfile)s; '''
    P.run()

############################################################
@transform( getNonCGIIntervals, regex(r"replicated_intervals/(\S+).cap_only.bed"), r"replicated_intervals/\1.cap_only.bed.load")
def loadNonCGIIntervals(infile, outfile):
    '''Load intervals not overlapping CGI into database '''
    track = P.snip( os.path.basename( infile ), ".cap_only.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_cap_only_bed
                      --index=contig,start
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( replicatedIntervals, regex(r"replicated_intervals/(\S+).replicated.bed"), r"replicated_intervals/\1.cgi_only.bed")
def getPredictedCGIIntervals(infile, outfile):
    '''identify predicted CGI intervals not overlapping CAPseq intervals for each dataset'''

    CGI = PARAMS["bed_ucsc_cgi"]
    dataset_name =  P.snip( os.path.basename( CGI ), ".bed")
    statement = '''intersectBed -a %(CGI)s -b %(infile)s -v > %(outfile)s; '''
    P.run()

############################################################
@transform( getPredictedCGIIntervals, regex(r"replicated_intervals/(\S+).cgi_only.bed"), r"replicated_intervals/\1.cgi_only.bed.load")
def loadPredictedCGIIntervals(infile, outfile):
    '''Load predicted CGI intervals not overlapping CAP-seq intervals into database '''
    track = P.snip( os.path.basename( infile ), ".cgi_only.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop,source"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_cgi_only_bed
                      --index=contig,start
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( replicatedIntervals, regex(r"replicated_intervals/(\S+).replicated.bed"), r"replicated_intervals/\1.cgi_overlap")
def getCGIOverlapCount(infile, outfile):
    '''identify intervals overlapping CGI for each datasets'''

    CGI = P.asList(PARAMS["bed_cgi"])
    for island in CGI:
       dataset_name =  P.snip( os.path.basename( island ), ".bed")
       statement = '''echo %(dataset_name)s >> %(outfile)s; intersectBed -a %(infile)s -b %(island)s -u | wc -l >> %(outfile)s; '''
       P.run()
    statement = '''sed -i '{N;s/\\n/\\t/}' %(outfile)s; '''
    P.run()

############################################################
@transform( getCGIOverlapCount, regex(r"replicated_intervals/(\S+).cgi_overlap"), r"replicated_intervals/\1.cgi_overlap.load")
def loadCGIOverlapCount(infile, outfile):
    '''Load intervals overlapping CGI into database '''
    track = P.snip( os.path.basename( infile ), ".cgi_overlap" ).replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_cgi
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( replicatedIntervals, regex(r"replicated_intervals/(\S+).replicated.bed"), r"replicated_intervals/\1.chipseq")
def getChipseqOverlap(infile, outfile):
    '''identify intervals overlapping chipseq intervals for each datasets'''

    chipseq = P.asList(PARAMS["bed_chipseq"])
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
@transform( getChipseqOverlap, regex(r"replicated_intervals/(\S+).chipseq"), r"replicated_intervals/\1.chipseq.load")
def loadChipseqIntervals(infile, outfile):
    '''Load intervals overlapping chipseq into database '''
    track = P.snip( os.path.basename( infile ), ".chipseq" ).replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_chipseq
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( replicatedIntervals, regex(r"replicated_intervals/(\S+).replicated.bed"), r"replicated_intervals/\1.capseq")
def getCapseqOverlap(infile, outfile):
    '''identify intervals overlapping capseq intervals for each datasets'''

    capseq = P.asList(PARAMS["bed_capseq"])
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
@transform( getCapseqOverlap, regex(r"replicated_intervals/(\S+).capseq"), r"replicated_intervals/\1.capseq.load")
def loadCapseqIntervals(infile, outfile):
    '''Load intervals overlapping capseq into database '''
    track = P.snip( os.path.basename( infile ), ".capseq" ).replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_capseq
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( replicatedIntervals, regex(r"replicated_intervals/(\S+).replicated.bed"), r"replicated_intervals/\1.chromatin")
def getChromatinMarkOverlap(infile, outfile):
    '''identify intervals overlapping chromatin mark intervals for each datasets'''

    chromatin = P.asList(PARAMS["bed_chromatin"])
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
@transform( getChromatinMarkOverlap, regex(r"replicated_intervals/(\S+).chromatin"), r"replicated_intervals/\1.chromatin.load")
def loadChromatinMarkIntervals(infile, outfile):
    '''Load intervals overlapping chromatin marks into database '''
    track = P.snip( os.path.basename( infile ), ".chromatin" ).replace(".","_").replace("-","_")
    header = "track,overlap"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_chromatin
                      --header=%(header)s
                      --allow-empty
                   > %(outfile)s '''
    P.run()

############################################################
@merge( "gat/*.bed", "external_interval_sets.stats" )
def getExternalBedStats(infiles, outfile):
    '''Calculate statistics for external bed files '''
    chromatin = P.asList(PARAMS["bed_chromatin"])
    capseq = P.asList(PARAMS["bed_capseq"])
    chipseq = P.asList(PARAMS["bed_chipseq"])
    CGI = P.asList(PARAMS["bed_cgi"])
    extBed = chromatin + capseq + chipseq + CGI

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
## Compare bed files using GAT
@files(PARAMS["samtools_genome"]+".fai", "gat/"+PARAMS["genome"]+".bed.gz")
def buildGATWorkspace(infile, outfile ):
    '''Build genomic workspace file for GAT '''
    statement = '''cat %(infile)s | awk 'OFS="\\t" {print $1,0,$2,"workspace"}' | gzip > %(outfile)s '''
    P.run()

############################################################
@follows( buildGATWorkspace )
@follows( mkdir("gat") )
@merge( replicatedIntervals, "gat/external_dataset_gat.tsv" )
def runExternalDatasetGAT(infiles, outfile):
    '''Run genome association tester on bed files '''

    to_cluster = True
    segfiles = ""
    for x in infiles:
        track = P.snip(os.path.basename(x), ".replicated.bed")
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
@transform( runExternalDatasetGAT, regex(r"gat/(\S+).tsv"), r"\1.load" )
def loadExternalDatasetGAT(infile, outfile):
    '''Load genome association tester results into database '''
    statement = """cat %(infile)s | grep -v "^#" | python %(scriptsdir)s/csv2db.py 
                         --table=external_dataset_gat_results gat/external_dataset_gat.tsv
                    > %(outfile)s"""
    P.run()

############################################################
############################################################
############################################################
## Compare intervals from different peak callers
@transform( mergeIntervals, regex(r"(\S+)/macs/(\S+).merged.bed"), r"\1/macs/\2.macs.sicer.bed")
def getSicerOverlap(infile, outfile):
    '''identify intervals overlapping SICER intervals for each datasets'''
    track = P.snip( os.path.basename( infile ), ".merged.bed")
    macsdir = os.path.dirname(infile)
    sicer = track + ".sicer.bed"
    sicerdir = macsdir.replace("macs","sicer")
    statement = '''intersectBed -a %(infile)s -b %(sicerdir)s/%(sicer)s -u > %(outfile)s; '''
    P.run()

############################################################
@transform( getSicerOverlap, regex(r"(\S+)/macs/(\S+).macs.sicer.bed"), r"\1/macs/\2.macs.sicer.load")
def loadSicerIntervals(infile, outfile):
    '''Load intervals overlapping SICER intervals into database '''
    track = P.snip( os.path.basename( infile ), ".macs.sicer.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_macs_sicer_intervals_shared
                      --header=%(header)s
                      --ignore-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( mergeIntervals, regex(r"(\S+)/macs/(\S+).merged.bed"), r"\1/macs/\2.macs.zinba.bed")
def getZinbaOverlap(infile, outfile):
    '''identify intervals overlapping ZINBA intervals for each datasets'''
    track = P.snip( os.path.basename( infile ), ".merged.bed")
    macsdir = os.path.dirname(infile)
    zinba = track + ".zinba.bed"
    zinbadir = macsdir.replace("macs","zinba")
    statement = '''intersectBed -a %(infile)s -b %(zinbadir)s/%(zinba)s -u > %(outfile)s; '''
    P.run()

############################################################
@transform( getZinbaOverlap, regex(r"(\S+)/macs/(\S+).macs.zinba.bed"), r"\1/macs/\2.macs.zinba.load")
def loadZinbaIntervals(infile, outfile):
    '''Load intervals overlapping ZINBA intervals into database '''
    track = P.snip( os.path.basename( infile ), ".macs.zinba.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_macs_zinba_intervals_shared
                      --header=%(header)s
                      --ignore-empty
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## ANNOTATE INTERVALS
@transform( replicatedIntervals, suffix(".replicated.bed"), ".annotations" )
def annotateIntervals( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set. '''
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
@transform( annotateIntervals, suffix(".annotations"), ".annotations.load" )
def loadAnnotations( infile, outfile ):
    '''load interval annotations: genome architecture '''
    P.load( infile, outfile, "--index=gene_id" )

############################################################
@transform( replicatedIntervals, suffix(".replicated.bed"), ".tss" )
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
@transform( annotateTSS, suffix( ".tss"), ".tss.load" )
def loadTSS( infile, outfile ):
    '''load interval annotations: distance to transcription start sites '''
    P.load( infile, outfile, "--index=gene_id --index=closest_id --index=id5 --index=id3" )

############################################################
@transform( replicatedIntervals, suffix(".replicated.bed"), ".tts" )
def annotateTTS( infile, outfile ):
    '''compute distance to TTS'''

    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_tts_bed"] )

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
@transform( annotateTTS, suffix( ".tts"), ".tts.load" )
def loadTTS( infile, outfile ):
    '''load interval annotations: distance to transcription termination sites '''
    P.load( infile, outfile, "--index=gene_id --index=closest_id --index=id5 --index=id3" )

############################################################
@transform( replicatedIntervals, suffix(".replicated.bed"), ".repeats" )
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
@transform( annotateRepeats, suffix(".repeats"), ".repeats.load" )
def loadRepeats( infile, outfile ):
    '''load interval annotations: repeats'''
    P.load( infile, outfile, "--index=gene_id --allow-empty" )

############################################################
@transform( replicatedIntervals, suffix(".replicated.bed"), ".composition" )
def annotateComposition( infile, outfile ):
    '''Establish the nucleotide composition of intervals'''

    to_cluster = True

    statement = """cat %(infile)s 
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-na
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( replicatedIntervals, suffix(".replicated.bed"), ".control.composition" )
def annotateControlComposition( infile, outfile ):
    '''Establish the nucleotide composition of control intervals'''

    to_cluster = True
    track= P.snip( os.path.basename(infile), ".replicated.bed")
    dirname= os.path.dirname(infile)

    statement = """cat %(infile)s | python %(scriptsdir)s/bed2bed.py -m shift -g %(genome_dir)s/%(genome)s --offset=-10000 -S %(dirname)s/%(track)s.control.bed;
                   cat %(dirname)s/%(track)s.control.bed
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-na
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s.log
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( (annotateComposition, annotateControlComposition), suffix( ".composition"), ".composition.load" )
def loadComposition( infile, outfile ):
    '''Load the nucleotide composition of intervals'''

    track= P.snip( os.path.basename(infile), ".composition").replace(".","_").replace("-","_")

    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=%(track)s_composition 
                         --index=gene_id
                 > %(outfile)s; """
    P.run()

############################################################
@files( PARAMS["annotations_dir"] + "/tss.bed.gz", "tss.composition" )
def annotateTSSComposition( infile, outfile ):
    '''Establish the nucleotide composition of tss intervals'''

    to_cluster = True
    tss_extend = 500

    statement = """zcat %(infile)s 
                   | slopBed -i stdin -g %(samtools_genome)s.fai -b %(tss_extend)s
                   | python %(scriptsdir)s/bed2gff.py --as-gtf 
                   | python %(scriptsdir)s/gtf2table.py 
                      	 --counter=composition-na
                      	 --counter=composition-cpg
                      	 --log=%(outfile)s.log
                         --genome-file=%(genome_dir)s/%(genome)s
                   > %(outfile)s; """
    P.run()

############################################################
@transform( annotateTSSComposition, suffix( ".composition"), ".composition.load" )
def loadTSSComposition( infile, outfile ):
    '''Load the nucleotide composition of tss intervals'''

    statement = """cat %(infile)s | python ~/src/csv2db.py 
                         --table=tss_comp
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
                      	 --counter=composition-na
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
@files( PARAMS["bed_ucsc_cgi"], "cgi.annotations" )
def annotatePredictedCGIs( infile, outfile ):
    '''classify predicted CGI intervals according to their location with respect to the gene set. '''
    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_annotation_gff"] )

    # Give each row a unique identifier
    statement = """cat %(infile)s 
                   | awk '{print $1,$2,$3,$4-NR}'
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
@transform( annotatePredictedCGIs, suffix(".annotations"), ".annotations.load" )
def loadCGIAnnotations( infile, outfile ):
    '''load CGI annotations: genome architecture '''
    P.load( infile, outfile, "--index=gene_id" )

############################################################
@files( ( PARAMS["annotations_dir"] + "/tss.bed.gz", PARAMS["bed_ucsc_cgi"]), "overlap.tss.cgi" )
def getTSSCGIOverlap( infiles, outfile ):
    '''Establish overlap between cgi and tss intervals'''

    tss, cgi = infiles
    tss_extend = 500

    to_cluster = True

    outtemp1 = P.getTempFile()
    tmpfilename1 = outtemp1.name

    statement = """zcat %(tss)s | slopBed -i stdin -g %(samtools_genome)s.fai -b %(tss_extend)s > tss_extended.bed; """
    statement += """echo "shared" >> %(outfile)s; intersectBed -a tss_extended.bed -b %(cgi)s | wc -l >> %(outfile)s; """
    statement += """echo "tss" >> %(outfile)s; intersectBed -a tss_extended.bed -b %(cgi)s -v | wc -l >> %(outfile)s; """
    statement += """echo "cgi" >> %(outfile)s; intersectBed -a %(cgi)s -b tss_extended.bed -v | wc -l >> %(outfile)s; """
    statement += '''sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''
    P.run()

############################################################
@transform( getTSSCGIOverlap, regex(r"overlap.tss.cgi"), r"overlap.tss.cgi.load")
def loadTSSCGIOverlap(infile, outfile):
    '''load TSS CGI overlap into database'''

    header = "track,intervals"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                        --table=tss_cgi_venn
                        --header=%(header)s
                   > %(outfile)s '''
    P.run()

############################################################
############################################################
############################################################
## Compare intervals genomic features using GAT
@follows(buildGATWorkspace)
@merge( replicatedIntervals, "gat/genomic_features_gat.tsv" )
def runGenomicFeaturesGAT(infiles, outfile):
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
@transform( runGenomicFeaturesGAT, regex(r"gat/(\S+).tsv"), r"\1.load" )
def loadGenomicFeaturesGAT(infile, outfile):
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
                       --genome-file=%(genome_dir)s/%(genome)s 
                       --output-format=bigwig 
                       --output-filename=%(outfile)s 
                   %(infile)s > %(outfile)s.log '''
    P.run()

############################################################
@merge( exportBigwig, "bigwig.view")
def viewBigwig( infiles, outfile ):

    outs = open( outfile, "w" )
    outs.write( "# paste the following into the UCSC browser: \n" )

    try:
        os.makedirs( PARAMS["ucsc_dir"] )
    except OSError:
        pass
    
    for src in infiles:
        dest = os.path.join( PARAMS["ucsc_dir"], src ) 
        if not os.path.exists( dest ) or os.path.getmtime(src) > os.path.getmtime(dest):
            shutil.copyfile( src, dest )
        track = src[:-len(".bigwig")]
        url = PARAMS["ucsc_url"] % src 
        outs.write( '''track type=bigWig name="%(track)s" description="%(track)s" bigDataUrl=%(url)s\n''' % locals() )
    outs.close()

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
@transform( annotateTSS, suffix(".tss"), r".cgi_genes")
def exportCGIGeneLists( infile, outfile):
    '''Export list of genes associated with replicated CAPseq intervals '''

    max_gene_dist = PARAMS["intervals_max_gene_dist"]

    # Connect to DB
    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename( infile ), ".tss" ).replace("-","_")

    # Extract data from db
    cc = dbhandle.cursor()
    query = '''SELECT distinct id5, dist5 FROM %(track)s_tss where dist5 < %(max_gene_dist)s ORDER by id5;''' % locals()
    cc.execute( query )

    # Write to file
    outs = open( outfile, "w")
    for result in cc:
        gene_id, distance = result
        outs.write( "%s\t%i\n" % (gene_id, distance) )
    cc.close()
    outs.close()

############################################################
############################################################
############################################################
## Pipeline organisation
@follows( buildBAM, 
          dedup, loadPicardDuplicateStats, 
          buildPicardAlignStats, loadPicardAlignStats,
          buildBAMStats, loadBAMStats)
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

@follows( runSICER, loadSICER, 
          summarizeSICER, loadSICERSummary, 
          exportSicerAsBed )
def buildIntervalsSicer():
    '''Find peaks using SICER using control sample'''
    pass

@follows( runZinba, loadZinba, exportZinbaAsBed )
def buildIntervalsZinba():
    '''Find peaks using ZINBA using control sample'''
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
    '''Assess number of intervals above different fold chnage thresholds'''
    pass

@follows( pairwiseIntervals, loadPairwiseIntervals, 
          uniqueIntervals, loadUniqueIntervals, 
          sharedIntervals, loadSharedIntervals,
          replicatedIntervals, loadReplicatedIntervals)
def comparePeaks():
    '''Compare intervals across tracks'''
    pass

@follows( getCGIIntervals, loadCGIIntervals,
          getNonCGIIntervals, loadNonCGIIntervals,
          getPredictedCGIIntervals, loadPredictedCGIIntervals,
          getCGIOverlapCount, loadCGIOverlapCount,
          getChipseqOverlap, loadChipseqIntervals,
          getCapseqOverlap, loadCapseqIntervals,
          getChromatinMarkOverlap, loadChromatinMarkIntervals,
          getExternalBedStats, loadExternalBedStats,
          buildGATWorkspace, runExternalDatasetGAT, loadExternalDatasetGAT )
def compareExternal():
    '''Compare intervals external bed files'''
    pass

@follows( getSicerOverlap, loadSicerIntervals )
def compareCallers():
    '''Compare intervals from different peak callers'''
    pass

@follows( annotateIntervals, loadAnnotations, 
          annotateTSS, loadTSS, 
          annotateTTS, loadTTS,
          annotateRepeats, loadRepeats,
          annotateComposition, annotateControlComposition, loadComposition,
          annotateTSSComposition, loadTSSComposition,
          annotateCGIComposition, loadCGIComposition,
          annotatePredictedCGIs, loadCGIAnnotations,
          getTSSCGIOverlap, loadTSSCGIOverlap,
          runGenomicFeaturesGAT, loadGenomicFeaturesGAT )
def annotation():
    '''Annotate interval location and composition'''
    pass

@follows( exportBigwig, exportMacsModel )
def export():
    '''export files'''
    pass

@follows( mapReads,
          NormaliseBAMFiles,
          buildIntervalsMacs,
          buildIntervalsMacsNoControl,
          mergePeaks,
          background,
          FoldChangeThreshold,
          comparePeaks,
          compareExternal,
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



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

import Pipeline as P
P.getParameters(  ["%s.ini" % __file__[:-len(".py")],  "pipeline.ini" ] )
PARAMS = P.PARAMS

###################################################################
###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
# load all tracks - exclude input/control tracks
Sample = PipelineTracks.Sample3

#TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( [x for x in glob.glob( "*.fastq.gz" ) if PARAMS["tracks_control"] not in x], "(\S+).fastq.gz" )
TRACKS = PipelineTracks.Tracks( Sample ).loadFromDirectory( 
    [ x.replace("../","") for x in glob.glob( "*.export.txt.gz" ) if PARAMS["tracks_control"] not in x ],
      "(\S+).export.txt.gz" ) +\
      PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
          [ x.replace("../","") for x in glob.glob( "*.sra" ) if PARAMS["tracks_control"] not in x ], 
          "(\S+).sra" ) +\
          PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
              [x.replace("../","") for x in glob.glob( "*.fastq.gz" ) if PARAMS["tracks_control"] not in x], 
              "(\S+).fastq.gz" ) +\
              PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
                  [x.replace("../","") for x in glob.glob( "*.fastq.1.gz" ) if PARAMS["tracks_control"] not in x], 
                  "(\S+).fastq.1.gz" ) +\
                  PipelineTracks.Tracks( PipelineTracks.Sample3 ).loadFromDirectory( 
                      [ x.replace("../","") for x in glob.glob( "*.csfasta.gz" ) if PARAMS["track_control"] not in x], 
                        "(\S+).csfasta.gz" )
                                      
for X in TRACKS:
    print "TRACK=", X, "\n"

###################################################################
###################################################################
###################################################################
@transform( ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz" ),
            regex( r"(\S+).(export.txt.gz|fastq.1.gz|fastq.gz|sra|csfasta.gz)"), 
            r"\1.bam" )
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
@transform( buildBAM, suffix( ".bam"), ".dedup.bam")
def dedup(infiles, outfile):
        '''Remove duplicate alignments from BAM files.'''
        to_cluster = USECLUSTER
        track = P.snip( outfile, ".bam" )
        statement = '''MarkDuplicates INPUT=%(infiles)s  ASSUME_SORTED=true OUTPUT=%(outfile)s 
                       METRICS_FILE=%(track)s.dupstats REMOVE_DUPLICATES=true 
                       VALIDATION_STRINGENCY=SILENT; 
                       samtools index %(outfile)s; ''' % locals()
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
############################################################    
############################################################
@follows( dedup )
@files( [ ("%s.dedup.bam" % x, "%s.macs" % x ) for x in TRACKS ] )
def runMACSsolo( infile, outfile ):
    '''Run MACS for peak detection.'''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".dedup.bam" )
    statement = '''macs14 -t%(infile)s 
                          --name=%(track)s
                          --format=BAM
                          --wig -S
                          %(macs_options)s 
                   >& %(outfile)s;''' 
    P.run() 
    
############################################################
@transform( runMACSsolo, regex(r"(\S+).macs"),
            inputs( (r"\1.macs", r"\1.dedup.bam")), 
            r"\1.macs.load" )
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
@transform( loadMACSsolo, regex(r"(\S+).macs.load"), r"\1.macs.bed" )
def exportIntervalsAsBedsolo( infile, outfile ):
    '''Export list of intervals passing fold change threshold to file '''
    fc = PARAMS["intervals_min_fc"]
    PIntervals.exportMacsIntervalsAsBed( infile, outfile, fc )
    
############################################################
############################################################ 
############################################################ 
@transform( dedup, regex(r"(\S+).dedup.bam"), r"\1.bw" )
def bamToWig( infile, outfile ):
    ''' convert bam to bigwig '''
    statement = '''python %%(scriptsdir)s/bam2wiggle.py --output-format=bigwig %(infile)s %(outfile)s ''' % locals()
    P.run()

###################################################################
@transform( dedup, regex(r"(\S+).dedup.bam"), r"\1.bed" )
def bamToBed( infile, outfile ):
    ''' convert bam to bigwig '''
    statement = '''bamToBed -i %(infile)s > %(outfile)s; 
                   bgzip %(outfile)s;
                   tabix %(outfile)s.gz; ''' % locals()
    P.run()

###################################################################
@follows( buildBAM, buildPicardAlignStats,
          loadPicardAlignStats, dedup,
          loadPicardDuplicateStats )
def align():
    '''align fastq files to genome using Bowtie and convert BAM to bigwig'''
    pass

@follows( runMACSsolo, loadMACSsolo, summarizeMACSsolo, 
          loadMACSsoloSummary, exportIntervalsAsBedsolo )
def macs():
    '''align fastq files to genome using Bowtie and convert BAM to bigwig'''
    pass

@follows( bamToWig, bamToBed )
def convert():
    '''align fastq files to genome using Bowtie and convert BAM to bigwig'''
    pass
    
@follows( align, macs, convert )
def full():
    '''align fastq files to genome using Bowtie and convert BAM to bigwig'''
    pass



if __name__== "__main__":
    sys.exit( P.main(sys.argv) )
    

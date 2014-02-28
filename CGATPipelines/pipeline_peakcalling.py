"""=====================
Peak calling pipeline
=====================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The peak calling pipeline applies various peak calling algorithms
on mapped reads.

Overview
========

pipeline_peakcalling takes as input reads aligned to genomic sequence
as :term:`bam` formatted files and calls peaks. The pipeline
implements several peak callers:

macs_
   Model-based Analysis of ChIP-Seq (MACS), for identifying
   transcript factor binding sites. MACS captures the influence of
   genome complexity to evaluate the significance of enriched ChIP
   regions, and MACS improves the spatial resolution of binding sites
   through combining the information of both sequencing tag position
   and orientation. MACS can be easily used for ChIP-Seq data alone,
   or with control sample with the increase of specificity.

macs2_ 
   MACS 2 is the new release of the MACS peak caller. Among other
   improvements it adds support for handling paired end reads.

spp_ 
   SPP is a R package especially designed for the analysis of
   Chip-Seq data from Illummina platform. The package was developed by
   Peter Park's group from Harvard Medical School.

zinba_ 
   ZINBA (Zero Inflated Negative Binomial Algorithm) is a
   computational and statistical framework used to call regions of the
   genome enriched for sequencing reads originating from a diverse
   array of biological experiments. We collectively refer to the
   sequencing data derived from these experiments as DNA-seq,
   including FAIRE-seq, ChIP-seq, and DNAase-seq experiments

sicer_narrow 
    A clustering approach for identification of enriched
    domains from histone modification ChIP-Seq data.  The types of
    region called by the sicer alogrithm reflect the paramaters it is
    run with. In this pipeline sicer is run by default with two
    different parameter sets, to allow the simultaneous detection of
    narrower and broader regions of enrichmnet.

sicer_broad
    (See above)

peakranger_ranger 
    PeakRanger is a multi-purpose, ultrafast ChIP Seq
    peak caller. It is used in the modENCODE project and included in
    the iPlant pipeline system.  PeakRanger v1.02 was developed in
    Dr.Lincoln Stein's lab at OICR and is now in continual development
    at Dr.Helen Hobbs's lab of the McDermott Center of UT
    Southwestern.  PeakRanger can run two separate algorithms for peak
    detection, "ranger" for point-source binding event detection and
    "ccat" for the detection of broader regions of enrichment.  In
    this pipeline, both alorighms are presented as separate peak
    callers for convience.

peakranger_ccat 
    PeakRanger is here run using the CCAT alogorithm (See:
    Xu, H., L. Handoko, et al. (2010).A signal-noise model for
    significance analysis of ChIP-seq with negative
    control.Bioinformatics 26(9): 1199-1204)

scripture
    As of version 2, scripture has added a ChIP-Seq module. The current
    version (5/11/2013) is still marked beta. Scripture implements a
    continuous segmentation algorithm. It scans windows of a fixed size
    across the genome, identifies significant windows and merges/trims
    to obtain peaks. Significance is assessed using the scan statistic
    (Wallenstein and Neff, 1987).
    Scripture in the pipeline is used according to Garber et al. (2012)
    (PMID:22940246). Significant windows are filtered against input or
    other control data (Garber et al.: whole cell extract) by computing
    an enrichment scrore. Only significant windows above a certain 
    enrichment score are kept.

Peak callers have different strengths and weaknesses. Some might work
well on broad peaks such as some histone marks, others work better for
narrow, sharp peaks. Many callers these days attempt to call both
types of peaks.  This pipeline implements the following nomenclature
for peaks:

.. glossary::

   region
      A broad region defined by a peak caller. Regions are usually
      created in the first step of peak calling, before peak
      refinement or detection of subpeaks is performed.

   summit
      A narrow region defined by a caller. These are the output of any
      peak refinement or subpeak detection steps

   interval
      Generic term for an interval in which the read density is higher
      then expected. Both a :term:`region` or a :term:`summit` are an
      :term:`interval`.

   peak
      Within a :term:`region` or :term:`summit` the position with the
      highest base coverage.

The pipeline computes some basic measures to validate peak calling. In
order to fully annotate peaks, use :doc:`pipeline_intervals`.

.. note::

   The pipeline currently expects that mulit-mapping reads (reads
   mapping to multiple locations) have been removed.  

QC
---

The pipeline implements the following QC measures. See :pmid:`22955991`.

NSC
   normalized strand coefficient.

RSC
  relative strand correlacion.

The pipeline will also do and IDR analysis (see
https://sites.google.com/site/anshulkundaje/projects/idr) for spp
called peaks.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.

The sphinxreport report requires a :file:`conf.py` and
:file:`sphinxreport.ini` file (see :ref:`PipelineReporting`). To start
with, use the files supplied with the Example_ data.

Input
-----

Mapped reads
++++++++++++

The principal input of this pipeline is a collection of reads mapped to a reference genome.
Mapped reads are imported by placing files are linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.genome.bam

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. 

Please not the suffix ``genome.bam`` which is required to distinguish the input :term:`bam`
formatted files from those that are created in the pipeline after duplication removal (``.prep.bam``
and ``.call.bam``)

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
|spp_                |1.11               |R package                                       |
+--------------------+-------------------+------------------------------------------------+
|macs_               |>1.4               |                                                |
+--------------------+-------------------+------------------------------------------------+
|zinba_              |2.02.03            |R package                                       |
+--------------------+-------------------+------------------------------------------------+
|peakranger          |1.15               |                                                |
+--------------------+-------------------+------------------------------------------------+
|sicer               |1.1                |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`. For each peak caller there are tables
called:

<track>_<caller>_regions
<track>_<caller>_summits

Each of these tables contains the following columns:

+------------------+--------------------------------------------------------------+
|*Column*          |*Content*                                                     |
+------------------+--------------------------------------------------------------+
|avgval            |Average read depth in interval                                |
+------------------+--------------------------------------------------------------+
|contig            |Contig                                                        |
+------------------+--------------------------------------------------------------+
|control_avgval    |Average read depth in control within inter val                |
|                  |                                                              |
+------------------+--------------------------------------------------------------+
|control_length    |Interval length                                               |
+------------------+--------------------------------------------------------------+
|control_npeaks    |Number of peaks in control                                    |
+------------------+--------------------------------------------------------------+
|control_nreads    |Number of control reads in interval                           |
+------------------+--------------------------------------------------------------+
|control_peakcenter|Peak center of control                                        |
+------------------+--------------------------------------------------------------+
|control_peakval   |Number of reads at peak in control                            |
+------------------+--------------------------------------------------------------+
|end               |End coordinate of interval                                    |
+------------------+--------------------------------------------------------------+
|length            |Length of interval                                            |
+------------------+--------------------------------------------------------------+
|npeaks            |Number of peaks in interval                                   |
+------------------+--------------------------------------------------------------+
|nreads            |Number of reads in interval                                   |
+------------------+--------------------------------------------------------------+
|peakcenter        |Peak center in interval                                       |
+------------------+--------------------------------------------------------------+
|peakval           |Number of reads at peak                                       |
+------------------+--------------------------------------------------------------+
|start             |246251                                                        |
+------------------+--------------------------------------------------------------+

The unprocessed output files created by the peak callers are in individual subdirectories
for each caller (:file:`macs.dir`, :file:`zinba.dir`, etc.).

IDR analysis
------------

The output of the IDR analysis is in the :file:`idr.dir` directory.

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

.. _macs: http://liulab.dfci.harvard.edu/MACS/00README.html
.. _macs2: http://liulab.dfci.harvard.edu/MACS/00README.html
.. _spp: http://compbio.med.harvard.edu/Supplements/ChIP-seq/tutorial.html
.. _sicer: http://home.gwu.edu/~wpeng/Software.htm
.. _zinba: http://code.google.com/p/zinba/
.. _peakranger: http://ranger.sourceforge.net/
.. _scripture: http://www.broadinstitute.org/software/scripture/home

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
import CGAT.Bed as Bed
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta

import CGATPipelines.PipelinePeakcalling as PipelinePeakcalling
import CGATPipelines.PipelineMotifs as PipelineMotifs
import CGATPipelines.PipelineTracks as PipelineTracks
import CGATPipelines.PipelineMappingQC as PipelineMappingQC

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
Sample = PipelineTracks.AutoSample
Sample.attributes = ('experiment', 'condition', 'tissue')

# collect tracks, exclude any control tracks
TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    [x for x in glob.glob("*.genome.bam") if
     PARAMS["tracks_control"] not in x],
    "(\S+).genome.bam")

EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

##################################################################
def getControl( track, suffix='.genome.bam'):
    '''return appropriate control(s) for a track.
    '''
    prefix = PipelineTracks.FILE_SEPARATOR.join((track.experiment,
                                                 PARAMS["tracks_control"]))
    controlfiles = glob.glob(prefix + "*" + suffix)
    controls = [Sample(filename=x[:-len(suffix)]) for x in controlfiles]
    # if multiple, filter by number of parts
    if len(controls) > 1:
        controls = [x for x in controls if 
                    len(x.attributes) == len(track.attributes)]
    return controls

def getControlFile( controls, pattern ):
    if not controls:
        L.warn("controls for track '%s' not found " % (track))
        controlfile = None
    else:
        controlfile = pattern % controls[0].asFile()

    return controlfile

def getUnstimulated( track ):
    '''return unstimulated condition for a track
    '''
    n = track.clone()
    n.condition = PARAMS["tracks_unstimulated"]
    return n

def getSubtracted( track ):
    '''return "subtracted" condition for a track
    '''
    n = track.clone()
    n.condition += PARAMS["tracks_subtract"]
    return n

def getUnsubtracted( track ):
    '''return "unsubtracted" condition for a track
    '''
    n = track.clone()
    if PARAMS["tracks_subtract"]:
        if n.condition.endswith( PARAMS["tracks_subtract"] ):
            n.condition = n.condition[:-len(PARAMS["tracks_subtract"])]
    return n

def getBamFiles( infile, suffix ):
    '''return associated bamfiles with infiles.'''
    track = P.snip( os.path.basename(infile), suffix)
    bamfile = P.snip( os.path.basename(infile), suffix ) + ".call.bam"
    assert os.path.exists( bamfile ), "bamfile %s does not exist" % bamfile

    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )
    if not os.path.exists( controlfile ):
        L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        controlfile = None

    return bamfile, controlfile

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

############################################################
############################################################
############################################################
@merge(["%s.genome.bam" % x.asFile() for x in TRACKS],
       'check_input.dummy')
def checkInput(infiles, outfile):
    '''perform a check if every track has a control/input.'''

    c = E.Counter()
    for infile in infiles:
        c.input += 1
        track = P.snip(infile, ".genome.bam")
        control = getControl(Sample(track))
        if len(control) == 0:
            E.warn('no control file found for %s' % track)
            c.notfound += 1
        elif len(control) > 1:
            E.warn('multiple control files found for %s: %s' % \
                   (track, control))
            c.multiple += 1
        else:
            E.info('track->control: %s -> %s' % (track, control[0]))
            c.found += 1
    E.info(c)

############################################################
############################################################
############################################################
@files(input == None, "regions.mask")
def makeMask(infile,outfile):
    '''Make a mask for filtering reads if required.
    '''
    if PARAMS["calling_filter_exons"] or PARAMS["calling_filter_regions"]:
        regions_to_filter = []

        if PARAMS["calling_filter_exons"]:
            regions_to_filter += PipelinePeakcalling.getExonLocations(PARAMS["calling_filter_exons"]) 

        if PARAMS["calling_filter_regions"]:
            regions_to_filter += Bed.iterator( IOTools.openFile( PARAMS["calling_filter_regions"]) )
            
        fh = IOTools.openFile(outfile,"w")

        for bed in itertools.chain( regions_to_filter ):
            fh.write( "%s\n" % "\t".join(map(str, (bed.contig, bed.start, bed.end))) )

        fh.close()
    else:
        P.touch(outfile)

############################################################
############################################################
############################################################
@transform( "*.genome.bam", suffix(".genome.bam"), add_inputs( makeMask), ".prep.bam")
def prepareBAMForPeakCalling(infiles, outfile):
    '''Prepare BAM files for peak calling.

        - unmapped reads are removed.

        - if the option "calling_deduplicate" is Picard.MarkDuplicates is run 
            to remove duplicate reads

        - reads may be filtered by exon or location 

           - to remove reads by exon, the option "calling_filter_exons" should specify a file containing 
             a list of ensembl gene identifiers (one per line)
           - to remove reads by location, the option "calling_filter_regions" should specify a bed file''
        
        The resulting bam file has a .prep.bam extension. Merging infiles is currently untested and the 
        methods only consider single end reads.
    '''
    bam_file, mask_file = infiles

    if PARAMS["calling_filter_exons"] or PARAMS["calling_filter_regions"]:
        mask = mask_file
    else:
        mask = None

    PipelinePeakcalling.buildBAMforPeakCalling(bam_file,outfile,PARAMS["calling_deduplicate"], mask )

############################################################
############################################################
############################################################
@merge( prepareBAMForPeakCalling, "preparation_stats.load" )
def loadDuplicationStats( infiles, outfile ):
    '''load output from Picard Deduplication step.'''
    
    PipelineMappingQC.loadPicardMetrics( infiles, outfile, 
                                         pipeline_suffix = ".prep.bam",
                                         suffix = ".picard_metrics" )
    

############################################################
############################################################
############################################################
if PARAMS["calling_normalize"]==True:
    '''Normalise the number of reads in a set of prepared bam files.

    The minimum number of reads in a prepared bam file is calculated and this
    number of reads is used as a threshold to randomly sample from each bam file 
    in order to create a set of bam files with near identical numbers of reads.
  
    This may result in considerable data loss. 

    Per experimental contrast normalisation could be preferable.
   
    Potentially usefull if using a peak caller that does not correct for tag
    count between experimental and input samples.
    '''
    # First count the number of reads in each bam
    @transform(prepareBAMForPeakCalling,suffix("prep.bam"),"prep.count")
    def countReadsInBAM(infile,outfile):
        to_cluster = True
        statement= '''samtools idxstats %s | awk '{s+=$3} END {print s}' > %s ''' % ( infile,outfile )
        P.run()

    # Get the lowest number of reads
    @merge(countReadsInBAM,"minreads")
    def minReads(infiles,outfile):
        counts = []
        countfiles = glob.glob("*.prep.count")
        for cf in countfiles:
            fh = IOTools.openFile(cf,"r")
            counts.append(int(fh.read()))
            fh.close()
        fh = IOTools.openFile(outfile,"w")
        fh.write(str(min(counts)))
        fh.close

    # Make normalised files
    @follows(minReads)
    @transform(prepareBAMForPeakCalling,
               regex(r"(.*).prep.bam"),
               inputs( (r"\1.prep.bam",r"\1.prep.count") ), 
               r"\1.call.bam")
    def normalizeBAM( infiles, outfile ):
        '''build a normalized BAM file such that all
    files have approximately the same number of 
    reads.
    '''   
        fh = IOTools.openFile("minreads")
        minreads = int(fh.read())
        fh.close
        PipelinePeakcalling.buildSimpleNormalizedBAM( infiles, 
                                                      outfile,
                                                      minreads)
else:
    @transform(prepareBAMForPeakCalling,
               suffix(".prep.bam"),
               ".call.bam")
    def normalizeBAM( infile, outfile):
        P.clone( infile, outfile )
        P.clone( infile + ".bai", outfile + ".bai" )


######################################################################
######################################################################
##                                                                  ##
##                 Statistics and QC Functions                      ## 
##                                                                  ##
######################################################################
######################################################################
@follows( mkdir("readstats.dir") )
@transform( normalizeBAM,
            regex("(.*).bam"),
            r"readstats.dir/\1.readstats" )
def buildBAMStats( infile, outfile ):
    '''count number of reads mapped, duplicates, etc.
    '''

    to_cluster = True

    statement = '''python
    %(scriptsdir)s/bam2stats.py
         --force
         --output-filename-pattern=%(outfile)s.%%s
    < %(infile)s
    > %(outfile)s
    '''

    P.run()

####################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''import bam statisticis.'''
    
    PipelineMappingQC.loadBAMStats( infiles, outfile )


######################################################################
@follows( normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "%s.data_quality" % x.asFile() ) for x in TRACKS ] )
def checkDataQuality( infile, outfile ):
    '''uses peakranger to check data quality.

    nr is the signal/noise ratio.
    '''
    
    track = P.snip( infile, ".call.bam" )
    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )

    if not os.path.exists(controlfile):
        L.warn("controlfile '%s' for track '%s' not found " % (controlfile, track ))
        P.touch(outfile)
        return

    to_cluster = True
    statement = '''peakranger nr --format bam %(infile)s %(controlfile)s 
                   | awk -v FS=":" '/Estimated noise rate/ { printf("estimated_noise_rate\\n%%f\\n", $2) }' > %(outfile)s'''
    P.run()

####################################################################
@merge( checkDataQuality, "data_quality.load" )
def loadDataQuality( infiles, outfile ):
    '''load data quality information.'''

    P.concatenateAndLoad( infiles, outfile, 
                          regex_filename = "(.*).data_quality",
                          has_titles = False)

####################################################################
@follows( normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "%s.library_complexity" % x.asFile() ) for x in TRACKS ] )
def checkLibraryComplexity( infile, outfile ):
    '''uses peakranger to check library complexity.'''
    
    to_cluster = True
    statement = '''peakranger lc --format bam %(infile)s > %(outfile)s'''
    P.run()


######################################################################
######################################################################
##                                                                  ##
##                    **** Call Peaks ****                          ## 
##                                                                  ##
######################################################################
######################################################################


######################################################################
######################################################################
##                                                                  ##
##                            MACS1.4                               ## 
##                                                                  ##
######################################################################
######################################################################

@follows( mkdir("macs.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "macs.dir/%s.macs" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithMACS( infile, outfile ):
    '''run MACS for peak detection.

    output bed files are compressed and indexed.
    '''
    track = P.snip( infile, ".call.bam" )

    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )

    PipelinePeakcalling.runMACS( infile, outfile, controlfile)

############################################################
@transform( callPeaksWithMACS,
            regex(r"(.*).macs"),
            r"\1_macs.load" )
def loadMACS( infile, outfile ):
    '''load macs results.'''
    bamfile, controlfile = getBamFiles( infile, ".macs" )
    PipelinePeakcalling.loadMACS( infile, outfile, bamfile, controlfile )

############################################################
@follows( mkdir(os.path.join( PARAMS["exportdir"], "macs" ) ) )
@transform( callPeaksWithMACS,
            regex(r"(.*)/(.*).macs"),
            add_inputs( os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_contigs"] ) ),
            (os.path.join( PARAMS["exportdir"], "macs", r"\2.macs.treat.bw" ) ,
            os.path.join( PARAMS["exportdir"], "macs", r"\2.macs.control.bw" ) ) )
def cleanMACS( infiles, outfiles ):
    '''clean up MACS - build bigwig file and remove wig files.'''
  
    to_cluster = True
    infile, contigfile = infiles
    outdir = os.path.join( PARAMS["exportdir"], "macs" )

    for indir, outfile in zip( 
        ( os.path.join( infile + "_MACS_wiggle", "treat" ),
          os.path.join( infile + "_MACS_wiggle", "control" )),
        outfiles ):

        if os.path.exists( indir ):

            statement = '''
        zcat %(indir)s/*.wig.gz 
        | awk '/track/ { if (skip) {next} skip=1; } { print }'
        | python %(scriptsdir)s/wig2wig.py 
                --method=sanitize-genome
                --log=%(outfile)s.log
                --genome=%(genome_dir)s/%(genome)s
        | wigToBigWig /dev/stdin %(contigfile)s %(outfile)s'''
            
            P.run()
        
            shutil.rmtree( indir )
        
############################################################
@merge( callPeaksWithMACS, "macs.summary" )
def summarizeMACS( infiles, outfile ):
    '''summarize MACS results.''' 
    PipelinePeakcalling.summarizeMACS( infiles, outfile )

############################################################
@merge( callPeaksWithMACS, "macs_fdr.summary" )
def summarizeMACSFDR( infiles, outfile ):
    '''summarize MACS results.''' 
    PipelinePeakcalling.summarizeMACSFDR( infiles, outfile )

############################################################
@transform( summarizeMACS,
            suffix(".summary"),
            "_summary.load" )
def loadMACSSummary( infile, outfile ):
    '''load macs summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
@transform( summarizeMACSFDR,
            suffix("_fdr.summary"),
            "_fdr.load" )
def loadMACSSummaryFDR( infile, outfile ):
    '''load macs summary.'''
    P.load( infile, outfile, "--index=track", transpose="fdr" )


######################################################################
######################################################################
##                                                                  ##
##                          MACS version 2                          ## 
##                                                                  ##
######################################################################
######################################################################

@follows( mkdir("macs2.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "macs2.dir/%s.macs2" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithMACS2( infile, outfile ):
    '''run MACS 2 for peak detection.
    output bed files are compressed and indexed.
    '''
    track = P.snip( infile, ".call.bam" )
    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )
    PipelinePeakcalling.runMACS2( infile, outfile, controlfile)

############################################################
@transform( callPeaksWithMACS2,
            regex(r"(.*).macs2"),
            r"\1_macs2.load" )
def loadMACS2( infile, outfile ):
    '''load macs results.'''
    bamfile, controlfile = getBamFiles( infile, ".macs2" )
    PipelinePeakcalling.loadMACS2( infile, outfile, bamfile, controlfile )

############################################################
@merge( callPeaksWithMACS2, "macs2.summary" )
def summarizeMACS2( infiles, outfile ):
    '''summarize MACS results.''' 
    PipelinePeakcalling.summarizeMACS2( infiles, outfile )

############################################################
@merge( callPeaksWithMACS2, "macs2_fdr.summary" )
def summarizeMACS2FDR( infiles, outfile ):
    '''summarize MACS results.''' 
    PipelinePeakcalling.summarizeMACS2FDR( infiles, outfile )

############################################################
@transform( summarizeMACS2,
            suffix(".summary"),
            "_summary.load" )
def loadMACS2Summary( infile, outfile ):
    '''load macs2 summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
@transform( summarizeMACS2FDR,
            suffix("_fdr.summary"),
            "_fdr.load" )
def loadMACS2SummaryFDR( infile, outfile ):
    '''load macs2 summary.'''
    P.load( infile, outfile, "--index=track", transpose="fdr" )

######################################################################
######################################################################
##                                                                  ##
##                              Zinba                               ## 
##                                                                  ##
######################################################################
######################################################################

@follows( mkdir("zinba.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "zinba.dir/%s.zinba" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithZinba( infiles, outfile ):
    '''run Zinba for peak detection.'''
    infile, controlfile = infiles

    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )

    if os.path.exists( os.path.join( outfile + "_files" , outfile + ".model")):
        PipelinePeakcalling.runZinba( infile, 
                                  outfile, 
                                  controlfile, 
                                  action = "model" )

    elif os.path.exists( os.path.join( outfile + "_files" , outfile + ".list")):
        PipelinePeakcalling.runZinba( infile, 
                                  outfile, 
                                  controlfile, 
                                  action = "predict" )

    else:
        PipelinePeakcalling.runZinba( infile, outfile, controlfile, action = "full" )
        
############################################################
@transform( callPeaksWithZinba,
            suffix( ".zinba" ),
            "_zinba.load" )
def loadZinba( infile, outfile ):
    '''load Zinba results.'''

    bamfile, controlfile = getBamFiles( infile, ".zinba" )
    PipelinePeakcalling.loadZinba( infile, outfile, 
                                   bamfile, 
                                   controlfile = controlfile )


######################################################################
######################################################################
##                                                                  ##
##                       (SICER) Narrower                           ## 
##                                                                  ##
######################################################################
######################################################################
@follows( mkdir("sicer.narrow.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "sicer.narrow.dir/%s.narrow.sicer" % x.asFile() ) for x in TRACKS ] )
def callNarrowerPeaksWithSICER( infile, outfile ):
    '''run SICER for peak detection.'''
    track = P.snip( infile, ".call.bam" )

    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )
    PipelinePeakcalling.runSICER( infile, outfile, controlfile, "narrow")


######################################################################
######################################################################
##                                                                  ##
##                       (SICER) Broader                            ## 
##                                                                  ##
######################################################################
######################################################################
@follows( mkdir("sicer.broad.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "sicer.broad.dir/%s.broad.sicer" % x.asFile() ) for x in TRACKS ] )
def callBroaderPeaksWithSICER( infile, outfile ):
    '''run SICER for peak detection.'''
    track = P.snip( infile, ".call.bam" )

    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )
    PipelinePeakcalling.runSICER( infile, outfile, controlfile, "broad")

######################################################################
######################################################################
##                                                                  ##
##                  SICER Narrower + Broader                        ## 
##                                                                  ##
######################################################################
######################################################################
#@transform( [ callNarrowerPeaksWithSICER, callBroaderPeaksWithSICER ], suffix(".sicer"), "_sicer.load" )
@transform( [ callNarrowerPeaksWithSICER, callBroaderPeaksWithSICER ], regex(r"(sicer.)(.*)(.dir/)([^.]*).([^.]*).sicer"), r"\1\2\3\4_\5Sicer.load" )
def loadSICER(infile, outfile ):
    '''load sicer results.'''
    mode = infile.split(".")[1]
    bamfile, controlfile = getBamFiles( infile, "."+mode+".sicer" )
    PipelinePeakcalling.loadSICER( infile, outfile, bamfile, controlfile, mode)

############################################################
@merge( [ callNarrowerPeaksWithSICER, callBroaderPeaksWithSICER ], "sicer.summary" )
def summarizeSICER( infiles, outfile ):
    '''summarize SICER results.'''
    PipelinePeakcalling.summarizeSICER( infiles, outfile )

############################################################
@transform( summarizeSICER, regex(r"(sicer.)(.*)(.summary)"), r"\1_\2_summary.load" )
def loadSICERSummary( infile, outfile ):
    '''load sicer summary.'''
    P.load( infile, outfile, "--index=track" )

######################################################################
######################################################################
##                                                                  ##
##                       (PeakRanger) Ranger                        ## 
##                                                                  ##
######################################################################
######################################################################
@follows( mkdir("peakranger.ranger.dir/"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "peakranger.ranger.dir/%s.peakranger" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithPeakRanger( infile, outfile ):
    '''run PeakRanger Ranger for peak detection.'''
    track = P.snip( infile, ".call.bam" )
    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )
    PipelinePeakcalling.runPeakRanger( infile, outfile, controlfile)

@transform( callPeaksWithPeakRanger,
            regex(r"(.*).peakranger"),
            r"\1_peakranger.load" )
def loadPeakRanger( infile, outfile ):
    '''load macs results.''' 
    bamfile, controlfile = getBamFiles(infile, ".peakranger")
    PipelinePeakcalling.loadPeakRanger(infile, 
                                       outfile, 
                                       bamfile, 
                                       controlfile ,
                                       "peaks")

############################################################
@merge( callPeaksWithPeakRanger, "peakranger.ranger.summary" )
def summarizePeakRanger( infiles, outfiles ):
    '''summarize peakranger results.'''
    PipelinePeakcalling.summarizePeakRanger( infiles, outfiles )

############################################################
@transform( summarizePeakRanger, suffix(".summary"), "_summary.load" )
def loadPeakRangerSummary( infile, outfile ):
    '''load Peakranger summarys.'''
    P.load( infile, outfile, "--index=track" )

######################################################################
######################################################################
##                                                                  ##
##                       (PeakRanger) CCAT                          ## 
##                                                                  ##
######################################################################
######################################################################

@follows( mkdir("peakranger.ccat.dir/"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "peakranger.ccat.dir/%s.ccat" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithPeakRangerCCAT( infile, outfile ):

    '''run Peak Ranger CCAT for broad peak detection.'''
    track = P.snip( infile, ".call.bam" )
    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )
    PipelinePeakcalling.runPeakRangerCCAT( infile, outfile, controlfile)

##########################################################
@transform( callPeaksWithPeakRangerCCAT,
            regex(r"(.*).ccat"),
            r"\1_ccat.load" )
def loadPeakRangerCCAT( infile, outfile ):
    '''load macs results.''' 
    bamfile,controlfile = getBamFiles( infile,".ccat")
    PipelinePeakcalling.loadPeakRanger( infile, outfile, bamfile, controlfile , "regions")

############################################################
@merge( callPeaksWithPeakRangerCCAT, "peakranger.ccat.summary" )
def summarizePeakRangerCCAT( infiles, outfiles ):
    '''summarize peakranger results.'''
    PipelinePeakcalling.summarizePeakRanger( infiles, outfiles )

############################################################
@transform( summarizePeakRanger, suffix(".summary"), "_summary.load" )
def loadPeakRangerSummaryCCAT( infile, outfile ):
    '''load Peakranger summarys.'''
    P.load( infile, outfile, "--index=track" )

######################################################################
######################################################################
##                                                                  ##
##                           BroadPeak                              ## 
##                                                                  ##
######################################################################
######################################################################

@follows( mkdir( "broadpeak.dir" ), normalizeBAM ) 
@files( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] + ".fasta" ),
        "broadpeak.dir/genome_windows.bed.gz" )
def buildBroadPeakGenomeWindows( infile, outfile ):
    windows = PARAMS["broadpeak_genome_windows"]
    tmpf = P.getTempFilename( "." )
    PipelinePeakcalling.createGenomeWindows( infile, outfile, windows )


@transform( [ "%s.call.bam" % x.asFile() for x in TRACKS ], 
            regex( "(.+)/(.+).call.bam" ), 
            add_inputs( buildBroadPeakGenomeWindows ),
            r"./broadpeak.dir/\2.bedgraph" )
def buildBroadPeakBedgraphFiles( infiles, outfile ):
    logfile = outfile + ".log"
    overlap = str( float( PARAMS["broadpeak_read_length"] )
                   / float( 2*PARAMS["broadpeak_genome_windows"] ) )
    if PARAMS["broadpeak_remove_background"]:
        remove_background = "true"
        track = P.snip( infile, ".call.bam" )
        controls = getControl(Sample(track))
        controlfile = getControlFile( controls, "%s.call.bam" )
        infiles.append( controlfile )
    else: 
        remove_background = "false"


    params = [ overlap, remove_background ]

    P.info( "Creating BroadPeak bedgraph with the following options:\n"
            " subtract INPUT from SAMPLE = %s\n"
            " create .bdg file with window size of %s\n"
            " reads must overlap window by %s bases to be counted"  
            % ( params[1], 
                PARAMS["broadpeak_genome_windows"], 
                params[0] ) )

    P.submit( "/ifs/devel/projects/proj010/PipelineProj010", 
              "createBroadPeakBedgraphFile", 
              params, 
              infiles, 
              outfile, 
              toCluster = True, 
              logfile = logfile,
              jobOptions = "-l mem_free=20G" )

@transform( buildBroadPeakBedgraphFiles, 
            regex( "(.+)/(.+).bedgraph" ),
            r"\1/\2" )
def callPeaksWithBroadPeak( infile, outfile ):
    if isinstance(PARAMS["broadpeak_genome_size"], int):
        genome_size = PARAMS["broadpeak_genome_size"]
    elif PARAMS["broadpeak_genome_size"].startswith( "hg" ):
        genome_size = 3107677273
    elif PARAMS["broadpeak_genome_size"].startswith( "mm" ):
        genome_size = 2730871774
    else:
        ValueError( "Broadpeak genome size not recognised,"
                    " try entering a numeric value" )

    training_set = PARAMS["broadpeak_intervals"]

    logfile = P.snip( infile, ".bedgraph" ) + ".log"

    PipelinePeakcalling.runBroadPeak( infile, 
                                      outfile, 
                                      logfile, 
                                      genome_size, 
                                      training_set)

@merge( callPeaksWithBroadPeak, "./broadpeak.dir/broadpeak.stats" )
def summarizeBroadPeak( infiles, outfile):
    if PARAMS["broadpeak_intervals"]:
        PipelinePeakcalling.summarizeBroadPeak( infiles, outfile, intervals)
    else:
        PipelinePeakcalling.summarizeBroadPeak( infiles, outfile )

@transform( summarizeBroadPeak, suffix( ".stats" ), ".load" )
def loadBroadPeak( infile, outfile ):
    P.load( infile, outfile )

######################################################################
######################################################################
##                                                                  ##
##                              SPP                                 ## 
##                                                                  ##
######################################################################
######################################################################

@follows( mkdir("spp.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "spp.dir/%s.spp" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithSPP( infile, outfile ):
    '''run SICER for peak detection.'''
    track = P.snip( infile, ".call.bam" )
    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )
    PipelinePeakcalling.runSPP( infile, outfile, controlfile)

############################################################
@transform( callPeaksWithSPP,
            regex(r"(.*).spp"),
            r"\1_spp.load" )
def loadSPP( infile, outfile ):
    '''load spp results.'''
    bamfile, controlfile = getBamFiles( infile, ".spp" )
    PipelinePeakcalling.loadSPP( infile, outfile, bamfile, controlfile )

############################################################
@merge( callPeaksWithSPP, "spp.summary" )
def summarizeSPP( infiles, outfile ):
    '''summarize SPP results.'''
    PipelinePeakcalling.summarizeSPP( infiles, outfile )

############################################################
@transform( summarizeSPP, suffix(".summary"), "_summary.load" )
def loadSPPSummary( infile, outfile ):
    '''load sicer summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
############################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "quality" ) ),
          mkdir( "spp.dir" ),
          normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "spp.dir/%s.qual" % x.asFile() ) for x in TRACKS ] )
def estimateSPPQualityMetrics( infile, outfile ):
    '''estimate ChIP-Seq quality metrics using SPP'''

    to_cluster = True
    job_options= "-l mem_free=4G"
    track = P.snip(infile, ".call.bam" )
    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )
    if controlfile is None:
        raise ValueError( "idr analysis requires a control")
    
    executable = P.which( "run_spp.R" )
    if executable == None:
        raise ValueError( "could not find run_spp.R" )

    statement = '''
    Rscript %(executable)s -c=%(infile)s -i=%(controlfile)s -rf -savp -out=%(outfile)s
    >& %(outfile)s.log'''
    
    P.run()

    if os.path.exists( track + ".pdf" ):
        dest = os.path.join( PARAMS["exportdir"], "quality", track + ".pdf" )
        if os.path.exists( dest ): os.unlink( dest )
        shutil.move( track + ".pdf", dest )

############################################################
############################################################
############################################################
@merge( estimateSPPQualityMetrics, "spp_quality.load" )
def loadSPPQualityMetrics( infiles, outfile ):
    '''load spp quality metrics.'''
    P.concatenateAndLoad( infiles, outfile,
                          regex_filename = "spp.dir/(.*).qual",
                          has_titles = False,
                          header = "track,bamfile,mapped_reads,estFragLen,corr_estFragLen,phantomPeak,corr_phantomPeak,argmin_corr,min_corr,nsc,rsc,quality")

######################################################################
######################################################################
##                                                                  ##
##        IDR Analysis (using relaxed SPP peak calling)             ## 
##                                                                  ##
######################################################################
######################################################################
@follows( mkdir("idr.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "idr.dir/%s.spp" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithSPPForIDR( infile, outfile ):
    '''run SICER for peak detection.'''
    track = P.snip( infile, ".call.bam" )
    controls = getControl(Sample(track))
    controlfile = getControlFile( controls, "%s.call.bam" )

    job_options= "-l mem_free=4G"

    to_cluster = True

    if controlfile is None:
        raise ValueError( "idr analysis requires a control")

    executable = P.which( "run_spp.R" )
    if executable == None:
        raise ValueError( "could not find run_spp.R" )

    statement = '''
    Rscript %(executable)s -c=%(infile)s -i=%(controlfile)s -npeak=%(idr_npeaks)s 
            -odir=idr.dir -savr -savp -rf -out=%(outfile)s
    >& %(outfile)s.log'''
    
    P.run()

    track = P.snip(infile, ".bam" )

    if os.path.exists( track + ".pdf" ):
        shutil.move( infile + ".pdf", os.path.join( PARAMS["exportdir"], "idr" ))

############################################################
@collate( callPeaksWithSPPForIDR, 
          regex( r"idr.dir/(.+)-[^-]+.spp" ),
          r"idr.dir/\1.idr")
def applyIDR( infiles, outfile ):
    '''apply IDR analysis.'''

    job_options= "-l mem_free=4G"

    to_cluster = True
    chromosome_table = os.path.join(PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_contigs"])

    for infile1, infile2 in itertools.combinations( infiles, 2 ):
        E.info( "applyIDR: processing %s and %s" % (infile1,infile2))

        basename1 = os.path.basename( infile1 )
        basename2 = os.path.basename( infile2 )

        track1 = P.snip( basename1, ".spp" )
        control1 = getControl(Sample(track1)).asFile()
        track2 = P.snip( basename2, ".spp" )
        control2 = getControl(Sample(track2)).asFile()
        
        statement = '''
          python %(scriptsdir)s/WrapperIDR.py 
                 --action=run
                 --output-prefix=%(track1)s_vs_%(track2)s.idr
                 --chromosome-table=%(chromosome_table)s
                 idr.dir/%(track1)s.call_VS_%(control1)s.call.regionPeak.gz 
                 idr.dir/%(track2)s.call_VS_%(control2)s.call.regionPeak.gz 
          >> %(outfile)s'''

        P.run()

############################################################
@follows(mkdir( os.path.join( PARAMS["exportdir"], "idr" ) ) )
@transform( applyIDR, 
            suffix(".idr"),
            ".plot")
def plotIDR( infile, outfile ):
    '''plot IDR results.'''

    to_cluster = True

    track = P.snip( infile, ".idr")
    files = glob.glob( track + "*.idr-em.sav" )
    files = " ".join([ P.snip(x, "-em.sav" ) for x in files ])
    output_prefix = os.path.join( PARAMS["exportdir"], "idr", os.path.basename(track ) )
    statement = '''
    python %(scriptsdir)s/WrapperIDR.py
               --action=plot
               --output-prefix=%(output_prefix)s
               %(files)s
    > %(outfile)s'''

    P.run()

######################################################################
######################################################################
##                                                                  ##
##                           Scripture                              ## 
##                                                                  ##
######################################################################
######################################################################
@follows( mkdir("scripture.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "scripture.dir/%s.scripture" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithScripture( infile, outfile ):
    '''run SICER for peak detection.'''
    track = P.snip( infile, ".call.bam" )
    controls = getControl(Sample(track))
    controlfile = getControlFile(controls, "%s.call.bam")

    contig_sizes = os.path.join( PARAMS["annotations_dir"],
                                 PARAMS_ANNOTATIONS["interface_contigs" ] )

    PipelinePeakcalling.runScripture( infile, 
                                      outfile, 
                                      contig_sizes)

@transform( callPeaksWithScripture, suffix(".scripture"), ".load")
def loadScripture( infile, outfile ):
    '''load scripture data.'''
    bamfile, controlfile = getBamFiles( infile, ".scripture" )
    PipelinePeakcalling.loadScripture( infile, outfile, bamfile, controlfile )

######################################################################
######################################################################
##                                                                  ##
##            ___ Collate peak calling results ___                  ## 
##                                                                  ##
######################################################################
######################################################################

CALLINGTARGETS, SUMMARYTARGETS = [], []
mapToCallingTargets = { 'macs': loadMACS,
                        'macs2' : loadMACS2,
                        'zinba': loadZinba,
                        'sicer': loadSICER,
                        'peakranger': loadPeakRanger,
                        'ccat' : loadPeakRangerCCAT ,
                        'spp': loadSPP,
                        'scripture' : loadScripture,
                        }

mapToSummaryTargets = { 'macs': [loadMACSSummary, loadMACSSummaryFDR],
                        'macs2': [loadMACS2Summary, loadMACS2SummaryFDR],
                        'sicer': [loadSICERSummary],
                        'spp' : [loadSPPSummary],
                        'peakranger' : [loadPeakRangerSummary],
                        'ccat' : [loadPeakRangerSummaryCCAT],
                        }

for x in P.asList( PARAMS["peakcallers"]):
    CALLINGTARGETS.append( mapToCallingTargets[x] )
    if x in mapToSummaryTargets:
        SUMMARYTARGETS.extend( mapToSummaryTargets[x] )

###################################################################
###################################################################
###################################################################
@follows( *(CALLINGTARGETS + SUMMARYTARGETS) )
def calling(): pass

############################################################
############################################################
############################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "bedfiles" ) ) )
@transform( CALLINGTARGETS, regex("(.*)/(.*).load"), 
            (os.path.join( PARAMS["exportdir"], "bedfiles", r"\2.peaks.bed.gz" ),
             os.path.join( PARAMS["exportdir"], "bedfiles", r"\2.regions.bed.gz" ),
             os.path.join( PARAMS["exportdir"], "bedfiles", r"\2.summits.bed.gz" ) ) )
def exportIntervalsAsBed( infile, outfiles ):
    '''export all intervals as bed files.'''

    outfile_peaks, outfile_regions, outfile_summits = outfiles
    track = P.snip( os.path.basename(infile), ".load" ) 

    #PipelinePeakcalling.exportIntervalsAsBed( infile, outfile_regions, "%s_regions" % P.quote(track) )

    dbh = connect()
    tablename = "%s_peaks" % P.quote(track) 
    if tablename in Database.getTables( dbh ):
        PipelinePeakcalling.exportIntervalsAsBed( infile, outfile_peaks, tablename )
    else:
        E.warn( "no table %s - empty bed file output" % tablename )
        P.touch( outfile_peaks )

    dbh = connect()
    tablename = "%s_regions" % P.quote(track) 
    if tablename in Database.getTables( dbh ):
        PipelinePeakcalling.exportIntervalsAsBed( infile, outfile_regions, tablename )
    else:
        E.warn( "no table %s - empty bed file output" % tablename )
        P.touch( outfile_regions )

    dbh = connect()
    tablename = "%s_summits" % P.quote(track) 
    if tablename in Database.getTables( dbh ):
        PipelinePeakcalling.exportIntervalsAsBed( infile, outfile_summits, tablename )
    else:
        E.warn( "no table %s - empty bed file output" % tablename )
        P.touch( outfile_summits )

###################################################################
###################################################################
###################################################################
# Targets for the annotation of intervals.
###################################################################
@split( exportIntervalsAsBed, os.path.join( PARAMS["exportdir"], "bedfiles", "*.bed.gz" ) )
def flattenBedFiles( infile, outfile ):
    '''dummy target - merge all files in exportIntervalsAsBed'''


def getPeakShift( track, method ):
    '''return peak shift for track and method.'''
    dbh = connect()
    result = Database.executewait( dbh, "SELECT shift FROM %(method)s_summary where track = '%(track)s'" % locals())
    return result.fetchone()[0]

###################################################################
###################################################################
###################################################################
@follows( mkdir( "peakshapes.dir" ) )
@transform( flattenBedFiles,
            regex(".*/(.*).bed.gz"),
            r"peakshapes.dir/\1.peakshape.tsv.gz" )
def buildPeakShapeTable( infile, outfile ):
    '''build a table with peak shape parameters.'''
    
    to_cluster = True

    # compute suffix (includes method name)
    track, method, section = re.match( "(.*)_(.*)\.(.*).bed.gz", os.path.basename(infile) ).groups()
    
    suffix = "_%s.%s.bed.gz" % (method, section)
    bamfile, controlfile = getBamFiles( infile, suffix )
   
    shift = getPeakShift( track, method )

    if shift:
        E.info( "applying read shift %i for track %s" % (shift, track ) )

    options = []
    if controlfile:
        options.append( "--control-file=%s" % controlfile )
    options = " ".join( options )

    statement = '''python %(scriptsdir)s/bam2peakshape.py
                      --window-size=%(peakshape_window_size)i
                      --bin-size=%(peakshape_bin_size)i
                      --output-filename-pattern="%(outfile)s.%%s"
                      --force
                      --shift=%(shift)i
                      --sort=peak-height
                      --sort=peak-width
                      %(options)s
                      --log=%(outfile)s.log
                      %(bamfile)s %(infile)s
                   | gzip
                   > %(outfile)s
                '''
    P.run()

###################################################################
@transform( buildPeakShapeTable, suffix(".tsv.gz"), ".load" )
def loadPeakShapeTable( infile, outfile ):
    '''load peak shape information.'''
    P.load( infile, outfile, "--ignore-column=bins --ignore-column=counts --allow-empty" )

############################################################
############################################################
############################################################
## targets to do with the analysis of replicates
############################################################
# dummy task - flatten the nested list of result files created
# by exportIntervalsAsBed
@split( exportIntervalsAsBed, os.path.join( PARAMS["exportdir"], "bedfiles", "*.bed.gz" ) )
def allIntervalsAsBed( infile, outfile ): pass

############################################################
############################################################
############################################################
@follows( mkdir( "reproducibility.dir"))
@collate( allIntervalsAsBed,
          regex( os.path.join( PARAMS["exportdir"], "bedfiles", r"(.+)_(.+)\.(.+).bed.gz" )),
          r"reproducibility.dir/\1.\3.reproducibility")
def makeReproducibilityOfMethods( infiles, outfile ):
    '''compute overlap between intervals.

    Note the exon percentages are approximations assuming that there are
    not more than one intervals in one set overlapping one in the other set.
    '''
    PipelinePeakcalling.makeReproducibility( infiles, outfile )

############################################################
############################################################
############################################################
@follows( mkdir( "reproducibility.dir"))
@collate( allIntervalsAsBed, 
          regex( os.path.join( PARAMS["exportdir"], "bedfiles", r"(.+)-[^-]+_(.+)\.(.+).bed.gz" )),
          r"reproducibility.dir/\1-\2.\3.reproducibility")
def makeReproducibilityOfReplicates( infiles, outfile ):
    '''compute overlap between intervals.

    Note the exon percentages are approximations assuming that there are
    not more than one intervals in one set overlapping one in the other set.
    '''
    PipelinePeakcalling.makeReproducibility( infiles, outfile )

############################################################
############################################################
############################################################
@transform( (makeReproducibilityOfMethods, makeReproducibilityOfReplicates), suffix(".reproducibility"), "_reproducibility.load" )
def loadReproducibility( infile, outfile ):
    '''load Reproducibility results
    '''
    P.load( infile, outfile, options="--allow-empty" )

############################################################
############################################################
############################################################
@follows( loadReproducibility )
def reproducibility(): pass


###################################################################
@follows(loadBAMStats,
         loadDuplicationStats,
         loadSPPQualityMetrics)
def qc():
    pass

###################################################################
@follows( calling, exportIntervalsAsBed, qc )
def full():
    pass

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

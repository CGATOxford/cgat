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

pipeline_peakcalling takes as input reads aligned to genomic sequence as :term:`bam` formatted files
and calls peaks. The pipeline implements several peak callers:

macs_
   Model-based Analysis of ChIP-Seq (MACS), for identifying transcript factor binding sites. MACS captures the influence of genome complexity to evaluate the significance of enriched ChIP regions, and MACS improves the spatial resolution of binding sites through combining the information of both sequencing tag position and orientation. MACS can be easily used for ChIP-Seq data alone, or with control sample with the increase of specificity.

spp_
   SPP is a R package especially designed for the analysis of Chip-Seq data from Illummina platform. The package was developed by Peter Park's group from Harvard Medical School.

zinba_
   ZINBA (Zero Inflated Negative Binomial Algorithm) is a computational and statistical framework used to call regions of the genome enriched for sequencing reads originating from a diverse array of biological experiments. We collectively refer to the sequencing data derived from these experiments as DNA-seq, including FAIRE-seq, ChIP-seq, and DNAase-seq experiments

sicer_
    A clustering approach for identification of enriched domains from histone modification ChIP-Seq data.

peakranger_
    PeakRanger is a multi-purpose, ultrafast ChIP Seq peak caller. It is used in the modENCODE project and included in the iPlant pipeline system.
    PeakRanger v1.02 was developed in Dr.Lincoln Stein's lab at OICR and is now in continual development at Dr.Helen Hobbs's lab of the McDermott Center of UT Southwestern.


Peak callers have different strengths and weaknesses. Some might work well on broad peaks such as some histone
marks, others work better for narrow, sharp peaks. Many callers these days attempt to call both types of peaks.
This pipeline implements the following nomenclature for peaks:


.. glossary::

   region
      A broad region defined by a peak caller. Regions are usually created in the first step of
      peak calling, before peak refinement or detection of subpeaks is performed.

   summit
      A narrow region defined by a caller. These are the output of any peak refinement
      or subpeak detection steps

   interval
      Generic term for an interval in which the read density is higher then expected. Both
      a :term:`region` or a :term:`summit` are an :term:`interval`
      
   peak
      Within a :term:`region` or :term:`summit` the position with the highest base coverage.

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
.. _spp: http://compbio.med.harvard.edu/Supplements/ChIP-seq/tutorial.html
.. _sicer: http://home.gwu.edu/~wpeng/Software.htm
.. _zinba: http://code.google.com/p/zinba/
.. _peakranger: http://ranger.sourceforge.net/

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

import PipelinePeakcalling as PipelinePeakcalling
import PipelineMotifs as PipelineMotifs
import PipelineGeneset as PGeneset
import PipelineTracks
import PipelineMapping

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
import Pipeline as P
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
import PipelineTracks

Sample = PipelineTracks.Sample3

# collect tracks, exclude any control tracks
TRACKS = PipelineTracks.Tracks( Sample ).loadFromDirectory( 
    [ x for x in glob.glob( "*.genome.bam" ) if PARAMS["tracks_control"] not in x ],
    "(\S+).genome.bam" )

ALL = Sample()
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

###################################################################
def getControl( track ):
    '''return appropriate control for a track
    '''
    n = track.clone()
    n.condition = PARAMS["tracks_control"]
    n.replicate = "R1"
    return n

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
    assert os.path.exists( bamfile )

    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()
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
############o################################################
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


####################################################################
####################################################################
####################################################################
## Peak calling
####################################################################
@follows( mkdir("macs.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "macs.dir/%s.macs" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithMACS( infile, outfile ):
    '''run MACS for peak detection.

    output bed files are compressed and indexed.
    '''
    track = P.snip( infile, ".call.bam" )

    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        controlfile = None

    PipelinePeakcalling.runMACS( infile, outfile, controlfile)

############################################################
############################################################
############################################################
@transform( callPeaksWithMACS,
            regex(r"(.*).macs"),
            r"\1_macs.load" )
def loadMACS( infile, outfile ):
    '''load macs results.'''
    bamfile, controlfile = getBamFiles( infile, ".macs" )
    PipelinePeakcalling.loadMACS( infile, outfile, bamfile, controlfile )
        
############################################################
############################################################
############################################################
@merge( callPeaksWithMACS, "macs.summary" )
def summarizeMACS( infiles, outfile ):
    '''summarize MACS results.''' 
    PipelinePeakcalling.summarizeMACS( infiles, outfile )

############################################################
############################################################
############################################################
@merge( callPeaksWithMACS, "macs_fdr.summary" )
def summarizeMACSFDR( infiles, outfile ):
    '''summarize MACS results.''' 
    PipelinePeakcalling.summarizeMACSFDR( infiles, outfile )

############################################################
############################################################
############################################################
@transform( summarizeMACS,
            suffix(".summary"),
            "_summary.load" )
def loadMACSSummary( infile, outfile ):
    '''load macs summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
############################################################
############################################################
@transform( summarizeMACSFDR,
            suffix("_fdr.summary"),
            "_fdr.load" )
def loadMACSSummaryFDR( infile, outfile ):
    '''load macs summary.'''
    P.load( infile, outfile, "--index=track", transpose="fdr" )

############################################################
############################################################
############################################################
## Zinba
############################################################

############################################################
############################################################
############################################################
@follows( mkdir("zinba.dir"), normalizeBAM )
@files( [ (("%s.call.bam" % (x.asFile()), 
            "%s.call.bam" % (getControl(x).asFile())), 
           "zinba.dir/%s.zinba" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithZinba( infiles, outfile ):
    '''run Zinba for peak detection.'''
    infile, controlfile = infiles

    if os.path.exists( os.path.join( outfile + "_files" , outfile + ".model")):
        PipelinePeakcalling.runZinba( infile, outfile, controlfile, action = "predict" )
    elif os.path.exists( os.path.join( outfile + "_files" , outfile + ".list")):
        PipelinePeakcalling.runZinba( infile, outfile, controlfile, action = "model" )
    else:
        PipelinePeakcalling.runZinba( infile, outfile, controlfile, action = "full" )
        
############################################################
############################################################
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

############################################################
############################################################
############################################################
## SICER
############################################################
@follows( mkdir("sicer.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "sicer.dir/%s.sicer" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithSICER( infile, outfile ):
    '''run SICER for peak detection.'''
    track = P.snip( infile, ".call.bam" )

    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        controlfile = None

    PipelinePeakcalling.runSICER( infile, outfile, controlfile)

############################################################
############################################################
############################################################
@transform( callPeaksWithSICER, suffix(".sicer"), "_sicer.load" )
def loadSICER(infile, outfile ):
    '''load sicer results.'''
    bamfile, controlfile = getBamFiles( infile, ".sicer" )
    PipelinePeakcalling.loadSICER( infile, outfile, bamfile, controlfile )

############################################################
############################################################
############################################################
@merge( callPeaksWithSICER, "sicer.summary" )
def summarizeSICER( infiles, outfile ):
    '''summarize SICER results.'''
    PipelinePeakcalling.summarizeSICER( infiles, outfile )

############################################################
############################################################
############################################################
@transform( summarizeSICER, suffix(".summary"), "_summary.load" )
def loadSICERSummary( infile, outfile ):
    '''load sicer summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
############################################################
############################################################
## PeakRanger
############################################################
@follows( mkdir("peakranger.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "peakranger.dir/%s.peakranger" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithPeakRanger( infile, outfile ):
    '''run SICER for peak detection.'''
    track = P.snip( infile, ".call.bam" )
    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        controlfile = None

    PipelinePeakcalling.runPeakRanger( infile, outfile, controlfile)

############################################################
############################################################
############################################################
@transform( callPeaksWithPeakRanger,
            regex(r"(.*).peakranger"),
            r"\1_peakranger.load" )
def loadPeakRanger( infile, outfile ):
    '''load macs results.'''
    bamfile, controlfile = getBamFiles( infile, ".peakranger" )
    PipelinePeakcalling.loadPeakRanger( infile, outfile, bamfile, controlfile )

############################################################
############################################################
############################################################
## PeakRanger
############################################################
@follows( mkdir("spp.dir"), normalizeBAM )
@files( [ ("%s.call.bam" % (x.asFile()), 
           "spp.dir/%s.spp" % x.asFile() ) for x in TRACKS ] )
def callPeaksWithSPP( infile, outfile ):
    '''run SICER for peak detection.'''
    track = P.snip( infile, ".call.bam" )
    controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()

    if not os.path.exists( controlfile ):
        L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
        controlfile = None

    PipelinePeakcalling.runSPP( infile, outfile, controlfile)

############################################################
############################################################
############################################################
@transform( callPeaksWithSPP,
            regex(r"(.*).spp"),
            r"\1_spp.load" )
def loadSPP( infile, outfile ):
    '''load spp results.'''
    bamfile, controlfile = getBamFiles( infile, ".spp" )
    PipelinePeakcalling.loadSPP( infile, outfile, bamfile, controlfile )

############################################################
############################################################
############################################################
@merge( callPeaksWithSPP, "spp.summary" )
def summarizeSPP( infiles, outfile ):
    '''summarize SPP results.'''
    PipelinePeakcalling.summarizeSPP( infiles, outfile )

############################################################
############################################################
############################################################
@transform( summarizeSPP, suffix(".summary"), "_summary.load" )
def loadSPPSummary( infile, outfile ):
    '''load sicer summary.'''
    P.load( infile, outfile, "--index=track" )

############################################################
############################################################
############################################################
CALLINGTARGETS, SUMMARYTARGETS = [], []
mapToCallingTargets = { 'macs': loadMACS,
                        'zinba': loadZinba,
                        'sicer': loadSICER,
                        'peakranger': loadPeakRanger,
                        'spp': loadSPP,
                        }

mapToSummaryTargets = { 'macs': [loadMACSSummary, loadMACSSummaryFDR],
                        'sicer': [loadSICERSummary],
                        'spp' : [loadSPPSummary],
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
@transform( CALLINGTARGETS, regex(".load"), ".bed.gz" )
def exportIntervalsAsBed( infile, outfile ):
    '''export all intervals as bed files.'''
    PipelinePeakcalling.exportIntervalsAsBed( infile, outfile )

###################################################################
@follows( calling )
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

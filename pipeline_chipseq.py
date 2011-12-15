################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_chipseq.py 2900 2010-04-13 14:38:00Z andreas $
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
=================
ChIP-Seq pipeline
=================

:Author: Andreas Heger
:Release: $Id: pipeline_chipseq.py 2900 2010-04-13 14:38:00Z andreas $
:Date: |today|
:Tags: Python

The ChIP-Seq pipeline imports reads from one or more ChIP-Seq experiments and
performs the following tasks:

   * align reads to the genome
   * call peaks 
   * annotate intervals with respect to a reference gene set
   * describe de-novo motifs
   * find motifs within intervals

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. The pipeline looks for a configuration file in several places:

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

    python <codedir>/pipeline_chipseq.py config

The following sections and parameters probably should be changed from the default 
values:

.. todo::
   describe important parameters

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.


Input
-----

Reads
++++++

Input are :file:`.export.txt.gz`-formatted files from Illumina, :term:`fastq.gz` files,
or :term:`csfasta.gz` files.


The files should be labeled in the following way::

   sample-condition-replicate.<suffix>.gz

For example::
 
   GM00855-D3-R1.<suffix>.gz
   GM00855-D3-R2.<suffix>.gz
   GM00855-input-R1.<suffix>.gz
   GM00855-unstim-R1.<suffix>.gz
   GM00855-unstim-R2.<suffix>.gz

Note that neither ``sample``, ``condition`` or ``replicate`` should contain 
``_`` (underscore) and ``.`` (dot) characters as these are used by the pipeline
to delineate tasks.

Optional inputs
+++++++++++++++

Optinally, peaks can be supplied as :term:`bed` formatted files. These peak files
will then be processed in the same way as peaks called within the pipeline. Use the
option ``tracks_extra`` to declare any additional tracks.

Additional peak files can be associated with one of the :term:`bam` files created
by the pipeline. This permits counting the number of tags inside peaks, finding the
peak summit, etc. In order to associated a peak file with a :term:`bam` formatted
file, define a section in the pipeline.ini file. For example::


   [mycalls.bed]
   track=tissue-sample-agg

will process the file ``mycalls.bed`` exactly the same way as the track ``tissue-sample-agg``.
Replicates can be specified explicitely::

   [mycalls.bed]
   replicates=tissue1-sample1-R1,tissue1-sample1-R2

will associate the file ``mycalls.bed`` with the replicates ``tissue1-sample1-R1`` 
and ``tissue1-sample1-R2``. Note that globs don't work for this yet, all
replicates have to be specified explicitely.

Reference motifs
++++++++++++++++

Reference motifs are described by fasta-formatted multiple alignments, see for
example Jaspar download. The motifs are build by running MEME on the file.

Reference motifs should end in the suffix ".motif.fasta", for example,
:file:`rxrvdr.motif.fasta`.

Requirements
------------

The pipeline requires the information from the following pipelines:

:doc:`pipeline_annotations`
   set the configuration variable :py:data:`annotations_database` and 
   :py:data:`annotations_dir`.

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

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_chipseq.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_chipseq.tgz
   tar -xvzf pipeline_chipseq.tgz
   cd pipeline_chipseq
   python <srcdir>/pipeline_chipseq.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` as well.

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
import Database
from ruffus import *
import csv
import sqlite3
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import IOTools
import GTF, GFF, Bed
import pysam
import numpy
import gzip

import PipelineChipseq as PipelineChipseq
import PipelineMotifs as PipelineMotifs
import PipelineGeneset as PGeneset
import PipelineTracks
import PipelineMapping

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
import Pipeline as P
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ],
    defaults = {
        'annotations_dir' : "" })

PARAMS = P.PARAMS

PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

###################################################################
###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
# load all tracks - exclude input/control tracks
Sample = PipelineTracks.Sample3

suffixes = ["export.txt.gz",
            "sra",
            "fastq.gz",
            "cfastq.1.gz",
            "csfasta.gz" ]

TRACKS = sum( itertools.chain( [ PipelineTracks.Tracks( Sample ).loadFromDirectory( 
        [ x for x in glob.glob( "*.%s" % s ) if PARAMS["tracks_control"] not in x ],
        "(\S+).%s" % s ) for s in suffixes ] ), 
              PipelineTracks.Tracks( Sample ) )

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

# tracks for subtraction of unstim condition
if "tracks_subtract" in PARAMS and PARAMS["tracks_subtract"]:
    TOSUBTRACT = [ x for x in EXPERIMENTS \
                       if not x.condition == PARAMS["tracks_unstimulated"] and \
                       getUnstimulated( x ) in EXPERIMENTS ]
else:
    TOSUBTRACT = []

# compound targets : correlation between tracks
TRACKS_CORRELATION = TRACKS_MASTER + list(TRACKS) + [ getSubtracted( x ) for x in TOSUBTRACT ]

# print "EXP=", EXPERIMENTS
# print "COND=", CONDITIONS
# print "TISSUES=", TISSUES
# print "TOSUBTRACT=", TOSUBTRACT
# print "MASTER=", TRACKS_MASTER
# print "CORRELATION=", TRACKS_CORRELATION

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
## General preparation tasks
###################################################################

###################################################################
###################################################################
###################################################################
## load number of reads
###################################################################
@transform( ("*.export.txt.gz",
             "*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra",
             "*.csfasta.gz" ),
            regex( r"(\S+).(export.txt.gz|fastq.1.gz|fastq.gz|sra|csfasta.gz)"), 
            r"\1.nreads" )
def countReads( infile, outfile ):
    '''count number of reads in input files.'''
    to_cluster = True
    m = PipelineMapping.Counter()
    statement = m.build( (infile,), outfile ) 
    P.run()

@merge(countReads, "reads_summary.load" )
def loadReadCounts( infiles, outfile ):
    '''load read counts into database.'''

    outf = P.getTempFile()
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
## Mapping
###################################################################

if PARAMS["mapping_mapper"] == "bowtie":
    
    ############################################################
    ############################################################
    ############################################################
    @transform( ("*.export.txt.gz",
                 "*.fastq.1.gz", 
                 "*.fastq.gz",
                 "*.sra",
                 "*.csfasta.gz" ),
                regex( r"(\S+).(export.txt.gz|fastq.1.gz|fastq.gz|sra|csfasta.gz)"), 
                r"\1.genome.bam" )
    def buildBAM( infile, outfile ):
        '''re-map eland formatted reads with bowtie
        '''
        to_cluster = True

        job_options= "-pe dedicated %i -R y" % PARAMS["bowtie_threads"]
        m = PipelineMapping.Bowtie()
        reffile = PARAMS["samtools_genome"]
        statement = m.build( (infile,), outfile ) 
        P.run()

elif PARAMS["mapping_mapper"] == "bwa":
    
    ############################################################
    ############################################################
    ############################################################
    @transform( ("*.export.txt.gz",
                 "*.fastq.1.gz", 
                 "*.fastq.gz",
                 "*.sra",
                 "*.csfasta.gz" ),
                regex( r"(\S+).(export.txt.gz|fastq.1.gz|fastq.gz|sra|csfasta.gz)"), 
                r"\1.genome.bam" )
    def buildBAM( infile, outfile ):
        '''re-map eland formatted reads with bowtie
        '''
        to_cluster = True

        job_options= "-pe dedicated %i -R y" % PARAMS["bwa_threads"]
        m = PipelineMapping.BWA()
        statement = m.build( (infile,), outfile ) 
        P.run()

else:
    raise ValueError("unknown mapper %s" % PARAMS["mapping_mapper"] )

###################################################################
###################################################################
###################################################################
## Read procesing
###################################################################

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
            regions_to_filter += PipelineChipseq.getExonLocations(PARAMS["calling_filter_exons"]) 

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
@transform( buildBAM, suffix(".genome.bam"), add_inputs( makeMask), ".prep.bam")
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

    PipelineChipseq.buildBAMforPeakCalling(bam_file,outfile,PARAMS["calling_deduplicate"], mask )

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
        PipelineChipseq.buildSimpleNormalizedBAM( infiles, 
                                                  outfile,
                                                  minreads)
else:
    @transform(prepareBAMForPeakCalling,
               suffix(".prep.bam"),
               ".call.bam")
    def normalizeBAM( infile, outfile):
        P.clone( infile, outfile )
        P.clone( infile + ".bai", outfile + ".bai" )

############################################################
############################################################
############################################################
@follows( normalizeBAM )
@files( [( [ "%s.call.bam" % (y.asFile()) for y in EXPERIMENTS[x]], 
           "%s.readcorrelations.gz" % x.asFile()) 
         for x in EXPERIMENTS ] )
def makeReadCorrelation( infiles, outfile ):
    '''compute correlation between reads.
    '''

    # deactivated - solve memory problems first
    P.touch( outfile )
    return

    to_cluster = True
    # job_options = "-l mem_free=8000M"

    infiles = " ".join(infiles)

    statement = '''
    python %(scriptsdir)s/bam_correlation.py 
           --log=%(outfile)s.log 
           --genome=%(genome_dir)s/%(genome)s 
           %(infiles)s 
    | gzip > %(outfile)s
    ''' 
    
    P.run()

############################################################
############################################################
############################################################
@merge( makeReadCorrelation, "readcorrelation.table" )
def makeReadCorrelationTable( infiles, outfile ):
    '''compute correlation between reads.
    '''

    correlation_cutoff = 5
    outf = IOTools.openFile(outfile, "w" )
    outf.write("track\tcutoff\tpositions\tcoefficient\n" )

    if len(infiles) < 2:
        E.warn( "no replicates - no read correlation computed" )
        P.touch( outfile )

    def selector( infile ):
        correlation_cutoff = 5
        for line in infile:
            try:
                data = map(int, line[:-1].split("\t")[2:])
            except ValueError:
                continue
            
            # require at least two data points
            if len(data) < 2:
                break

            # only take position with minimum number of reads
            for x in range(len(data)):
                if data[x] < correlation_cutoff: break
            else:
                yield data

    for infile in infiles:
        mat = numpy.array( [ x for x in selector( IOTools.openFile(infile, "r") )] )
        if len(mat) == 0: continue

        corr = numpy.corrcoef(mat[:,0], mat[:,1])
        outf.write( "%s\t%i\t%i\t%f\n" % (
            infile,
            correlation_cutoff,
            len(mat),
            corr[0,1] ) )
        outf.flush()
        
    outf.close()

############################################################
############################################################
############################################################
@transform( (buildBAM, 
             prepareBAMForPeakCalling, 
             normalizeBAM),
            suffix(".bam"),
            ".readstats" )
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

#########################################################################
#########################################################################
#########################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''import bam statistics.'''

    header = ",".join( [P.snip( x, ".readstats") for x in infiles] )
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

####################################################################
####################################################################
####################################################################
## Peak calling
####################################################################
if PARAMS["calling_caller"] == "macs":

    ############################################################
    ############################################################
    ############################################################
    @follows( normalizeBAM )
    @files( [ ("%s.call.bam" % (x.asFile()), 
               "%s.macs" % x.asFile() ) for x in TRACKS ] )
    def runMACS( infile, outfile ):
        '''run MACS for peak detection.'''
        track = P.snip( infile, ".call.bam" )
        controlfile = "%s.call.bam" % getControl(Sample(track)).asFile()
        if not os.path.exists( controlfile ):
            L.warn( "no controlfile '%s' for track '%s' not found " % (controlfile, track ) )
            controlfile = None

        PipelineChipseq.runMACS( infile, outfile, controlfile)

    ############################################################
    ############################################################
    ############################################################
    @transform( runMACS,
                regex(r"(.*).macs"),
                inputs( (r"\1.macs", r"\1.call.bam") ),
                r"\1_macs.load" )
    def loadMACS( infiles, outfile ):
        infile, bamfile = infiles
        PipelineChipseq.loadMACS( infile, outfile, bamfile )
        
    ############################################################
    ############################################################
    ############################################################
    @merge( runMACS, "macs.summary" )
    def summarizeMACS( infiles, outfile ):
        '''summarize MACS results.''' 
        PipelineChipseq.summarizeMACS( infiles, outfile )

    @merge( runMACS, "macs_fdr.summary" )
    def summarizeMACSFDR( infiles, outfile ):
        '''summarize MACS results.''' 
        PipelineChipseq.summarizeMACSFDR( infiles, outfile )

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
    @follows( loadMACSSummary, loadMACSSummaryFDR )
    @transform( loadMACS, suffix("_macs.load"), ".bed.gz" )
    def exportIntervalsAsBed( infile, outfile ):
        PipelineChipseq.exportIntervalsAsBed( infile, outfile )

elif PARAMS["calling_caller"] == "zinba":

    ############################################################
    ############################################################
    ############################################################
    @follows( normalizeBAM )
    @files( [ (("%s.call.bam" % (x.asFile()), 
                "%s.call.bam" % (getControl(x).asFile())), 
               "%s.zinba" % x.asFile() ) for x in TRACKS ] )
    def runZinba( infiles, outfile ):
        '''run Zinba for peak detection.'''
        infile, controlfile = infiles
        PipelineChipseq.runZinba( infile, outfile, controlfile )
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runZinba,
                regex(r"(.*).zinba"),
                inputs( (r"\1.zinba", r"\1.call.bam" ) ),
                r"\1_zinba.load" )
    def loadZinba( infiles, outfile ):
        infile, bamfile = infiles
        track = P.snip( infile, ".zinba" )
        control = "%s.call.bam" % (getControl(Sample(track)).asFile())
        PipelineChipseq.loadZinba( infile, outfile, 
                                   bamfile, 
                                   controlfile = control )
    
    ############################################################
    ############################################################
    ############################################################
    @transform( loadZinba, suffix("_zinba.load"), ".bed.gz" )
    def exportIntervalsAsBed( infile, outfile ):
        PipelineChipseq.exportIntervalsAsBed( infile, outfile )

else:
    raise ValueError("unknown peak caller %s" % PARAMS["calling_caller"] )

############################################################
############################################################
############################################################
@follows( exportIntervalsAsBed )
@files( [( [ "%s.bed.gz" % y.asFile() for y in EXPERIMENTS[x]], 
           "%s.bed.gz" % x.asFile()) 
         for x in EXPERIMENTS ] )
def combineExperiment( infiles, outfile ):
    '''combine replicates between experiments.

    The replicates are combined using intersection.
    '''
    PipelineChipseq.intersectBedFiles( infiles, outfile )

############################################################
############################################################
############################################################
@follows( exportIntervalsAsBed )    
@files( [( [ "%s.bed.gz" % y.asFile() for y in CONDITIONS[x]], 
           "%s.bed.gz" % x.asFile()) 
         for x in CONDITIONS ] )
def combineCondition( infiles, outfile ):
    '''combine conditions between cell lines. 

    The conditions are merged via intersection.
    '''
    PipelineChipseq.intersectBedFiles( infiles, outfile )

############################################################
############################################################
############################################################
@follows( exportIntervalsAsBed )    
@files( [( [ "%s.bed.gz" % y.asFile() for y in TISSUES[x]], 
           "%s.bed.gz" % x.asFile()) 
         for x in TISSUES ] )
def combineTissue( infiles, outfile ):
    '''combine conditions between cell lines. 

    The conditions are merged via intersection.
    '''
    PipelineChipseq.intersectBedFiles( infiles, outfile )

############################################################
############################################################
############################################################
@follows( combineExperiment, combineCondition )
@files( [ ( ("%s.bed.gz" % x.asFile(), 
             "%s.bed.gz" % getUnstimulated(x).asFile()),
            "%s.bed.gz" % getSubtracted(x).asFile() ) for x in TOSUBTRACT ] )
def subtractUnstimulated( infiles, outfile ):
    '''remove the unstimulated data sets from the individual tracks. 
    '''

    infile, subtract = infiles
    PipelineChipseq.subtractBedFiles( infile, subtract, outfile )

############################################################
############################################################
############################################################
@follows(normalizeBAM)
@transform( ( combineExperiment,
              combineCondition, 
              subtractUnstimulated ) +
            P.asTuple(PARAMS["tracks_extra"]),
            suffix(".bed.gz"),
            "_bed.load" )
def loadIntervalsFromBed( infile, outfile ):
    '''load intervals from :term:`bed` formatted files into database.
    
    These :term:`bed` formatted intervals are derived by 
    merging/intersecting the various tracks or have been
    placed explicitely into the directory.

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

    tmpfile = P.getTempFile()

    headers = ("AvgVal","DisttoStart","GeneList","Length","PeakCenter","PeakVal","Position","interval_id","nCpGs","nGenes","nPeaks","nProbes","nPromoters", "contig","start","end" )

    tmpfile.write( "\t".join(headers) + "\n" )

    avgval,contig,disttostart,end,genelist,length,peakcenter,peakval,position,start,interval_id,ncpgs,ngenes,npeaks,nprobes,npromoters = \
        0,"",0,0,"",0,0,0,0,0,0,0,0,0,0,0,

    samfiles, offsets = [], []

    track = Sample( filename = P.snip( infile, ".bed.gz") )
    associated_track = track

    E.info( "loading data for track %s" % track )
    
    # get associated bam files
    if "%s_track" % track in PARAMS:
        associated_track = PARAMS[ "%s_track" % track ]
        E.info("using %s as associated track for %s" % (associated_track, track))
        
    if "%s_replicates" % track in PARAMS:
        replicates = P.asList( PARAMS["%s_replicates" % replicates] )
    else:
        # get replicates / aggregated tracks associated with track
        # remove subtraction as not relevant for tag counting
        unsubtracted_track = getUnsubtracted ( associated_track )
        
        try:
            replicates = PipelineTracks.getSamplesInTrack( unsubtracted_track, TRACKS )
        except KeyError:
            replicates = []

    if len(replicates) == 0:
        E.warn( "no replicates associated with track %s" % track )
        
    # setup files
    for t in replicates:
        fn = "%s.call.bam" % (t.asFile())
        assert os.path.exists( fn ), "could not find bamfile %s for track %s" % ( fn, str(t))
        samfiles.append( pysam.Samfile( fn,  "rb" ) )
        fn = "%s.macs" % t.asFile()
        if os.path.exists( fn ):
            if PARAMS["calling_caller"] == "macs":
                offsets.append( PipelineChipseq.getPeakShiftFromMacs( fn ) )
            elif PARAMS["calling_caller"] == "zinba":
                offsets.append( PipelineChipseq.getPeakShiftFromZinba( fn ) )
            else: 
                raise ValueError("unknown peak caller %s" % PARAMS["calling_caller"] )

    mlength = int(PARAMS["calling_merge_min_interval_length"])

    c = E.Counter()

    # count tags
    for bed in Bed.iterator( IOTools.openFile(infile, "r") ): 

        c.input += 1

        if "name" not in bed:
            bed.name = c.input
        
        # remove very short intervals
        if bed.end - bed.start < mlength: 
            c.skipped_length += 1
            continue

        if replicates:
            npeaks, peakcenter, length, avgval, peakval, nprobes = \
                PipelineChipseq.countPeaks( bed.contig, bed.start, bed.end, samfiles, offsets )

            # nreads can be 0 if the intervals overlap only slightly
            # and due to the binning, no reads are actually in the overlap region.
            # However, most of these intervals should be small and have already be deleted via 
            # the merge_min_interval_length cutoff.
            # do not output intervals without reads.
            if nprobes == 0:
                c.skipped_reads += 1

        else:
            npeaks, peakcenter, length, avgval, peakval, nprobes = ( 1, 
                                                                     bed.start + (bed.end - bed.start) // 2, 
                                                                     bed.end - bed.start, 
                                                                     1, 
                                                                     1,
                                                                     1 )
            
        c.output += 1
        tmpfile.write( "\t".join( map( str, (avgval,disttostart,genelist,length,
                                             peakcenter,peakval,position, bed.name,
                                             ncpgs,ngenes,npeaks,nprobes,npromoters, 
                                             bed.contig,bed.start,bed.end) )) + "\n" )

    if c.output == 0:
        E.warn( "%s - no aggregate intervals" )
        
 
    tmpfile.close()

    tmpfilename = tmpfile.name
    tablename = "%s_intervals" % track.asTable()
    
    statement = '''
    python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=interval_id 
              --table=%(tablename)s
    < %(tmpfilename)s 
    > %(outfile)s
    '''

    P.run()
    os.unlink( tmpfile.name )

    L.info( "%s\n" % str(c) )

############################################################
############################################################
############################################################
@merge( exportIntervalsAsBed, "merged.bed.gz" )
def buildMergedIntervals( infiles, outfile ):
    '''combine all experiments.

    The replicates are combined using a merge.
    '''
    PipelineChipseq.mergeBedFiles( infiles, outfile )

############################################################
############################################################
############################################################
## master target for this section
############################################################
@follows( loadIntervalsFromBed, 
          buildMergedIntervals )
def buildIntervals():
    pass

def getBamPeakPairParams():
    '''retrieves filenames for all vs. all comparsion (bed vs bam)'''

    bamfiles = glob.glob("*.bam")
    peakfiles = glob.glob("*.bed.gz")
    
    for bam in bamfiles:
        for peak in peakfiles:
            result = "%s_vs_%s.coverage.gz" % (P.snip(bam,".bam"), 
                                               P.snip(peak, ".bed.gz"))
            yield (bam,peak), result

###################################################################
###################################################################
###################################################################
@transform( normalizeBAM, suffix(".call.bam"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"] )),
            ".readprofile.transcripts.tsv.gz" )
def buildReadProfileOfTranscripts( infiles, outfile ):
    '''build a table with peak shape parameters.'''
    
    to_cluster = True

    bamfile, gtffile = infiles
    track = P.snip( bamfile, '.call.bam' )
    
    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                      --output-filename-pattern="%(outfile)s.%%s"
                      --f
                      --reporter=transcript
                      --method=geneprofile 
                      --method=tssprofile 
                      --normalization=total-sum
                      %(bamfile)s %(gtffile)s
                   > %(outfile)s
                '''
    P.run()

###################################################################
###################################################################
###################################################################
@transform( exportIntervalsAsBed, suffix(".bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"] )),
            ".intervalprofile.transcripts.tsv.gz" )
def buildIntervalProfileOfTranscripts( infiles, outfile ):
    '''build a table with peak shape parameters.'''
    
    to_cluster = True

    bedfile, gtffile = infiles
    track = P.snip( bedfile, '.bed.gz' )
    
    statement = '''python %(scriptsdir)s/bam2geneprofile.py
                      --output-filename-pattern="%(outfile)s.%%s"
                      --force
                      --reporter=transcript
                      --method=geneprofile 
                      --method=tssprofile 
                      %(bedfile)s %(gtffile)s
                   > %(outfile)s
                '''
    P.run()


###################################################################
###################################################################
###################################################################
@transform( exportIntervalsAsBed,
            suffix(".bed.gz"),
            ".peakshape.tsv.gz" )
def buildPeakShapeTable( infile, outfile ):
    '''build a table with peak shape parameters.'''
    
    to_cluster = True

    track = TRACKS.factory( filename = infile[:-len(".bed.gz")] )

    fg_replicates = PipelineTracks.Aggregate( TRACKS, track = track )[track]
    if len(fg_replicates) > 1:
        raise NotImplementedError( "peakshape of aggregate tracks not implemented" )

    bamfile = "%s.call.bam" % (fg_replicates[0])

    # get peak shift
    shift = PipelineChipseq.getPeakShift( track )
    if shift == None:
        E.warn ("could not get peak shift for track %s" % track )
        return

    E.info( "applying shift %i for track %s" % (shift, track ) )

    statement = '''python %(scriptsdir)s/bam2peakshape.py
                      --window-size=%(calling_peakshape_window_size)i
                      --bin-size=%(calling_peakshape_bin_size)i
                      --output-filename-pattern="%(outfile)s.%%s"
                      --force
                      --shift=%(shift)i
                      --sort=peak-height
                      --sort=peak-width
                      %(bamfile)s %(infile)s
                   > %(outfile)s
                '''
    P.run()

###################################################################
###################################################################
###################################################################
@follows( buildIntervals, mkdir("coverage.dir") )
@files( [ ( ("%s.prep.bam" % x, "%s.bed.gz" % y), 
            "coverage.dir/%s_vs_%s.coverage.gz" % (x,y)) for (x,y) in itertools.product( TRACKS, repeat = 2 ) ] )
def peakCoverage(infiles, outfile):
    '''uses coverageBed to count the total number of reads under peaks'''

    to_cluster = True
    bamfile, bedfile = infiles
    if P.isEmpty( bedfile ):
        P.touch( outfile )
    else:
        statement = '''coverageBed -abam %(bamfile)s -b %(bedfile)s | gzip > %(outfile)s ''' 
        P.run()                  

###################################################################
###################################################################
###################################################################
@merge( peakCoverage, "reads_under_peaks.tsv")
def buildReadCoverageTable(infiles, outfile):
    '''Counts reads from .coverage files and writes out a matrix of the results from all .bed vs. all .bam files'''
    
    out = IOTools.openFile(outfile, "w")

    bams = []
    beds = []
    data = {}

    for coverageFile in infiles:
        name = P.snip(os.path.basename( coverageFile ), ".coverage.gz")
        name = name.split("_")
        bam = name[0]
        bed = name[2] 
        value = 0
        for line in (IOTools.openFile(coverageFile).readlines()):
            value += int(line.split("\t")[5])
        
        if bam not in bams:
            bams.append(bam)
        if bed not in beds:
            beds.append(bed)
        
        data[(bam, bed)] = value

    table = numpy.empty((len(bams), len(beds)))
    for i in range(len(bams)):
        for j in range(len(beds)):
            table[i][j] = data[bams[i], beds[j]]    

    table = zip(bams, table)
        
    out.write ("track" + "\t" + "\t".join(beds) + "\n")
    for i in range(len(bams)):
        out.write(table[i][0] + "\t" + "\t".join((map(str, table[i][1]))) + "\n")


   # table = numpy.empty((len(beds), len(bams)), dtype=numpy.int )
    #for i in range(len(beds)):
     #   for j in range(len(bams)):
      #      table[i][j] = data[bams[j], beds[i]]    

   # table = zip(beds, table)
        
   # out.write ("track" + "\t" + "\t".join(bams) + "\n")
   # for i in range(len(beds)):
    #    out.write(table[i][0] + "\t" + "\t".join((map(str, table[i][1]))) + "\n")
    
###################################################################
###################################################################
###################################################################
@transform( buildReadCoverageTable, suffix(".tsv"), ".load")
def loadReadCoverageTable( infile, outfile ):
    '''load read coverage table.'''
    P.load( infile, outfile )

###################################################################
###################################################################
###################################################################
## general import
###################################################################

############################################################
############################################################
############################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "wig") ) )
@follows( exportIntervalsAsBed )
@transform( (buildBAM, normalizeBAM),
            regex("(.*).bam"),
            r"%s/\1.bigwig" % os.path.join( PARAMS["exportdir"], "wig") )
def exportBigwig( infile, outfile ):
    '''convert BAM to bigwig file.'''

    to_cluster = True

    track = re.sub( ".(call|genome).bam", "", infile)

    shift = PipelineChipseq.getPeakShift( track ) 
    if shift == None: 
        E.warn( "no shift found for %s - using 0" % track )
        shift = 0
        extend = 0
    else: 
        shift /= 2
        extend = shift

    statement = '''python %(scriptsdir)s/bam2wiggle.py 
                      --output-format=bigwig 
                      --output-filename=%(outfile)s 
                      --shift=%(shift)i
                      --extend=%(extend)i
                      --log=%(outfile)s.log
                %(infile)s '''
     
    P.run()

############################################################
############################################################
############################################################
## 
############################################################
@merge( exportBigwig, "bigwiginfo.tsv" )
def buildBigwigInfo( infiles, outfile ):
    '''collate summary of bigwig files.'''

    pattern = '''<( bigWigInfo %s | perl -p -e "s/: /\\t/" )'''
    
    p = " ".join( [pattern % infile for infile in infiles ] )
    headers = ",".join( [ P.snip( os.path.basename( x ), ".bigwig") for x in infiles ] )

    statement = '''python %(scriptsdir)s/combine_tables.py
                   --headers=%(headers)s
                         %(p)s
                > %(outfile)s
                '''
    P.run()

############################################################
############################################################
############################################################
# fix: does not work due to exportIntervalsAsBed
# (loadIntervalsFromBed, exportIntervalsAsBed ),
@follows( buildIntervals )
@transform( [ "%s.bed.gz" % x.asFile() for x in TRACKS_MASTER ] +\
                P.asList( PARAMS["tracks_extra"] ),
            suffix(".bed.gz"),
            ".fasta" )
def exportMotifSequences( infile, outfile ):
    '''export sequences for all intervals.
    '''
    to_cluster = True
    statement = '''zcat %(infile)s 
    | python %(scriptsdir)s/bed2fasta.py 
              --genome-file=%(genome_dir)s/%(genome)s 
              --masker=%(motifs_masker)s
              --mode=intervals
              --log=%(outfile)s.log
    > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
# (loadIntervalsFromBed, exportIntervalsAsBed ),
# TODO: fix, causes a problem due to exportIntervalsAsBed
@follows( buildIntervals )
@transform( [ "%s.bed.gz" % x.asFile() for x in TRACKS_MASTER] +\
                P.asList( PARAMS["tracks_extra"] ),
            suffix(".bed.gz"),
            ".controlfasta" )
def exportMotifControlSequences( infile, outfile ):
    '''for each interval, export the left and right 
    sequence segment of the same size.
    '''

    to_cluster = True
    statement = '''zcat %(infile)s 
    | python %(scriptsdir)s/bed2fasta.py 
              --genome-file=%(genome_dir)s/%(genome)s 
              --masker=%(motifs_masker)s
              --mode=leftright
              --log=%(outfile)s.log
    > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
@merge( (exportIntervalsAsBed, 
         combineExperiment, 
         combineCondition,
         ),
        "intervals.overlap" )
def buildOverlap( infiles, outfile ):
    '''compute overlap between intervals.
    '''

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
############################################################
############################################################
@transform( buildOverlap, suffix(".overlap"), "_overlap.load" )
def loadOverlap( infile, outfile ):
    '''load overlap results.
    '''

    tablename = "overlap"

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=set1 
              --index=set2 
              --table=%(tablename)s 
    < %(infile)s > %(outfile)s
    '''

    P.run()

if 0:
    ############################################################
    ############################################################
    ############################################################
    @merge( (exportUCSCEncodeTracks, 
             combineConditions,
             combineUnstim,
             combineExperiment ),
            "ucsc.overlap" )
    def buildUCSCOverlap( infiles, outfile ):
        '''compute overlap between intervals and the ucsc tracks.
        '''

        if os.path.exists(outfile): 
            os.rename( outfile, outfile + ".orig" )
            options = "--update=%s.orig" % outfile
        else:
            options = ""

        to_cluster = True

        infiles = " ".join(infiles)
        statement = '''
            python %(scriptsdir)s/diff_bed.py --tracks %(options)s %(infiles)s > %(outfile)s
            '''

        P.run()

    ############################################################
    ############################################################
    ############################################################
    @transform( buildUCSCOverlap, suffix(".overlap"), "_overlap.import" )
    def importUCSCOverlap( infile, outfile ):
        '''import overlap results.
        '''

        tablename = "ucsc_overlap"

        statement = '''
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
                  --index=set1 \
                  --index=set2 \
                  --table=%(tablename)s \
        < %(infile)s > %(outfile)s
        '''

        P.run()

############################################################
############################################################
############################################################
## targets to do with the analysis of replicates
############################################################
@follows( buildIntervals )
@files( [ ([ "%s.bed.gz" % y.asFile() for y in EXPERIMENTS[x]], 
           "%s.reproducibility" % x.asFile()) for x in EXPERIMENTS ] )
def makeReproducibility( infiles, outfile ):
    '''compute overlap between intervals.

    Note the exon percentages are approximations assuming that there are
    not more than one intervals in one set overlapping one in the other set.
    '''

    dbhandle = connect()

    if len(infiles) < 2:
        P.touch(outfile)
        return

    data = []
    for replicate in infiles:
        cc = dbhandle.cursor()
        tablename = "%s_intervals" % P.quote( P.snip(replicate, ".bed.gz") )
        statement = "SELECT contig, start, end, peakval FROM %(tablename)s" % locals()
        cc.execute( statement )
        data.append( cc.fetchall() )

    if len(data) == 0:
        E.warn( "no data for %s" % outfile )
        P.touch( outfile )
        return
    
    try:
        ma = int(max( [ x[3] for x in itertools.chain( *data ) ] ) + 1)
    except ValueError:
        E.warn( "no data for %s" % outfile )
        P.touch( outfile )
        return
    
    nexons1, nexons2, nexons_overlapping = [0] * ma, [0] * ma, [0] * ma
    nbases1, nbases2, nbases_overlapping = [0] * ma, [0] * ma, [0] * ma

    idx = IndexedGenome.IndexedGenome()
    for contig, start, end, peakval in data[0]:
        peakval = int(peakval)
        for x in range( 0, peakval):
            nexons1[x] += 1
            nbases1[x] += end - start
        idx.add( contig, start, end, peakval )

    for contig, start, end, peakval in data[1]:
        peakval = int(peakval)
        for x in range( 0, peakval):
            nexons2[x] += 1
            nbases2[x] += end - start

        try:
            intervals = list(idx.get( contig, start, end ))
        except KeyError:
            continue

        if len(intervals) == 0: continue

        for other_start,other_end,other_peakval in intervals:
            ovl = min(end, other_end ) - max(start, other_start )
            for x in range( 0, min( other_peakval, peakval) ):
                nexons_overlapping[x] += 1
                nbases_overlapping[x] += ovl

    outs = IOTools.openFile( outfile, "w" )
    outs.write("peakval\tnexons1\tnexons2\tnexons_union\tnexons_ovl\tpexons_ovl\tpexons_union\tnbases1\tnbases2\tnbases_union\tnbases_ovl\tpbases_ovl\tpbases_union\n" )
    total_bases = nbases1[0] + nbases2[0] - nbases_overlapping[0]
    for x in range(0, ma):
        bases_union = nbases1[x] + nbases2[x] - nbases_overlapping[x]
        exons_union = nexons1[x] + nexons2[x] - nexons_overlapping[x]
        outs.write( "\t".join( (str(x),
                                str(nexons1[x]),
                                str(nexons2[x]),
                                str(exons_union),
                                str(nexons_overlapping[x]),
                                IOTools.prettyPercent( nbases_overlapping[x], bases_union),
                                IOTools.prettyPercent( bases_union, total_bases),
                                str(nbases1[x]),
                                str(nbases2[x]),
                                str(bases_union),
                                str(nbases_overlapping[x]),
                                IOTools.prettyPercent( nbases_overlapping[x], bases_union),
                                IOTools.prettyPercent( bases_union, total_bases),
                                )) + "\n" )

    outs.close()

@transform( makeReproducibility, suffix(".reproducibility"), "_reproducibility.load" )
def loadReproducibility( infile, outfile ):
    '''load Reproducibility results
    '''
    P.load( infile, outfile, options="--allow-empty" )

@follows( loadReproducibility )
def reproducibility(): pass
    
############################################################
############################################################
############################################################
##
############################################################
@follows( buildIntervals )
@files_re( ["%s.bed.gz" % x.asFile() for x in TRACKS_CORRELATION],
           combine("(.*).bed.gz"),
           "peakval.correlation" )
def makePeakvalCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PipelineChipseq.makeIntervalCorrelation( infiles, outfile, "peakval", "merged.bed.gz" )

@follows( buildIntervals )
@files_re( ["%s.bed.gz" % x.asFile() for x in TRACKS_CORRELATION],
           combine("(.*).bed.gz"),
           "avgval.correlation" )
def makeAvgvalCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PipelineChipseq.makeIntervalCorrelation( infiles, outfile, "avgval", "merged.bed.gz" )

@follows( buildIntervals )
@files_re( ["%s.bed.gz" % x.asFile() for x in TRACKS_CORRELATION],
           combine("(.*).bed.gz"),
           "length.correlation" )
def makeLengthCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PipelineChipseq.makeIntervalCorrelation( infiles, outfile, "length", "merged.bed.gz" )

@transform( ( makePeakvalCorrelation, makeAvgvalCorrelation, makeLengthCorrelation ),
            suffix(".correlation"),
            "_correlation.load")
def loadCorrelation( infile, outfile ):
    '''load correlation data.'''
    P.load( infile, outfile, "--index=id --map=default:float --map=id:int" )

############################################################
############################################################
############################################################
@transform( "*.motif.fasta",
            suffix(".motif.fasta"),
            ".motif" )
def buildReferenceMotifs( infile, outfile ):
    '''build motif from jasper sites by running MEME.
    '''

    tmpdir = tempfile.mkdtemp()
    tmpname = os.path.join( tmpdir, "tmpfile") 
    tmpfile = IOTools.openFile( tmpname, "w")
    for line in IOTools.openFile(infile,"r"):
        tmpfile.write( re.sub( "[ \t]+", "_", line) )
    tmpfile.close()

    statement = '''
    meme %(tmpname)s -mod oops -dna -revcomp -nmotifs 1 -text -oc %(tmpdir)s > %(outfile)s
    '''
    P.run()

    shutil.rmtree( tmpdir )

############################################################
############################################################
############################################################
@transform( exportMotifSequences, suffix(".fasta"), ".meme")
def runMEME( infile, outfile ):
    '''run MEME to find motifs.

    In order to increase the signal/noise ratio,
    MEME is not run on all intervals but only the 
    top 10% of intervals (peakval) are used. 
    Also, only the segment of 200 bp around the peak
    is used and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker
    '''
    dbhandle = connect()

    track = infile[:-len(".fasta")]

    PipelineMotifs.runMEME( track, outfile, dbhandle )

############################################################
############################################################
############################################################
##
############################################################
@transform( exportMotifSequences, suffix(".fasta"), ".glam2")
def runGLAM2( infile, outfile ):
    '''run glam2 on all intervals and motifs.

    In order to increase the signal/noise ratio,
    MEME is not run on all intervals but only the 
    top 10% of intervals (peakval) are used. 
    Also, only the segment of 200 bp around the peak
    is used and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker
    '''
    PipelineMotifs.runGLAM2( infile, outfile )

############################################################
############################################################
############################################################
## selecting which motifs to run for motif search
############################################################
if PARAMS["tomtom_master_motif"] != "":
    ############################################################
    ############################################################
    ## filter motifs by a master motif
    ############################################################
    ############################################################
    ############################################################
    @follows( buildReferenceMotifs )
    @transform( runMEME, suffix(".meme"), ".tomtom" )
    def runTomTom( infile, outfile ):
        '''compare ab-initio motifs against tomtom.'''
        to_cluster = True
        statement = '''
           tomtom -text %(tomtom_master_motif)s %(infile)s > %(outfile)s
           '''
        P.run()

    ############################################################
    ############################################################
    ############################################################
    @transform( runTomTom, suffix(".tomtom"), "_tomtom.load" )
    def loadTomTom( infile, outfile ):
        '''compare ab-initio motifs against tomtom.'''

        tmpfile = P.getTempFile()

        tmpfile.write( "\t".join( \
            ("query_id", "target_id",
             "optimal_offset", "pvalue", "evalue", "qvalue", 
             "overlap", "query_consensus",
             "target_consensus", "orientation") ) + "\n" )

        for line in IOTools.openFile( infile, "r" ):
            if line.startswith( "#Query" ): continue
            data = line[:-1].split("\t")
            (query_id, 
             target_id, 
             optimal_offset, 
             pvalue, 
             evalue,
             qvalue, 
             overlap, 
             query_consensus,
             target_consensus, 
             orientation) = data
            tmpfile.write( line )

        tmpfile.close()

        tmpname = tmpfile.name
        tablename = P.toTable( outfile )

        statement = '''
        python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                  --allow-empty 
                  --table=%(tablename)s 
        < %(tmpname)s 
        > %(outfile)s
        '''

        P.run()
        os.unlink( tmpname )

    @transform( runTomTom, suffix(".tomtom"), ".motif" )
    def filterMotifs( infile, outfile ):
        '''scan motifs found with MEME against master motif using tomtom.
        '''

        infile_meme = infile[:-len(".tomtom")] + ".meme"

        selected = []
        max_pvalue = float(PARAMS["tomtom_filter_pvalue"])
        for line in IOTools.openFile( infile, "r" ):
            if line.startswith( "#Query" ): continue
            (query_id, target_id, 
             optimal_offset, pvalue, evalue, qvalue, overlap, query_consensus,
             target_consensus, orientation) = line[:-1].split("\t")
            if float(pvalue) <= max_pvalue:
                selected.append( target_id )

        L.info( "%s: keeping %i motifs" % (infile, len(selected) ) )

        PipelineMotifs.filterMotifsFromMEME( infile_meme, outfile, selected )
        
else:
    ############################################################
    ############################################################
    ## no master motif
    ############################################################
    ############################################################
    ############################################################
    @follows(runMEME)
    def runTomTom(): pass

    @follows(runTomTom)
    def loadTomTom(): pass
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runMEME, suffix(".meme"), ".motif" )
    def filterMotifs( infile, outfile ):
        '''take the top scoring motif from meme runs
        '''

        PipelineMotifs.filterMotifsFromMEME( infile, outfile, ["1"] )

@follows( filterMotifs, buildReferenceMotifs )
@transform( "*.motif", suffix(".motif"), ".motif")
def collectMotifs( infile, outfile ):
    '''dummy target - merge motifs.'''

############################################################
############################################################
############################################################
@merge( (collectMotifs), "motif_info.load" )
def loadMotifInformation( infiles, outfile ):
    '''load information about motifs into database.'''
    
    outf = P.getTempFile()

    outf.write("motif\n" )

    for infile in infiles:
        if IOTools.isEmpty( infile ): continue
        motif = P.snip( infile, ".motif" )
        outf.write( "%s\n" % motif )

    outf.close()

    P.load( outf.name, outfile )
    
    os.unlink( outf.name )

############################################################
############################################################
############################################################
@transform( exportMotifSequences,
            suffix(".fasta"),
            ".motifseq_stats.load" )
def loadMotifSequenceComposition( infile, outfile ):
    '''compute sequence composition of sequences used for ab-initio search.'''

    to_cluster = True

    tablename = P.toTable( outfile )

    statement = '''
    python %(scriptsdir)s/fasta2table.py 
        --section=na
        --log=%(outfile)s.log
    < %(infile)s
    | python %(scriptsdir)s/csv2db.py
        %(csv2db_options)s
        --table=%(tablename)s
    > %(outfile)s'''
    
    P.run()

############################################################
############################################################
############################################################
# todo: fix, causes a problem: exportSequences and exportMASTControlSequences,
@follows( collectMotifs )
@files_re( (exportMotifSequences, exportMotifControlSequences),
           "(\S+).controlfasta",
           [ r"\1.controlfasta", r"\1.fasta",  glob.glob("*.motif")],
           r"\1.mast.gz" )
def runMAST( infiles, outfile ):
    '''run mast on all intervals and motifs.

    Collect all results for an E-value up to 10000 so that
    all sequences are output and MAST curves can be computed. 

    10000 is a heuristic.
    '''
    PipelineMotifs.runMAST( infiles, outfile )

############################################################
############################################################
############################################################
@transform( exportMotifSequences,
            suffix(".fasta"),
            ".bioprospector")
def runBioProspector( infiles, outfile ):
    '''run bioprospector for motif discovery.

    Bioprospector is run on only the top 10% of peaks.
    '''
    dbhandle = connect()

    PipelineMotifs.runBioProspector( infiles, outfile, dbhandle )

############################################################
############################################################
############################################################
##
############################################################
@transform( runBioProspector, suffix(".bioprospector"), "_bioprospector.load")
def loadBioProspector( infile, outfile ):
    '''load results from bioprospector.'''
    PipelineMotifs.loadBioProspector( infile, outfile )

############################################################
############################################################
############################################################
@transform( runMAST,
            suffix(".mast.gz"),
            "_mast.load" )
def loadMAST( infile, outfile ):
    '''parse mast file and load into database.

    Parse several motif runs and add them to the same
    table.

    Add columns for the control data as well.
    '''
    PipelineMotifs.loadMAST( infile, outfile )

############################################################
############################################################
############################################################
##
############################################################
@follows( collectMotifs )
@files_re( (exportMotifSequences, exportMotifControlSequences),
           "(\S+).controlfasta",
           [ r"\1.controlfasta", r"\1.fasta",  glob.glob("*.glam2")],
           r"\1.glam2scan" )
def runGLAM2SCAN( infiles, outfile ):
    '''run glam2scan on all intervals and motifs.
    '''
    dbhandle = connect()

    PipelineMotifs.runGLAM2SCAN( infiles, outfile, dbhandle )

############################################################
############################################################
############################################################
@transform( runGLAM2SCAN,
            suffix(".glam2scan"),
            "_glam.load" )
def loadGLAM2SCAN( infile, outfile ):
    '''parse mast file and load into database.

    Parse several motif runs and add them to the same
    table.
    '''
    PipelineMotifs.loadGLAM2SCAN( infile, outfile )

############################################################
############################################################
############################################################
@follows( buildIntervals )
@files( [ ("%s.bed.gz" % x.asFile(), 
           "%s.annotations" % x.asFile() ) for x in TRACKS_MASTER ] )
def annotateIntervals( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set.
    '''
    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_geneset_all_gtf"] )

    statement = """
    zcat < %(infile)s 
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
############################################################
############################################################
@follows( buildIntervals )
@files( [ ("%s.bed.gz" % x.asFile(), 
           "%s.tss" % x.asFile() ) for x in TRACKS_MASTER ] )
def annotateTSS( infile, outfile ):
    '''compute distance to TSS'''

    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_tss_bed"] )

    statement = """
    zcat < %(infile)s 
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
############################################################
############################################################
@follows( buildIntervals )
@files( [ ("%s.bed.gz" % x.asFile(), 
           "%s.repeats" % x.asFile() ) for x in TRACKS_MASTER ] )
def annotateRepeats( infile, outfile ):
    '''count the overlap between intervals and repeats.'''

    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_repeats_gff"] )

    statement = """
    zcat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python %(scriptsdir)s/gtf2table.py \
		--counter=overlap \
		--log=%(outfile)s.log \
		--filename-gff=%(annotation_file)s \
		--genome-file=%(genome_dir)s/%(genome)s
	> %(outfile)s"""

    P.run()

############################################################
@transform( annotateIntervals, suffix(".annotations"), "_annotations.load" )
def loadAnnotations( infile, outfile ):
    '''load interval annotations: genome architecture
    '''
    P.load( infile, outfile, "--index=gene_id --allow-empty" )

############################################################
@transform( annotateTSS, suffix( ".tss"), "_tss.load" )
def loadTSS( infile, outfile ):
    '''load interval annotations: distance to transcription start sites
    '''
    P.load( infile, outfile, "--index=gene_id --index=closest_id --index=id5 --index=id3 --allow-empty" )

############################################################
@transform( annotateRepeats, suffix(".repeats"), "_repeats.load" )
def loadRepeats( infile, outfile ):
    '''load interval annotations: repeats
    '''
    P.load( infile, outfile, "--index=gene_id --allow-empty" )

############################################################
############################################################
############################################################
## count coverage within intervals for each track against
## the unstimulated tracks
############################################################
@follows( subtractUnstimulated )
@files( [ ("%s.bed.gz" % x.asFile(), "%s.readcounts" % x.asFile() ) for x in TOSUBTRACT ] )
def buildIntervalCounts( infile, outfile ):
    '''count read density in bed files comparing stimulated versus unstimulated binding.
    '''
    track = TRACKS.factory( filename = outfile[:-len(".readcounts")] )
    unstim = getUnstimulated( track )

    fg_replicates = PipelineTracks.Aggregate( TRACKS, track = track )[track]
    bg_replicates = PipelineTracks.Aggregate( TRACKS, track = unstim )[unstim]

    PipelineChipseq.buildIntervalCounts( infile, outfile, track, fg_replicates, bg_replicates )

############################################################
############################################################
############################################################
@transform( buildIntervalCounts, 
            suffix(".readcounts"), 
            "_readcounts.load" )
def loadIntervalCounts( infile, outfile ):
    '''load interval counts.'''
    P.load( infile, outfile, "--index=gene_id --allow-empty" )

############################################################
############################################################
############################################################
## export section start
############################################################

############################################################
############################################################
############################################################
## export bigwig tracks
############################################################
@files_re(exportBigwig, combine("(.*).bigwig"), "bigwig.view")
def viewBigwig( infiles, outfile ):

    outs = IOTools.openFile( outfile, "w" )
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
############################################################
############################################################
## export intervals
############################################################
@follows ( mkdir( 'export' ))
@merge( "*.bed.gz", ("export/intervals_%s.bed.gz" % PARAMS["version"], "intervals.view") )
def viewIntervals( infiles, outfiles ):

    outfile_bed, outfile_code = outfiles
    outs = IOTools.openFile( outfile_bed, "w" )
    version = PARAMS["version"]
    for infile in infiles:

        track = infile[:-len(".bed")]
        
        outs.write( '''track name="interval_%(track)s_%(version)s" description="Intervals in %(track)s - version %(version)s" visibility=2\n''' % locals() )

        with IOTools.openFile(infile,"r") as f:
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

    outs = IOTools.openFile( outfile_code, "w" )
    outs.write( "#paste the following into the UCSC browser:\n" )
    outs.write( "http://wwwfgu.anat.ox.ac.uk/~andreas/ucsc_tracks/%(filename)s\n" % locals())
    outs.close()

############################################################
############################################################
############################################################
## export section end
############################################################

############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
##
############################################################
@follows( makeReadCorrelationTable,
          loadBAMStats )
def mapping():
    '''map reads'''
    pass

@follows( mapping,
          buildIntervals, 
          loadReadCoverageTable,
          buildPeakShapeTable,
          buildReadProfileOfTranscripts,
          buildIntervalProfileOfTranscripts,
          buildBigwigInfo)
def intervals():
    '''compute binding intervals.'''
    pass

@follows( exportBigwig, viewIntervals, viewBigwig )
def export():
    '''export files.'''
    pass

@follows( buildReferenceMotifs,
          exportMotifSequences,
          runMEME,
          runTomTom, loadTomTom )
#          runBioProspector )
#          runGLAM2, 
def discover_motifs():
    '''run motif discovery.'''
    pass

@follows( filterMotifs,
          exportMotifControlSequences,
          loadMotifSequenceComposition,
          loadMotifInformation,
          runMAST, loadMAST )
#          runGLAM2SCAN, loadGLAM2SCAN )
def detect_motifs():
    '''run motif detection.'''
    pass

@follows( loadCorrelation, 
          loadOverlap,
          reproducibility)
def correlation():
    '''run the correlation targets.'''
    pass

@follows( annotateIntervals, loadAnnotations, 
          annotateTSS, loadTSS, 
          annotateRepeats, loadRepeats,
          #annotateTSSIntervalAssociations, loadTSSIntervalAssociations,
          #annotateTSSIntervalDistance, loadTSSIntervalDistance, 
          buildIntervalCounts, loadIntervalCounts )
def annotation():
    '''run the annotation targets.'''
    pass

###################################################################
###################################################################
###################################################################
## export targets
###################################################################
@merge( intervals,  "view_mapping.load" )
def createViewMapping( infile, outfile ):
    '''create view in database for alignment stats.

    This view aggregates all information on a per-track basis.

    The table is built from the following tracks:
    
    bam_stats: .call
    '''

    tablename = P.toTable( outfile )

    # can not create views across multiple database, so use table
    view_type = "TABLE"
    
    dbhandle = connect()
    Database.executewait( dbhandle, "DROP %(view_type)s IF EXISTS %(tablename)s" % locals() )

    statement = '''
    CREATE %(view_type)s %(tablename)s AS
    SELECT SUBSTR( b.track, 1, LENGTH(b.track) - LENGTH( '.genome')) AS track, *
    FROM bam_stats AS b
    WHERE b.track LIKE "%%.genome" 
    ''' % locals()

    Database.executewait( dbhandle, statement )

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
@follows( intervals, 
          discover_motifs, 
          detect_motifs,
          correlation, 
          annotation,
          views)
def full():
    '''run the full pipeline.'''
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
          mkdir("%s/intervals" % PARAMS["web_dir"]),
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
        "bamfiles" : glob.glob( "*.bam" ) + glob.glob( "*.bam.bai" ),
        # "genesets": [ "lincrna.gtf.gz", "abinitio.gtf.gz" ],
        "intervals" : glob.glob("*.bed"),
        # "classification": glob.glob("*.class.tsv.gz") ,
        #"differential_expression" : glob.glob( "*.cuffdiff.dir" ),
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

    # print( "# tracks found: %s" % TRACKS_ALL )
    # print( "# tracks by experiment: %s" % TRACKS_EXPERIMENTS )
    # print( "# tracks by condition: %s" % TRACKS_CONDITIONS )
    # print( "# tracks by tissue: %s" % TRACKS_TISSUES )

    # fails for config command
    # check compatibility
    # assert PARAMS["genome"] == PARAMS_ANNOTATIONS["genome"] 
    # sanity checks
    # assert len(TRACKS_ALL) > 0

    sys.exit( P.main(sys.argv) )


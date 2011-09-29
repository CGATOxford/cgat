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
import MAST, GTF, GFF, Bed
import pysam
import numpy
import gzip
import Masker
import Glam2Scan
import Motifs
import Bioprospector

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
            regions_to_filter += PipelineChipseq.getBedLocations(PARAMS["calling_filter_regions"])
        fh = open(outfile,"w")
        for r in regions_to_filter:
            fh.write("\t".join([r[0],str(r[1]),str(r[2])]))
            fh.write("\n")
        fh.close()
    else:
        P.touch(outfile)

############################################################
############o################################################
############################################################
@transform( buildBAM, suffix(".genome.bam"), add_inputs( makeMask), ".prep.bam")
def prepBAMForPeakCalling(infiles, outfile):
    '''Prepare BAM files for peak calling.

        - unmapped reads are removed.

        - if the option "calling_deduplicate" is true duplicate reads are simplistically removed 
            (using samtools rmdup may be superior)

        - reads may be filtered by exon or location 

           - to remove reads by exon, the option "calling_filter_exons" should specify a file containing 
             a list of ensembl gene identifiers (one per line)
           - to remove reads by location, the option "calling_filter_regions" should specify a bed file''
        
        The resulting bam file has a .prep.bam extension. Merging infiles is currently untested and the 
        methods only consider single end reads.
    '''
    bam_file, mask_file = infiles

    if PARAMS["calling_filter_exons"] or PARAMS["calling_filter_regions"]:
        mask = PipelineChipseq.buildQuicksectMask(mask_file)
    else:
        mask = None

    PipelineChipseq.buildBAMforPeakCalling(bam_file,outfile,PARAMS["calling_deduplicate"],mask)

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
    @transform(prepBAMForPeakCalling,suffix("prep.bam"),"prep.count")
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
            fh = open(cf,"r")
            counts.append(int(fh.read()))
            fh.close()
        fh = open(outfile,"w")
        fh.write(str(min(counts)))
        fh.close

    # Make normalised files
    @follows(minReads)
    @transform(prepBAMForPeakCalling,
               regex(r"(.*).prep.bam"),
               inputs( (r"\1.prep.bam",r"\1.prep.count") ), 
               r"\1.call.bam")
    def normalizeBAM( infiles, outfile ):
        '''build a normalized BAM file such that all
    files have approximately the same number of 
    reads.
  '''   
        fh = open("minreads")
        minreads = int(fh.read())
        fh.close
        PipelineChipseq.buildSimpleNormalizedBAM( infiles, 
                                             outfile,
                                             minreads)
else:
    @transform(prepBAMForPeakCalling,
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
    outf = open(outfile, "w" )
    outf.write("track\tcutoff\tpositions\tcoefficient\n" )

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
        mat = numpy.array( [ x for x in selector( gzip.open(infile, "r") )] )
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
             prepBAMForPeakCalling, 
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
    @files( [ (("%s.call.bam" % (x.asFile()), 
                "%s.call.bam" % (getControl(x).asFile())), 
               "%s.macs" % x.asFile() ) for x in TRACKS ] )
    def runMACS( infiles, outfile ):
        '''run MACS for peak detection.'''
        infile, controlfile = infiles

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
    @transform( loadMACS, suffix("_macs.load"), ".bed" )
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
        '''run MACS for peak detection.'''
        infile, controlfile = infiles

        # to_cluster = False
        
        job_options= "-pe dedicated %i -R y" % PARAMS["zinba_threads"]

        mappability_dir = os.path.join( PARAMS["zinba_mappability_dir"], 
                                 PARAMS["genome"],
                                 "%i" % PARAMS["zinba_read_length"],
                                 "%i" % PARAMS["zinba_alignability_threshold"],
                                 "%i" % PARAMS["zinba_fragment_size"])

        if not os.path.exists( mappability_dir ):
            raise OSError("mappability not found, expected to be at %s" % mappability_dir )

        bit_file = os.path.join( PARAMS["zinba_index_dir"], 
                                 PARAMS["genome"] ) + ".2bit"
        if not os.path.exists( bit_file):
            raise OSError("2bit file not found, expected to be at %s" % bit_file )

        options = []
        if controlfile:
            options.append( "--control-filename=%(controlfile)s" % locals() )
            
        options = " ".join(options)
            
        statement = '''
        python %(scriptsdir)s/WrapperZinba.py
               --input-format=bam
               --fdr-threshold=%(zinba_fdr_threshold)f
               --fragment-size=%(zinba_fragment_size)s
               --threads=%(zinba_threads)i
               --bit-file=%(bit_file)s
               --mappability-dir=%(mappability_dir)s
               %(options)s
        %(infile)s %(outfile)s
        >& %(outfile)s.log
        '''
        
        P.run()
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runZinba,
                regex(r"(.*).zinba"),
                inputs( (r"\1.zinba", r"\1.call.bam") ),
                r"\1_zinba.load" )
    def loadZinba( infiles, outfile ):
        infile, bamfile = infiles
        PipelineChipseq.loadZinba( infile, outfile, bamfile )
    
    ############################################################
    ############################################################
    ############################################################
    @transform( loadZinba, suffix("_zinba.load"), ".bed" )
    def exportIntervalsAsBed( infile, outfile ):
        PipelineChipseq.exportIntervalsAsBed( infile, outfile )

else:
    raise ValueError("unknown peak caller %s" % PARAMS["calling_caller"] )

############################################################
############################################################
############################################################
@follows( exportIntervalsAsBed )
@files( [( [ "%s.bed" % y.asFile() for y in EXPERIMENTS[x]], 
           "%s.bed" % x.asFile()) 
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
@files( [( [ "%s.bed" % y.asFile() for y in CONDITIONS[x]], 
           "%s.bed" % x.asFile()) 
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
@files( [( [ "%s.bed" % y.asFile() for y in TISSUES[x]], 
           "%s.bed" % x.asFile()) 
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
@files( [ ( ("%s.bed" % x.asFile(), 
             "%s.bed" % getUnstimulated(x).asFile()),
            "%s.bed" % getSubtracted(x).asFile() ) for x in TOSUBTRACT ] )
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
            suffix(".bed"),
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

    track = Sample( filename = P.snip( infile, ".bed") )
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
    for line in open(infile, "r"):
        c.input += 1
        contig, start, end, interval_id = line[:-1].split()[:4]

        start, end = int(start), int(end)

        # remove very short intervals
        if end-start < mlength: 
            c.skipped_length += 1
            continue

        if replicates:
            npeaks, peakcenter, length, avgval, peakval, nprobes = \
                PipelineChipseq.countPeaks( contig, start, end, samfiles, offsets )

            # nreads can be 0 if the intervals overlap only slightly
            # and due to the binning, no reads are actually in the overlap region.
            # However, most of these intervals should be small and have already be deleted via 
            # the merge_min_interval_length cutoff.
            # do not output intervals without reads.
            if nprobes == 0:
                c.skipped_reads += 1

        else:
            npeaks, peakcenter, length, avgval, peakval, nprobes = ( 1, 
                                                                     start + (end - start) // 2, 
                                                                     end - start, 
                                                                     1, 
                                                                     1,
                                                                     1 )
            
        c.output += 1
        tmpfile.write( "\t".join( map( str, (avgval,disttostart,genelist,length,
                                             peakcenter,peakval,position,interval_id,
                                             ncpgs,ngenes,npeaks,nprobes,npromoters, 
                                             contig,start,end) )) + "\n" )

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
@merge( exportIntervalsAsBed, "merged.bed" )
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
    peakfiles = glob.glob("*.bed")
    
    for bam in bamfiles:
        for peak in peakfiles:
            result = "%s_vs_%s.coverage.gz" % (P.snip(bam,".bam"), 
                                               P.snip(peak, ".bed"))
            yield (bam,peak), result

###################################################################
###################################################################
###################################################################
@follows( buildIntervals )
@files( [ ( ("%s.prep.bam" % x, "%s.bed" % y), 
            "%s_vs_%s.coverage.gz" % (x,y)) for (x,y) in itertools.product( TRACKS, repeat = 2 ) ] )
def peakCoverage(infile, outfile):
    '''uses coverageBed to count the total number of reads under peaks'''

    tracks = infile
    to_cluster = True
    statement = '''coverageBed -abam %s -b %s | gzip > %%(outfile)s '''  % (infile[0], infile[1]) 
    P.run()                  

###################################################################
###################################################################
###################################################################
@merge( peakCoverage, "reads_under_peaks.tsv")
def buildReadCoverageTable(infiles, outfile):
    '''Counts reads from .coverage files and writes out a matrix of the results from all .bed vs. all .bam files'''
    
    out = open(outfile, "w")

    bams = []
    beds = []
    data = {}

    for coverageFile in infiles:
        name = P.snip(coverageFile, ".coverage.gz")
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
@files_re( (buildBAM, normalizeBAM),
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
############################################################
############################################################
# fix: does not work due to exportIntervalsAsBed
# (loadIntervalsFromBed, exportIntervalsAsBed ),
@follows( buildIntervals )
@transform( [ "%s.bed" % x.asFile() for x in TRACKS_MASTER ] +\
                P.asList( PARAMS["tracks_extra"] ),
            suffix(".bed"),
            ".fasta" )
def exportMotifSequences( infile, outfile ):
    '''export sequences for all intervals.
    '''
    
    track = P.snip( infile, ".bed" )

    dbhandle = connect()

    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] ) )

    cc = dbhandle.cursor()
    statement = "SELECT interval_id, contig, start, end FROM %s_intervals" % P.quote( track )
    cc.execute( statement )

    outs = open( outfile, "w")

    for result in cc:
        interval, contig, start, end = result
        id = "%s_%s %s:%i..%i" % (track, str(interval), contig, start, end)
        
        seq = fasta.getSequence( contig, "+", start, end)
        outs.write( ">%s\n%s\n" % (id, seq))
    
    cc.close()
    outs.close()

############################################################
############################################################
############################################################
# (loadIntervalsFromBed, exportIntervalsAsBed ),
# TODO: fix, causes a problem due to exportIntervalsAsBed
@follows( buildIntervals )
@transform( [ "%s.bed" % x.asFile() for x in TRACKS_MASTER] +\
                P.asList( PARAMS["tracks_extra"] ),
            suffix(".bed"),
            ".controlfasta" )
def exportMotifControlSequences( infile, outfile ):
    '''for each interval, export the left and right 
    sequence segment of the same size.
    '''
    
    track = P.snip( infile, ".bed")

    dbhandle = connect()

    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] ) )

    table = "%s_intervals" % (P.quote( track ) )
    cc = dbhandle.cursor()
    statement = "SELECT interval_id, contig, start, end FROM %s" % table
    cc.execute( statement )

    outs = open( outfile, "w")

    for result in cc:
        interval, contig, xstart, xend = result
        l = xend - xstart
        lcontig = fasta.getLength( contig )
        start, end = max(0,xstart-l), xend-l
        id = "%s_%s_l %s:%i..%i" % (track, str(interval), contig, start, end)
        seq = fasta.getSequence( contig, "+", start, end)
        outs.write( ">%s\n%s\n" % (id, seq))

        start, end = xstart+l, min(lcontig,xend+l)
        id = "%s_%s_r %s:%i..%i" % (track, str(interval), contig, start, end)
        seq = fasta.getSequence( contig, "+", start, end)
        outs.write( ">%s\n%s\n" % (id, seq))
    
    cc.close()
    outs.close()

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
@files( [ ([ "%s.bed" % y.asFile() for y in EXPERIMENTS[x]], 
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
        tablename = "%s_intervals" % P.quote( P.snip(replicate, ".bed") )
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

    outs = open( outfile, "w" )
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
@files_re( ["%s.bed" % x.asFile() for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "peakval.correlation" )
def makePeakvalCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PipelineChipseq.makeIntervalCorrelation( infiles, outfile, "peakval", "merged.bed" )

@follows( buildIntervals )
@files_re( ["%s.bed" % x.asFile() for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "avgval.correlation" )
def makeAvgvalCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PipelineChipseq.makeIntervalCorrelation( infiles, outfile, "avgval", "merged.bed" )

@follows( buildIntervals )
@files_re( ["%s.bed" % x.asFile() for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "length.correlation" )
def makeLengthCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PipelineChipseq.makeIntervalCorrelation( infiles, outfile, "length", "merged.bed" )

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
    tmpfile = open( tmpname, "w")
    for line in open(infile,"r"):
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
    '''run MEME on all intervals and motifs.

    In order to increase the signal/noise ratio,
    MEME is not run on all intervals but only the 
    top 10% of intervals (peakval) are used. 
    Also, only the segment of 200 bp around the peak
    is used and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker
    '''
    to_cluster = True
    # job_options = "-l mem_free=8000M"

    target_path = os.path.join( os.path.abspath(PARAMS["exportdir"]), "meme", outfile )

    dbhandle = connect()
    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] ) )

    track = P.snip(infile, ".fasta")

    cc = dbhandle.cursor()
    statement = "SELECT peakcenter, interval_id, contig FROM %s_intervals ORDER BY peakval DESC" % P.quote(track)
    cc.execute( statement )
    data = cc.fetchall()
    cc.close()
    
    cutoff = max( int( len(data) * PARAMS["meme_proportion"] ) + 1, int( PARAMS["meme_min_sequences"] ) )
    E.info( "using %i out of %i peaks for motif detection" % (cutoff, len(data)))

    # maximum size of data set (in characters)
    maxsize = int(PARAMS["meme_max_size"])

    L.info( "runMeme: %s: using at most %i sequences for pattern finding" % (track, cutoff) )

    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] ) )

    tmpdir = tempfile.mkdtemp( dir = "." )
    tmpfasta =  os.path.join( tmpdir, "in.fa")
    outs = open( tmpfasta , "w" )

    masker = Masker.MaskerDustMasker()

    hw = int(PARAMS["meme_halfwidth"])
    current_size, nseq = 0, 0
    for peakcenter, id, contig in data[:cutoff]:
        start, end = peakcenter - hw, peakcenter + hw
        id = "%s_%s %s:%i..%i" % (track, str(id), contig, start, end)
        seq = fasta.getSequence( contig, "+", start, end )

        if PARAMS["motifs_masker"] == "repeatmasker":
            # the genome sequence is repeat masked
            masked_seq = seq
        elif PARAMS["motifs_masker"] == "dustmasker":
            masked_seq = masker( seq.upper() )

        # hard mask softmasked characters
        masked_seq = re.sub( "[a-z]","N", masked_seq )
        current_size += len(masked_seq)
        if current_size >= maxsize: 
            L.info( "runMEME: %s: maximum size (%i) reached - only %i sequences output" % (track, maxsize, nseq))
            break
        nseq += 1
        outs.write( ">%s\n%s\n" % (id, masked_seq))

    outs.close()

    if nseq == 0:
        E.warn( "%s: no sequences - meme skipped" % infile )
        P.touch( outfile )
    else:
        statement = '''
        meme %(tmpfasta)s -dna -revcomp -mod %(meme_model)s -nmotifs %(meme_nmotifs)s -oc %(tmpdir)s -maxsize %(maxsize)s %(meme_options)s > %(outfile)s.log
        '''
        P.run()

        # copy over results
        try:
            os.makedirs( os.path.dirname( target_path ) )
        except OSError: 
            # ignore "file exists" exception
            pass

        if os.path.exists( target_path ): shutil.rmtree( target_path )
        shutil.move( tmpdir, target_path )

        shutil.copyfile( os.path.join(target_path, "meme.txt"), outfile)

############################################################
############################################################
############################################################
def writeSequencesForIntervals( track, filename,
                                full = False,
                                halfwidth = None,
                                maxsize = None,
                                proportion = None,
                                masker = None,
                                offset = 0,
                                shuffled = False ):
    '''build a sequence set for motif discovery.

    If num_shuffles is set, shuffled copies are created as well with
    the shuffled number appended to the filename.

    The sequences are masked before shuffling (is this appropriate?)

    If *full* is set, the whole intervals will be output, otherwise
    only the region around the peak given by *halfwidth*

    If *maxsize* is set, the output is truncated at *maxsize* characters.

    If proportion is set, only the top *proportion* intervals are output
    (sorted by peakval).
    '''

    dbhandle = connect()

    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] ) )

    cc = dbhandle.cursor()
        
    tablename = "%s_intervals" % P.quote( track )
    if full:
        statement = '''SELECT start, end, interval_id, contig 
                       FROM %(tablename)s 
                       ORDER BY peakval DESC''' % locals()
    elif halfwidth:
        statement = '''SELECT peakcenter - %(halfwidth)s, peakcenter + %(halfwidth)s,
                       interval_id, contig 
                       FROM %(tablename)s 
                       ORDER BY peakval DESC''' % locals()
    else:
        raise ValueError("either specify full or halfwidth" )
    
    cc.execute( statement )
    data = cc.fetchall()
    cc.close()

    if proportion:
        cutoff = int(len(data) * proportion) + 1
    else:
        cutoff = len(data)
        
        L.info( "writeSequencesForIntervals %s: using at most %i sequences for pattern finding" % (track, cutoff) )
    L.info( "writeSequencesForIntervals %s: masker=%s" % (track,masker))

    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"]) )
    masker_object = Masker.MaskerDustMasker()

    sequences = []
    current_size, nseq = 0, 0
    new_data = []
    for start, end, id, contig in data[:cutoff]:
        lcontig = fasta.getLength( contig )
        start, end = max(0, start + offset), min(end + offset, lcontig)
        if start >= end:
            L.info( "writeSequencesForIntervals %s: sequence %s is empty: start=%i, end=%i, offset=%i - ignored" % (track, id, start, end, offset))
            continue
        seq = fasta.getSequence( contig, "+", start, end )
        sequences.append( seq )
        new_data.append( (start, end, id, contig) )
        current_size += len(seq)
        if maxsize and current_size >= maxsize: 
            L.info( "writeSequencesForIntervals %s: maximum size (%i) reached - only %i sequences output (%i ignored)" % (track, maxsize, nseq,
                                                                                                                          len(data) - nseq ))
            break
        nseq += 1
        
    data = new_data
            
    def maskSequences( sequences, masker ):
        
        if masker == "repeatmasker":
            # the genome sequence is repeat masked
            masked_seq = sequences
        elif masker == "dust":
            masked_seq = [ masker_object( x.upper() ) for x in sequences ]
        else:
            masked_seq = [x.upper() for x in sequences ]

        # hard mask softmasked characters
        masked_seq = [re.sub( "[a-z]","N", x) for x in masked_seq ]
        return masked_seq

    if shuffled:
        # note that shuffling is done on the unmasked sequences
        # Otherwise N's would be interspersed with real sequence
        # messing up motif finding unfairly. Instead, masking is
        # done on the shuffled sequence.
        sequences = [ list(x) for x in sequences ]
        for sequence in sequences: random.shuffle(sequence)
        sequences = maskSequences( ["".join(x) for x in sequences ], masker )
        
    outs = open(filename, "w" )
    for sequence, d in zip( maskSequences( sequences, masker ), data ):
        start, end, id, contig = d
        id = "%s_%s %s:%i..%i" % (track, str(id), contig, start, end)
        outs.write( ">%s\n%s\n" % (id, sequence ) )
    outs.close()
    
    return len(sequences)

# writeSequencesForIntervals( "runEstSub",
#                             "test1.fasta",
#                             full = True,
#                             offset = 0 )

# writeSequencesForIntervals( "runEstSub",
#                             "test2.fasta",
#                             full = True,
#                             offset = -10000)


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
    to_cluster = True

    target_path = os.path.join( os.path.abspath( PARAMS["exportdir"] ), "glam2", outfile )
    track = infile[:-len(".fasta")]

    tmpdir = tempfile.mkdtemp()
    tmpfasta =  os.path.join( tmpdir, "in.fa")
    
    nseq = writeSequencesForIntervals( track, tmpfasta,
                                       full = False,
                                       halfwidth = int(PARAMS["meme_halfwidth"]),
                                       maxsize = int(PARAMS["meme_max_size"]),
                                       proportion = PARAMS["meme_proportion"] )

    min_sequences = int(nseq / 10.0)
    statement = '''
    %(execglam2)s -2 -O %(tmpdir)s %(glam2_options)s -z %(min_sequences)i n %(tmpfasta)s > %(outfile)s.log
    '''
    P.run()

    # copy over results
    try:
        os.makedirs( os.path.dirname( target_path ) )
    except OSError: 
        # ignore "file exists" exception
        pass

    if os.path.exists( target_path ): shutil.rmtree( target_path )
    shutil.move( tmpdir, target_path )

    shutil.copyfile( os.path.join(target_path, "glam2.txt"), outfile)


############################################################
############################################################
############################################################
## selecting which motifs to run
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

        for line in open( infile, "r" ):
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
        for line in open( infile, "r" ):
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

############################################################
############################################################
############################################################
@merge( (buildReferenceMotifs, filterMotifs), "motif_info.load" )
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
@follows( buildReferenceMotifs, filterMotifs )
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
    to_cluster = True

    # job_options = "-l mem_free=8000M"

    controlfile, dbfile, motiffiles  = infiles
    controlfile = dbfile[:-len(".fasta")] + ".controlfasta"

    if IOTools.isEmpty( dbfile ):
        P.touch( outfile )
        return

    if not os.path.exists( controlfile ):
        raise P.PipelineError( "control file %s for %s does not exist" % (controlfile, dbfile))

    # remove previous results
    if os.path.exists(outfile): os.remove( outfile )
    
    tmpdir = P.getTempDir( "." )
    tmpfile = P.getTempFilename( "." )

    for motiffile in motiffiles:
        if IOTools.isEmpty( motiffile ):
            L.info( "skipping empty motif file %s" % motiffile )
            continue
        
        of = open(tmpfile, "a")
        motif, x = os.path.splitext( motiffile )
        of.write(":: motif = %s ::\n" % motif )
        of.close()
        
        statement = '''
        cat %(dbfile)s %(controlfile)s 
        | mast %(motiffile)s - -nohtml -oc %(tmpdir)s -ev %(mast_evalue)f %(mast_options)s >> %(outfile)s.log;
        cat %(tmpdir)s/mast.txt >> %(tmpfile)s
        '''
        P.run()

    statement = "gzip < %(tmpfile)s > %(outfile)s" 
    P.run()

    shutil.rmtree( tmpdir )
    os.unlink( tmpfile )

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

    # bioprospector currently not working on the nodes
    to_cluster = False

    # only use new nodes, as /bin/csh is not installed
    # on the old ones.
    # job_options = "-l mem_free=8000M"

    tmpfasta = P.getTempFilename( "." )
    track = outfile[:-len(".bioprospector")]
    nseq = writeSequencesForIntervals( track,
                                       tmpfasta,
                                       full = True,
                                       masker = "dust",
                                       proportion = PARAMS["bioprospector_proportion"] )

    if nseq == 0:
        E.warn( "%s: no sequences - bioprospector skipped" % track )
        P.touch( outfile )
    else:
        statement = '''
    BioProspector -i %(tmpfasta)s %(bioprospector_options)s -o %(outfile)s > %(outfile)s.log
    '''
        P.run()

    os.unlink( tmpfasta )

############################################################
############################################################
############################################################
##
############################################################
@transform( runBioProspector, suffix(".bioprospector"), "_bioprospector.load")
def loadBioProspector( infile, outfile ):
    '''load results from bioprospector.'''

    tablename = outfile[:-len(".load")]
    target_path = os.path.join( os.path.abspath( PARAMS["exportdir"] ), "bioprospector" )

    try:
        os.makedirs( target_path )
    except OSError: 
        pass

    track = infile[:-len(".bioprospector")]

    results = Bioprospector.parse( open(infile, "r") )

    tmpfile = P.getTempFile()
    tmpfile.write( "id\tmotif\tstart\tend\tstrand\tarrangement\n" )
    
    for x, motifs in enumerate( results ):
        outname = os.path.join( target_path, "%s_%02i.png" % (track, x ) )
        Bioprospector.build_logo( [y.sequence for y in motifs.matches],
                                  outname )

        for match in motifs.matches:

            distance = abs(match.start + match.width1 - (match.end - match.width2))

            if match.strand in ( "+-", "-+"):
                arrangement = "ER" 
            elif match.strand in ("++", "--"):
                arrangement = "DR"
            else:
                arrangement = "SM"
                distance = 0
                
            arrangement += "%i" % distance
            strand = match.strand[0]

            id = re.sub( ".*_", "", match.id )
            tmpfile.write( "%s\t%i\t%i\t%i\t%s\t%s\n" % \
                           (id,
                            x,
                            match.start,
                            match.end,
                            strand,
                            arrangement) )
    tmpfile.close()
    tmpfilename = tmpfile.name

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
    --allow-empty \
    -b sqlite \
    --index=id \
    --index=motif \
    --index=id,motif \
    --table=%(tablename)s \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
# todo: fix, causes a problem: exportSequences and exportMASTControlSequences,
@files_re( (exportMotifSequences, exportMotifControlSequences),
           "(\S+).controlfasta",
           [ r"\1.controlfasta", r"\1.fasta"],
           r"\1.regexmotif" )
def runRegexMotifSearch( infiles, outfile ):
    '''run a regular expression search on sequences.
    compute counts.
    '''

    motif = "[AG]G[GT]T[CG]A"
    reverse_motif = "T[GC]A[CA]C[TC]"

    controlfile, dbfile  = infiles
    controlfile = dbfile[:-len(".fasta")] + ".controlfasta"
    if not os.path.exists( controlfile ):
        raise P.PipelineError( "control file %s for %s does not exist" % (controlfile, dbfile))

    motifs = []
    for x in range( 0, 15):
        motifs.append( ("DR%i" % x, re.compile( motif + "." * x + motif, re.IGNORECASE ) ) )
    for x in range( 0, 15):
        motifs.append( ("ER%i" % x, re.compile( motif + "." * x + reverse_motif, re.IGNORECASE ) ) )

    db_positions = Motifs.countMotifs( open( dbfile, "r" ), motifs )
    control_positions = Motifs.countMotifs( open( controlfile, "r" ), motifs )

    db_counts, control_counts = Motifs.getCounts( db_positions), Motifs.getCounts( control_positions )
    db_seqcounts, control_seqcounts = Motifs.getOccurances( db_positions), Motifs.getCounts( control_positions )
    
    ndb, ncontrol = len( db_positions), len( control_positions )
    outf = open( outfile, "w" )
    outf.write( "motif\tmotifs_db\tmotifs_control\tseq_db\tseq_db_percent\tseq_control\tseq_control_percent\tfold\n" )
    for motif, pattern in motifs:
        try:
            fold = float(db_seqcounts[motif]) * ncontrol / (ndb * control_seqcounts[motif])
        except ZeroDivisionError:
            fold = 0
            
        outf.write( "%s\t%i\t%i\t%i\t%s\t%i\t%s\t%5.2f\n" % \
                    (motif,
                     db_counts[motif],
                     control_counts[motif],
                     db_seqcounts[motif],
                     IOTools.prettyPercent( db_seqcounts[motif], ndb),
                     control_seqcounts[motif],
                     IOTools.prettyPercent( control_seqcounts[motif], ncontrol),
                     fold) )

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

    tablename = P.toTable( outfile )

    tmpfile = tempfile.NamedTemporaryFile(delete=False)

    tmpfile.write( MAST.Match().header +\
                   "\tmotif\tcontig" \
                   "\tl_evalue\tl_pvalue\tl_nmatches\tl_length\tl_start\tl_end" \
                   "\tr_evalue\tr_pvalue\tr_nmatches\tr_length\tr_start\tr_end" \
                   "\tmin_evalue\tmin_pvalue\tmax_nmatches" + "\n" )

    lines = IOTools.openFile(infile).readlines()
    chunks = [x for x in range(len(lines)) if lines[x].startswith("::") ]
    chunks.append( len(lines) )

    for chunk in range(len(chunks)-1):
        # use real file, as MAST parser can not deal with a
        # list of lines
        tmpfile2 = tempfile.NamedTemporaryFile(delete=False)        
        try:
            motif = re.match( ":: motif = (\S+) ::", lines[chunks[chunk]]).groups()[0]
        except AttributeError:
            raise P.PipelineError("parsing error in line '%s'" % lines[chunks[chunk]])

        tmpfile2.write( "".join( lines[chunks[chunk]+1:chunks[chunk+1]]) )
        tmpfile2.close()
        mast = MAST.parse( open(tmpfile2.name, "r") )
        os.unlink( tmpfile2.name )        

        # collect control data
        full_matches = []
        controls = collections.defaultdict( dict )
        for match in mast.matches:
            m = match.id.split("_")
            if len(m) == 2:
                full_matches.append( match )
                continue

            track, id, pos = m
            controls[id][pos] = (match.evalue, match.pvalue, match.nmotifs, match.length, match.start, match.end )

        for match in full_matches:
            match.id = match.id.split("_")[1]
            # move to genomic coordinates
            contig, start, end = re.match( "(\S+):(\d+)..(\d+)", match.description).groups()
            if match.nmotifs > 0:
                start, end = int(start), int(end)
                match.start += start
                match.end += start
                match.positions = [ x + start for x in match.positions ]

            id = match.id
            if id not in controls:
                P.warn( "no controls for %s - increase MAST evalue" % id )

            if "l" not in controls[id]: controls[id]["l"] = (float(PARAMS["mast_evalue"]), 1, 0, 0, 0, 0)
            if "r" not in controls[id]: controls[id]["r"] = (float(PARAMS["mast_evalue"]), 1, 0, 0, 0, 0)

            min_evalue = min( controls[id]["l"][0], controls[id]["r"][0])
            min_pvalue = min( controls[id]["l"][1], controls[id]["r"][1])
            max_nmatches = max( controls[id]["l"][2], controls[id]["r"][2])

            tmpfile.write( str(match) + "\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
                               (motif,contig,
                                "\t".join( map(str, controls[id]["l"] )),
                                "\t".join( map(str, controls[id]["r"] )),
                                str(min_evalue),
                                str(min_pvalue),
                                str(max_nmatches),
                                ) + "\n" )
            
    tmpfile.close()
    tmpfilename = tmpfile.name

    statement = '''
    python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              -b sqlite 
              --index=id 
              --index=motif 
              --index=id,motif 
              --table=%(tablename)s 
              --allow-empty
              --map=base_qualities:text 
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run()
    os.unlink( tmpfile.name )

############################################################
############################################################
############################################################
##
############################################################
#    @files_re( "run*.controlfasta",
#           "(\S+).controlfasta",
#           [ r"\1.controlfasta", r"\1.fasta", glob.glob("*.glam2")],
#           r"\1.glam2scan" )

@follows( buildReferenceMotifs, runGLAM2 )
@files_re( (exportMotifSequences, exportMotifControlSequences),
           "(\S+).controlfasta",
           [ r"\1.controlfasta", r"\1.fasta",  glob.glob("*.glam2")],
           r"\1.glam2scan" )
def runGLAM2SCAN( infiles, outfile ):
    '''run glam2scan on all intervals and motifs.
    '''

    to_cluster = True
    # only use new nodes, as /bin/csh is not installed
    # on the old ones.
    # job_options = "-l mem_free=8000M"

    controlfile, dbfile, motiffiles  = infiles
    controlfile = dbfile[:-len(".fasta")] + ".controlfasta"
    if not os.path.exists( controlfile ):
        raise P.PipelineError( "control file %s for %s does not exist" % (controlfile, dbfile))

    if os.path.exists(outfile): os.remove( outfile )
    
    for motiffile in motiffiles:
        of = open(outfile, "a")
        motif, x = os.path.splitext( motiffile )
        of.write(":: motif = %s ::\n" % motif )
        of.close()
        
        statement = '''
        cat %(dbfile)s %(controlfile)s | %(execglam2scan)s -2 -n %(glam2scan_results)i n %(motiffile)s - >> %(outfile)s
        '''
        P.run( **dict( locals().items() + PARAMS.items() ) )
        
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

    tablename = outfile[:-len(".load")]
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    tmpfile.write( "motif\tid\tnmatches\tscore\tscores\tncontrols\tmax_controls\n" )

    lines = open(infile).readlines()
    chunks = [x for x in range(len(lines)) if lines[x].startswith("::") ]
    chunks.append( len(lines) )

    for chunk in range(len(chunks)-1):

        # use real file, as parser can not deal with a
        # list of lines
        
        try:
            motif = re.match( ":: motif = (\S+) ::", lines[chunks[chunk]]).groups()[0]
        except AttributeError:
            raise P.PipelineError("parsing error in line '%s'" % lines[chunks[chunk]])

        if chunks[chunk]+1 == chunks[chunk+1]:
            L.warn( "no results for motif %s - ignored" % motif )
            continue
        
        tmpfile2 = tempfile.NamedTemporaryFile(delete=False)        
        tmpfile2.write( "".join( lines[chunks[chunk]+1:chunks[chunk+1]]) )
        tmpfile2.close()
        glam = Glam2Scan.parse( open(tmpfile2.name, "r") )
                    
        os.unlink( tmpfile2.name )        

        # collect control data
        full_matches = collections.defaultdict( list )
        controls = collections.defaultdict( list )
        for match in glam.matches:
            m = match.id.split("_")
            track, id = m[:2]
            if len(m) == 2:
                full_matches[id].append( match )
            else:
                controls[id].append( match.score )

        for id, matches in full_matches.iteritems():
            
            nmatches = len(matches)
            scores = [x.score for x in matches ]
            score = max(scores)
            # move to genomic coordinates
            #contig, start, end = re.match( "(\S+):(\d+)..(\d+)", match.id).groups()
            #start, end = int(start), int(end)
            #match.start += start
            #match.end += start
            contig = ""

            if id not in controls:
                P.warn( "no controls for %s - increase evalue?" % id )
            
            c = controls[id]
            if len(c) == 0: mmax = "" 
            else: mmax = max(c)

            tmpfile.write( "\t".join( map(str,
                    (motif, id, 
                     nmatches,
                     score,
                     ",".join(map(str,scores)),
                     len(c),
                     mmax))) + "\n" )

    tmpfile.close()
    tmpfilename = tmpfile.name

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              -b sqlite \
              --index=id \
              --index=motif \
              --index=id,motif \
              --table=%(tablename)s \
              --map=base_qualities:text \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run()
    os.unlink( tmpfile.name )

############################################################
############################################################
############################################################
@follows( buildIntervals )
@files( [ ("%s.bed" % x.asFile(), 
           "%s.annotations" % x.asFile() ) for x in TRACKS_MASTER ] )
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
############################################################
############################################################
@follows( buildIntervals )
@files( [ ("%s.bed" % x.asFile(), 
           "%s.tss" % x.asFile() ) for x in TRACKS_MASTER ] )
def annotateTSS( infile, outfile ):
    '''compute distance to TSS'''

    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_tss_bed"] )

    statement = """
    cat < %(infile)s 
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
@files( [ ("%s.bed" % x.asFile(), 
           "%s.repeats" % x.asFile() ) for x in TRACKS_MASTER ] )
def annotateRepeats( infile, outfile ):
    '''count the overlap between intervals and repeats.'''

    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_repeats_gff"] )

    statement = """
    cat < %(infile)s |\
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
@files( [ ("%s.bed" % x.asFile(), "%s.readcounts" % x.asFile() ) for x in TOSUBTRACT ] )
def buildIntervalCounts( infile, outfile ):
    '''count read density in bed files comparing stimulated versus
    unstimulated binding.
    '''
    track = TRACKS.factory( filename = outfile[:-len(".readcounts")] )

    samfiles_fg, samfiles_bg = [], []

    # collect foreground and background bam files
    fg_replicates = PipelineTracks.Aggregate( TRACKS, track = track )[track]
    for replicate in fg_replicates:
        samfiles_fg.append( "%s.call.bam" % replicate.asFile() )

    unstim = getUnstimulated( track )
    bg_replicates = PipelineTracks.Aggregate( TRACKS, track = unstim )[unstim]

    for replicate in bg_replicates:
        samfiles_bg.append( "%s.call.bam" % replicate.asFile())
        
    samfiles_fg = [ x for x in samfiles_fg if os.path.exists( x ) ]
    samfiles_bg = [ x for x in samfiles_bg if os.path.exists( x ) ]

    samfiles_fg = ",".join(samfiles_fg)
    samfiles_bg = ",".join(samfiles_bg)

    tmpfile1 = P.getTempFilename( os.getcwd() ) + ".fg"
    tmpfile2 = P.getTempFilename( os.getcwd() ) + ".bg"

    # start counting
    to_cluster = True

    statement = """
    cat < %(infile)s 
    | python %(scriptsdir)s/bed2gff.py --as-gtf 
    | python %(scriptsdir)s/gtf2table.py 
                --counter=read-coverage 
                --log=%(outfile)s.log 
                --bam-file=%(samfiles_fg)s 
    > %(tmpfile1)s"""
    P.run()

    if samfiles_bg:
        statement = """
        cat < %(infile)s 
        | python %(scriptsdir)s/bed2gff.py --as-gtf 
        | python %(scriptsdir)s/gtf2table.py 
                    --counter=read-coverage 
                    --log=%(outfile)s.log 
                    --bam-file=%(samfiles_bg)s 
        > %(tmpfile2)s"""
        P.run()

        statement = '''
        python %(toolsdir)s/combine_tables.py 
               --add-file-prefix 
               --regex-filename="[.](\S+)$" 
        %(tmpfile1)s %(tmpfile2)s > %(outfile)s
        '''

        P.run()

        os.unlink( tmpfile2 )
        
    else:
        statement = '''
        python %(toolsdir)s/combine_tables.py 
               --add-file-prefix 
               --regex-filename="[.](\S+)$" 
        %(tmpfile1)s > %(outfile)s
        '''

        P.run()

    os.unlink( tmpfile1 )

############################################################
@transform( buildIntervalCounts, 
            suffix(".readcounts"), 
            "_readcounts.load" )
def loadIntervalCounts( infile, outfile ):
    '''load interval counts.'''
    P.load( infile, outfile, "--index=gene_id" )

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
############################################################
############################################################
## export intervals
############################################################
@follows ( mkdir( 'export' ))
@merge( "*.bed", ("export/intervals_%s.bed" % PARAMS["version"], "intervals.view") )
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

@follows( buildIntervals, 
          makeReadCorrelationTable,
          loadBAMStats)
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
          runTomTom, loadTomTom,
          runBioProspector )
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

@follows( intervals, 
          discover_motifs, 
          detect_motifs,
          correlation, 
          annotation )
def full():
    '''run the full pipeline.'''
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

###################################################################
###################################################################
###################################################################
@follows( createViewMapping )
def views():
    pass

###################################################################
###################################################################
###################################################################
@follows( views, mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )

###################################################################
###################################################################
###################################################################
@follows( views, mkdir( "report" ) )
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


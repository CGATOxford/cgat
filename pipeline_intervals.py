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
Interval pipeline
=================

:Author: Andreas Heger
:Release: $Id: pipeline_chipseq.py 2900 2010-04-13 14:38:00Z andreas $
:Date: |today|
:Tags: Python

The interval pipeline takes a several set of :term:`bed` formatted genomic intervals 
and annotates them.

It performs the following analyses:
   * Peak location
      * requires a bam-file to be associated with each bed-file.
   * Motif discovery using MEME
   * Motif detection using MAST
      * requires a set of known motifs
   * Overlap with genomic context.
   * Aggregate transcript/gene profiles
   * 

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
Sample = PipelineTracks.Sample

TRACKS = PipelineTracks.Tracks( Sample ).loadFromDirectory( glob.glob("*.bed.gz"),
                                                            "(\S+).bed.gz" )

TRACKS_BEDFILES = ["%s.bed.gz" % x for x in TRACKS]

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
# 
###################################################################
def getAssociatedBAMFiles( track ):
    '''return a list of BAM files associated with a track.

    By default, this method searches for ``track.bam`` 
    file in the current directory and returns an offset of 0.

    Associations can be defined in the .ini file in the section
    [bams]. For example, the following snippet associates track
    track1 with the bamfiles :file:`track1.bam` and :file:`track2.bam`::
    
       [bams]
       track1=track1.bam,track2.bam

    Glob expressions are permitted.

    Offsets are used to shift tags in ChIP experiments. Offsets
    need to be defined in the [offsets] sections. If no offsets
    are defined, the method returns a list of 0 offsets.

    Offsets need to be defined in the same order as the bam files::
    
       [offsets]
       track1=120,200

    returns a list of BAM files and offsets.
    '''
    fn = track.asFile()
    bamfiles = []
    if "bams_%s" % fn.lower() in PARAMS:
        for ff in P.asList( PARAMS["bams_%s" % fn.lower() ] ):
            bamfiles.extend( glob.glob( ff ) )
    else:
        bamfiles = glob.glob( "%s.bam" % fn )
        
        
    if "offsets_%s" % fn.lower() in PARAMS:
        offsets = map(int, P.asList( PARAMS["offsets_%s" % fn.lower() ] ))
    else:
        offsets = [0] * len(bamfiles)

    return bamfiles, offsets

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

############################################################
############################################################
############################################################
@transform( TRACKS_BEDFILES,
            suffix(".bed.gz"),
            "_intervals.load" )
def loadIntervals( infile, outfile ):
    '''load intervals from :term:`bed` formatted files into database.
    
    These :term:`bed` formatted intervals are derived by 
    merging/intersecting the various tracks or have been
    placed explicitely into the directory.

    Also, if a :term:`bam` file is associated with a :term:`bed`
    file, re-evaluate the intervals by counting reads within
    the interval. In contrast to the initial pipeline, the
    genome is not binned. 

       nprobes: number of reads in interval
       peakcenter: position with maximum number of reads in interval
       avgval: average coverage within interval

    If *replicates* is true, only replicates will be considered
    for the counting. Otherwise the counts aggregate both replicates
    and conditions.
    '''

    bedfile = infile

    tmpfile = P.getTempFile()

    headers = ("avgval","disttostart","genelist","length","peakcenter","peakval","position","interval_id","npeaks","nprobes", "contig","start","end" )

    tmpfile.write( "\t".join(headers) + "\n" )

    avgval,contig,disttostart,end,genelist,length,peakcenter,peakval,position,start,interval_id,npeaks,nprobes = \
        0,"",0,0,"",0,0,0,0,0,0,0,0

    track = Sample( filename = P.snip( infile, ".bed.gz") )

    bamfiles, offsets = getAssociatedBAMFiles( track )

    if bamfiles:
        E.info( "%s: associated bamfiles = %s" % (track, bamfiles))
    else:
        E.info( "%s: no bamfiles associated" % (track))

    # open all bamfiles
    samfiles = [ pysam.Samfile( fn,  "rb" ) for fn in bamfiles ]

    c = E.Counter()

    # count tags
    for bed in Bed.iterator( IOTools.openFile(infile, "r") ): 

        c.input += 1

        if "name" not in bed:
            bed.name = c.input
        
        if samfiles:
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
                                             npeaks,nprobes, 
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
@follows( mkdir( os.path.join( PARAMS["exportdir"], "bed" ) ) )
@transform( TRACKS_BEDFILES,
            regex(r"(.*).bed.gz"),
            os.path.join( PARAMS["exportdir"], "bed", r"\1.bed.gz") )
def indexIntervals( infile, outfile ):
    '''index intervals.
    '''

    statement = '''zcat %(infile)s | bgzip > %(outfile)s; tabix -p bed %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################
@transform( TRACKS_BEDFILES,
            suffix(".bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_genomic_context_bed"] ) ),
            ".contextstats" )
def buildContextStats( infiles, outfile ):
    '''build mapping context stats.

    Examines the genomic context to where reads align.

    A read is assigned to the genomic context that it
    overlaps by at least 50%. Thus some reads mapping
    several contexts might be dropped.
    '''

    infile, reffile = infiles

    min_overlap = 0.5

    to_cluster = True
    statement = '''
       python %(scriptsdir)s/rnaseq_bam_vs_bed.py
              --min-overlap=%(min_overlap)f
              --log=%(outfile)s.log
              %(infile)s %(reffile)s
       > %(outfile)s
       '''

    P.run()

############################################################
############################################################
############################################################
@merge( buildContextStats, "context_stats.load" )
def loadContextStats( infiles, outfile ):
    """load context mapping statistics."""

    header = ",".join( [P.snip( os.path.basename(x), ".contextstats") for x in infiles] )
    filenames = " ".join( infiles  )
    tablename = P.toTable( outfile )

    statement = """python %(scriptsdir)s/combine_tables.py
                      --headers=%(header)s
                      --missing=0
                      --skip-titles
                   %(filenames)s
                | perl -p -e "s/bin/track/; s/\?/Q/g"
                | python %(scriptsdir)s/table2table.py --transpose
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
                """
    P.run()

############################################################
############################################################
############################################################
@transform( TRACKS_BEDFILES,
            suffix(".bed.gz"),
            ".annotations" )
def annotateIntervals( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set.
    '''
    to_cluster = True

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_annotation_gff"] )

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
@transform( TRACKS_BEDFILES,
            suffix(".bed.gz"),
            ".binding.tsv.gz" )
def annotateBinding( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set.
    '''
    to_cluster = True

    geneset = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_geneset_all_gtf"] )

    statement = """
    zcat < %(geneset)s
    | python %(scriptsdir)s/gtf2table.py 
		--counter=position 
		--counter=binding-pattern
		--log=%(outfile)s.log 
		--filename-gff=%(infile)s
		--genome-file=%(genome_dir)s/%(genome)s
    | gzip
    > %(outfile)s"""
    
    P.run()

############################################################
############################################################
############################################################
@transform( TRACKS_BEDFILES,
            suffix(".bed.gz"),
            ".tss" )
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
@transform( TRACKS_BEDFILES,
            suffix(".bed.gz"),
            ".repeats" )
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
@transform( annotateBinding, suffix(".binding.tsv.gz"), "_binding.load" )
def loadBinding( infile, outfile ):
    '''load interval binding: genome architecture
    '''
    P.load( infile, outfile, "--index=gene_id --allow-empty --map=pattern:str" )

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

###################################################################
###################################################################
###################################################################
@follows( mkdir( "transcriptprofiles" ) )
@transform( indexIntervals,
            regex(r".*/([^/].*)\.bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_geneset_all_gtf"] )),
            r"transcriptprofiles/\1.transcriptprofile.tsv.gz" )
def buildIntervalProfileOfTranscripts( infiles, outfile ):
    '''build a table with the overlap profile of transcripts.'''
    
    to_cluster = True

    bedfile, gtffile = infiles
    
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
    
############################################################
############################################################
############################################################
## master target for this section
############################################################
@follows( loadIntervals ) 
def buildIntervals():
    pass

###################################################################
###################################################################
###################################################################
@follows( mkdir( "peakshapes" ) )
@transform( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            r"peakshapes/\1.peakshape.tsv.gz" )
def buildPeakShapeTable( infile, outfile ):
    '''build a table with peak shape parameters.'''
    
    to_cluster = True

    track = TRACKS.factory( filename = infile[:-len(".bed.gz")] )
    
    track = Sample( filename = P.snip( infile, ".bed.gz") )

    bamfiles, offsets = getAssociatedBAMFiles( track )

    if bamfiles:
        E.info( "%s: associated bamfiles = %s" % (track, bamfiles))
    else:
        E.warn( "%s: no bamfiles associated - target skipped" % (track))
        P.touch( outfile )
        return

    if len(bamfiles) > 1:
        raise NotImplementedError( "peakshape with multiple bamfiles not implement" )

    shift = offsets[0]
    bamfile = bamfiles[0]

    if shift:
        E.info( "applying read shift %i for track %s" % (shift, track ) )

    statement = '''python %(scriptsdir)s/bam2peakshape.py
                      --window-size=%(peakshape_window_size)i
                      --bin-size=%(peakshape_bin_size)i
                      --output-filename-pattern="%(outfile)s.%%s"
                      --force
                      --shift=%(shift)i
                      --sort=peak-height
                      --sort=peak-width
                      --log=%(outfile)s.log
                      %(bamfile)s %(infile)s
                   | gzip
                   > %(outfile)s
                '''
    P.run()

@transform( buildPeakShapeTable, suffix(".tsv.gz"), ".load" )
def loadPeakShapeTable( infile, outfile ):
    '''load peak shape information.'''
    P.load( infile, outfile, "--ignore-column=bins --ignore-column=counts --allow-empty" )

###################################################################
###################################################################
###################################################################
## general import
###################################################################

############################################################
############################################################
############################################################
@merge( TRACKS_BEDFILES,
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

############################################################
############################################################
############################################################
@transform( loadIntervals,
            suffix("_intervals.load"),
            ".motifs.fasta" )
def exportMotifSequences( infile, outfile ):
    '''export sequences for motif discovery.

    This method requires the _interval tables.

    For motif discovery, only the sequences with the highest S/N ratio are supplied.
    
    1. The top *motifs_proportion* intervals sorted by peakval
    2. Only a region +/- *motifs_halfwidth* around the peak 
    3. At least *motifs_min_sequences*. If there are not enough sequences
          to start with, all will be used.
    4. At most *motifs_max_size* sequences will be output.
    '''
    track = P.snip( infile, "_intervals.load" )
    dbhandle = connect()
        
    nseq = PipelineMotifs.writeSequencesForIntervals( track, 
                                                      outfile,
                                                      dbhandle,
                                                      full = False,
                                                      masker = PARAMS['motifs_masker'],
                                                      halfwidth = int(PARAMS["motifs_halfwidth"]),
                                                      maxsize = int(PARAMS["motifs_max_size"]),
                                                      proportion = PARAMS["motifs_proportion"],
                                                      min_sequences = PARAMS["motifs_min_sequences"] )

    if nseq == 0:
        E.warn( "%s: no sequences - meme skipped" % outfile)
        P.touch( outfile )


############################################################
############################################################
############################################################
@follows( mkdir( "motifs" ) )
@transform( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            r"motifs/\1.foreground.fasta" )
def exportMotifIntervalSequences( infile, outfile ):
    '''export sequences for motif detection.

    This method requires the _interval tables.
    '''
    PipelineMotifs.exportSequencesFromBedFile( infile, outfile,
                                               masker = PARAMS['motifs_masker'])

@follows( mkdir( "motifs" ) )
@transform( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            r"motifs/\1.control.fasta" )
def exportMotifControlSequences( infile, outfile ):
    '''for each interval, export the left and right 
    sequence segment of the same size.
    '''
    PipelineMotifs.exportSequencesFromBedFile( infile, outfile,
                                               masker = PARAMS['motifs_masker'],
                                               mode = "leftright" )


############################################################
############################################################
############################################################
@transform( exportMotifSequences, suffix(".motifs.fasta"), ".meme")
def runMeme( infile, outfile ):
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

    track = P.snip( infile, ".motifs.fasta" )

    PipelineMotifs.runMEMEOnSequences( infile, outfile )

############################################################
############################################################
############################################################
@merge( runMeme, "meme_summary.load" )
def loadMemeSummary( infiles, outfile ):
    '''load information about motifs into database.'''
    
    outf = P.getTempFile()

    outf.write("track\n" )

    for infile in infiles:
        if IOTools.isEmpty( infile ): continue
        motif = P.snip( infile, ".meme" )
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
        --log=%(outfile)s
    < %(infile)s
    | python %(scriptsdir)s/csv2db.py
        %(csv2db_options)s
        --table=%(tablename)s
    > %(outfile)s'''
    
    P.run()

############################################################
############################################################
############################################################
@merge( "*.motif", "motif_info.load" )
def loadMotifInformation( infiles, outfile ):
    '''load information about motifs into database.'''
    
    outf = P.getTempFile()

    outf.write("motif\n" )

    for infile in infiles:
        if IOTools.isEmpty( infile ): continue
        motif = P.snip( infile, ".motif" )
        outf.write( "%s\n" % motif )

    outf.close()

    P.load( outf.name, outfile, "--allow-empty" )
    
    os.unlink( outf.name )

############################################################
############################################################
############################################################
@files_re( (exportMotifIntervalSequences, exportMotifControlSequences),
           "(\S+).control.fasta",
           [ r"\1.control.fasta", r"\1.foreground.fasta",  glob.glob("*.motif")],
           r"\1.mast.gz" )
def runMast( infiles, outfile ):
    '''run mast on all intervals and motifs.

    Collect all results for an E-value up to 10000 so that
    all sequences are output and MAST curves can be computed. 

    10000 is a heuristic.
    '''
    PipelineMotifs.runMAST( infiles, outfile )

############################################################
############################################################
############################################################
@transform( runMast,
            suffix(".mast.gz"),
            "_mast.load" )
def loadMast( infile, outfile ):
    '''parse mast file and load into database.

    Parse several motif runs and add them to the same
    table.

    Add columns for the control data as well.
    '''
    PipelineMotifs.loadMAST( infile, outfile )

############################################################
############################################################
############################################################
## compute overlap with genomic features
############################################################
@follows( mkdir("context_gat.dir") )
@transform( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_genomic_context_bed"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_mapability_bed"] % PARAMS["gat_mapability"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_gc_profile_bed"] ),
                        ),
            r"context_gat.dir/\1.gat.tsv.gz" )
def runGATOnGenomicContext( infiles, outfile ):
    '''run gat agains genomic context.

    The workspace is composed of all mapable regions.
    Enrichment is controlled by isochores.

    To be rigorous, FDR should be re-computed after merging all
    analyses.
    '''
    
    bedfile, annofile, workspacefile, isochorefile = infiles

    to_cluster = True
    outdir = "context_gat.dir"

    statement = '''gatrun.py
         --segments=%(bedfile)s
         --annotations=%(annofile)s
         --workspace=%(workspacefile)s
         --num-samples=%(gat_num_samples)i
         --isochores=%(isochorefile)s
         --force
         --ignore-segment-tracks
         --output-filename-pattern=%(outfile)s.%%s
         -v 5
         --log=%(outfile)s.log
         | gzip
         > %(outfile)s'''

    P.run()

@follows( mkdir("annotations_gat.dir") )
@transform( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_annotation_gff"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_mapability_bed"] % PARAMS["gat_mapability"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_gc_profile_bed"] ), 
                        ),
            r"annotations_gat.dir/\1.gat.tsv.gz" )
def runGATOnGenomicAnnotations( infiles, outfile ):
    '''run gat agains genomic context.

    The workspace is composed of all mapable regions.
    Enrichment is controlled by isochores.

    To be rigorous, FDR should be re-computed after merging all
    analyses.
    '''

    bedfile, annofile, workspacefile, isochorefile = infiles    

    to_cluster = True
    outdir = "annotations_gat.dir"

    statement = '''gatrun.py
         --segments=%(bedfile)s
         --annotations=<(zcat %(annofile)s | awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\n",$1,$4,$5,$3);}')
         --workspace=%(workspacefile)s
         --isochores=%(isochorefile)s
         --num-samples=%(gat_num_samples)i
         --force
         --ignore-segment-tracks
         --output-filename-pattern=%(outfile)s.%%s
         -v 5
         --log=%(outfile)s.log
         | gzip
         > %(outfile)s'''

    P.run()

@transform( (runGATOnGenomicContext,
             runGATOnGenomicAnnotations ),
            regex("(.*)_gat.dir/(.*).gat.tsv.gz" ),
            r"\1_gat.dir/gat_\1_\2.load" )
def loadGat( infile, outfile ):
    P.load( infile, outfile )

@follows( loadGat )
def gat(): pass


############################################################
############################################################
############################################################
## export section start
############################################################

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
@follows( loadPeakShapeTable )
def annotate_withreads():
    pass

@follows( loadAnnotations, loadTSS, loadRepeats, loadContextStats, loadBinding )
def annotate_intervals(): pass

# @follows( mapping,
#           buildIntervals, 
#           loadReadCoverageTable,
#           buildPeakShapeTable,
#           buildReadProfileOfTranscripts,
#           buildIntervalProfileOfTranscripts )
# def intervals():
#     '''compute binding intervals.'''
#     pass

# @follows( exportBigwig, viewIntervals, viewBigwig )
# def export():
#     '''export files.'''
#     pass

# @follows( buildReferenceMotifs,
#           exportMotifSequences,
#           runMEME,
#           runTomTom, loadTomTom,
#           runBioProspector )
# #          runGLAM2, 
# def discover_motifs():
#     '''run motif discovery.'''
#     pass

# @follows( filterMotifs,
#           exportMotifControlSequences,
#           loadMotifSequenceComposition,
#           loadMotifInformation,
#           runMAST, loadMAST )
# #          runGLAM2SCAN, loadGLAM2SCAN )
# def detect_motifs():
#     '''run motif detection.'''
#     pass

# @follows( loadCorrelation, 
#           loadOverlap,
#           reproducibility)
# def correlation():
#     '''run the correlation targets.'''
#     pass

# @follows( annotateIntervals, loadAnnotations, 
#           annotateTSS, loadTSS, 
#           annotateRepeats, loadRepeats,
#           #annotateTSSIntervalAssociations, loadTSSIntervalAssociations,
#           #annotateTSSIntervalDistance, loadTSSIntervalDistance, 
#           buildIntervalCounts, loadIntervalCounts )
# def annotation():
#     '''run the annotation targets.'''
#     pass

# ###################################################################
# ###################################################################
# ###################################################################
# ## export targets
# ###################################################################
# @merge( intervals,  "view_mapping.load" )
# def createViewMapping( infile, outfile ):
#     '''create view in database for alignment stats.

#     This view aggregates all information on a per-track basis.

#     The table is built from the following tracks:
    
#     bam_stats: .call
#     '''

#     tablename = P.toTable( outfile )

#     # can not create views across multiple database, so use table
#     view_type = "TABLE"
    
#     dbhandle = connect()
#     Database.executewait( dbhandle, "DROP %(view_type)s IF EXISTS %(tablename)s" % locals() )

#     statement = '''
#     CREATE %(view_type)s %(tablename)s AS
#     SELECT SUBSTR( b.track, 1, LENGTH(b.track) - LENGTH( '.genome')) AS track, *
#     FROM bam_stats AS b
#     WHERE b.track LIKE "%%.genome" 
#     ''' % locals()

#     Database.executewait( dbhandle, statement )

#     P.touch( outfile )

# ###################################################################
# ###################################################################
# ###################################################################
# @follows( createViewMapping )
# def views():
#     pass

###################################################################
###################################################################
###################################################################
@follows( annotate_intervals, 
          annotate_withreads,
          runMeme,
          loadMemeSummary,
          loadMast,
          loadMotifInformation,
          loadMotifSequenceComposition )
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
@follows( mkdir( "%s/bedfiles" % PARAMS["web_dir"]), 
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
        "intervals" : glob.glob(os.path.join( PARAMS["exportdir"], "bed", "*.bed.gz" ))+\
            glob.glob(os.path.join( PARAMS["exportdir"], "bed", "*.bed.gz.tbi" )),
        }

    bams = []
    
    for targetdir, filenames in exportfiles.iteritems():
        if len(filenames) == 0:
            E.warn( "no files for target '%s'" % targetdir)
        for src in filenames:
            dest = "%s/%s/%s" % (web_dir, targetdir, os.path.basename(src))
            if dest.endswith( ".bam"): bams.append( dest )
            dest = os.path.abspath( dest )
            destdir = os.path.dirname( dest )
            if not os.path.exists( destdir ):
                os.makedirs( destdir)
                
            if not os.path.exists( dest ):
                E.debug( "creating symlink from %s to %s" % (src, dest))
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


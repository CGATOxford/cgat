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
(see :ref:`PipelineReporting`). To start with, use the files supplied with the
Example_ data.


Input
-----

Intervals
+++++++++

Input are :term:`bed`-formatted files of intervals. Intervals should be at least 
bed4 formatted, i.e., each interval should be labelled (uniquely).

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

import CGAT.Experiment as E
import logging as L
import CGAT.Database as Database
from ruffus import *
import csv
import sqlite3
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGAT.MatrixTools as MatrixTools
import pysam
import numpy
import gzip
import xml.etree.ElementTree

import PipelineChipseq as PipelineChipseq
import PipelineMotifs as PipelineMotifs
import PipelineGeneset as PGeneset
import CGATPipelines.PipelineTracks as PipelineTracks
import CGATPipelines.PipelineMapping as PipelineMapping

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
import CGAT.Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
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

    Default tracks and offsets can be specified using a placeholder ``%``. The
    following will associate all tracks with the same bam file::

        [bams]
        %=all.bam


    '''
    fn = track.asFile()
    bamfiles = glob.glob( "%s.bam" % fn )

    if bamfiles == []:
        if "bams_%s" % fn.lower() in PARAMS:
            for ff in P.asList( PARAMS["bams_%s" % fn.lower() ] ):
                bamfiles.extend( glob.glob( ff ) )
        else:
            for pattern, value in P.CONFIG.items( "bams" ):
                if "%" in pattern:
                    p = re.sub( "%", "\S+", pattern )
                    if re.search( p, fn ):
                        bamfiles.extend( glob.glob( value ) )

    offsets = []
    if "offsets_%s" % fn.lower() in PARAMS:
        offsets = map(int, P.asList( PARAMS["offsets_%s" % fn.lower() ] ))
    else:
        for pattern, value in P.CONFIG.items( "offsets" ):
            if "%" in pattern:
                p = re.sub( "%", "\S+", pattern )
                if re.search( p, fn ):
                    offsets.extend( map( int, value.split(",") ) )

    if offsets == []:
        offsets = [0] * len(bamfiles)

    if len(bamfiles) != len(offsets):
        raise ValueError("number of BAM files %s is not the same as number of offsets: %s" % (str(bamfiles), str(offsets)))


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

    headers = ("avgval","disttostart","genelist","length","peakcenter","peakval","position","interval_id","npeaks","nprobes", "contig","start","end","score" )

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

        # The fifth field of a bed file can be used to supply a score. Our iterator returns 
        # the optional fields as a "fields array". The first of these is the interval name, 
        # and the second the score. The score may be more is better or less is better.
        if len(bed.fields)>1:
            value = bed.fields[1]
            if value != "":
                score = value
            else:
                score = 1
        else:
            score = 1
            
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
                                             bed.contig,bed.start,bed.end,score) )) + "\n" )

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
    statement = '''zcat %(infile)s | sort -k1,1 -k2,2n | bgzip > %(outfile)s; tabix -p bed %(outfile)s'''
    P.run()

############################################################
############################################################
############################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "peaks" ) ) )
@transform( loadIntervals,
            regex(r"(.*)_intervals.load"),
            os.path.join( PARAMS["exportdir"], "peaks", r"\1.peak.bed.gz") )
def exportPeakLocations( infile, outfile ):
    '''export peak locations
    '''

    dbh = connect()
    outf = IOTools.openFile( outfile, "w" )
    cc = dbh.cursor()
    table = P.toTable(infile) 
    for x in cc.execute( """SELECT contig, peakcenter, peakcenter+1, interval_id, peakval 
                                   FROM %(table)s """ % locals() ):
        outf.write( "\t".join( map(str, x) ) + "\n" )
    outf.close()


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

    job_options = "-l mem_free=4G"

    to_cluster = True
    statement = '''
       python %(scriptsdir)s/bam_vs_bed.py
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
def annotateIntervalsFull( infile, outfile ):
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
@transform( exportPeakLocations,
            suffix(".bed.gz"),
            ".annotations" )
def annotateIntervalsPeak( infile, outfile ):
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
def annotateBindingFull( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set.
    
    The reference gene set is 

    Binding is counted for the full intervals.
    '''
    to_cluster = True

    geneset = os.path.join( PARAMS["annotations_dir"],
                            PARAMS_ANNOTATIONS[PARAMS["geneset_binding"]] )

    statement = """
    zcat < %(geneset)s
    | awk '$2 == "protein_coding"'
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
@transform( exportPeakLocations,
            suffix(".bed.gz"),
            ".binding.tsv.gz" )
def annotateBindingPeak( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set.

    Binding is counted for peaks.
    '''
    to_cluster = True

    geneset = os.path.join( PARAMS["annotations_dir"],
                            PARAMS_ANNOTATIONS[PARAMS["geneset_binding"]] )

    statement = """
    zcat < %(geneset)s
    | awk '$2 == "protein_coding"'
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
@transform( exportPeakLocations,
            suffix(".bed.gz"),
            ".tss" )
def annotateTSSPeaks( infile, outfile ):
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
############################################################
############################################################
@transform( TRACKS_BEDFILES,
            suffix(".bed.gz"),
            ".nuc" )
def annotateNucleotides( infile, outfile ):
    '''get the nucleotide composition of the intervals'''

    to_cluster = True

#    statement = '''zcat %s | cut -f1,2,3 | python %s/bed2fasta.py -g %s/%s 
#                   | sed 's/[0-9]*\s\(chr[^:]*\):\([0-9]*\)..\([0-9]*\)\s(+)/\\1|\\2|\\3/g' 
#                   | python %s/fasta2table.py -s na | sed 's/id/contig|start|end/g' | tr '|' '\\t' > %s''' \
#        % (infile,PARAMS["scriptsdir"],PARAMS["genome_dir"],PARAMS["genome"],PARAMS["scriptsdir"],outfile)

    #The bed file is cut to ensure each entry is assigned a unique name from bed2gff - possibly it would be better to validate interval files at the start of the pipeline and assign unique identifiers.
    statement = '''zcat %(infile)s | cut -f1,2,3
                   | python %(scriptsdir)s/bed2gff.py --as-gtf
                   | python %(scriptsdir)s/gtf2table.py --counter=position --counter=composition-na --counter=composition-cpg \
                   --genome-file=%(genome_dir)s/%(genome)s > %(outfile)s
                   '''          
    P.run()


############################################################
@transform( (annotateIntervalsFull, annotateIntervalsPeak), suffix(".annotations"), "_annotations.load" )
def loadAnnotations( infile, outfile ):
    '''load interval annotations: genome architecture
    '''
    P.load( infile, outfile, "--index=gene_id --allow-empty" )

############################################################
@transform( (annotateBindingFull, annotateBindingPeak), suffix(".binding.tsv.gz"), "_binding.load" )
def loadBinding( infile, outfile ):
    '''load interval binding: genome architecture
    '''
    P.load( infile, outfile, "--index=gene_id --allow-empty --map=pattern:str" )

############################################################
@transform( (annotateTSS, annotateTSSPeaks), suffix( ".tss"), "_tss.load" )
def loadTSS( infile, outfile ):
    '''load interval annotations: distance to transcription start sites
    '''
    P.load( infile, outfile, "--index=gene_id --index=closest_id --index=upstream_id --index=downstream_id --allow-empty" )

############################################################
@transform( annotateRepeats, suffix(".repeats"), "_repeats.load" )
def loadRepeats( infile, outfile ):
    '''load interval annotations: repeats
    '''
    P.load( infile, outfile, "--index=gene_id --allow-empty" )

############################################################
@transform( annotateNucleotides, suffix(".nuc"), "_nuc.load" )
def loadNucleotides( infile, outfile ):
    '''load interval annotations: nucleotide composition
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
                      --normalize-profile=none
                      --normalize-profile=area
                      --normalize-profile=counts
                      %(bedfile)s %(gtffile)s
                   > %(outfile)s
                '''
    P.run()


###################################################################
###################################################################
###################################################################
@follows( mkdir( "transcriptprofiles" ) )
@split( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            [r"transcriptprofiles/\1.withoverlap.gtf.gz",r"transcriptprofiles/\1.woutoverlap.gtf.gz",r"\1.tss.withoverlap.gtf.gz",r"\1.tss.woutoverlap.gtf.gz"] )

def prepareGTFsByOverlapWithIntervals( infile, outfiles ):
    '''Prepare GTF file of overlapping and non-overlapping genes for profile plots'''
    
    to_cluster = True

    track = TRACKS.factory( filename = infile[:-len(".bed.gz")] )
    track = Sample( filename = P.snip( infile, ".bed.gz") )

    out1, out2, out3, out4 = outfiles
    geneset = PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]

    statement = '''
                  intersectBed -u -a %(annotations_dir)s/%(geneset)s -b %(track)s.bed.gz  | 
                  python %(scriptsdir)s/gff2bed.py --is-gtf -v 0 | cut -f4 | sort | uniq > %(track)s_overlapping_genes; 
                  checkpoint;
                  zgrep -f %(track)s_overlapping_genes %(annotations_dir)s/%(geneset)s | gzip > %(out1)s; checkpoint;
                  zgrep -v -f %(track)s_overlapping_genes %(annotations_dir)s/%(geneset)s | gzip > %(out2)s; checkpoint;
                  zgrep -f %(track)s_overlapping_genes %(annotations_dir)s/tss.gene.gtf | gzip > %(out3)s; checkpoint;
                  zgrep -v -f %(track)s_overlapping_genes %(annotations_dir)s/tss.gene.gtf | gzip > %(out4)s; checkpoint;
                '''
    P.run()

############################################################
############################################################
############################################################
@transform( prepareGTFsByOverlapWithIntervals,
            regex("(.*tss.*).gtf.gz"),
            r"\1.nuc" )
def annotateTSSNucleotides( infile, outfile ):
    '''get the nucleotide composition of the intervals'''
    to_cluster = True
    #statement = '''zcat %s | awk '{OFS="\\t"}{print $1,$4-50,$5+50}' 
    #               | python %s/bed2fasta.py -g %s/%s 
    #               | sed 's/[0-9]*\s\(chr[^:]*\):\([0-9]*\)..\([0-9]*\)\s(+)/\\1|\\2|\\3/g' 
    #               | python %s/fasta2table.py -s na | sed 's/id/contig|start|end/g' 
    #               | tr '|' '\\t' > %s''' % (infile,PARAMS["scriptsdir"],PARAMS["genome_dir"],PARAMS["genome"],PARAMS["scriptsdir"],outfile)

    statement = '''zcat %(infile)s | slopBed -b 50 -g %(genome_dir)s/%(genome)s.fasta.fai
                   | python %(scriptsdir)s/gtf2table.py --counter=position --counter=composition-na --counter=composition-cpg \
                   --genome-file=%(genome_dir)s/%(genome)s > %(outfile)s
                   '''          
    P.run()

############################################################
@transform( annotateTSSNucleotides, suffix(".nuc"), "_nuc.load" )
def loadTSSNucleotides( infile, outfile ):
    '''load interval annotations: nucleotide composition
    '''
    P.load( infile, outfile, "--index=gene_id --allow-empty" )

###################################################################
###################################################################
###################################################################
@transform( prepareGTFsByOverlapWithIntervals,
            regex("(transcriptprofiles.*).gtf.gz"),
            r"\1.tsv.gz" )

def buildGenesByIntervalsProfiles( infile, outfile ):
    '''Make gene profile plots.'''
    
    to_cluster = True

    track = TRACKS.factory( filename = infile[len("transcriptprofiles/"):-len(".withoverlap.tsv.gz")] )   

    bamfiles, offsets = getAssociatedBAMFiles( track )
    geneset = PARAMS_ANNOTATIONS["interface_geneset_all_gtf"]

    if bamfiles:
        E.info( "%s: associated bamfiles = %s" % (track, bamfiles))
    else:
        E.warn( "%s: no bamfiles associated - target skipped" % (track))
        P.touch( outfile )
        P.touch( outfile[:-len(".tsv.gz")]+".geneprofile.counts.tsv.gz" )
        return

    if len(bamfiles) > 1:
        raise NotImplementedError( "peakshape with multiple bamfiles not implement" )
    bamfile=bamfiles[0]

    outpat = outfile[:-len(".tsv.gz")]
    statement = '''

    python %(scriptsdir)s/bam2geneprofile.py
                      --output-filename-pattern="%(outpat)s.%%s"
                      --force
                      --reporter=transcript
                      --method=geneprofile 
                      --method=tssprofile 
                      --normalize-profile=none
                      --normalize-profile=area
                      --normalize-profile=counts
                      %(bamfile)s <(zcat %(infile)s)
                   > %(outfile)s ;
                '''
    P.run()

############################################################
@transform( buildGenesByIntervalsProfiles,
            suffix(".tsv.gz"),
            r"\1.geneprofile.counts.load" )

def loadByIntervalProfiles( infile, outfile ):
    '''load interval annotations: nucleotide composition
    '''
    countsfile = infile[:-len(".tsv.gz")]+".geneprofile.counts.tsv.gz"
    P.load( countsfile, outfile, "--index=gene_id --allow-empty" )

    
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
        
    p = P.substituteParameters( **locals() )
    nseq = PipelineMotifs.writeSequencesForIntervals( track, 
                                                      outfile,
                                                      dbhandle,
                                                      full = False,
                                                      masker = P.asList(p['motifs_masker']),
                                                      halfwidth = int(p["motifs_halfwidth"]),
                                                      maxsize = int(p["motifs_max_size"]),
                                                      proportion = p["motifs_proportion"],
                                                      min_sequences = p["motifs_min_sequences"] )

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
## run against database of known motifs
############################################################
@transform( runMeme, suffix(".meme"), ".tomtom" )
def runTomTom( infile, outfile ):
    '''compare ab-initio motifs against tomtom.'''
    PipelineMotifs.runTomTom( infile, outfile )

@transform( runTomTom, suffix(".tomtom"), "_tomtom.load" )
def loadTomTom( infile, outfile ):
    '''load tomtom results'''

    tablename = P.toTable( outfile )

    resultsdir = os.path.join( os.path.abspath(PARAMS["exportdir"]), "tomtom", infile )
    xml_file = os.path.join( resultsdir, "tomtom.xml" ) 

    if not os.path.exists( xml_file ):
        E.warn( "no tomtom output - skipped loading " )
        P.touch( outfile )
        return

    # get the motif name from the xml file

    tree = xml.etree.ElementTree.ElementTree()
    tree.parse( xml_file )
    motifs =  tree.find( "targets" )
    name2alt = {}
    for motif in motifs.getiterator( "motif" ):
        name = motif.get("name")
        alt = motif.get("alt")
        name2alt[name] = alt

    tmpfile = P.getTempFile(".")
    
    # parse the text file
    for line in IOTools.openFile( infile ):
        if line.startswith( "#Query"):
            tmpfile.write( "target_name\tquery_id\ttarget_id\toptimal_offset\tpvalue\tevalue\tqvalue\tOverlap\tquery_consensus\ttarget_consensus\torientation\n" )
            continue
        data = line[:-1].split("\t" )
        target_name = name2alt[data[1]]
        tmpfile.write( "%s\t%s" % (target_name, line ) )
    tmpfile.close()
    
    P.load( tmpfile.name, outfile )
    
    os.unlink( tmpfile.name )

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
@follows( loadMotifInformation, mkdir( os.path.join( PARAMS["exportdir"], "motifs" ) ) )
@merge( loadMast, "motifs.export" )
def exportMotifLocations( infiles, outfile ):
    '''export motif locations. There will be a bed-file per motif.

    Overlapping motif matches in different tracks will be merged.
    '''

    dbh = connect()
    cc = dbh.cursor()

    motifs = [ x[0] for x in cc.execute( "SELECT motif FROM motif_info" ).fetchall()]

    
    for motif in motifs:

        tmpf = P.getTempFile()
        
        for infile in infiles:
            table = P.toTable(infile) 
            track = P.snip( table, "_mast" )
            for x in cc.execute( """SELECT contig, start, end, '%(track)s', evalue
                                   FROM %(table)s WHERE motif = '%(motif)s' AND start IS NOT NULL""" % locals() ):
                tmpf.write( "\t".join( map(str, x) ) + "\n" )
        tmpf.close()

        outfile = os.path.join( PARAMS["exportdir"], "motifs", "%s.bed.gz" % motif )
        tmpfname = tmpf.name 

        statement = '''mergeBed -i %(tmpfname)s -nms | gzip > %(outfile)s'''
        P.run()

        os.unlink( tmpf.name )

############################################################
############################################################
############################################################
## compute overlap with genomic features
############################################################
@follows( mkdir("gat_context.dir") )
@transform( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_genomic_context_bed"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_mapability_filtered_bed"] % PARAMS["gat_mapability"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_gc_profile_bed"] ),
                        ),
            r"gat_context.dir/\1.gat.tsv.gz" )
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

    statement = '''gat-run.py
         --segments=%(bedfile)s
         --annotations=%(annofile)s
         --workspace=%(workspacefile)s
         --num-samples=%(gat_num_samples)i
         --isochores=%(isochorefile)s
         --force
         --ignore-segment-tracks
         --output-filename-pattern=%(outfile)s.%%s
         --output-counts-pattern=%(outfile)s.%%s.counts.gz
         -v 5
         --log=%(outfile)s.log
         | gzip
         > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
## compute overlap with genomic annotations
############################################################
@follows( mkdir("gat_annotations.dir") )
@transform( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_annotation_gff"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_mapability_filtered_bed"] % PARAMS["gat_mapability"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_gc_profile_bed"] ), 
                        ),
            r"gat_annotations.dir/\1.gat.tsv.gz" )
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

    statement = '''gat-run.py
         --segments=%(bedfile)s
         --annotations=<(zcat %(annofile)s | awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\n",$1,$4,$5,$3);}')
         --workspace=%(workspacefile)s
         --isochores=%(isochorefile)s
         --num-samples=%(gat_num_samples)i
         --force
         --ignore-segment-tracks
         --output-filename-pattern=%(outfile)s.%%s
         --output-counts-pattern=%(outfile)s.%%s.counts.gz
         -v 5
         --log=%(outfile)s.log
         | gzip
         > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
## compute overlap with gene structure
############################################################
@follows( mkdir("gat_genestructure.dir") )
@transform( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_genestructure_gff"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_mapability_filtered_bed"] % PARAMS["gat_mapability"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_gc_profile_bed"] ), 
                        ),
            r"gat_genestructure.dir/\1.gat.tsv.gz" )
def runGATOnGeneStructure( infiles, outfile ):
    '''run gat on gene structures

    The purpose of this gat run is to test for differential location
    of intervals in parts of certain gene structures.

    The workspace is composed of all mapable regions. The workspace is 
    restricted to annotations in order to reduce the effect of
    intergenic depletion. Furthermore, the workspace is restricted
    to those parts that contain segments and annotations in order
    to avoid a gene bias (only genes of a certain structure contain
    segments).

    Enrichment is controlled by isochores.

    To be rigorous, FDR should be re-computed after merging all
    analyses.
    '''

    bedfile, annofile, workspacefile, isochorefile = infiles    

    to_cluster = True
    outdir = "gat_genestructure.dir"

    statement = '''gat-run.py
         --segments=%(bedfile)s
         --annotations=<(zcat %(annofile)s | awk '{printf("%%s\\t%%i\\t%%i\\t%%s\\n",$1,$4,$5,$3);}')
         --workspace=%(workspacefile)s
         --isochores=%(isochorefile)s
         --num-samples=%(gat_num_samples)i
         --counter=segment-midoverlap
         --truncate-workspace-to-annotations
         --restrict-workspace
         --force
         --ignore-segment-tracks
         --output-filename-pattern=%(outfile)s.%%s
         --output-counts-pattern=%(outfile)s.%%s.counts.gz
         -v 5
         --log=%(outfile)s.log
         | gzip
         > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
## compute overlap with gene annotations
############################################################
@follows( mkdir("gat_functions.dir") )
@transform( TRACKS_BEDFILES,
            regex("(.*).bed.gz"),
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_genomic_function_bed"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_genomic_function_tsv"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_territories_gff"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_mapability_filtered_bed"] % PARAMS["gat_mapability"] ),
                        os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_gc_profile_bed"] ), 
                        ),
            r"gat_functions.dir/\1.gat.tsv.gz" )
def runGATOnGeneAnnotations( infiles, outfile ):
    '''run gat against genes and their annotations.

    The workspace is composed of all mapable regions within gene territories.
    Enrichment is controlled by isochores.

    To be rigorous, FDR should be re-computed after merging all
    analyses.
    '''

    # requires a large amount of memory
    job_options = "-l mem_free=20G"
    
    bedfile, annofile, descriptionfile, workspacefile, mapabilityfile, isochorefile = infiles    

    to_cluster = True

    statement = '''gat-run.py
         --segments=%(bedfile)s
         --annotations=%(annofile)s
         --workspace=<(zcat %(workspacefile)s | awk '{printf("%%s\\t%%i\\t%%i\\n",$1,$4,$5);}')
         --workspace=%(mapabilityfile)s
         --isochores=%(isochorefile)s
         --num-samples=%(gat_num_samples)i
         --descriptions=%(descriptionfile)s
         --force
         --ignore-segment-tracks
         --output-filename-pattern=%(outfile)s.%%s
         --output-counts-pattern=%(outfile)s.%%s.counts.gz
         -v 5
         --log=%(outfile)s.log
         | gzip
         > %(outfile)s'''

    P.run()

@transform( (runGATOnGenomicContext,
              runGATOnGenomicAnnotations,
              runGATOnGeneAnnotations,
              runGATOnGeneStructure),
            regex("gat_(.*).dir/(.*).gat.tsv.gz" ),
            r"gat_\1.dir/gat_\1_\2.load" )
def loadGat( infile, outfile ):
    '''load individual gat results.'''

    P.load( infile, outfile )

@collate( (runGATOnGenomicContext,
           runGATOnGenomicAnnotations,
           runGATOnGeneAnnotations),
          regex("gat_(.*).dir/.*.gat.tsv.gz" ),
          r"gat_\1.summary.tsv.gz" )
def summarizeGAT( infiles, outfile ):
    '''summarize GAT results.

    outputs a log2fold and pvalue table for results.
    
    The results are filtered. Remove all rows that 
       * have no significant results (fdr)
       * expected overlap less than 1kb.
    '''

    col_headers = [ P.snip( os.path.basename(x), ".gat.tsv.gz") for x in infiles]

    # get qvalues
    qval_matrix, qval_row_headers = MatrixTools.buildMatrixFromTables( infiles, 
                                                                       column = "qvalue",
                                                                       column_header = "annotation",
                                                                       default = 1.0 )
    # output qvalues
    column = "qvalue"
    with IOTools.openFile( outfile + "." + column + ".gz", "w") as outf:
        IOTools.writeMatrix( outf, qval_matrix, qval_row_headers, col_headers )

    ncols = len(infiles)
    min_qvalue = PARAMS["gat_fdr"]
    min_expected = PARAMS["gat_min_expected"]

    E.info( "read matrix with %i rows and %i columns" % qval_matrix.shape )

    # set all values to 1.0 that are above fdr threshold
    qval_matrix[numpy.where( qval_matrix > min_qvalue )] = 1.0
    s = numpy.sum( qval_matrix, 1 )
    take_qvalue = numpy.where( s != ncols )[0]
    E.info( "taking %i rows after qvalue filtering" % len(take_qvalue))

    # remove small results
    expected_matrix, expected_row_headers = MatrixTools.buildMatrixFromTables( infiles, 
                                                                               column = "expected",
                                                                               column_header = "annotation",
                                                                               default = 1.0 )

    expected_matrix[numpy.where( expected_matrix < min_expected )] = 0
    s = numpy.sum( expected_matrix, 1 )
    take_expected = numpy.where( s != 0 )[0]
    E.info( "taking %i rows after min_expected < % i filtering" % ( len(take_expected), min_expected))

    take = sorted( list( set(take_qvalue).intersection( set(take_expected) )) )

    E.info( "taking %i rows after filtering" % ( len(take)))

    for (column, default) in ( ("fold", 1.0), ( "pvalue", 1.0) ):
        # output single tables for fold and pvalue
        matrix, row_headers = MatrixTools.buildMatrixFromTables( infiles, 
                                                                 column = column,
                                                                 column_header = "annotation",
                                                                 default = default )

        assert qval_row_headers == row_headers

        # set those items with qvalue > fdr to default
        matrix[numpy.where( qval_matrix > min_qvalue )] = default
        
        # remove all rows that have no significant results
        row_headers = numpy.take( row_headers, take )
        matrix = numpy.take( matrix, take, 0 )

        if column == "fold": 
            matrix[numpy.where( matrix < 0.00001 )] = 0.00001
            matrix = numpy.log2( matrix )

        col_headers = [ P.snip( os.path.basename(x), ".gat.tsv.gz") for x in infiles]
        with IOTools.openFile( outfile + "." + column + ".gz", "w") as outf:
            IOTools.writeMatrix( outf, matrix, row_headers, col_headers )

    P.touch( outfile )
    
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
@follows( loadPeakShapeTable,
          buildIntervalProfileOfTranscripts )
def annotate_withreads():
    pass

@follows( loadAnnotations, loadTSS, loadRepeats, loadContextStats, loadBinding, loadNucleotides, loadTSSNucleotides )
def annotate_intervals(): pass

# @follows( mapping,
#           buildIntervals, 
#           loadReadCoverageTable,
#           buildReadProfileOfTranscripts)
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
          loadByIntervalProfiles,
          runMeme,
          loadMemeSummary,
          loadTomTom,
          loadMotifInformation,
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


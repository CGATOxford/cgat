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

:Author: Andreas Heger
:Release: $Id: pipeline_chipseq.py 2900 2010-04-13 14:38:00Z andreas $
:Date: |today|
:Tags: Python

=================
ChIP-Seq pipeline
=================

The ChIP-Seq pipeline imports reads from one or more ChIP-Seq experiments and
performs the following tasks:

   * align reads to the genome
   * call peaks using MACS (TODO: at QuEST)
   * annotate intervals with respect to a reference gene set
   * describe de-novo motifs using MEME and BioProspector
   * find canonical motifs using MAST

Maybe:
   * perform enrichment analysis against UCSC tracks
   * perform GO analysis of peaks close to protein coding genes.

Configuration
=============

Input
=====

Reads
-----

Input are :file:`_export.txt.gz`-formatted files from Illumina. The files should be
labeled in the following way::

   sample-condition-replicate_export.txt.gz

For example::

   GM00855-D3-R1_export.txt.gz
   GM00855-D3-R2_export.gz
   GM00855-input-R1_export.gz
   GM00855-unstim-R1_export.txt.gz
   GM00855-unstim-R2_export.txt.gz

Note that neither ``sample``, ``condition`` or ``replicate`` should contain 
``_`` (underscore) and ``.`` (dot) characters as these are used by the pipeline.

Output
======

The results of the computation are all stored in an sqlite relational
database file.

Usage
=====

Type::

   python <script_name>.py --help

for command line help.

Code
----

TODO: currently the bed files and the intervals are inconsistent 
    (due to filtering, there are more intervals in the bed files than
     in the table. The ids do correspond).

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections

import Experiment as E
import Pipeline as P
from ruffus import *
import csv
import sqlite3
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import IOTools
import MAST, GTF, GFF, Bed, Stats
import cStringIO
import pysam
import numpy
import gzip
import Masker
import Glam2Scan
import fileinput
import Motifs
import gff2annotator
import Bioprospector

import pipeline_chipseq_intervals as PIntervals

import pipeline_vitaminD_annotator as PAnnotator

import pipeline_vitaminD_motifs as PMotifs
import PipelineGeneset as PGeneset

BAM_SUFFIX = ".norm.bam"

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################
P.getParameters( 
    ["%s.ini" % __file__[:-len(".py")],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
def splitTrack( track ):
    '''split track into tuple tissue, condition, replicate.'''

    d = track.split("-")
    tissue, condition, replicate = None, None, None
    if len(d) == 3: tissue, condition, replicate = d
    elif len(d) == 2: tissue, condition = d
    elif len(d) == 1: tissue = d[0]
    return tissue, condition, replicate

def getControl( track ):
    '''return appropriate control for a track
    
    Assumption: a control has no replicate.
    '''
    tissue, condition, replicate = splitTrack(track)
    return "-".join( (tissue, PARAMS["tracks_control"] ) )

def isReplicate( track ):
    '''return true if track is a replicate.'''
    return splitTrack(track)[2] != None

## all tracks
TRACKS_ALL=[ x[:-len( "_export.txt.gz")] for x in glob.glob( "*_export.txt.gz" ) if PARAMS["tracks_control"] not in x]


## all experiments (combinations of sample + condition)
TRACKS_EXPERIMENTS = list(set([ "-".join( splitTrack(x)[0:2]) for x in TRACKS_ALL ] ) )


## aggregate by condition
TRACKS_CONDITIONS = list(set([ "-".join( ("agg", splitTrack(x)[1])) for x in TRACKS_ALL ] ) )

## aggregate by tissue
TRACKS_TISSUES = list(set([ "-".join( (splitTrack(x)[0], "agg")) for x in TRACKS_ALL ] ) )

## todo: discover these automatically
REPLICATES=PARAMS["tracks_replicates"]
CONDITIONS=PARAMS["tracks_conditions"]
TISSUES=PARAMS["tracks_tissues"]

TRACKS_NOT_UNSTIM = [ "%s" % x for x in TRACKS_CONDITIONS + TRACKS_EXPERIMENTS if not re.search("Unstim",x)]
TRACKS_SUBTRACT = [ "%sSub" % x for x in TRACKS_NOT_UNSTIM ]
TRACKS_MASTER = TRACKS_EXPERIMENTS + TRACKS_CONDITIONS + TRACKS_SUBTRACT
TRACKS_CORRELATION = TRACKS_MASTER + TRACKS_ALL

RUNS=[ "%02i" % x for x in range(1,17) ]

# tracks with reference ChIPSeq Intervals
TRACKS_REFERENCE = ["refGR"]

# tracks to use for running annotator
TRACKS_ANNOTATE = TRACKS_MASTER + TRACKS_REFERENCE

# tracks to use for running annotator
TRACKS_ANNOTATOR = TRACKS_MASTER + TRACKS_REFERENCE

# tracks to submit for motif analysi
TRACKS_MOTIFS = TRACKS_MASTER + TRACKS_REFERENCE
 
TRACKS_ROI = ( "pritchard_ASIA", "pritchard_CEU", "pritchard_YRI" )

TRACKS_REGIONS = ["roi", "gwas", "gwas_merged", "selection" ]

FILEPATTERN="run%(tissue)s%(condition)s%(replicate)s"

###################################################################
###################################################################

if os.path.exists("conf.py"):
    E.info( "reading additional configuration from conf.py" )
    execfile("conf.py")

###################################################################
###################################################################
###################################################################
## General preparation tasks
###################################################################

############################################################
############################################################
############################################################
## 
############################################################
@files( "genome.fasta", "genome_gc.bed" )
def buildGenomeGCSegmentation( infile, outfile ):
    '''segment the genome into windows according to G+C content.'''

    to_cluster = True
    
    old_statement = '''
    python %(scriptsdir)s/fasta2bed.py \
        --method=GCProfile \
        --segmentation-threshold=300 \
        --gap-size=10000 \
        --log=%(outfile)s.log \
    < %(genome)s.fasta > %(outfile)s'''

    statement = '''
    python %(scriptsdir)s/fasta2bed.py \
        --method=fixed-width-windows --window-size=1000 \
        --log=%(outfile)s.log \
    < %(genome)s.fasta > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
## 
############################################################
@files( buildGenomeGCSegmentation, "annotator_gc.bed" )
def buildAnnotatorGC( infile, outfile ):
    '''compute G+C regions.'''

    to_cluster = True
    statement = '''
    python %(scriptsdir)s/bed2bed.py \
        --method=bins \
        --num-bins=%(annotator_gc_bins)s \
        --binning-method=%(annotator_gc_method)s \
        --log=%(outfile)s.log \
    < %(infile)s > %(outfile)s'''

    P.run()

###################################################################
###################################################################
###################################################################
## Version 1: import from given bed files
###################################################################
if PARAMS["mapping_mapper"] == "bowtie":
    
    ############################################################
    ############################################################
    ############################################################
    @transform( "*_export.txt.gz", 
                regex("(\S+)_export.txt.gz"), 
                r"\1.bam" )
    def buildBAM( infile, outfile ):
        '''re-map eland formatted reads with bowtie
        '''
        to_cluster = True

        # require 4Gb of free memory
        job_options = "-l mem_free=4000M"

        job_options += " -q server_jobs.q"

        tmpfilename = P.getTempFilename()

        prefix = outfile[:-len(".bam")]

        statement = '''
        gunzip < %(infile)s |\
        awk '$11 != "QC" || $10 ~ /(\d+):(\d+):(\d+)/ \
             { if ($1 != "") { readname=sprintf( "%%s_%%s:%%s:%%s:%%s:%%s", $1,$2,$3,$4,$5,$6);}
              else { readname=sprintf( "%%s:%%s:%%s:%%s:%%s", $1,$3,$4,$5,$6); }
              printf("@%%s\\n%%s\\n+\\n%%s\\n",readname,$9,$10);}' |\
        bowtie --sam %(bowtie_options)s %(bowtie_index)s - 2>%(outfile)s.log |\
        samtools import %(genome)s - %(tmpfilename)s >& %(outfile)s.log;
        samtools sort %(tmpfilename)s %(prefix)s;
        samtools index %(outfile)s;
        rm -f %(tmpfilename)s
        '''

        P.run()

        if os.path.exists( tmpfilename ):
            os.unlink( tmpfilename )
else:
    raise ValueError("unknown mapper %s" % PARAMS["mapping_mapper"] )

###################################################################
###################################################################
###################################################################
## process reads
###################################################################

############################################################
############################################################
############################################################
@transform( buildBAM,
            suffix(".bam"),
            ".readstats" )
def buildBAMStats( infile, outfile ):
    '''count number of reads mapped, duplicates, ....'''
    PIntervals.buildBAMStats( infile, outfile )

############################################################
############################################################
############################################################
@transform( buildBAMStats,
            regex(r"(.*).readstats"),
            inputs( (r"\1.bam", r"\1.readstats") ),
            r"\1.norm.bam" )
def normalizeBAMPerReplicate( infiles, outfile ):
    '''build a normalized BAM file such that all
    files have approximately the same number of 
    reads.

    Duplicated reads are removed at the same time.
    '''
    PIntervals.buildNormalizedBAM( (infiles,), 
                                   outfile,
                                   PARAMS["calling_normalize"])

if PARAMS["calling_caller"] == "macs":

    ############################################################
    ############################################################
    ############################################################
    @follows( normalizeBAMPerReplicate )
    @files( [ (("%s.norm.bam" % x, "%s.norm.bam" % getControl(x)), "%s.macs" % x ) for x in TRACKS_ALL ] )
    def runMACS( infiles, outfile ):
        '''run MACS for peak detection.'''
        infile, controlfile = infiles

        statement = '''
            macs14 -t %(infile)s 
                   -c %(controlfile)s
                   --diag
                   --name=%(outfile)s
                   --format=BAM
                   %(macs_options)s 
            >& %(outfile)s''' 
    
        P.run() 
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runMACS,
                regex(r"(.*).macs"),
                inputs( (r"\1.macs", r"\1.norm.bam") ),
                r"\1_macs.import" )
    def importMACS( infiles, outfile ):
        infile, bamfile = infiles
        PIntervals.importMACS( infile, outfile, bamfile )
        
    ############################################################
    ############################################################
    ############################################################
    @merge( runMACS, "macs.summary" )
    def summarizeMACS( infiles, outfile ):
        '''run MACS for peak detection.'''

        PIntervals.summarizeMACS( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( summarizeMACS,
                suffix(".summary"),
                "_summary.import" )
    def importMACSSummary( infile, outfile ):
        '''import macs summary.'''
        PIntervals.importMACSSummary( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( importMACS, suffix("_macs.import"), ".bed" )
    def exportIntervalsAsBed( infile, outfile ):
        PIntervals.exportIntervalsAsBed( infile, outfile )

else:
    raise ValueError("unknown peak caller %s" % PARAMS["calling_caller"] )

############################################################
############################################################
############################################################
@follows( exportIntervalsAsBed )
@files( [ ( ["%s%s.bed" % (x,y) for y in REPLICATES],
            "%s.bed" % x, 
            x) for x in TRACKS_EXPERIMENTS ] )
def combineReplicates( infile, outfile, track ):
    '''combine replicates between experiments.

    The replicates are combined using intersection.
    '''

    replicates = [ "%s%s.bed" % (track, x) for x in REPLICATES  if os.path.exists( "%s%s.bed" % (track, x) ) ]
    PIntervals.intersectBedFiles( replicates, outfile )

############################################################
############################################################
############################################################
@follows( exportIntervalsAsBed, combineReplicates )
@files( [ ( ["run%s%s.bed" % (y,x) for y in TISSUES],
            "run%s.bed" % x, 
            x) for x in CONDITIONS ] )
def combineConditions( infile, outfile, track ):
    '''combine conditions between cell lines. 

    The conditions are merged via intersection.
    '''

    conditions = [ "run%s%s.bed" % (x,track) for x in TISSUES ]
    PIntervals.intersectBedFiles( conditions, outfile )

############################################################
############################################################
############################################################
@follows( combineConditions )
@files( [ ( ("%s.bed" % x, "runUnstim.bed"),
            "%sSub.bed" % x, 
            x) for x in TRACKS_EXPERIMENTS + TRACKS_CONDITIONS if not re.search("Unstim",x) ] )
def combineUnstim( infiles, outfile, track ):
    '''remove the unstimulated data sets from the individual tracks. 
    '''

    infile, subtract = infiles
    PIntervals.subtractBedFiles( infile, subtract, outfile )

############################################################
############################################################
############################################################
@transform( (combineReplicates,
             combineConditions, 
             combineUnstim, 
             normalizeBAMPerReplicate),
            suffix(".bed"),
            "_bed.import" )
def importCombinedIntervals( infiles, outfile ):
    PIntervals.importCombinedIntervals( infiles, outfile )

############################################################
############################################################
############################################################
## master target for this section
############################################################
@follows( importCombinedIntervals, importMACSSummary )
def buildIntervals():
    pass

###################################################################
###################################################################
###################################################################
## general import
###################################################################

############################################################
############################################################
############################################################
@follows( exportIntervalsAsBed )
@files_re( ["%s.bed" % x for x in TRACKS_ALL],
          combine( "(.*).bed" ),
          "merged.bed" )
def makeMerged( infiles, outfile ):
    '''combine all experiments.

    The replicates are combined using a merge.
    '''
    PIntervals.mergeBedFiles( infiles, outfile )

        

############################################################
############################################################
############################################################
@files( [ ("%s.idx" % PARAMS["genome"], "genome.gff"), ] )
def buildGenome( infile, outfile ):

    fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )

    entry = GFF.Entry()
    entry.start = 0
    entry.feature = "contig"
    entry.source = "genome"
    outs = open( outfile, "w" )

    for contig, size in fasta.getContigSizes( with_synonyms = False ).iteritems():
        entry.contig = contig
        entry.end = int(size)
        outs.write( "%s\n" % str(entry) )
    outs.close()


############################################################
############################################################
############################################################
@files_re( (buildBAM, normalizeBAMPerReplicate),
           "(.*).bam",
           r"\1.bigwig" )
def buildBigwig( infile, outfile ):
    '''convert BAM to bigwig file.'''

    # no bedToBigBed on the 32 bit cluster
    to_cluster = True
    job_queue = "server_jobs.q"

    statement = '''python %(scriptsdir)s/bam2wiggle.py \
                --genome-file=%(genome)s \
                --output-format=bigwig \
                --output-filename=%(outfile)s \
                %(infile)s \
                > %(outfile)s.log'''
     
    P.run()

############################################################
############################################################
############################################################
# fix: does not work due to exportIntervalsAsBed
# (importCombinedIntervals, exportIntervalsAsBed ),
@transform( [ "%s.bed" % x for x in TRACKS_MASTER ],
            suffix(".bed"),
            ".fasta" )
def exportMotifSequences( infile, outfile ):
    '''export sequences for all intervals.
    '''
    
    track = infile[:-len(".bed")]
    dbhandle = sqlite3.connect( PARAMS["database"] )

    fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )

    cc = dbhandle.cursor()
    statement = "SELECT interval_id, contig, start, end FROM %s_intervals" % track
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
# (importCombinedIntervals, exportIntervalsAsBed ),
# TODO: fix, causes a problem due to exportIntervalsAsBed
@transform( [ "%s.bed" % x for x in TRACKS_MASTER],
            suffix(".bed"),
            ".controlfasta" )
def exportMotifControlSequences( infile, outfile ):
    '''for each interval, export the left and right 
    sequence segment of the same size.
    '''
    
    track = infile[:-len(".bed")]
    dbhandle = sqlite3.connect( PARAMS["database"] )

    fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )

    cc = dbhandle.cursor()
    statement = "SELECT interval_id, contig, start, end FROM %s_intervals" % track
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

if 0:
    ############################################################
    ############################################################
    ############################################################
    @merge( (PARAMS["overlap"], exportIntervalsAsBed ),
            "intervals.overlap" )
    def buildOverlap( infile, outfile ):
        '''compute overlap between intervals.
        '''

        if os.path.exists(outfile): 
            os.rename( outfile, outfile + ".orig" )
            options = "--update=%s.orig" % outfile
        else:
            options = ""

        statement = '''
            python %(scriptsdir)s/diff_bed.py %(options)s *.bed > %(outfile)s
            '''

        P.run()

    ############################################################
    ############################################################
    ############################################################
    @transform( buildOverlap, suffix(".overlap"), "_overlap.import" )
    def importOverlap( infile, outfile ):
        '''import overlap results.
        '''

        tablename = "overlap"

        statement = '''
        csv2db.py %(csv2db_options)s 
                  --index=set1 
                  --index=set2 
                  --table=%(tablename)s 
        < %(infile)s > %(outfile)s
        '''

        P.run()

    ############################################################
    ############################################################
    ############################################################
    @merge( (exportUCSCEncodeTracks, 
             combineConditions,
             combineUnstim,
             combineReplicates ),
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

        P.run( **dict( locals().items() + PARAMS.items() ) )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildUCSCOverlap, suffix(".overlap"), "_overlap.import" )
    def importUCSCOverlap( infile, outfile ):
        '''import overlap results.
        '''

        tablename = "ucsc_overlap"

        statement = '''
        csv2db.py %(csv2db_options)s \
                  --index=set1 \
                  --index=set2 \
                  --table=%(tablename)s \
        < %(infile)s > %(outfile)s
        '''

        P.run( **dict( locals().items() + PARAMS.items() ) )


############################################################
############################################################
############################################################
## targets to do with the analysis of replicates
############################################################
if REPLICATES:

    @files( [( ["%s%s.bed" % (x,y) for y in REPLICATES], "%s.reproducibility" % x, x) for x in TRACKS_EXPERIMENTS ] )
    def makeReproducibility( infile, outfile, track ):
        '''compute overlap between intervals.

        Note the exon percentages are approxmations assuming that there are
        not more than one intervals in one set overlapping one in the other set.
        '''

        dbhandle = sqlite3.connect( PARAMS["database"] )

        data = []
        for replicate in REPLICATES:
            cc = dbhandle.cursor()
            statement = "SELECT contig, start, end, peakval FROM %(track)s%(replicate)s_intervals" % locals()
            cc.execute( statement )
            data.append( cc.fetchall() )

        ma = int(max( [ x[3] for x in itertools.chain( *data ) ] ) + 1)

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

    @files_re( makeReproducibility, "(.*).reproducibility", r"\1_reproducibility.import" )
    def importReproducibility( infile, outfile ):
        '''import Reproducibility results
        '''

        tablename = outfile[:-len(".import")] 

        statement = '''
        csv2db.py %(csv2db_options)s \
                  --table=%(tablename)s \
        < %(infile)s > %(outfile)s
        '''

        P.run( **dict( locals().items() + PARAMS.items() ) )

    @follows( importReproducibility )
    def reproducibility(): pass

else:

    @follows( buildGeneSet )
    def reproducibility(): pass

############################################################
############################################################
############################################################
@files( [( ["%s%s.norm.bam" % (x,y) for y in REPLICATES], "%s.readcorrelations.gz" % x) for x in TRACKS_EXPERIMENTS ] )
def makeReadCorrelation( infiles, outfile ):
    '''compute correlation between reads.
    '''

    to_cluster = True
    job_options = "-l mem_free=8000M"

    infiles = " ".join(infiles)

    statement = '''
    python %(scriptsdir)s/bam_correlation.py \
    --log=%(outfile)s.log \
    --genome=%(genome)s \
    %(infiles)s | gzip > %(outfile)s
    ''' 
    
    P.run( **dict( locals().items() + PARAMS.items() ) )

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
            
            for x in range(len(data)):
                if data[x] < correlation_cutoff: break
            else:
                yield data

    for infile in infiles:
        mat = numpy.array( [ x for x in selector( gzip.open(infile, "r") )] )
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
##
############################################################
@follows( makeMerged )
@files_re( ["%s.bed" % x for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "peakval.correlation" )
def makePeakvalCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PIntervals.makeIntervalCorrelation( infiles, outfile, "peakval" )

@follows( makeMerged )
@files_re( ["%s.bed" % x for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "avgval.correlation" )
def makeAvgvalCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PIntervals.makeIntervalCorrelation( infiles, outfile, "avgval" )

@follows( makeMerged )
@files_re( ["%s.bed" % x for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "length.correlation" )
def makeLengthCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PIntervals.makeIntervalCorrelation( infiles, outfile, "length" )

@transform( ( makePeakvalCorrelation, makeAvgvalCorrelation, makeLengthCorrelation ),
            suffix(".correlation"),
            "_correlation.import")
def importCorrelation( infile, outfile ):

    tablename = outfile[:-len(".import")] 

    statement = '''
    csv2db.py %(csv2db_options)s \
              --index=id \
              --map=default:float \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )




############################################################
############################################################
############################################################
@files_re( "*.sites",
           "(.*).sites",
           r"\1.motif" )
def makeMotifs( infile, outfile ):
    '''build motifs from jasper sites by running MEME.'''

    tmpdir = tempfile.mkdtemp()
    tmpname = os.path.join( tmpdir, "tmpfile") 
    tmpfile = open( tmpname, "w")
    for line in open(infile,"r"):
        tmpfile.write( re.sub( "[ \t]+", "_", line) )
    tmpfile.close()

    statement = '''
    %(execmeme)s %(tmpname)s -mod oops -dna -revcomp -nmotifs 1 -text -oc %(tmpdir)s > %(outfile)s
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

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
    job_queue = "medium_jobs.q"
    job_options = "-l mem_free=8000M"

    target_path = os.path.join( os.path.abspath(PARAMS["exportdir"]), "meme", outfile )

    dbhandle = sqlite3.connect( PARAMS["database"] )
    fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )

    track = infile[:-len(".fasta")]
    cc = dbhandle.cursor()
    statement = "SELECT peakcenter, interval_id, contig FROM %s_intervals ORDER BY peakval DESC" % track
    cc.execute( statement )
    data = cc.fetchall()
    cc.close()

    cutoff = len(data) // 10
    # maximum size of data set (in characters)
    maxsize = int(PARAMS["meme_max_size"])

    E.info( "runMeme: %s: using at most %i sequences for pattern finding" % (track, cutoff) )

    fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )

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

        if PARAMS["meme_masker"] == "repeatmasker":
            # the genome sequence is repeat masked
            masked_seq = seq
        elif PARAMS["meme_masker"] == "dustmasker":
            masked_seq = masker( seq.upper() )

        # hard mask softmasked characters
        masked_seq = re.sub( "[a-z]","N", masked_seq )
        current_size += len(masked_seq)
        if current_size >= maxsize: 
            E.info( "runMEME: %s: maximum size (%i) reached - only %i sequences output" % (track, maxsize, nseq))
            break
        nseq += 1
        outs.write( ">%s\n%s\n" % (id, masked_seq))

    outs.close()

    statement = '''
    %(execmeme)s %(tmpfasta)s -dna -revcomp -mod %(meme_model)s -nmotifs %(meme_nmotifs)s -oc %(tmpdir)s -maxsize %(maxsize)s %(meme_options)s > %(outfile)s.log
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

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

    dbhandle = sqlite3.connect( PARAMS["database"] )

    fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )

    cc = dbhandle.cursor()
        
    if full:
        statement = '''SELECT start, end, interval_id, contig FROM
        %(track)s_intervals ORDER BY peakval DESC''' % locals()
    elif halfwidth:
        statement = '''SELECT peakcenter - %(halfwidth)s, peakcenter + %(halfwidth)s,
                       interval_id, contig FROM %(track)s_intervals ORDER BY peakval DESC''' % locals()
    else:
        raise ValueError("either specify full or halfwidth" )
    
    cc.execute( statement )
    data = cc.fetchall()
    cc.close()

    if proportion:
        cutoff = int(len(data) * proportion)
    else:
        cutoff = len(data)
        
        E.info( "writeSequencesForIntervals %s: using at most %i sequences for pattern finding" % (track, cutoff) )
    E.info( "writeSequencesForIntervals %s: masker=%s" % (track,masker))

    fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )
    masker_object = Masker.MaskerDustMasker()

    sequences = []
    current_size, nseq = 0, 0
    new_data = []
    for start, end, id, contig in data[:cutoff]:
        lcontig = fasta.getLength( contig )
        start, end = max(0, start + offset), min(end + offset, lcontig)
        if start >= end:
            E.info( "writeSequencesForIntervals %s: sequence %s is empty: start=%i, end=%i, offset=%i - ignored" % (track, id, start, end, offset))
            continue
        seq = fasta.getSequence( contig, "+", start, end )
        sequences.append( seq )
        new_data.append( (start, end, id, contig) )
        current_size += len(seq)
        if maxsize and current_size >= maxsize: 
            E.info( "writeSequencesForIntervals %s: maximum size (%i) reached - only %i sequences output (%i ignored)" % (track, maxsize, nseq,
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
    job_queue = "server_jobs.q"

    target_path = os.path.join( os.path.abspath( PARAMS["exportdir"] ), "glam2", outfile )
    track = infile[:-len(".fasta")]

    tmpdir = tempfile.mkdtemp()
    tmpfasta =  os.path.join( tmpdir, "in.fa")
    
    nseq = writeSequencesForIntervals( track, tmpfasta,
                                       full = False,
                                       halfwidth = int(PARAMS["meme_halfwidth"]),
                                       maxsize = int(PARAMS["meme_max_size"]),
                                       proportion = 0.1 )

    min_sequences = int(nseq / 10.0)
    statement = '''
    %(execglam2)s -2 -O %(tmpdir)s %(glam2_options)s -z %(min_sequences)i n %(tmpfasta)s > %(outfile)s.log
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

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

if PARAMS["motifs_tomtom_master_motif"] != "":
    ############################################################
    ############################################################
    ## filter motifs by a master motif
    ############################################################
    ############################################################
    ############################################################
    @follows( makeMotifs )
    @transform( runMEME, suffix(".meme"), ".tomtom" )
    def runTomTom( infile, outfile ):
        '''compare ab-initio motifs against tomtom.'''
        to_cluster = True
        statement = '''
           tomtom -text -query %(tomtom_master_motif)s -target %(infile)s > %(outfile)s
           '''
        P.run()

    ############################################################
    ############################################################
    ############################################################
    @transform( runTomTom, suffix(".tomtom"), "_tomtom.import" )
    def importTomTom( infile, outfile ):
        '''compare ab-initio motifs against tomtom.'''

        tmpfile = P.getTempFile()

        tmpfile.write( "\t".join( \
            ("query_id", "target_id",
             "optimal_offset", "pvalue", "qvalue", "overlap", "query_consensus",
             "target_consensus", "orientation") ) + "\n" )

        for line in open( infile, "r" ):
            if line.startswith( "#Query" ): continue
            (query_id, target_id, 
             optimal_offset, pvalue, qvalue, overlap, query_consensus,
             target_consensus, orientation) = line[:-1].split("\t")
            tmpfile.write( line )

        tmpfile.close()

        tmpname = tmpfile.name
        tablename = outfile[:-len(".import")]

        statement = '''
        csv2db.py %(csv2db_options)s \
                  --allow-empty \
                  --table=%(tablename)s \
        < %(tmpname)s > %(outfile)s
        '''

        P.run( **dict( locals().items() + PARAMS.items() ) )
        os.unlink( tmpname )

    @transform( runTomTom, suffix(".tomtom"), ".motif" )
    def filterMotifs( infile, outfile ):
        '''scan motifs found with MEME against master motif using tomtom.
        '''

        infile_meme = infile[:-len(".tomtom")] + ".meme"

        selected = []
        max_pvalue = float(PARAMS["motifs_tomtom_filter_pvalue"])
        for line in open( infile, "r" ):
            if line.startswith( "#Query" ): continue
            (query_id, target_id, 
             optimal_offset, pvalue, qvalue, overlap, query_consensus,
             target_consensus, orientation) = line[:-1].split("\t")
            if float(pvalue) <= max_pvalue:
                selected.append( target_id )

        E.info( "%s: keeping %i motifs" % (infile, len(selected) ) )

        PMotifs.filterMotifsFromMEME( infile_meme, outfile, selected )
        
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
    def importTomTom(): pass
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runMEME, suffix(".meme"), ".motif" )
    def filterMotifs( infile, outfile ):
        '''take the top scoring motif from meme runs
        '''

        PMotifs.filterMotifsFromMEME( infile, outfile, ["1"] )

############################################################
############################################################
############################################################
# todo: fix, causes a problem: exportSequences and exportMASTControlSequences,
@follows( makeMotifs, filterMotifs )
@files_re( (exportMotifSequences, exportMotifControlSequences),
           "(\S+).controlfasta",
           [ r"\1.controlfasta", r"\1.fasta",  glob.glob("*.motif")],
           r"\1.mast" )
def runMAST( infiles, outfile ):
    '''run mast on all intervals and motifs.

    Collect all results for an E-value up to 10000 so that
    all sequences are output and MAST curves can be computed. 

    10000 is a heuristic.
    '''
    to_cluster = True
    job_queue = "medium_jobs.q"
    # only use new nodes, as /bin/csh is not installed
    # on the old ones.
    job_options = "-l mem_free=8000M"

    controlfile, dbfile, motiffiles  = infiles
    controlfile = dbfile[:-len(".fasta")] + ".controlfasta"
    if not os.path.exists( controlfile ):
        raise P.PipelineError( "control file %s for %s does not exist" % (controlfile, dbfile))

    if os.path.exists(outfile): os.remove( outfile )
    
    for motiffile in motiffiles:
        if IOTools.isEmpty( motiffile ):
            E.info( "skipping empty motif file %s" % motiffile )
            continue
        
        of = open(outfile, "a")
        motif, x = os.path.splitext( motiffile )
        of.write(":: motif = %s ::\n" % motif )
        of.close()
        
        statement = '''
        cat %(dbfile)s %(controlfile)s 
        mast %(motiffile)s -stdin -stdout -text -ev %(motif_mast_evalue)f >> %(outfile)s
        '''
        P.run( **dict( locals().items() + PARAMS.items() ) )

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
    to_cluster = True

    # only use new nodes, as /bin/csh is not installed
    # on the old ones.
    job_options = "-l mem_free=8000M"

    tmpfasta = P.getTempFilename( "." )
    track = outfile[:-len(".bioprospector")]
    nseq = writeSequencesForIntervals( track,
                                       tmpfasta,
                                       full = True,
                                       masker = "dust",
                                       proportion = 0.10 )

    statement = '''
    bioprospector -i %(tmpfasta)s %(bioprospector_options)s -o %(outfile)s > %(outfile)s.log
    '''
    P.run()

    os.unlink( tmpfasta )

############################################################
############################################################
############################################################
##
############################################################
@transform( runBioProspector, suffix(".bioprospector"), "_bioprospector.import")
def importBioProspector( infile, outfile ):
    '''import results from bioprospector.'''

    tablename = outfile[:-len(".import")]
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
    csv2db.py %(csv2db_options)s \
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
            suffix(".mast"),
            "_mast.import" )
def importMAST( infile, outfile ):
    '''parse mast file and import into database.

    Parse several motif runs and add them to the same
    table.

    Add columns for the control data as well.
    '''

    tablename = os.path.basename( infile )
    tmpfile = tempfile.NamedTemporaryFile(delete=False)

    tmpfile.write( MAST.Match().header +\
                   "\tmotif\tcontig" \
                   "\tl_evalue\tl_pvalue\tl_nmatches\tl_length\tl_start\tl_end" \
                   "\tr_evalue\tr_pvalue\tr_nmatches\tr_length\tr_start\tr_end" \
                   "\tmin_evalue\tmin_pvalue\tmax_nmatches" + "\n" )

    lines = open(infile).readlines()
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
    csv2db.py %(csv2db_options)s \
              -b sqlite \
              --index=id \
              --index=motif \
              --index=id,motif \
              --table=%(tablename)s \
              --map=base_qualities:text \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )
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

@follows( makeMotifs, runGLAM2 )
@files_re( (exportMotifSequences, exportMotifControlSequences),
           "(\S+).controlfasta",
           [ r"\1.controlfasta", r"\1.fasta",  glob.glob("*.glam2")],
           r"\1.glam2scan" )
def runGLAM2SCAN( infiles, outfile ):
    '''run glam2scan on all intervals and motifs.
    '''

    to_cluster = True
    job_queue = "medium_jobs.q"
    # only use new nodes, as /bin/csh is not installed
    # on the old ones.
    job_options = "-l mem_free=8000M"

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
            "_glam.import" )
def importGLAM2SCAN( infile, outfile ):
    '''parse mast file and import into database.

    Parse several motif runs and add them to the same
    table.
    '''

    tablename = outfile[:-len(".import")]
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
            E.warn( "no results for motif %s - ignored" % motif )
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
    csv2db.py %(csv2db_options)s \
              -b sqlite \
              --index=id \
              --index=motif \
              --index=id,motif \
              --table=%(tablename)s \
              --map=base_qualities:text \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )
    os.unlink( tmpfile.name )

############################################################
############################################################
############################################################
@files( [ ("%s.bed" % x, "%s.annotations" % x ) for x in TRACKS_MASTER ] )
def annotateIntervals( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set.
    '''
    to_cluster = True
    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python /home/andreas/gpipe/gtf2table.py \
		--counter=position \
		--counter=classifier-chipseq \
		--section=exons \
		--counter=length \
		--log=%(outfile)s.log \
		--filename-gff=%(annotation)s \
		--genome-file=%(genome)s
	> %(outfile)s"""

    P.run()

def straightImport( infile, outfile, options = "" ):
    '''straight import from gtf2table results.

    The table name is given by outfile without the
    ".import" suffix.
    '''

    tablename = outfile[:-len(".import")]

    statement = '''
    csv2db.py %(csv2db_options)s \
              --index=gene_id \
              %(options)s \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@follows( annotateIntervals )
@files_re( "*.annotations", "(.*).annotations", r"\1_annotations.import" )
def importAnnotations( infile, outfile ):
    '''import interval annotations: genome architecture
    '''
    straightImport( infile, outfile )

############################################################
############################################################
############################################################
@files( [ ("%s.bed" % x, "%s.tss" % x ) for x in TRACKS_MASTER ] )
def annotateTSS( infile, outfile ):

    to_cluster = True
    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python /home/andreas/gpipe/gtf2table.py \
		--counter=distance-tss \
		--log=%(outfile)s.log \
		--filename-gff=%(transcripts_tss)s \
		--genome-file=%(genome)s
	> %(outfile)s"""

    P.run()

############################################################
############################################################
############################################################
@files( [ ("%s.bed" % x, "%s.repeats" % x ) for x in TRACKS_MASTER ] )
def annotateRepeats( infile, outfile ):
    '''count the overlap between intervals and repeats.'''

    to_cluster = True
    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python /home/andreas/gpipe/gtf2table.py \
		--counter=overlap \
		--log=%(outfile)s.log \
		--filename-gff=%(repeats)s \
		--genome-file=%(genome)s
	> %(outfile)s"""

    P.run()

############################################################
############################################################
############################################################
@files_re( annotateTSS, "(.*).tss", r"\1_tss.import" )
def importTSS( infile, outfile ):
    '''import interval annotations: distance to transcription start sites
    '''
    straightImport( infile, outfile, options = "--index=closest_id --index=id5 --index=id3" )

############################################################
############################################################
############################################################
@transform( annotateRepeats, suffix(".repeats"), "_repeats.import" )
def importRepeats( infile, outfile ):
    '''import interval annotations: repeats
    '''
    straightImport( infile, outfile )

############################################################
############################################################
############################################################
@files( [ ("%s.bed" % x, "%s.probetss" % x ) for x in TRACKS_MASTER ] )
def annotateProbesetTSS( infile, outfile ):

    to_cluster = True
    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python /home/andreas/gpipe/gtf2table.py \
		--counter=distance-tss \
		--log=%(outfile)s.log \
		--filename-gff=%(probeset)s \
		--genome-file=%(genome)s
	> %(outfile)s"""

    P.run()

############################################################
############################################################
############################################################
## count coverage within intervals for each track agains all
## the Unstim tracks
############################################################
@files( [ ("%s.bed" % x, "%s.readcounts" % x ) for x in TRACKS_NOT_UNSTIM ] )
def buildIntervalCounts( infile, outfile ):
    '''count read density in bed files.
    '''
    track = outfile[:-len(".readcounts")]

    tissue, condition, replicate = splitTrack( track )

    samfiles_fg, samfiles_bg = [], []

    # save original, because condition is later assigned to
    one = tissue and condition
    two = condition

    if REPLICATES:
        if one:
            for replicate in REPLICATES:
                samfiles_fg.append( FILEPATTERN % locals() + BAM_SUFFIX )
            condition = "Unstim"
            for replicate in REPLICATES:
                samfiles_bg.append( FILEPATTERN % locals() + BAM_SUFFIX )
        elif two:
            for tissue in TISSUES:
                for replicate in REPLICATES:
                    samfiles_fg.append( FILEPATTERN % locals() + BAM_SUFFIX )
            condition = "Unstim"
            for tissue in TISSUES:
                for replicate in REPLICATES:
                    samfiles_bg.append( FILEPATTERN % locals() + BAM_SUFFIX )

    else:
        samfiles_fg.append( FILEPATTERN % locals() + BAM_SUFFIX )
        condition = "Unstim"
        if os.path.exists(FILEPATTERN % locals() + BAM_SUFFIX ):
            samfiles_bg.append( FILEPATTERN % locals() + BAM_SUFFIX )
        else:
            samfiles_bg = samfiles_fg
            
    samfiles_fg = ",".join(samfiles_fg)
    samfiles_bg = ",".join(samfiles_bg)

    tmpfile1 = P.getTempFilename( os.getcwd() ) + ".fg"
    tmpfile2 = P.getTempFilename( os.getcwd() ) + ".bg"

    to_cluster = True

    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python /home/andreas/gpipe/gtf2table.py \
		--counter=read-coverage \
		--log=%(outfile)s.log \
		--bam-file=%(samfiles_fg)s \
	> %(tmpfile1)s"""

    P.run( **dict( locals().items() + PARAMS.items() ) )

    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python /home/andreas/gpipe/gtf2table.py \
		--counter=read-coverage \
		--log=%(outfile)s.log \
		--bam-file=%(samfiles_bg)s \
	> %(tmpfile2)s"""

    P.run( **dict( locals().items() + PARAMS.items() ) )

    statement = '''
    python %(toolsdir)s/combine_tables.py \
           --add-file-prefix \
           --regex-filename="[.](\S+)$" \
    %(tmpfile1)s %(tmpfile2)s > %(outfile)s
    '''

    P.run()
    os.unlink( tmpfile1 )
    os.unlink( tmpfile2 )


@files_re( buildIntervalCounts, "(.*).readcounts", r"\1_readcounts.import" )
def importIntervalCounts( infile, outfile ):
    '''import interval counts.'''
    tablename = outfile[:-len(".import")]

    statement = '''
    csv2db.py %(csv2db_options)s \
              --index=gene_id \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

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
@files_re(buildBigwig,combine("(.*).bigwig"), "bigwig.view")
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

@follows( )
def prepare():
    pass

@follows( buildIntervals )
def intervals():
    pass

@follows( buildBigwig, viewIntervals, viewBigwig )
def export():
    pass

@follows( makeMotifs,
          exportMotifSequences,
          runMEME,
          runTomTom, importTomTom,
          filterMotifs,
          runGLAM2, 
          exportMotifControlSequences,
          runBioProspector )
def build_motifs():
    pass

@follows( runMAST, importMAST,
          runGLAM2SCAN, importGLAM2SCAN )
def find_motifs():
    pass

@follows( importCorrelation, 
          reproducibility)
def correlation():
    pass

@follows( annotateIntervals, importAnnotations, 
          annotateTSS, importTSS, 
          annotateRepeats, importRepeats,
          #annotateTSSIntervalAssociations, importTSSIntervalAssociations,
          #annotateTSSIntervalDistance, importTSSIntervalDistance, 
          buildIntervalCounts, importIntervalCounts )
def annotation():
    pass

@follows( prepare, intervals, 
          build_motifs, find_motifs,
          correlation, annotation )
def full():
    pass

if __name__== "__main__":

    print( "tracks found: %s" % TRACKS_ALL )
    print( "tracks by experiment: %s" % TRACKS_EXPERIMENTS )
    print( "tracks by condition: %s" % TRACKS_CONDITIONS )
    print( "tracks by tissue: %s" % TRACKS_TISSUES )

    sys.exit( P.main(sys.argv) )


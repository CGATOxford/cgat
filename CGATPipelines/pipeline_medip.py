################################################################################
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
====================
MeDIP pipeline
====================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

   1. Align to genome using gapped alignment (BWA)
   2. Calculate alignment and coverage statistics (BAMStats)
   3. Identify differentially methylated regions (DMRs)
   4. Filter DMRs
   5. Calculate DMR statistics
   6. Produce report (SphinxReport)

Methods
=======

Read processing
---------------

For medip-seq analysis, the following filtering steps are typically applied to the mapped data:

   1. removing duplicate reads
   2. removing reads with a mapping quality of less than 10

Medip analysis
--------------

The medip analysis makes use of the MEDIPS R package by `Chavez et al. <http://medips.molgen.mpg.de/>`_ 
(see :pmid:`PMID: 20802089`).

Briefly, the data is processed in the following way:

1. Quality control
   1. Saturation analysis
   2. Computing CpG coverage
   3. Computing CpG enrichment
   
2. Normalization
   1. Output data normalized by total read depth (rpm - reads per million)
   2. Output normalized relative methylation scores (rms)
   3. Output normalized absolute methylation scores (ams)


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

Input
-----

Reads are imported by placing files or linking to files in the :term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while ``replicate`` denotes
the :term:`replicate` within an :term:`experiment`. The ``suffix`` determines the file type.
The following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format. The two fastq files must be sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input files. Thus it might be
   difficult to mix different formats.

Requirements
------------

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|Stampy              |>=0.9.0            |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|BWA                 |                   |read mapping                                    |
+--------------------+-------------------+------------------------------------------------+
|SAMtools            |                   |filtering, SNV / indel calling                  |
+--------------------+-------------------+------------------------------------------------+
|BEDTools            |                   |filtering, SNV / indel calling                  |
+--------------------+-------------------+------------------------------------------------+
|sra-tools           |                   |extracting reads from .sra files                |
+--------------------+-------------------+------------------------------------------------+
|picard              |>=1.38             |bam/sam files. The .jar files need to be in your|
|                    |                   | CLASSPATH environment variable.                |
+--------------------+-------------------+------------------------------------------------+
|vcf-tools           |                   |                                                |
+--------------------+-------------------+------------------------------------------------+
|BAMStats            |                   |                                                |
+--------------------+-------------------+------------------------------------------------+

Pipeline output
===============

?

Example
=======

ToDo: make exome sequencing example


Code
====

"""

# load modules
from ruffus import *

import Experiment as E
import logging as L
import Database
import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random, csv
import numpy, sqlite3
import GTF, IOTools, IndexedFasta
import PipelineGeneset
import PipelineMapping
import Stats
import PipelineTracks
import PipelineMappingQC
import PipelineMedip
import Pipeline as P
import Expression

from rpy2.robjects import r as R
import rpy2.robjects as ro

#########################################################################
#########################################################################
#########################################################################
# load options from the config file
P.getParameters( ["%s/pipeline.ini" % os.path.splitext(__file__)[0], 
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
## TRIM READS
## this should go elsewhere, the readqc pipeline?
@follows(mkdir("trim"))
@transform( "*.gz", regex( r"(\S+).gz"), r"trim/\1.gz" )
def trimReads( infile, outfile ):
    '''trim reads with FastX'''
    to_cluster = True
    first_base = PARAMS["trim_first_base"]
    last_base = PARAMS["trim_last_base"]

    statement = """zcat %(infile)s | fastx_trimmer -f %(first_base)s -l %(last_base)s -z -o %(outfile)s """
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( ("*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra"),
             regex( r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
             r"\1.dir")
def makeTrackDirectories( infile, outfile ):
    '''make track directories.'''
    os.mkdir( outfile )
    os.mkdir( os.path.join( outfile, "bam" ) )

#########################################################################
#########################################################################
#########################################################################
## Map reads to genome using BWA
@follows( makeTrackDirectories, mkdir( PARAMS["exportdir"] ) )
@transform( ("*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra"),
             regex( r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
             r"\1.dir/\1.genome.bam")
def mapReads(infiles, outfile):
    '''Map reads to the genome using BWA '''
    to_cluster = True
    job_options= "-pe dedicated %i -R y" % PARAMS["bwa_threads"]
    m = PipelineMapping.BWA()
    statement = m.build((infiles,), outfile) 
    P.run()

#########################################################################
#########################################################################
#########################################################################
#@transform( "*CD4*/bam/*.genome.bam",
@transform( mapReads,
            suffix(".genome.bam"),
            ".prep.bam" )
def prepareBAMs( infile, outfile ):
    '''filter bam files for medip-seq analysis.

    Optional steps include:

    * deduplication - remove duplicate reads
    * quality score filtering - remove reads below a certain quality score.

    '''
    to_cluster = True
    track = P.snip( outfile, ".bam" )

    tmpdir = P.getTempFilename()
    
    current_file = infile

    nfiles = 0
    statement = [ "mkdir %(tmpdir)s" ]

    if "filtering_quality" in PARAMS and PARAMS["filtering_quality"] > 0:
        next_file = "%(tmpdir)s/bam_%(nfiles)i.bam" % locals()
        statement.append( '''samtools view -q %%(filtering_quality)i -b 
                             %(current_file)s 
                             2>> %%(outfile)s.log 
                             > %(next_file)s ''' % locals())
        nfiles += 1
        current_file = next_file

    if "filtering_dedup" in PARAMS and PARAMS["filtering_dedup"]:
        # Picard's MarkDuplicates requries an explicit bam file.
        next_file = "%(tmpdir)s/bam_%(nfiles)i.bam" % locals()

        dedup_method = PARAMS["filtering_dedup_method"]
        
        if dedup_method == 'samtools':
            statement.append( '''samtools rmdup - - ''' )

        elif dedup_method == 'picard':
            statement.append('''MarkDuplicates INPUT=%(current_file)s
                                               OUTPUT=%(next_file)s
                                               ASSUME_SORTED=true 
                                               METRICS_FILE=%(outfile)s.duplicate_metrics
                                               REMOVE_DUPLICATES=TRUE 
                                               VALIDATION_STRINGENCY=SILENT
                                               2>> %%(outfile)s.log ''' % locals() )
        nfiles += 1
        current_file = next_file

    statement.append( "mv %%(current_file)s %(outfile)s" % locals() )
    statement.append( "rm -rf %(tmpdir)s" )
    statement.append( "samtools index %(outfile)s" )

    statement = " ; ".join( statement )

    P.run()

    os.unlink( tmpdir )

#########################################################################
#########################################################################
#########################################################################
@merge( prepareBAMs, "picard_duplicates.load" )
def loadPicardDuplicateStats( infiles, outfile ):
    '''Merge Picard duplicate stats into single table and load into SQLite.
    '''
    PipelineMappingQC.loadPicardDuplicateStats( infiles, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( (mapReads, prepareBAMs), 
            suffix(".bam"),
            ".picard_stats")
def buildPicardAlignmentStats( infile, outfile ):
    '''Gather BAM file alignment statistics using Picard '''

    PipelineMappingQC.buildPicardAlignmentStats( infile, outfile, 
                                                 os.path.join( PARAMS["bwa_index_dir"],
                                                               PARAMS["genome"] + ".fa" ) )

############################################################
############################################################
############################################################
@merge( buildPicardAlignmentStats, "picard_stats.load" )
def loadPicardAlignmentStats( infiles, outfile ):
    '''Merge Picard alignment stats into single table and load into SQLite.'''

    PipelineMappingQC.loadPicardAlignmentStats( infiles, outfile )
    
#########################################################################
#########################################################################
#########################################################################
@transform( (mapReads, prepareBAMs), 
            suffix(".bam"),
            ".gcstats" )
def buildPicardGCStats( infile, outfile ):
    '''Gather BAM file GC bias stats using Picard '''
    PipelineMappingQC.buildPicardGCStats( infile, outfile, 
                                                 os.path.join( PARAMS["bwa_index_dir"],
                                                               PARAMS["genome"] + ".fa" ) )


#########################################################################
#########################################################################
#########################################################################
@merge( buildPicardGCStats, "picard_gcbias_stats.load" )
def loadPicardGCStats( infiles, outfile ):
    '''Merge Picard insert size stats into single table and load into SQLite.'''
    
    tablename = P.toTable( outfile )
    outf = P.getTempFile()

    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".gcstats" )
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
                      %(csv2db_options)s
                      --index=track
                      --table=%(tablename)s 
                   > %(outfile)s '''
    P.run()

    os.unlink( tmpfilename )

#########################################################################
#########################################################################
#########################################################################
@transform( (mapReads, prepareBAMs),
            suffix(".bam"),
            ".readstats" )
def buildBAMStats( infile, outfile ):
    '''Count number of reads mapped, duplicates, etc. '''
    PipelineMappingQC.buildBAMStats( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''Import bam statistics into SQLite'''
    PipelineMappingQC.loadBAMStats( infiles, outfile )

#########################################################################
#########################################################################
#########################################################################
## Methylation analysis
@transform( prepareBAMs, suffix(".bam"), ".medips")
def runMEDIPS( infile, outfile ):
    '''run MEDIPS analysis - 
    outputs methylation profiles.
    '''

    to_cluster = True

    job_options = "-l mem_free=23G"

    statement = '''
    cat %(infile)s 
    | python %(scriptsdir)s/bam2bed.py
          --merge-pairs
          --min-insert-size=%(medips_min_insert_size)i
          --max-insert-size=%(medips_max_insert_size)i
          --log=%(outfile)s.log
          -
    | python %(scriptsdir)s/WrapperMEDIPS.py
         --ucsc-genome=%(genome)s
         --genome-file=%(genome_dir)s/%(genome)s
         --bigwig
         --input-format=bed 
         --extension=%(medips_extension)i
         --fragment-length=%(medips_fragment_length)i
         --force
         --bin-size=%(medips_bin_size)i
         --output-filename-pattern="%(outfile)s_%%s"
         -
    >& %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( runMEDIPS, suffix(".medips"), "_medips.load")
def loadMEDIPS( infile, outfile ):
    '''load medips results'''

    table_prefix = re.sub( "_prep", "", P.toTable( outfile ))
    
    table = table_prefix + "_coveredpos" 

    statement = """
    cat %(infile)s_saturation_coveredpos.csv
    | tail -n 3 
    | perl -p -e 's/\\"//g; s/[,;]/\\t/g; '
    | python %(scriptsdir)s/table2table.py --transpose
    | python %(scriptsdir)s/csv2db.py 
           %(csv2db_options)s
           --table=%(table)s
           --replace-header
           --header=coverage,ncovered,pcovered
    >> %(outfile)s
    """
    
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( prepareBAMs, suffix(".bam"), ".covered.bed.gz" )
def buildCoverageBed( infile, outfile ):
    '''build bed file with regions covered by reads.

    Intervals containing only few reads (tiling_min_reads) are removed.
    '''
    
    to_cluster = True

    statement = '''
    cat %(infile)s 
    | python %(scriptsdir)s/bam2bed.py
          --merge-pairs
          --min-insert-size=%(medips_min_insert_size)i
          --max-insert-size=%(medips_max_insert_size)i
          --log=%(outfile)s.log
          -
    | sort -k1,1 -k2,2n
    | cut -f 1,2,3
    | python %(scriptsdir)s/bed2bed.py
          --method=sanitize-genome
          --genome-file=%(genome_dir)s/%(genome)s
          --log=%(outfile)s.log
    | python %(scriptsdir)s/bed2bed.py
          --method=merge
          --merge-distance=%(medips_extension)i
          --log=%(outfile)s.log
          --merge-min-intervals=%(tiling_min_reads)i
    | gzip
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( prepareBAMs, 
            suffix(".bam"), 
            add_inputs( os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_cpg_bed"] ) ),
            ".cpg_coverage.gz" )
def buildCpGCoverage( infiles, outfile ):
    '''count number times certain CpG are covered by reads.

    Reads are processed in the same way as by buildCoverageBed.
    '''
    
    infile, cpg_file = infiles
    to_cluster = True

    job_options = "-l mem_free=10G"

    statement = '''
    cat %(infile)s 
    | python %(scriptsdir)s/bam2bed.py
          --merge-pairs
          --min-insert-size=%(medips_min_insert_size)i
          --max-insert-size=%(medips_max_insert_size)i
          --log=%(outfile)s.log
          -
    | coverageBed -a stdin -b %(cpg_file)s -counts
    | cut -f 6
    | python %(scriptsdir)s/data2histogram.py
    | gzip
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( buildCpGCoverage, "cpg_coverage.load" )
def loadCpGCoverage( infiles, outfile ):
    '''load cpg coverag data.'''
    P.mergeAndLoad( infiles, outfile, regex="/(.*).prep.cpg_coverage.gz",
                    row_wise = False )

#########################################################################
#########################################################################
#########################################################################
@merge( buildCoverageBed, "tiles_variable_width.bed.gz" )
def buildVariableWidthTiles( infiles, outfile ):
    '''bed file with intervals that are covered by reads in any of the experiments.
    '''
    
    infiles = " ".join( infiles )
    to_cluster = True

    statement = '''
    zcat %(infiles)s 
    | sort -k1,1 -k2,2n
    | python %(scriptsdir)s/bed2bed.py
          --method=merge
          --merge-distance=0
          --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
## Run DESeq to identify differentially methylated regions
@files( (( None, "tiles_fixednonovl.bed.gz"),) )
def buildNonoverlappingFixedWidthTiles( infile, outfile ):
    '''Build bed file segmenting entire genome using window x and shift y'''

    shift = PARAMS["tiling_nonoverlapping_window"]
    statement = '''python %(scriptsdir)s/genome_bed.py
                      -g %(genome_dir)s/%(genome)s
                      --window=%(tiling_nonoverlapping_window)i
                      --shift=%(shift)i
                      --log=%(outfile)s.log
                | awk '$1 !~ /%(tiling_remove_contigs)s/'
                | gzip
                > %(outfile)s'''
    P.run()

@files( (( None, "tiles_fixedovl.bed.gz"),) )
def buildOverlappingFixedWidthTiles( infile, outfile ):
    '''Build bed file segmenting entire genome using window x and shift y'''
    assert PARAMS["tiling_overlapping_window"] % 2 == 0
    shift = PARAMS["tiling_overlapping_window"] // 2

    statement = '''python %(scriptsdir)s/genome_bed.py
                      -g %(genome_dir)s/%(genome)s
                      --window=%(tiling_overlapping_window)i
                      --shift=%(shift)i
                      --log=%(outfile)s.log
                | awk '$1 !~ /%(tiling_remove_contigs)s/'
                | gzip
                > %(outfile)s'''
    P.run()

@transform( (buildNonoverlappingFixedWidthTiles,
             buildOverlappingFixedWidthTiles,
             buildVariableWidthTiles ), 
            suffix(".bed.gz"), 
            ".bed.gz")
def buildTiles( infile, outfile ):
    pass

#########################################################################
#########################################################################
#########################################################################
@transform( buildTiles,
            suffix(".bed.gz"),
            ".stats")
def buildTileStats( infile, outfile ):
    '''compute tiling window size statistics from bed file.'''

    use_cluster = True

    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/gff2histogram.py 
                   --force
                   --format=bed 
                   --data=size
                   --method=hist
                   --method=stats
                   --output-filename-pattern=%(outfile)s.%%s.tsv
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( buildTiles,
            suffix(".bed.gz"), 
            ".bigbed")
def buildBigBed( infile, outfile ):
    '''bed file with intervals that are covered by reads in any of the experiments.
    '''
    
    to_cluster = True
    to_cluster = False

    tmpfile = P.getTempFilename()

    contig_sizes = os.path.join( PARAMS["annotations_dir"], PARAMS_ANNOTATIONS["interface_contigs"] )

    statement = '''
    zcat %(infile)s > %(tmpfile)s;
    bedToBigBed %(tmpfile)s %(contig_sizes)s %(outfile)s;
    rm -f %(tmpfile)s
    '''
    P.run()

    try: os.unlink( tmpfile )
    except OSError: pass

#########################################################################
# add_inputs( buildGenomeTilingBed ),
def buildTiledReadCounts( infiles, outfile ):
    '''compute coverage of genome with reads.'''

    to_cluster = True

    infile, tiles = infiles

    # note: needs to set flags appropriately for
    # single-end/paired-end data sets
    # set filter options
    # for example, only properly paired reads
    paired = True
    if paired:
        flag_filter = "-f 0x2"
    else:
        flag_filter = ""

    statement = '''
    samtools view -b %(flag_filter)s -q %(deseq_min_mapping_quality)s %(infile)s 
    | python %(scriptsdir)s/bam2bed.py
          --merge-pairs
          --min-insert-size=%(medips_min_insert_size)i
          --max-insert-size=%(medips_max_insert_size)i
           --log=%(outfile)s.log
           - 
    | coverageBed -a stdin -b %(tiles)s 
    | sort -k1,1 -k2,2n
    | gzip > %(outfile)s '''
    
    P.run()

@transform( prepareBAMs,
            suffix(".bam"), 
            add_inputs( buildVariableWidthTiles ),
            r".variablewidth.tilecounts.bed.gz" )
def buildTiledReadCountsVariableWidth(infiles, outfile ):
    '''build read counds for variable width windows.'''
    buildTiledReadCounts( infiles, outfile )

@transform( prepareBAMs,
            suffix(".bam"), 
            add_inputs( buildNonoverlappingFixedWidthTiles ),
            r".fixedwidthnoovl.tilecounts.bed.gz" )
def buildTiledReadCountsFixedWidthNoOverlap(infiles, outfile ): 
    '''build read counds for fixed width windows.'''
    buildTiledReadCounts( infiles, outfile )

@transform( prepareBAMs,
            suffix(".bam"), 
            add_inputs( buildOverlappingFixedWidthTiles ),
            r".fixedwidthovl.tilecounts.bed.gz" )
def buildTiledReadCountsFixedWidthOverlap(infiles, outfile ):
    '''build read counds for fixed width windows.'''
    buildTiledReadCounts( infiles, outfile )

@transform( (buildTiledReadCountsVariableWidth,
             buildTiledReadCountsFixedWidthNoOverlap,
             buildTiledReadCountsFixedWidthOverlap),
            suffix(".bed.gz"),
            ".bed.gz" )
def buildAllTiledReadCounts( infile, outfile ):
    pass

#########################################################################
@follows( mkdir( "diff_methylation" ) )
@collate( buildAllTiledReadCounts,
          regex( ".*\.([^.]+).tilecounts.bed.gz"),
          r"diff_methylation/\1.counts.tsv.gz")
def aggregateTiledReadCounts( infiles, outfile ):
    '''aggregate tag counts for each window.

    coverageBed outputs the following columns:
    1) Contig
    2) Start
    3) Stop
    4) Name
    5) The number of features in A that overlapped (by at least one base pair) the B interval.
    6) The number of bases in B that had non-zero coverage from features in A.
    7) The length of the entry in B.
    8) The fraction of bases in B that had non-zero coverage from features in A.

    For bed: use column 5
    For bed6: use column 7
    For bed12: use column 13

    This method uses the maximum number of reads found in any interval as the tag count.

    Tiles with no counts will not be output.
    '''
    
    to_cluster = True

    src = " ".join( [ '''<( zcat %s | awk '{printf("%%s:%%i-%%i\\t%%i\\n", $1,$2,$3,$4 );}' ) ''' % x for x in infiles] )
    tmpfile = P.getTempFilename( "." )
    statement = '''paste %(src)s > %(tmpfile)s'''
    P.run()
    
    tracks = [ re.sub( "\..*", '', os.path.basename(x) ) for x in infiles ]

    outf = IOTools.openFile( outfile, "w")
    outf.write( "interval_id\t%s\n" % "\t".join( tracks ) )
    
    for line in open( tmpfile, "r" ):
        data = line[:-1].split("\t")
        genes = list(set([ data[x] for x in range(0,len(data), 2 ) ]))
        values = [ int(data[x]) for x in range(1,len(data), 2 ) ]
        if sum(values) == 0: continue
        assert len(genes) == 1, "paste command failed, wrong number of genes per line: '%s'" % line
        outf.write( "%s\t%s\n" % (genes[0], "\t".join(map(str, values) ) ) )
    
    outf.close()

    os.unlink(tmpfile)

#########################################################################
def loadMethylationData( infile, design_file ):
    '''load methylation data for deseq/edger analysis.
    
    This method creates various R objects:

    countsTable : data frame with counts. 
    groups : vector with groups

    '''
    
    E.info( "reading data")
    R( '''counts_table = read.delim( '%(infile)s', header = TRUE, 
                                                   row.names = 1, 
                                                   stringsAsFactors = TRUE )''' % locals() )

    E.info( "read data: %i observations for %i samples" % tuple(R('''dim(counts_table)''')))

    # Load comparisons from file
    R('''pheno = read.delim( '%(design_file)s', header = TRUE, stringsAsFactors = TRUE )''' % locals() )

    # Make sample names R-like - substitute - for . and add the .prep suffix
    R('''pheno[,1] = gsub('-', '.', pheno[,1]) ''')
    
    # Ensure pheno rows match count columns 
    R('''pheno2 = pheno[match(colnames(counts_table),pheno[,1]),,drop=FALSE]''' ) 

    # Subset data & set conditions
    R('''includedSamples <- pheno2$include == '1' ''')
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''conds <- pheno2$group[ includedSamples ]''')

    # Subset data & set conditions
    R('''includedSamples <- pheno2$include == '1' ''')
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''groups <- factor(pheno2$group[ includedSamples ])''')
    R('''pairs = factor(pheno2$pair[ includedSamples ])''')
    
    groups = R('''levels(groups)''')
    pairs = R('''levels(pairs)''')

    E.info( "filtered data: %i observations for %i samples" % tuple( R('''dim(countsTable)''') ) )
    
    return groups, pairs

#########################################################################
#########################################################################
#########################################################################
def runDE( infiles, outfile, method ):
    '''run DESeq or EdgeR.

    The job is split into smaller sections. The order of the input 
    data is randomized in order to avoid any biases due to chromosomes.
    '''

    to_cluster = True
    
    infile, design_file = infiles
    design = P.snip( os.path.basename(design_file), ".tsv")
    tiling = P.snip( os.path.basename( infile ), ".counts.tsv.gz" )

    outdir = os.path.join( PARAMS["exportdir"], "diff_methylation", "%s_%s_%s_" % (tiling, design, method ) )

    statement = '''zcat %(infile)s 
              | perl %(scriptsdir)s/randomize_lines.pl -h
              | %(cmd-farm)s
                  --input-header 
                  --output-header 
                  --split-at-lines=100000 
                  --cluster-options="-l mem_free=8G"
                  --log=%(outfile)s.log
                  --output-pattern=%(outdir)s%%s
                  --subdirs
              "python %(scriptsdir)s/Expression.py
              --method=%(method)s
              --filename-tags=-
              --filename-design=%(design_file)s
              --output-filename-pattern=%%DIR%%/
              --deseq-fit-type=%(deseq_fit_type)s
              --deseq-dispersion-method=%(deseq_dispersion_method)s
              --log=%(outfile)s.log
              --fdr=%(edger_fdr)f"
              | grep -v "warnings"
              | gzip
              > %(outfile)s '''

    P.run()


@follows( aggregateTiledReadCounts, mkdir( os.path.join( PARAMS["exportdir"], "diff_methylation")) )
@files( [ ( (data, design), 
            "diff_methylation/%s_%s.deseq.gz" % (P.snip(os.path.basename(data),".counts.tsv.gz"),
                                   P.snip(os.path.basename(design),".tsv" ) ) ) \
              for data, design in itertools.product( 
                                               glob.glob("diff_methylation/*.counts.tsv.gz"),
                                               P.asList(PARAMS["deseq_designs"]) ) ] )
def runDESeq( infiles, outfile ):
    '''estimate differential expression using DESeq.

    The final output is a table. It is slightly edited such that
    it contains a similar output and similar fdr compared to cuffdiff.
    '''

    runDE( infiles, outfile, "deseq" )

#########################################################################
#########################################################################
#########################################################################
@follows( aggregateTiledReadCounts, mkdir( os.path.join( PARAMS["exportdir"], "diff_methylation")) )
@files( [ ( (data, design), 
            "diff_methylation/%s_%s.edger.gz" % (P.snip(os.path.basename(data),".counts.tsv.gz"),
                                                 P.snip(os.path.basename(design),".tsv" ) ) ) \
              for data, design in itertools.product( 
                                               glob.glob("diff_methylation/*.counts.tsv.gz"),
                                               P.asList(PARAMS["deseq_designs"]) ) ] )
def runEdgeR( infiles, outfile ):
    '''estimate differential methylation using EdgeR
    
    This method applies a paired test. The analysis follows
    the example in chapter 11 of the EdgeR manual.
    '''

    runDE( infiles, outfile, "edger" )

#########################################################################
@transform( (runDESeq, runEdgeR), suffix(".gz"), ".merged.gz" )
def mergeDMRWindows( infile, outfile ):
    '''merge overlapping windows.'''

    to_cluster = True

    statement = '''
    zcat %(infile)s
    | python %(scriptsdir)s/medip_merge_intervals.py
          --log=%(outfile)s.log
          --output-filename-pattern=%(outfile)s.%%s.bed.gz
    | gzip
    > %(outfile)s
    '''

    P.run()

#########################################################################
@jobs_limit(1)
@transform( mergeDMRWindows, suffix(".merged.gz"), ".load" )
def loadDMRWindows( infile, outfile ):
    '''merge overlapping windows.'''
    P.load( infile, outfile, options = "--quick" )

#########################################################################
@collate( loadDMRWindows, regex( "(\S+)[.](\S+).load" ), r"\2_stats.tsv" )
def buildDMRStats( infiles, outfile ):
    '''compute differential methylation stats.'''
    tablenames = [P.toTable( x ) for x in infiles ] 
    method = P.snip( outfile, "_stats.tsv" )
    PipelineMedip.buildDMRStats( tablenames, method, outfile )

#########################################################################
@transform( buildDMRStats, suffix(".tsv"), ".load" )
def loadDMRStats( infile, outfile ):
    '''load DMR stats into table.'''
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( mergeDMRWindows,
            suffix(".merged.gz"),
            ".stats")
def buildDMRWindowStats( infile, outfile ):
    '''compute tiling window size statistics from bed file.'''

    to_cluster = True

    statement = '''
    zcat %(infile)s
    | grep -v 'contig'
    | python %(scriptsdir)s/gff2histogram.py 
                   --force
                   --format=bed 
                   --data=size
                   --method=hist
                   --method=stats
                   --output-filename-pattern=%(outfile)s.%%s.tsv
    > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( (buildTileStats, buildDMRWindowStats),
        "tileinfo.load" )
def loadTileStats( infiles, outfile ):
    '''load tiling stats into database.'''
    prefix = P.snip(outfile, ".load")

    files = " ".join( [ "%s.stats.tsv" % x for x in infiles ] )

    tablename = P.snip( outfile, ".load" ) + "_stats" 

    statement = """
    python %(scriptsdir)s/combine_tables.py 
           --cat=track 
           --regex-filename="(.*).stats.stats.tsv" 
           %(files)s
    | python %(scriptsdir)s/csv2db.py 
           %(csv2db_options)s
           --index=track
           --table=%(tablename)s 
    > %(outfile)s"""
    P.run()
   
    files = " ".join( [ "%s.hist.tsv" % x for x in infiles ] )

    tablename = P.snip( outfile, ".load" ) + "_hist" 
    
    statement = """
    python %(scriptsdir)s/combine_tables.py 
           --regex-filename="(.*).stats.hist.tsv" 
           --sort-keys=numeric
           --use-file-prefix
           %(files)s
    | python %(scriptsdir)s/csv2db.py 
           %(csv2db_options)s
           --index=track
           --table=%(tablename)s 
    >> %(outfile)s"""

    P.run()


#########################################################################
#########################################################################
#########################################################################
@transform( mergeDMRWindows, regex(  "(.*)\.(.*).merged.gz"), r"\1_\2.dmr.bed.gz" )
def buildDMRBed( infile, outfile ):
    '''output bed6 file with differentially methylated regions.

    Overlapping/book-ended entries are merged.

    The score is the average log fold change.
    '''
    
    to_cluster = True

    statement = '''zcat %(infile)s
    | python %(scriptsdir)s/csv_cut.py contig start end l2fold significant
    | awk '$5 == "1" {printf("%%s\\t%%i\\t%%i\\t%%i\\t%%f\\n", $1,$2,$3,++a,$4)}'
    | gzip > %(outfile)s'''

#    | mergeBed -i stdin -scores mean 

    P.run()

@merge( buildDMRBed, "dmr_overlap.tsv.gz" )
def computeDMROverlap( infiles, outfile ):
    '''compute overlap between bed sets.'''
    
    to_cluster = True

    if os.path.exists(outfile): 
        # note: update does not work due to quoting
        os.rename( outfile, "orig." + outfile )
        options = "--update=orig.%s" % outfile
    else:
        options = ""
    
    infiles = " ".join( infiles )

    # note: need to quote track names
    statement = '''
        python %(scriptsdir)s/diff_bed.py 
              --pattern-id=".*/(.*).dmr.bed.gz"
              --log=%(outfile)s.log 
              %(options)s %(infiles)s 
        | awk -v OFS="\\t" '!/^#/ { gsub( /-/,"_", $1); gsub(/-/,"_",$2); } {print}'
        | gzip
        > %(outfile)s
        '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( mergeDMRWindows, regex(  "(.*)\.(.*).merged.gz"), r"\1_\2.bed.gz" )
def buildMRBed( infile, outfile ):
    '''output bed6 file with methylated regions.

    All regions are output, even the insignificant ones.

    The score is the log fold change.
    '''
    
    outf = IOTools.openFile( outfile, "w" )
    c = E.Counter()
    for row in csv.DictReader( IOTools.openFile( infile ),
                               dialect = "excel-tab" ):
        c.input += 1

        contig, start, end = re.match("(.*):(\d+)-(\d+)", row["interval_id"] ).groups()
        c.output += 1
        outf.write( "\t".join( (contig, start, end, str(c.input), row["lfold"] ) ) + "\n" )
        
    outf.close()
    
    E.info( "%s" % str(c) )

#########################################################################
#########################################################################
#########################################################################
@follows( loadPicardDuplicateStats,
          loadPicardAlignmentStats,
          loadPicardGCStats,
          loadBAMStats )
def mapping(): pass

@follows( aggregateTiledReadCounts,
          loadDMRStats,
          buildDMRBed,
          computeDMROverlap)
def callDMRs(): pass

@follows( mapping, callDMRs) 
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
          mkdir("%s/medips" % PARAMS["web_dir"]),
          )
def publish():
    '''publish files.'''

    # publish web pages
    P.publish_report()

    # publish additional data
    web_dir = PARAMS["web_dir"]
    project_id = P.getProjectId()

    ucsc_urls = {
        "bam": 
        """track type=bam name="%(track)s" bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/%(dirname)s/%(filename)s""" ,
        "bigwig":
        """track type=bigWig name="%(track)s" bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/%(dirname)s/%(filename)s""" ,
        }
        
    # directory, files
    exportfiles = (
        ( "bamfiles", glob.glob( "*/*.genome.bam" ) + glob.glob( "*/*.genome.bam.bai" ), "bam" ),
        ( "bamfiles", glob.glob( "*/*.prep.bam" ) + glob.glob( "*/*.prep.bam.bai" ), "bam" ),
        ( "medips", glob.glob( "*/*.bigwig" ), "bigwig"),
        )
    
    ucsc_files = []

    for targetdir, filenames, datatype in exportfiles:
        for src in filenames:
            filename = os.path.basename(src)
            dest = "%s/%s/%s" % (web_dir, targetdir, filename)
            suffix = os.path.splitext( src )
            if suffix in ucsc_urls: ucsc_files.append( ( datatype, targetdir, filename ) )
            dest = os.path.abspath( dest )
            if not os.path.exists( dest ):
                os.symlink( os.path.abspath(src), dest )

    # output ucsc links
    for ucsctype, dirname, filename in ucsc_files:
        filename = os.path.basename( filename )
        track = P.snip( filename, ucsctype )
        print ucsc_urls[ucsctype] % locals()

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )



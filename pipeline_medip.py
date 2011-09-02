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
import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random
import numpy, sqlite3
import GTF, IOTools, IndexedFasta
import PipelineGeneset
import PipelineMapping
import Stats
import PipelineTracks
import PipelineMappingQC
import Pipeline as P

from rpy2.robjects import r as R
import rpy2.robjects as ro

USECLUSTER = True

#########################################################################
#########################################################################
#########################################################################
# load options from the config file
P.getParameters( ["%s.ini" % __file__[:-len(".py")], 
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
    to_cluster = USECLUSTER
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
    to_cluster = USECLUSTER
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
            ".alignstats")
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

    to_cluster = USECLUSTER

    job_options = "-l mem_free=32G"

    statement = '''
    cat %(infile)s 
    | python %(scriptsdir)s/bam2bed.py
          --merge-pairs=%(medips_fragment_length)i
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
@transform( prepareBAMs, suffix(".bam"), ".covered.bed.gz" )
def buildCoverageBed( infile, outfile ):
    '''build bed file with regions covered by reads.

    Intervals containing only a single read are removed.
    '''
    
    to_cluster = USECLUSTER

    statement = '''
    cat %(infile)s 
    | python %(scriptsdir)s/bam2bed.py
          --merge-pairs=%(medips_fragment_length)i
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
@merge( buildCoverageBed, "tiles_variable_width.bed.gz" )
def buildVariableWidthTiles( infiles, outfile ):
    '''bed file with intervals that are covered by reads in any of the experiments.
    '''
    
    infiles = " ".join( infiles )
    to_cluster = USECLUSTER

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
@files(PARAMS["samtools_genome"]+".fai", "tiles_fixed_width.bed.gz" )
def buildFixedWidthTiles( infile, outfile ):
    '''Build bed file segmenting entire genome using window x and shift y'''

    statement = '''python %(scriptsdir)s/genome_bed.py
                      -g %(infile)s
                      -w %(deseq_window)s
                      -s %(deseq_shift)s
                | gzip
                > %(outfile)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( (buildVariableWidthTiles, buildFixedWidthTiles), 
            suffix(".bed.gz"), 
            ".bigbed")
def buildBigBed( infile, outfile ):
    '''bed file with intervals that are covered by reads in any of the experiments.
    '''
    
    to_cluster = USECLUSTER
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

    to_cluster = USECLUSTER

    infile, tiles = infiles

    E.warn( "truncated to chr22!!!")

    # note: needs to set flags appropriately for
    # single-end/paired-end data sets
    # set filter options
    # for example, only properly paired reads
    paired = True
    if paired:
        flag_filter = "-f 0x2"
    else:
        flag_filter = ""

    statement = '''samtools view -b %(flag_filter)s -q %(deseq_min_mapping_quality)s %(infile)s chr22
                   | coverageBed -abam stdin -b %(tiles)s 
                   | sort -k1,1 -k2,2n
                   | grep chr22
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
            add_inputs( buildFixedWidthTiles ),
            r".fixedwidth.tilecounts.bed.gz" )
def buildTiledReadCountsFixedWidth(infiles, outfile ):
    '''build read counds for fixed width windows.'''
    buildTiledReadCounts( infiles, outfile )

#########################################################################
@follows( mkdir( "deseq" ) )
@collate( (buildTiledReadCountsVariableWidth,
           buildTiledReadCountsFixedWidth), 
          regex( ".*\.([^.]+)width.tilecounts.bed.gz"),
          r"deseq/\1.counts.tsv.gz")
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

    This method uses the maximum number of reads
    found in any interval as the tag count.
    '''
    
    to_cluster = USECLUSTER

    src = " ".join( [ '''<( zcat %s | awk '{printf("%%s:%%i-%%i\\t%%i\\n", $1,$2,$3,$4 );}' ) ''' % x for x in infiles] )
    tmpfile = P.getTempFilename( "." )
    statement = '''paste %(src)s > %(tmpfile)s'''
    P.run()
    
    tracks = [ re.sub( "\..*", '', os.path.basename(x) ) for x in infiles ]

    outf = IOTools.openFile( outfile, "w")
    outf.write( "interval\t%s\n" % "\t".join( tracks ) )
    
    for line in open( tmpfile, "r" ):
        data = line[:-1].split("\t")
        genes = list(set([ data[x] for x in range(0,len(data), 2 ) ]))
        values = [ data[x] for x in range(1,len(data), 2 ) ]
        assert len(genes) == 1, "paste command failed, wrong number of genes per line: '%s'" % line
        outf.write( "%s\t%s\n" % (genes[0], "\t".join(map(str, values) ) ) )
    
    outf.close()

    os.unlink(tmpfile)

#########################################################################
@follows( aggregateTiledReadCounts)
@files( [ ( (data, design), 
            "deseq/%s_%s.deseq" % (P.snip(os.path.basename(data),".counts.tsv.gz"),
                                   P.snip(os.path.basename(design),".tsv" ) ) ) \
              for data, design in itertools.product( 
                                               glob.glob("deseq/*.counts.tsv.gz"),
                                               P.asList(PARAMS["deseq_designs"]) ) ] )
def runDESeq( infiles, outfile ):
    '''estimate differential expression using DESeq.

    The final output is a table. It is slightly edited such that
    it contains a similar output and similar fdr compared to cuffdiff.
    '''

    infile, comp = infiles
    to_cluster = USECLUSTER
    outdir = "deseq"
    deseq_fdr = PARAMS["deseq_fdr"]

    # load library 
    R('''suppressMessages(library('DESeq'))''')

    # load count data
    E.info( "loading data" )

    R( '''counts_table = read.delim( '%(infile)s', header = TRUE, row.names = 1, stringsAsFactors = TRUE )''' % locals() )
    E.info( "read data: %i observations for %i samples" % tuple(R('''dim(counts_table)''')))

    # Remove windows with no data
    R( '''max_counts = apply(counts_table,1,max)''' )
    R( '''counts_table = counts_table[max_counts>0,]''')
    E.info( "removed %i empty columns" % tuple( R('''sum(max_counts == 0)''') ) )
    E.info( "trimmed data: %i observations for %i samples" % tuple( R('''dim(counts_table)''') ) )

    # Load comparisons from file
    comp_name = P.snip( os.path.basename(comp), ".tsv")
    R('''pheno = read.delim( '%(comp)s', header = TRUE, stringsAsFactors = TRUE )''' % locals() )

    # Make sample names R-like - substitute - for . and add the .prep suffix
    R('''pheno[,1] = gsub('-', '.', pheno[,1]) ''')
    
    # Ensure pheno rows match count columns 
    R('''pheno2 = pheno[match(colnames(counts_table),pheno[,1]),,drop=FALSE]''' ) 

    # Subset data & set conditions
    R('''includedSamples <- pheno2$include == '1' ''')
    R('''countsTable <- counts_table[ , includedSamples ]''')
    R('''conds <- pheno2$group[ includedSamples ]''')

    # Test if replicates exist
    min_reps = R('''min(table(conds)) ''')
    no_replicates = False
    if min_reps < 2:
        no_replicates = True

    ######## Run DESeq
    # Create Count data object
    E.info( "running DESeq" )
    R('''cds <-newCountDataSet( countsTable, conds) ''')

    # Estimate size factors
    R('''cds <- estimateSizeFactors( cds )''')

    # Estimate variance
    if no_replicates:
        R('''cds <- estimateVarianceFunctions( cds, method="blind" )''')
    else:
        R('''cds <- estimateVarianceFunctions( cds )''')

    # Plot scvplot
    size_factors = R('''sizeFactors( cds )''')
    R.png( '''%(outdir)s/%(comp_name)s_scvplot.png''' % locals() )
    R('''scvPlot( cds, ylim = c(0,3))''')
    R['dev.off']()

    # Generate heatmap of variance stabilised data
    R('''vsd <- getVarianceStabilizedData( cds )''' )
    R('''dists <- dist( t( vsd ) )''')
    R.png( '''%(outdir)s/%(comp_name)s_heatmap.png''' % locals() )
    R('''heatmap( as.matrix( dists ), symm=TRUE )''' )
    R['dev.off']()

    conditions = R('''levels(conds)''')
    for condition in conditions:
        if not no_replicates:
            R.png( '''%(outdir)s/%(comp_name)s_%(condition)s_fit.png''' % locals() )
            R('''diagForT <- varianceFitDiagnostics( cds, "%s" )''' % condition )
            R('''smoothScatter( log10(diagForT$baseMean), log10(diagForT$baseVar) )''')
            R('''lines( log10(fittedBaseVar) ~ log10(baseMean), diagForT[ order(diagForT$baseMean), ], col="red" )''')
            R['dev.off']()
            R.png( '''%(outdir)s/%(comp_name)s_%(condition)s_residuals.png''' % locals()  )
            R('''residualsEcdfPlot( cds, "%s" )''' % condition )
            R['dev.off']()

    # Differential expression
    L.info("calling differential expression")
    R('''res <- nbinomTest( cds, '%s', '%s' )''' % (conditions[0],conditions[1]))

    # Plot significance
    R.png( '''%(outdir)s/%(comp_name)s_significance.png''' % locals() )
    R('''plot( res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res$padj < %(deseq_fdr)s, "red", "black" ) )''' % locals() )
    R['dev.off']()

    outf = IOTools.openFile( outfile, "w" )
    isna = R["is.na"]

    L.info("Generating output")
    # Get column names from output and edit
    names = None
    if not names:
        names = list(R['res'].names)
        m = dict( [ (x,x) for x in names ])
        m.update( dict(
                pval = "pvalue", 
                baseMeanA = "value1", 
                baseMeanB = "value2",
                id = "test_id", 
                log2FoldChange = "lfold") )
        
        header = [ m[x] for x in names ] 
        outf.write( "Group1\tGroup2\t%s\tstatus\tsignificant\n" % "\t".join(header))
    else:
        if names != list(R['res'].names):
            raise ValueError( "different column headers in DESeq output: %s vs %s" % (names, list(R['res'].names)))

    # Parse results and parse to file
    rtype = collections.namedtuple( "rtype", names )
    for data in zip( *R['res']) :
        d = rtype._make( data )
        outf.write( "%s\t%s\t" % (conditions[0],conditions[1]))
        # set significant flag
        if d.padj <= deseq_fdr: signif = 1
        else: signif = 0

        # set lfold change to 0 if both are not expressed
        if d.baseMeanA == 0.0 and d.baseMeanB == 0.0:
            d = d._replace( foldChange = 0, log2FoldChange = 0 )

        if isna( d.pval ): status = "OK"
        else: status = "FAIL"

        outf.write( "\t".join( map(str, d) ))
        outf.write("\t%s\t%s\n" % (status, str(signif)))
            
    outf.close()

#########################################################################
@transform( runDESeq, suffix(".deseq"), "_deseq.load" )
def loadDESeq( infile, outfile ):
    '''load differential expression results.'''
    P.load( infile, outfile, "--index=group1 --index=group2 --index=test_id --allow-empty" )

#########################################################################
#########################################################################
#########################################################################
def buildExpressionStats( tables, method, outfile ):
    '''build expression summary statistics.
    
    Creates some diagnostic plots in

    <exportdir>/<method> directory.
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    def togeneset( tablename ):
        return re.match("([^_]+)_", tablename ).groups()[0]

    keys_status = "OK", "NOTEST", "FAIL", "NOCALL"

    outf = IOTools.openFile( outfile, "w" )
    outf.write( "\t".join( ("geneset", "level", "track1", "track2", "tested",
                            "\t".join( [ "status_%s" % x for x in keys_status ] ),
                            "significant",
                            "twofold" ) ) + "\n" )

    all_tables = set(Database.getTables( dbhandle ))
    outdir = os.path.join( PARAMS["exportdir"], method )

    for level in CUFFDIFF_LEVELS:

        for tablename in tables:

            tablename_diff = "%s_%s_diff" % (tablename, level)
            tablename_levels = "%s_%s_diff" % (tablename, level )
            geneset = togeneset( tablename_diff )
            if tablename_diff not in all_tables: continue

            def toDict( vals, l = 2 ):
                return collections.defaultdict( int, [ (tuple( x[:l]), x[l]) for x in vals ] )
            
            tested = toDict( Database.executewait( dbhandle,
                                               """SELECT track1, track2, COUNT(*) FROM %(tablename_diff)s 
                                    GROUP BY track1,track2""" % locals() ).fetchall() )
            status = toDict( Database.executewait( dbhandle,
                                                   """SELECT track1, track2, status, COUNT(*) FROM %(tablename_diff)s 
                                    GROUP BY track1,track2,status""" % locals() ).fetchall(), 3 )
            signif = toDict( Database.executewait( dbhandle,
                                                   """SELECT track1, track2, COUNT(*) FROM %(tablename_diff)s 
                                    WHERE significant
                                    GROUP BY track1,track2""" % locals() ).fetchall() )
            fold2 = toDict( Database.executewait( dbhandle,
                    """SELECT track1, track2, COUNT(*) FROM %(tablename_diff)s 
                                    WHERE (lfold >= 1 or lfold <= -1) AND significant
                                    GROUP BY track1,track2,significant""" % locals() ).fetchall() )
            
            for track1, track2 in itertools.combinations( EXPERIMENTS, 2 ):
                outf.write( "\t".join(map(str, (
                                geneset,
                                level,
                                track1,
                                track2,
                                tested[(track1,track2)],
                                "\t".join( [ str(status[(track1,track2,x)]) for x in keys_status]),
                                signif[(track1,track2)],
                                fold2[(track1,track2)] ) ) ) + "\n" )
                
            ###########################################
            ###########################################
            ###########################################
            # plot length versus P-Value
            data = Database.executewait( dbhandle, 
                                         '''SELECT i.sum, pvalue 
                                 FROM %(tablename_diff)s, 
                                 %(geneset)s_geneinfo as i 
                                 WHERE i.gene_id = test_id AND significant'''% locals() ).fetchall()

            # require at least 10 datapoints - otherwise smooth scatter fails
            if len(data) > 10:
                data = zip(*data)

                pngfile = "%(outdir)s/%(geneset)s_%(method)s_%(level)s_pvalue_vs_length.png" % locals()
                R.png( pngfile )
                R.smoothScatter( R.log10( ro.FloatVector(data[0]) ),
                                 R.log10( ro.FloatVector(data[1]) ),
                                 xlab = 'log10( length )',
                                 ylab = 'log10( pvalue )',
                                 log="x", pch=20, cex=.1 )

                R['dev.off']()

    outf.close()

#########################################################################
@merge( loadDESeq, "deseq_stats.tsv" )
def buildDESeqStats( infiles, outfile ):
    tablenames = [P.toTable( x ) for x in infiles ] 
    buildExpressionStats( tablenames, "deseq", outfile )

#########################################################################
@transform( buildDESeqStats, suffix(".tsv"), ".load" )
def loadDESeqStats( infile, outfile ):
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@follows( loadPicardDuplicateStats,
          loadPicardAlignmentStats,
          loadPicardGCStats,
          loadBAMStats )
def mapping(): pass

@follows( aggregateTiledReadCounts,
          runDESeq,
          loadDESeq,
          buildDESeqStats,
          loadDESeqStats )
def callDMRs(): pass

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



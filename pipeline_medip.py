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

The exome pipeline imports unmapped reads from one or more 
fastq or sra files and aligns them to the genome, then filters calls variants (SNVs and indels) 
and filters them by both depth/rate and regions of interest.

   1. Align to genome using gapped alignment (BWA)
   2. Calculate alignment and coverage statistics (BAMStats)
   3. Identify differentially methylated regions (DMRs)
   4. Filter DMRs
   5. Calculate DMR statistics
   6. Produce report (SphinxReport)

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
from rpy2.robjects import r as R

import Experiment as E
import logging as L
import Database
import sys, os, re, shutil, itertools, math, glob, time, gzip, collections, random
import numpy, sqlite3
import GTF, IOTools, IndexedFasta
import Tophat
import rpy2.robjects as ro
import PipelineGeneset
import PipelineMapping
import Stats
import PipelineTracks
import Pipeline as P

USECLUSTER = True

#########################################################################
#########################################################################
#########################################################################
# load options from the config file
P.getParameters( ["%s.ini" % __file__[:-len(".py")], "../medip.ini", "medip.ini" ] )
PARAMS = P.PARAMS

#########################################################################
#########################################################################
#########################################################################
@transform( ("*.fastq.1.gz", 
             "*.fastq.gz",
             "*.sra"),
             regex( r"(\S+).(fastq.1.gz|fastq.gz|sra)"),
             r"\1/bam/\1.bam")
def mapReads(infiles, outfile):
    '''Map reads to the genome using BWA '''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(outfile), ".bam" )
    try: os.mkdir( track )
    except OSError: pass
    try: os.mkdir( '''%(track)s/bam''' % locals() )
    except OSError: pass
    m = PipelineMapping.bwa()
    statement = m.build((infiles,), outfile) 
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( mapReads,
            regex( r"(\S+)/bam/(\S+).bam"),
            r"\1/bam/\2.dedup.bam")
def dedup(infiles, outfile):
        '''Remove duplicate alignments from BAM files.'''
        to_cluster = USECLUSTER
        track = P.snip( outfile, ".bam" )
        dedup_method = PARAMS["dedup_method"]
        if dedup_method == 'samtools':
            statement = '''samtools rmdup %(infiles)s %(outfile)s; ''' % locals()    
        elif dedup_method == 'picard':
            statement = '''MarkDuplicates INPUT=%(infiles)s  ASSUME_SORTED=true OUTPUT=%(outfile)s METRICS_FILE=%(track)s.dupstats VALIDATION_STRINGENCY=SILENT; ''' % locals()
        statement += '''samtools index %(outfile)s; ''' % locals()
        #print statement
        P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( dedup, "picard_duplicate_stats.load" )
def loadPicardDuplicateStats( infiles, outfile ):
    '''Merge Picard duplicate stats into single table and load into SQLite.'''

    tablename = P.toTable( outfile )

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

    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
               '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( dedup, 
            regex( r"(\S+)/bam/(\S+).bam"),
            r"\1/bam/\2.alignstats" )
def buildPicardAlignStats( infile, outfile ):
    '''Gather BAM file alignment statistics using Picard '''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".bam" )
    statement = '''CollectAlignmentSummaryMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(bwa_index_dir)s/%%(genome)s.fa ASSUME_SORTED=true OUTPUT=%(outfile)s VALIDATION_STRINGENCY=SILENT ''' % locals()
    P.run()

############################################################
############################################################
############################################################
@merge( buildPicardAlignStats, "picard_align_stats.load" )
def loadPicardAlignStats( infiles, outfile ):
    '''Merge Picard alignment stats into single table and load into SQLite.'''

    tablename = P.toTable( outfile )

    outf = P.getTempFile()

    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".dedup.alignstats" )
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

    statement = '''cat %(tmpfilename)s
                | python %(scriptsdir)s/csv2db.py
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s
               '''
    P.run()

    os.unlink( tmpfilename )

#########################################################################
#########################################################################
#########################################################################
@transform( dedup, 
            regex( r"(\S+)/bam/(\S+).bam"),
            r"\1/bam/\2.isizestats" )
def buildPicardInsertSizeStats( infile, outfile ):
    '''Gather BAM file insert size statistics using Picard '''
    to_cluster = USECLUSTER
    track = P.snip( os.path.basename(infile), ".bam" )
    statement = '''CollectInsertSizeMetrics INPUT=%(infile)s REFERENCE_SEQUENCE=%%(bwa_index_dir)s/%%(genome)s.fa ASSUME_SORTED=true OUTPUT=%(outfile)s HISTOGRAM_FILE=%(outfile)s.pdf VALIDATION_STRINGENCY=SILENT ''' % locals()
    P.run()

############################################################
############################################################
############################################################
@merge( buildPicardInsertSizeStats, "picard_isize_stats.load" )
def loadPicardInsertSizeStats( infiles, outfile ):
    '''Merge Picard insert size stats into single table and load into SQLite.'''

    tablename = P.toTable( outfile )
    outf = P.getTempFile()

    first = True
    for f in infiles:
        track = P.snip( os.path.basename(f), ".dedup.isizestats" )
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
                > %(outfile)s
               '''
    P.run()

    os.unlink( tmpfilename )

#########################################################################
#########################################################################
#########################################################################
@transform( dedup, 
            regex(r"(\S+)/bam/(\S+).bam"),
            r"\1/bam/\2.readstats" )
def buildBAMStats( infile, outfile ):
    '''Count number of reads mapped, duplicates, etc. '''
    to_cluster = USECLUSTER
    scriptsdir = PARAMS["general_scriptsdir"]
    statement = '''python %(scriptsdir)s/bam2stats.py --force 
                   --output-filename-pattern=%(outfile)s.%%s < %(infile)s > %(outfile)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( buildBAMStats, "bam_stats.load" )
def loadBAMStats( infiles, outfile ):
    '''Import bam statistics into SQLite'''

    scriptsdir = PARAMS["general_scriptsdir"]
    header = ",".join( [P.snip( os.path.basename(x), ".dedup.readstats") for x in infiles] )
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
                      --allow-empty
                      --index=track
                      --table=%(tablename)s 
                > %(outfile)s"""
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
                      --table=%(tname)s 
                      --allow-empty
                >> %(outfile)s """
        P.run()

#########################################################################
#########################################################################
#########################################################################
#@files(r"%s/%s.fa.fai" % (PARAMS("genome_dir"), PARAMS("genome")), "%(genome)s_tiling.bed")
def buildTiledGenomeBed( infile, outfile ):
    '''Build bed file segmenting entire genome using window x and shift y'''

    statement = '''python %(scriptsdir)s/genome_bed.py
                      -g %(genome_dir)s/%(genome)s.fa.fai
                      -w %%(deseq_window)s
                      -s %%(deseq_shift)s
                      -o %(outfile)s '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@follows(buildTiledGenomeBed)
@transform( dedup, regex(r"(\S+)/bam/(\S+).bam"), r"\1/bam/\2.counts.bed.gz" )
def buildTiledReadCounts( infiles, outfile ):
    '''compute coverage of genome with reads.'''

    infile, genome_bed = infiles

    to_cluster = USECLUSTER

    # note: needs to set flags appropriately for
    # single-end/paired-end data sets
    # set filter options
    # for example, only properly paired reads
    paired = False
    if paired:
        flag_filter = "-f 0x2"
    else:
        flag_filter = ""

    # note: the -split option only concerns the stream in A - multiple
    # segments in B are not split. Hence counting has to proceed via
    # single exons - this can lead to double counting if exon counts
    # are later aggregated.

    statement = '''samtools view -b %(flag_filter)s -q %(deseq_min_mapping_quality)s %(infile)s
                   | coverageBed -abam stdin -b %(genome_bed)s -split
                   | sort -k1,1 -k2,2n
                   | gzip
                   > %(outfile)s '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@collate(buildTiledReadCounts,
         regex(r"(\S+)/bam/(\S+).counts.bed.gz"),  
         r"deseq/\2.counts.tsv.gz")
def aggregateTiledReadCounts( infiles, outfile ):
    '''aggregate tag counts for each window.

    coverageBed adds the following four columns:

    1) The number of features in A that overlapped (by at least one base pair) the B interval.
    2) The number of bases in B that had non-zero coverage from features in A.
    3) The length of the entry in B.
    4) The fraction of bases in B that had non-zero coverage from features in A.

    For bed6: use column 7
    For bed12: use column 13

    This method uses the maximum number of reads
    found in any exon as the tag count.
    '''
    
    to_cluster = USECLUSTER

    # aggregate not necessary for bed12 files, but kept in
    src = " ".join( [ "<( zcat %s | sort -k4,4 | groupBy -i stdin -g 4 -c 7 -o max | sort -k1,1)" % x for x in infiles ] )

    tmpfile = P.getTempFilename( "." )
    
    statement = '''paste %(src)s > %(tmpfile)s'''
    P.run()

    tracks = [P.snip(x, ".counts.bed.gz" ) for x in infiles ]
    #tracks = [re.match( "exon_counts.dir/(\S+)_vs.*", x).groups()[0] for x in tracks ]

    outf = IOTools.openFile( outfile, "w")
    outf.write( "gene_id\t%s\n" % "\t".join( tracks ) )
    
    for line in open( tmpfile, "r" ):
        data = line[:-1].split("\t")
        genes = list(set([ data[x] for x in range(0,len(data), 2 ) ]))
        values = [ data[x] for x in range(1,len(data), 2 ) ]
        assert len(genes) == 1, "paste command failed, wrong number of genes per line"
        outf.write( "%s\t%s\n" % (genes[0], "\t".join(map(str, values) ) ) )
    
    outf.close()

    os.unlink( tmpfile )

#########################################################################
#########################################################################
#########################################################################
@follows( mkdir( os.path.join( PARAMS["exportdir"], "deseq" ) ) )
@transform( buildTiledReadCounts, suffix(".genome.tsv.gz"), ".deseq")
def runDESeq( infile, outfile ):
    '''estimate differential expression using DESeq.

    The final output is a table. It is slightly edited such that
    it contains a similar output and similar fdr compared to cuffdiff.

    Plots are:

    <geneset>_<method>_<level>_<track1>_vs_<track2>_significance.png
        fold change against expression level
    '''
    
    to_cluster = USECLUSTER

    outdir = os.path.join( PARAMS["exportdir"], "deseq" )

    geneset, method = outfile.split(".")
    level = "gene"


    # load data 
    R('''suppressMessages(library('DESeq'))''')
    R( '''counts_table <- read.delim( '%s', header = TRUE, row.names = 1, stringsAsFactors = TRUE )''' % infile )

    # get conditions to test
    # note that tracks in R use a '.' as separator
    tracks = R('''colnames(counts_table)''')
    map_track2column = dict( [ (y,x) for x,y in enumerate( tracks ) ] )
    
    sample2condition = [None] * len(tracks)
    conditions = []
    no_replicates = False
    for group, replicates in EXPERIMENTS.iteritems():
        if len(replicates) == 1:
            E.warn( "only one replicate in %s - replicates will be ignored in ALL data sets for variance estimation" )
            no_replicates = True

        for r in replicates:
            sample2condition[map_track2column[r.asR()]] = group.asR()
        conditions.append( group )

    ro.globalenv['conds'] = ro.StrVector(sample2condition)

    def build_filename2( **kwargs ):
        return "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(track1)s_vs_%(track2)s_%(section)s.png" % kwargs
    def build_filename1( **kwargs ):
        return "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(section)s_%(track)s.png" % kwargs
    def build_filename0( **kwargs ):
        return "%(outdir)s/%(geneset)s_%(method)s_%(level)s_%(section)s.png" % kwargs

    # this analysis follows the 'Analysing RNA-Seq data with the "DESeq" package'
    # tutorial 
    R('''cds <-newCountDataSet( counts_table, conds) ''')
    R('''cds <- estimateSizeFactors( cds )''')

    if no_replicates:
        R('''cds <- estimateVarianceFunctions( cds, method="blind" )''')
    else:
        R('''cds <- estimateVarianceFunctions( cds )''')

    L.info("creating diagnostic plots" ) 
    size_factors = R('''sizeFactors( cds )''')
    R.png( build_filename0( section = "scvplot", **locals() ) )
    R('''scvPlot( cds, ylim = c(0,3))''')
    R['dev.off']()

    R('''vsd <- getVarianceStabilizedData( cds )''' )
    R('''dists <- dist( t( vsd ) )''')
    R.png( build_filename0( section = "heatmap", **locals() ) )
    R('''heatmap( as.matrix( dists ), symm=TRUE )''' )
    R['dev.off']()

    for track in conditions:
        condition = track.asR()
        R.png( build_filename1( section = "fit", **locals() ) )
        R('''diagForT <- varianceFitDiagnostics( cds, "%s" )''' % condition )
        if not no_replicates:
            R('''smoothScatter( log10(diagForT$baseMean), log10(diagForT$baseVar) )''')
            R('''lines( log10(fittedBaseVar) ~ log10(baseMean), diagForT[ order(diagForT$baseMean), ], col="red" )''')
            R['dev.off']()
            R.png( build_filename1( section = "residuals", **locals() ) )
            R('''residualsEcdfPlot( cds, "%s" )''' % condition )
            R['dev.off']()

    L.info("calling differential expression")

    outf = IOTools.openFile( outfile, "w" )
    names = None
    fdr = PARAMS["cuffdiff_fdr"]
    isna = R["is.na"]

    for track1, track2 in itertools.combinations( conditions, 2 ):
        R('''res <- nbinomTest( cds, '%s', '%s' )''' % (track1.asR(),track2.asR()))

        R.png( build_filename2( section = "significance", **locals() ) )
        R('''plot( res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1,
                   col = ifelse( res$padj < %(cuffdiff_fdr)s, "red", "black" ) )''' % PARAMS )
        R['dev.off']()
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
            outf.write( "track1\ttrack2\t%s\tstatus\tsignificant\n" % "\t".join(header))
        else:
            if names != list(R['res'].names):
                raise ValueError( "different column headers in DESeq output: %s vs %s" % (names, list(R['res'].names)))

        rtype = collections.namedtuple( "rtype", names )
        
        for data in zip( *R['res']) :
            d = rtype._make( data )
            outf.write( "%s\t%s\t" % (track1,track2))
            # set significant flag
            if d.padj <= fdr: signif = 1
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
#########################################################################
#########################################################################
@transform( runDESeq, suffix(".deseq"), "_deseq.load" )
def loadDESeq( infile, outfile ):
    '''load differential expression results.'''
    # add gene level follow convention "<level>_diff"
    
    # if one expression value is 0, the log fc is inf or -inf.
    # substitute with 10

    tablename = P.snip( outfile, ".load") + "_gene_diff" 
    statement = '''cat %(infile)s
            | python %(scriptsdir)s/csv2db.py %(csv2db_options)s
              --allow-empty
              --index=track1
              --index=track2
              --index=test_id
              --table=%(tablename)s 
            > %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge( loadDESeq, "deseq_stats.tsv" )
def buildDESeqStats( infiles, outfile ):
    tablenames = [P.toTable( x ) for x in infiles ] 
    buildExpressionStats( tablenames, "deseq", outfile )

#########################################################################
#########################################################################
#########################################################################
@transform( buildDESeqStats, suffix(".tsv"), ".load" )
def loadDESeqStats( infile, outfile ):
    P.load( infile, outfile )

#########################################################################
#########################################################################
#########################################################################
@follows( mapReads,
          dedup,
          loadPicardDuplicateStats,
          buildPicardAlignStats,
          loadPicardAlignStats,
          buildPicardInsertSizeStats,
          loadPicardInsertSizeStats,
          buildBAMStats,
          loadBAMStats )
def full(): pass


@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "Starting documentation build process from scratch" )
    dirname, basenme = os.path.split( os.path.abspath( __file__ ) )
    docdir = os.path.join( dirname, "pipeline_docs", P.snip( basenme, ".py" ) )

    # requires libtk, which is not present on the nodes
    to_cluster = True
    job_options= "-pe dedicated %i -R y" % PARAMS["report_threads"]
    statement = '''rm -rf report _cache _static;
                   sphinxreport-build 
                       --num-jobs=%(report_threads)s
                   sphinx-build 
                       -b html 
                       -d %(report_doctrees)s
                       -c . 
                   %(docdir)s %(report_html)s
                   > report.log '''
    P.run()

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )



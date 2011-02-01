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

:Author: Andreas Heger, Tildon Grant Belgard
:Release: $Id$
:Date: |today|
:Tags: Python

====================
RNA-Seq pipeline
====================

The RNA-Seq pipeline imports unmapped reads from one or more
RNA-Seq experiments. It assumes the data derive from multiple 
tissues/conditions (:term:`experiment`) with one or more 
biological and/or technical replicates (:term:`replicate`).

The pipeline performs the following tasks:

   * analyse each experiment:
      * for each replicate
          * map reads using tophat for each term:`replicate` separately. 
          * predict splice isoforms and expression levels with :term:`cufflinks`.
          * estimate expression levels of reference gene set with :term:`cufflinks`.
          * annotate isoforms in replicates with genomic annotations
      * compare isoforms in replicates within each :term:`experiment` (:term:`cuffcompare`) 
          and to reference gene set.
      * TODO: summary statistics on reproducibily within each experiment
   * build a combined gene set including the reference gene set and isoforms predicted by :term:`cufflinks`.
      * compare all isoforms in all experiments+isoforms (:term:`cuffcompare`) to each other 
         and the reference gene set
      * TODO: summary statistics on isoforms with respect to gene set
   * estimate differential expression levels of transcripts
      * different gene sets
         * reference gene set
         * combined gene set
      * different methods
         * :term:`DESeg` (tag counting)
            * using FPKM
            * using tag-counting in union-intersection genes
         * :term:`cuffdiff``
      * TODO: summary statistics on differential expression
   

Usage
-----

Configuration
=============

Input
=====

Reads
-----

Reads are sorted and indexed :term:`BAM` formatted files labeled in the following way::

   sample-condition-replicate_export.txt.gz

If there are several lanes for a single replicate, merge the input files first.

Note that neither ``sample``, ``condition`` or ``replicate`` should contain 
``_`` (underscore) and ``.`` (dot) characters as these are used by the pipeline
to delineate tasks.

The pipeline expects paired end results.

Optional inputs
---------------

Notes
-----


Overview
--------

The pipeline has the following major compononts:

1. Compute mean insert size distribution from data.

2. Discover splice junctions (combineJunctions)

4. build and compare gene models

5. calculate RPKM values

6. Analyze samples
   cluster by expression similarity
   define groups of co-regulated genes

Code
----

"""

# load modules
from ruffus import *

import Experiment as E

import sys, os, re, shutil, itertools, math, glob, logging, time

import numpy, sqlite3
import GTF, IOTools
import Tophat
from rpy2.robjects import r as R
import rpy2.robjects as ro

import pipeline_chipseq_intervals as PIntervals

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file
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
## Helper functions mapping tracks to conditions, etc
###################################################################
import PipelineTracks

TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
            glob.glob( "*.sra" ), "(\S+).sra" )

ALL = PipelineTracks.Aggregate( TRACKS )
EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", "tissue" ) )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, labels = ("condition", ) )
TISSUES = PipelineTracks.Aggregate( TRACKS, labels = ("tissue", ) )

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("conf.py"): 
    E.info( "reading additional configuration from conf.py" )
    execfile("conf.py")

USECLUSTER=True

#########################################################################
#########################################################################
#########################################################################
##
#########################################################################
@transform( "*.sra", suffix(".sra"), ".bam" )
def mapReadsFromSRAWithTophat( infile, outfile ):
    '''map reads from short read archive sequences
    using tophat
    
    not paired-ended.

    TODO: unpack reads on scratch disk

    There is a bug in tophat 1.1.4, see:
    http://seqanswers.com/forums/showthread.php?t=8103
    and 
    http://seqanswers.com/forums/showthread.php?t=8817

    Applying the suggested fix:
    
    cd xxx.dir

    # Fix the SO:sorted header
    sed s/sorted/unsorted/ tmp/accepted_hits.sam > fixed.sam

    # Use picard to fix clipping
    java -jar /ifs/apps/bio/picard-tools-1.38/CleanSam.jar INPUT=fixed.sam OUTPUT=fixed2.sam

    # Convert to bam
    java -jar /ifs/apps/bio/picard-tools-1.38/SortSam.jar INPUT=fixed2.sam OUTPUT=accepted_hits.bam SORT_ORDER=coordinate 
    rm -f fixed.sam fixed2.sam
    '''

    tmpfilename = outfile + ".dir" #P.getTempFilename()
    
    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )

    # job_options= "-pe dedicated 4-8 -l mem_free=3G -R y"

    prefix = P.snip( infile, ".sra" )

    # tophat does a seek operation on the fq files, hence they
    # need to be unpacked into uncompressed real files

    # todo: use scratch dir for fastq files

    statement = '''
    fastq-dump %(infile)s ;
    tophat --output-dir %(tmpfilename)s
           --quiet
           %(tophat_options)s
           %(bowtie_index_dir)s/%(genome)s
           %(prefix)s.fastq
    >& %(outfile)s.log;
    mv %(tmpfilename)s/accepted_hits.bam %(outfile)s; 
    mv %(tmpfilename)s/tmp/junctions.bed %(prefix)s.junctions.bed; 
    rm -f %(prefix)s.fastq; 
    samtools index %(outfile)s
    '''

    # rm -rf %(tmpfilename)s

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform( "*.sra", suffix(".sra"), ".bam" )
def mapReadsFromSRAWithBowtie( infile, outfile ):
    '''map reads from short read archive sequences
    using bowtie
    
    TODO: unpack reads on scratch disk

    '''

    to_cluster = USECLUSTER

    tmpfilename = outfile + ".dir" #P.getTempFilename()
    
    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )

    # job_options= "-pe dedicated 4-8 -l mem_free=3G -R y"

    prefix = P.snip( infile, ".sra" )

    # tophat does a seek operation on the fq files, hence they
    # need to be unpacked into uncompressed real files

    # todo: use scratch dir for fastq files
    tmpfilename = P.getTempFilename()

    statement = '''
    fastq-dump %(infile)s ;
    bowtie --quiet --sam
           %(bowtie_options)s 
           %(bowtie_index_dir)s/%(genome)s 
           %(prefix)s.fastq 2>%(outfile)s.log
    | samtools import %(cufflinks_genome_dir)s/%(genome)s.fa - %(tmpfilename)s 1>&2 2>> %(outfile)s.log;
    samtools sort %(tmpfilename)s %(prefix)s;
    samtools index %(outfile)s;
    rm -f %(tmpfilename)s;
    rm -f %(prefix)s.fastq; 
    samtools index %(outfile)s
    '''

    P.run()
    
@transform( mapReadsFromSRAWithTophat, suffix(".bam"), ".logs.gz" )
def collectLogFiles( infile, outfile ):

    tmpfilename = infile + ".dir" 
    
    outf = IOTools.openFile( outfile, "w" )
    for f in glob.glob( os.path.join( tmpfilename, "logs", "*.log" ) ):
        outf.write( "##> %s\n" % f )
        with open( f, "r" ) as inf: 
            for line in inf:
                # segment_juncs.log contains a lot of putative junctions
                # remove these from the log file
                if line.startswith( "chr"): continue
                outf.write( line )

    outf.close()

############################################################
############################################################
############################################################
@transform( "*.bam",
            suffix(".bam"),
            ".readstats" )
def buildBAMStats( infile, outfile ):
    '''count number of reads mapped, duplicates, etc.'''
    PIntervals.buildBAMStats( infile, outfile )
    
#########################################################################
#########################################################################
#########################################################################
@merge( PARAMS_ANNOTATIONS["ensembl_filename_gtf"],
        "transcripts.gtf.gz" )
def buildTranscripts( infile, outfile ):
    '''sanitize transcripts file for cufflinks analysis.

    Merge exons separated by small introns.

    Transcript will be ignored with:
       * very long introns (max_intron_size) (otherwise, cufflinks complains)
       * known rRNA genes
    '''
    max_intron_size =  PARAMS["max_intron_size"]

    c = E.Counter()

    outf = IOTools.openFile( outfile , "w" )

    for all_exons in GTF.transcript_iterator( GTF.iterator( IOTools.openFile( infile )) ):

        c.info += 1

        is_ok = True

        # merge exons and cds
        all_exons.sort( key = lambda x: x.feature )
        new_exons = []

        for feature, exons in itertools.groupby( all_exons, lambda x: x.feature ):

            tmp = sorted( list( exons ), key = lambda x: x.start )
            
            l, n = tmp[0], []

            for e in tmp[1:]:
                d = e.start - l.end
                if d > max_intron_size:
                    is_ok = False
                    break
                elif d < 5:
                    l.end = max(e.end, l.end)
                    c.merged += 1
                    continue

                n.append( l )
                l = e

            n.append( l )
            new_exons.extend( n )

            if not is_ok: break

        if not is_ok: 
            E.info( "removing transcript %s" % all_exons[0].transcript_id )
            c.skipped += 1
            continue

        new_exons.sort( key = lambda x: x.start )

        for e in new_exons:
            # add chr prefix 
            outf.write( "chr%s\n" % str(e) )
            c.exons += 1

        c.output += 1

    E.info( "%s\n" % str(c) )

#########################################################################
#########################################################################
#########################################################################
@merge( buildTranscripts, 
        "unique_exons.bed.gz" )
def buildUniqueExons( infile, outfile ):
    '''build union/intersection genes according to Bullard et al. (2010) BMC Bioinformatics.
    '''

    to_cluster = USECLUSTER
    statement = '''
    gunzip < %(infile)s
    | python %(scriptsdir)s/gtf2gtf.py --intersect-transcripts --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2gff.py --is-gtf --crop-unique  --log=%(outfile)s.log
    | python %(scriptsdir)s/gff2bed.py --is-gtf --log=%(outfile)s.log
    | gzip 
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@transform("*.bam", 
           suffix(".bam"), 
           add_inputs( buildUniqueExons),
           ".coverage.bed.gz")
def buildExonCoverage( infiles, outfile ):
    '''compute coverage of exons with reads.
    '''

    infile, exons = infiles

    to_cluster = USECLUSTER

    # set filter options
    # for example, only properly paired reads
    flag_filter = "-f 0x2"
    flag_filter = ""

    statement = '''
    samtools view -b %(flag_filter)s %(infile)s
    | bamToBed -i stdin 
    | coverageBed -a stdin -b %(exons)s 
    | gzip
    > %(outfile)s
    '''

    P.run()


#########################################################################
#########################################################################
#########################################################################
@merge(buildExonCoverage,
       "genelevel_raw_tagcounts.tsv.gz" )
def buildRawGeneLevelTagCounts( infiles, outfile ):
    '''aggregate exon level tag counts for each gene.

    coverageBed adds the following four columns:

    1) The number of features in A that overlapped (by at least one base pair) the B interval.
    2) The number of bases in B that had non-zero coverage from features in A.
    3) The length of the entry in B.
    4) The fraction of bases in B that had non-zero coverage from features in A.
    '''

    to_cluster = USECLUSTER
    
    src = " ".join( [ "<( zcat %s | groupBy -i stdin -g 4 -c 7 -o sum | sort -k1,1)" % x for x in infiles ] )
    
    tmpfile = P.getTempFilename( "." )
    
    statement = '''paste %(src)s 
                > %(tmpfile)s'''
    
    P.run()

    outf = IOTools.openFile( outfile, "w")
    tracks = [P.snip(x, ".coverage.bed.gz" ) for x in infiles ]
    outf.write( "gene_id\t%s\n" % "\t".join( tracks ) )
    
    for line in open( tmpfile, "r" ):
        data = line[:-1].split("\t")
        genes = list(set([ data[x] for x in range(0,len(data), 2 ) ]))
        values = [ data[x] for x in range(1,len(data), 2 ) ]
        assert len(genes) == 1, "paste command failed, wrong number of genes per line"
        outf.write( "%s\t%s\n" % (genes[0], "\t".join(map(str, values) ) ) )
    
    outf.close()

#########################################################################
#########################################################################
#########################################################################
@transform("*.bam", suffix(".bam"), ".gtf.gz")
def buildGeneModels(infile, outfile):
    '''build transcript models - run cufflinks on each region seperately
    '''

    to_cluster = USECLUSTER    

    track = os.path.basename( outfile[:-len(".gtf")] )

    tmpfilename = P.getTempFilename( "." )

    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )
    
    infile = os.path.abspath( infile )
    outfile = os.path.abspath( outfile )

    statement = '''mkdir %(tmpfilename)s; 
    cd %(tmpfilename)s; 
    cufflinks --label %(track)s           
              --reference %(cufflinks_genome_dir)s/%(genome)s
              %(cufflinks_options)s
               %(infile)s 
    >& %(outfile)s.log;
    gzip < transcripts.gtf > %(outfile)s;
    mv genes.expr %(outfile)s.genes.expr;
    mv transcripts.expr %(outfile)s.transcripts.expr
    '''

    P.run()

    shutil.rmtree( tmpfilename )

#########################################################################
#########################################################################
#########################################################################
@transform( buildGeneModels, 
            suffix(".gtf.gz"),
            ".annotations.gz" )
def annotateTranscripts( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set.
    '''
    to_cluster = USECLUSTER

    annotation_file = os.path.join( PARAMS["annotations_dir"],
                                    PARAMS_ANNOTATIONS["interface_annotation"] )

    statement = """
    zcat < %(infile)s 
    | python %(scriptsdir)s/gtf2table.py 
                --reporter=transcripts
		--counter=position 
		--counter=classifier
		--section=exons 
		--counter=length 
		--log=%(outfile)s.log 
		--filename-gff=%(annotation_file)s 
		--genome-file=%(genome_dir)s/%(genome)s
    | gzip
    > %(outfile)s"""
    
    P.run()

############################################################
@transform( annotateTranscripts, suffix(".annotations"), "_annotations.load" )
def loadAnnotations( infile, outfile ):
    '''load interval annotations: genome architecture
    '''
    P.load( infile, outfile, "--index=gene_id" )

#########################################################################
#########################################################################
#########################################################################
@transform("*.bam", 
           suffix(".bam"), 
           add_inputs( buildTranscripts),
           ".ref.gtf.gz")
def estimateExpressionLevels(infiles, outfile):
    '''estimate expression levels against a set
    of reference gene models.
    '''

    to_cluster = USECLUSTER    

    track = os.path.basename( outfile[:-len(".gtf")] )

    tmpfilename = P.getTempFilename( "." )

    if os.path.exists( tmpfilename ):
        os.unlink( tmpfilename )
    
    bamfile, gtffile = infiles

    gtffile = os.path.abspath( gtffile )
    bamfile = os.path.abspath( bamfile )
    outfile = os.path.abspath( outfile )

    # increase max-bundle-length to 4.5Mb due to Galnt-2 in mm9 with a 4.3Mb intron.
    statement = '''mkdir %(tmpfilename)s; 
    cd %(tmpfilename)s; 
    cufflinks --label %(track)s      
              --GTF=<(gunzip < %(gtffile)s)
              --reference %(cufflinks_genome_dir)s/%(genome)s.fa
              %(cufflinks_options)s
              %(bamfile)s 
    >& %(outfile)s.log;
    gzip < transcripts.gtf > %(outfile)s;
    mv -f genes.expr %(outfile)s.genes.expr;
    mv -f transcripts.expr %(outfile)s.transcripts.expr
    '''

    P.run()

    shutil.rmtree( tmpfilename )

#########################################################################
#########################################################################
#########################################################################
@transform( (estimateExpressionLevels, buildGeneModels), 
            suffix(".gtf.gz"), 
            "_gene_expression.load")
def loadExpressionLevels( infile, outfile ):
    '''load expression level measurements.'''

    track = P.snip( outfile, "_gene_expression.load" )
    P.load( infile + ".genes.expr",
            outfile = track + "_gene_expression.load",
            options = "--index=gene_id" )

    tablename = track + "_transcript_expression"
    infile2 = infile + ".transcripts.expr"

    statement = '''cat %(infile2)s
    | perl -p -e "s/trans_id/transcript_id/"
    | csv2db.py %(csv2db_options)s 
              --index=transcript_id 
              --table=%(tablename)s 
    > %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
def runCuffCompare( infiles, outfile, reffile ):
    '''run cuffcompare.

    Will create a .tmap and .refmap input file for each input file.
    '''

    to_cluster = USECLUSTER

    tmpdir = P.getTempDir( "." )
    
    cmd_extract = "; ".join( [ "gunzip < %s > %s/%s" % (x,tmpdir,x) for x in infiles ] )

    inf = " ".join( ["%s/%s" % (tmpdir,x) for x in infiles] )

    statement = '''
    %(cmd_extract)s;
    cuffcompare -o %(outfile)s
                -s %(cufflinks_genome_dir)s/%(genome)s.fa
                -r <( gunzip < %(reffile)s)
                %(inf)s
    > %(outfile)s.log;
    gzip -f %(outfile)s.{combined.gtf,tracking,loci};
    '''

    P.run()

    shutil.rmtree( tmpdir )

#########################################################################
#########################################################################
#########################################################################
@follows( buildGeneModels )
@files( [( ([ "%s.gtf.gz" % y.asFile() for y in EXPERIMENTS.getTracks(x)], buildTranscripts), 
           "%s.cuffcompare" % x.asFile()) 
         for x in EXPERIMENTS ] )
def compareTranscriptsPerExperiment( infiles, outfile ):
    '''compare transcript models between replicates within each experiment.'''
    infiles, reffile = infiles
    runCuffCompare( infiles, outfile, reffile )

#########################################################################
#########################################################################
#########################################################################
@merge( buildGeneModels, "%s.cuffcompare" % ALL.keys()[0].asFile() )
def compareTranscriptsBetweenExperiments( infiles, outfile ):
    '''compare transcript models between replicates in all experiments.'''
    # needs to be paramterized, unfortunately @merge has no add_inputs
    reffile = "transcripts.gtf.gz"
    runCuffCompare( infiles, outfile, reffile )

#########################################################################
#########################################################################
#########################################################################
@transform( (compareTranscriptsBetweenExperiments, 
             compareTranscriptsPerExperiment),
            suffix(".cuffcompare"), 
            "_cuffcompare.load" )
def loadTranscriptComparison( infile, outfile ):
    '''load data from transcript comparison.'''

    tmpfile = P.getTempFilename()
    
    outf = open( tmpfile, "w")
    outf.write( "track\tcontig\t%s\n" % "\t".join( Tophat.CuffCompareResult.getHeaders() ) )
    result = Tophat.parseTranscriptComparison( IOTools.openFile( infile ))
    tracks = []
    for track, vv in result.iteritems():
        track = P.snip( os.path.basename(track), ".gtf.gz" )
        tracks.append( track )
        for contig, v in vv.iteritems():
            if v.is_empty: continue
            outf.write( "%s\t%s\t%s\n" % (P.quote( track ), contig, str(v) ) )
    outf.close()

    tablename = P.toTable( outfile ) + "_benchmark"

    statement = '''cat %(tmpfile)s
    | csv2db.py %(csv2db_options)s 
              --index=track
              --index=contig
              --table=%(tablename)s 
    > %(outfile)s
    '''

    P.run()
    
    os.unlink( tmpfile )

    tablename = P.toTable( outfile ) + "_tracking"
    
    headers = ",".join( ( "transfrag_id",
                          "locus_id",
                          "transcript_id",
                          "code", ) + tuple( tracks ) )

    statement = '''zcat %(infile)s.tracking.gz 
    | csv2db.py %(csv2db_options)s
              --header=%(headers)s
              --index=locus_id
              --index=transfrag_id
              --table=%(tablename)s 
    >> %(outfile)s
    '''

    P.run()

#########################################################################
#########################################################################
#########################################################################
@files( [ ( ([ "%s.bam" % xx.asFile() for xx in EXPERIMENTS[x] ], 
             [ "%s.bam" % yy.asFile() for yy in EXPERIMENTS[y] ]),
            "%s_vs_%s.expression" % (x.asFile(),y.asFile()) )
              for x,y in itertools.combinations( EXPERIMENTS, 2) ] )
def estimateDifferentialExpressionPairwise( infiles, outfile ):
    '''estimate differential expression using cuffdiff.

    Replicates are grouped.
    '''
    
    to_cluster = USECLUSTER
    reffile = "transcripts.gtf.gz"        

    outdir = outfile + ".dir" 
    try: os.mkdir( outdir )
    except OSError: pass

    reps = "%s    %s" % (",".join( infiles[0]),
                         ",".join( infiles[1]) )
    
    statement = '''
    cuffdiff -o %(outdir)s
             --verbose
             -r %(cufflinks_genome_dir)s/%(genome)s.fa
             -p 20
             <(gunzip < %(reffile)s)
             %(reps)s
    >& %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
@merge("*.bam", 
       "expression.compare")
def estimateDifferentialExpression( infiles, outfile ):
    '''estimate differential expression using cuffdiff.

    Replicates are grouped.
    '''
    
    to_cluster = USECLUSTER
    reffile = "transcripts.gtf.gz"        

    outdir = outfile + ".dir" 
    try: os.mkdir( outdir )
    except OSError: pass

    # replicates are separated by ","
    reps = []
    for group, replicates in EXPERIMENTS.iteritems():
        reps.append( ",".join( [ "%s.bam" % r.asFile() for r in replicates] ) )
    reps = "   ".join( reps )

    statement = '''
    cuffdiff -o %(outdir)s
             --verbose
             -r %(cufflinks_genome_dir)s/%(genome)s.fa
             -p 20
             <(gunzip < %(reffile)s)
             %(reps)s
    >& %(outfile)s
    '''
    P.run()

#########################################################################
#########################################################################
#########################################################################
def getLibrarySizes( infiles ):
    
    vals = []
    for infile in infiles:
        assert infile.endswith( ".readstats")
        val, cont = [ x[:-1].split("\t") for x in open(infile).readlines() if re.search( "\tmapped", x ) ][0]
        vals.append(int(val))
        
    return vals

#########################################################################
#########################################################################
#########################################################################
@merge( loadExpressionLevels,
        "genelevel_fpkm_tagcounts.tsv.gz")
def buildFPKMGeneLevelTagCounts( infiles, outfile ):
    '''build tag counts using normalized counts from tophat.

    These are gene-length normalized count levels.

    They are obtained by multiplying the FPKM value
    by the median library size.
    '''
    infiles = [ x for x in infiles if x.endswith( ".ref_gene_expression.load" ) ]

    tracks = [ P.snip( x,".ref_gene_expression.load" ) for x in infiles ]

    # get normalization values
    library_sizes = getLibrarySizes( [ "%s.readstats" % x for x in tracks ] )
    if len(library_sizes) == 0: raise ValueError("could not get library sizes" )

    median_library_size = numpy.median( library_sizes )

    # dbhandle = sqlite3.connect( os.path.join( PARAMS["annotations_dir"],
    #                                           PARAMS_ANNOTATIONS["interface_database"] ) )
    # cc = dbhandle.cursor()    
    # median_gene_length = numpy.median( [ x for x in cc.execute( "SELECT sum FROM gene_stats") ] )

    scale = median_library_size / 1000000.0

    E.info( "normalization: median library size=%i, factor=1.0 / %f" % \
                (median_library_size, scale) )

    # normalize
    results = []
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()    

    for track in tracks:
        table = "%s_ref_gene_expression" % P.quote(track)
        statement = "SELECT gene_id, FPKM / %(scale)f FROM %(table)s" % locals()
        results.append( dict( cc.execute( statement ).fetchall() ) )
    
    outf = IOTools.openFile( outfile, "w" )
    gene_ids = set()
    for x in results: gene_ids.update( x.keys() )

    outf.write( "gene_id\t%s\n" % "\t".join( tracks ) )
    for gene_id in gene_ids:
        outf.write( "%s\t%s\n" % ( gene_id, "\t".join( [str(int(x[gene_id])) for x in results ] ) ) )
    outf.close()
            

#########################################################################
#########################################################################
#########################################################################
@transform( (buildRawGeneLevelTagCounts,
             buildFPKMGeneLevelTagCounts),
            suffix(".tsv.gz"),
            ".diff")
def estimateDifferentialExpressionDESeq( infile, outfile ):
    '''estimate differential expression using DESeq
    '''
    
    to_cluster = USECLUSTER

    R('''library('DESeq')''')

    # load data 
    R( '''counts_table <- read.delim( '%s', header = TRUE, row.names = 1, stringsAsFactors = TRUE )''' % infile )

    # get conditions to test
    # note that tracks in R use a '.' as separator
    tracks = R('''colnames(counts_table)''')

    map_track2column = dict( [ (y,x) for x,y in enumerate( tracks ) ] )
    
    conds = [None] * len(tracks)
    for group, replicates in EXPERIMENTS.iteritems():
        for r in replicates:
            conds[map_track2column[r.asR()]] = group.asR()

    ro.globalenv['conds'] = ro.StrVector(conds)

    # this analysis follows the 'Analysing RNA-Seq data with the "DESeq" package'
    # tutorial 
    R('''cds <-newCountDataSet( counts_table, conds) ''')
    R('''cds <- estimateSizeFactors( cds )''')
    R('''cds <- estimateVarianceFunctions( cds )''')

    E.info("creating diagnostic plots" ) 
    size_factors = R('''sizeFactors( cds )''')
    R.png( "%s_scvplot.png" % outfile )
    R('''scvPlot( cds, ylim = c(0,3))''')
    R['dev.off']()

    R('''vsd <- getVarianceStabilizedData( cds )''' )
    R('''dists <- dist( t( vsd ) )''')
    R.png( "%s_heatmap.png" % outfile )
    R('''heatmap( as.matrix( dists ), symm=TRUE )''' )
    R['dev.off']()

    for track in conds:
        R.png( "%s_fit_%s.png" % (outfile, track) )
        R('''diagForT <- varianceFitDiagnostics( cds, "%s" )''' % track )
        R('''smoothScatter( log10(diagForT$baseMean), log10(diagForT$baseVar) )''')
        R('''lines( log10(fittedBaseVar) ~ log10(baseMean), diagForT[ order(diagForT$baseMean), ], col="red" )''')
        R['dev.off']()
        R.png( "%s_residuals_%s.png" % (outfile, track) )
        R('''residualsEcdfPlot( cds, "%s" )''' % track )
        R['dev.off']()

    E.info("calling differential expression")

    for x,y in itertools.combinations( conds, 2 ):
        R('''res <- nbinomTest( cds, '%s', '%s' )''' % (x,y))
        R.png( "%s_diff_%s_vs_%s.png" % (outfile, x, y) )
        R('''plot( res$baseMean, res$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res$padj < .1, "red", "black" ) )''')
        R['dev.off']()


#########################################################################
#########################################################################
#########################################################################
@follows(mapReadsFromSRAWithTophat )
def mapReads(): pass

@follows( mapReads,
          buildGeneModels,
          loadTranscriptComparison,
          estimateExpressionLevels,
          estimateDifferentialExpression )
def full(): pass

def export(): pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

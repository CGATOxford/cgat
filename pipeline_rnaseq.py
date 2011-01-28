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

The RNA-Seq pipeline imports mapped reads from one or more
RNA-Seq experiments and performs the following tasks:

   * predict transcripts
   * estimate expression levels of transcripts
   * predict differentially expressed transcripts/genes

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

EXPERIMENTS = PipelineTracks.Aggregate( TRACKS, "condition", "tissue" )
CONDITIONS = PipelineTracks.Aggregate( TRACKS, "condition" )
TISSUES = PipelineTracks.Aggregate( TRACKS, "tissue" )

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

def getInsertSizes( reads ):
    return 200, 50

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

@transform("*.bam", suffix(".bam"), ".gtf.gz")
def buildGeneModels(infile, outfile):
    '''build transcript models - run cufflinks on each region seperately'''

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

@transform("*.bam", 
           suffix(".bam"), 
           add_inputs( buildTranscripts),
           ".ref.gtf.gz")
def estimateExpressionLevels(infiles, outfile):
    '''estimate expression levels.'''

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
    mv genes.expr %(outfile)s.genes.expr;
    mv transcripts.expr %(outfile)s.transcripts.expr
    '''

    P.run()

    shutil.rmtree( tmpfilename )

@transform( (estimateExpressionLevels,
             buildGeneModels), 
            suffix(".gtf.gz"), 
            ".gtf.load")
def loadGeneModels( infile, outfile ):
    
    track = P.snip( outfile, ".gtf.load" )
    P.load( infile + ".genes.expr",
            outfile = track + ".genes.load",
            options = "--index=gene_id" )

    tablename = track + "_transcripts_expression"
    infile2 = infile + ".transcripts.expr"

    statement = '''cat %(infile2)s
    | perl -p -e "s/trans_id/transcript_id/"
    | csv2db.py %(csv2db_options)s 
              --index=gene_id 
              --table=%(tablename)s 
    > %(outfile)s
    '''

    P.run()

@merge( buildGeneModels, "transcripts.compare" )
# suffix(".gtf.gz"), 
# add_inputs( buildTranscripts ),
# ".cuffcompare" )
def compareTranscripts( infiles, outfile ):
    '''run cuffcompare.

    Will create a .tmap and .refmap input file for each input file.
    '''

    to_cluster = USECLUSTER
    reffile = "transcripts.gtf.gz"

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
    gzip %(outfile)s.combined.gtf
    '''

    P.run()

    shutil.rmtree( tmpdir )

@transform( compareTranscripts, suffix(".compare"), "_compare.load" )
def loadTranscriptComparison( infile, outfile ):
    '''load data from transcript comparison.'''


    tmpfile = P.getTempFilename()
    
    outf = open( tmpfile, "w")
    outf.write( "track\tcontig\t%s\n" % "\t".join( Tophat.CuffCompareResult.getHeaders() ) )
    result = Tophat.parseTranscriptComparison( IOTools.openFile( infile ))
    for track, vv in result.iteritems():
        track = P.snip( os.path.basename(track), ".gtf.gz" )
        for contig, v in vv.iteritems():
            if v.is_empty: continue
            outf.write( "%s\t%s\t%s\n" % (P.quote( track ), contig, str(v) ) )
    outf.close()

    tablename = P.toTable( outfile )

    statement = '''cat %(tmpfile)s
    | csv2db.py %(csv2db_options)s 
              --index=track
              --index=contig
              --table=%(tablename)s 
    > %(outfile)s
    '''

    P.run()
    
    os.unlink( tmpfile )

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
    replicates = getReplicateGroups( [ P.snip( x ) for x in infiles ] )
    reps = []
    for track, r in replicates.iteritems():
        reps.append( ",".join( [ "%s.bam" % buildTrack(*rr) for rr in r] ) )
    reps = " ".join( reps )

    statement = '''
    cuffdiff -o %(outdir)s
             --verbose
             -r %(cufflinks_genome_dir)s/%(genome)s.fa
             -p 20
             <(gunzip < %(reffile)s | grep "chr20" | grep "protein_coding") 
             %(reps)s
    >& %(outfile)s
    '''
    P.run()

@merge("*.genes.load",
       "gene_tagcounts.tsv.gz")
def buildExpressionTable( infiles, outfile ):
    '''estimate differential expression using DESeq
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    cc = dbhandle.cursor()    

    results = []
    tracks = [ P.snip( x,".genes.load" ) for x in infiles ]
    for track in tracks:
        table = "%s_ref_genes" % P.quote(track)
        statement = "SELECT gene_id, FPKM FROM %(table)s" % locals()
        results.append( dict( cc.execute( statement ).fetchall() ) )
    
    outf = IOTools.openFile( outfile, "w" )
    gene_ids = set()
    for x in results: gene_ids.update( x.keys() )

    outf.write( "gene_id\t%s\n" % "\t".join( tracks ) )
    for gene_id in gene_ids:
        outf.write( "%s\t%s\n" % ( gene_id, "\t".join( [str(int(x[gene_id])) for x in results ] ) ) )
    outf.close()
            

@transform( buildExpressionTable,
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
    tracks = [ re.sub( "\.", "-", x) for x in R('''colnames(counts_table)''') ]

    map_track2column = dict( [ (y,x) for x,y in enumerate( tracks ) ] )
    
    conds = [None] * len(tracks)
    for group, replicates in getReplicateGroups( tracks ).iteritems():
        for r in replicates:
            conds[map_track2column[buildTrack(*r)]] = ".".join(group)
    print conds
    
    # R(''' cds <- newCountDataSet( counts_table, conds )''')

@follows(mapReadsFromSRAWithTophat )
def mapReads(): pass

@follows( mapReads,
          buildGeneModels,
          estimateExpressionLevels,
          compareTranscripts,
          estimateDifferentialExpression )
def full(): pass

def export(): pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

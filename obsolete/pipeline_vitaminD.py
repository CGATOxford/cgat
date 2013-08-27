################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_vitaminD.py 2870 2010-03-03 10:20:29Z andreas $
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
pipeline_vitaminD.py - vitaminD project pipeline
================================================

:Author: Andreas Heger
:Release: $Id: pipeline_vitaminD.py 2870 2010-03-03 10:20:29Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

vitamin D project pipeline

Input

   experimental tracks should be called run<cellline><condition><replicate>
   control tracks should be called control<cellline>

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

TODO: currently the bed files and the intervals are inconsistent 
    (due to filtering, there are more intervals in the bed files than
     in the table. The ids do correspond).

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
import CGAT.Pipeline as P
from ruffus import *
import csv
import sqlite3
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.MAST as MAST
import CGAT.GTF as GTF
import CGAT.GFF as GFF
import CGAT.Bed as Bed
import CGAT.Stats as Stats
import cStringIO
import pysam
import numpy
import gzip
import CGAT.Expression as Expression
import CGAT.Masker as Masker
import CGAT.Glam2Scan as Glam2Scan
import fileinput
import CGAT.Motifs as Motifs

try: import nubiscan
except ImportError: pass

import gff2annotator
import CGAT.Bioprospector as Bioprospector

import pipeline_vitaminD_annotator as PAnnotator
import pipeline_vitaminD_expression as PExpression
import PipelineChipseq as PipelineChipseq
import PipelineMotifs as PipelineMotifs
import PipelineGeneset as PGeneset
import time

if os.path.exists("conf.py"):
    execfile("conf.py")

TARGET_ANNOTATION= 'ensembl_regions.gff'
TARGET_GENESET= 'ensembl.gtf'
TARGET_PROMOTORS = 'promotors.gtf'
TARGET_TSS = 'tss.gtf'
TARGET_REPEATS = 'repeats.gff'
TARGET_TRANSCRIPTS = 'transcripts.gtf.gz'
TARGET_PROBESET = 'probeset.gtf'
TARGET_TRANSCRIPTS_TSS = 'transcripts_tss.gtf'
TARGET_TRANSCRIPTS_PROMOTORS = 'transcripts_promotors.gtf'
TARGET_ANNOTATOR_GENETERRITORIES='annotator_geneterritories.gff'
TARGET_MAPPABILITY='mappability.bed'
BAM_SUFFIX = ".norm.bam"

PARAMS = P.getParameters()

###################################################################
###################################################################
###################################################################
## General preparation tasks
###################################################################

############################################################
############################################################
############################################################
@files( "genome.fasta", "genome.fa" )
def indexGenome( infile, outfile ):
    '''index the genome for samtools.

    Samtools does not like long lines, so create a new file
    with split lines (what a waste).
    '''

    statement = '''fold %(infile)s > %(outfile)s'''
    P.run()
    
    pysam.faidx( outfile )

############################################################
############################################################
############################################################
## get UCSC tables
############################################################
def getUCSCTracks(infile = PARAMS["filename_ucsc_encode"]):
    '''return a list of UCSC tracks from infile.'''
    tables = []
    with open(infile) as f:
        for line in f:
            if line.startswith("#"): continue
            tablename = line[:-1].strip()
            if tablename == "": continue
            tables.append( tablename )
    return tables

############################################################
############################################################
############################################################
## 
############################################################
@files( "genome.fasta", "genome_gc.bed" )
def buildGenomeGCSegmentation( infile, outfile ):
    '''segment the genome into windows according to G+C content.'''

    to_cluster = True
    
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

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
## import UCSC encode tracks
############################################################
@posttask( touch_file("ucsc_encode.import") )
@files( PARAMS["filename_ucsc_encode"], "ucsc_encode.import") 
def importUCSCEncodeTracks( infile, outfile ):
    
    statement = '''
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -B -e "SELECT * FROM %(tablename)s" %(ucsc_database)s |\
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --table=%(tablename)s \
    >> %(outfile)s

    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    cc = dbhandle.cursor()
    tables = set( [ x[0] for x in cc.execute( "SELECT name FROM sqlite_master WHERE type='table'") ] )
    cc.close()
    
    for tablename in getUCSCTracks( infile ):
        if tablename in tables:
            E.info( "skipping %(tablename)s - already exists" % locals())
            continue            
            
        E.info( "importing %(tablename)s" % locals() )
        P.run()

############################################################
############################################################
############################################################
## import UCSC encode tracks
############################################################
@files( PARAMS["filename_mappability"], PARAMS["annotator_mappability"] ) 
def exportUCSCMappabilityTrackToBed( infile, outfile ):
    '''convert wiggle track with mappability information to a bed-formatted file
    with only the mappable regions in the genome.'''
    
    infile = gzip.open( infile, "r" )

    outf = open( outfile, "w" )
    
    for line in infile:
        if line.startswith("fixedStep" ):
            contig, start, step = re.match( "fixedStep chrom=(\S+) start=(\d+) step=(\d+)", line ).groups()
            start, step = int(start)-1, int(step)
            end = start + step
            last_val = None
        else:
            val = int(line)
            if last_val != val:
                if last_val == 1:
                    outf.write( "\t".join( (contig, str(start), str(end ) ) ) + "\n" )
                start = end
                end = start + step
            else:
                end += step
            last_val = val
    outf.close()
    
############################################################
############################################################
############################################################
## export UCSC encode tracks as bed
############################################################
@transform( importUCSCEncodeTracks, suffix(".import"), ".bed")
def exportUCSCEncodeTracks( infile, outfile ):

    dbhandle = sqlite3.connect( PARAMS["database"] )

    outs = open(outfile, "w")
    for tablename in getUCSCTracks():
        outs.write( "track name=%s\n" % tablename )
        
        cc = dbhandle.cursor()
        statement = "SELECT chrom, chrostart, chroend FROM %s ORDER by chrom, chrostart" % (tablename)
        cc.execute( statement )
        for contig, start, end in cc:
            outs.write("%s\t%i\t%i\n" % (contig, start, end) )
    outs.close()

############################################################
############################################################
############################################################
@files( PARAMS["filename_ensembl_geneset"], "transcript_info.import" )
def importTranscriptInformation( infile, outfile ):
    '''import the transrcipt set.'''
    PGeneset.importTranscriptInformation( infile, outfile )

############################################################
############################################################
############################################################
@files( PARAMS["filename_ensembl_geneset"], "gene_info.import" )
def importGeneInformation( infile, outfile ):
    '''import the transrcipt set.'''
    PGeneset.importGeneInformation( infile, outfile )

###################################################################
###################################################################
###################################################################
##
###################################################################
@files( ( ( PARAMS["filename_regions_of_interest"], "roi.import" ),
          ( PARAMS["filename_selection"], "selection.import" ) ) )
def importRegionsOfInterest( infile, outfile ):
    '''import regions of interest.'''
    
    table = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
               --dialect=excel \
               --index=class \
               --index=roi_id \
               --map=roi_id:str \
               --map=start:int \
               --map=end:int \
               --table=%(table)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

###################################################################
###################################################################
###################################################################
##
###################################################################
@merge( PARAMS["filename_regions_of_interest"], "snps.import" )
def importSNPs( infile, outfile ):
    '''import snps from :term:`roi` file.'''

    reader = csv.DictReader( open(infile,"rU") )
    ids = {}
    for row in reader:
        
        roi_id, contig, start, end, snps =\
                row["roi_id"], row["contig"], row["start"], row["end"], row["snp"]
        snps = snps.split(",")
        for snp in snps:
            s = snp.strip()
            if s == "na" or s == "": continue
            try:
                ids[snp.strip()] = (roi_id, contig, int(start), int(end))
            except ValueError:
                # ignore empty roi
                continue
            
    tmpf = P.getTempFile()
    tmpf.write("roi_id\tsnp\tcontig\tpos\n" )

    counter = E.Counter()
    
    found = set()
    inf = gzip.open( PARAMS["filename_dbsnp"], "r" )
    for line in inf:
        data = line[:-1].split("\t")
        contig, snp, pos = data[1], data[4], data[2]
        if snp not in ids: continue

        roi_id, roi_contig, roi_start, roi_end = ids[snp]

        if contig != roi_contig or not roi_start <= int(pos) < roi_end:
            counter.out_of_range += 1
            continue

        counter.accepted += 1

        tmpf.write( "%s\t%s\t%s\t%s\n" % (roi_id, snp, contig, pos ) )
        found.add( snp )

    inf.close()

    table = outfile[:-len(".import")]
    
    tmpfilename = tmpf.name
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
               --index=roi_id \
               --map=roi_id:str \
               --table=%(table)s \
    < %(tmpfilename)s > %(outfile)s
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )
    
    os.unlink( tmpfilename )

    outfile = open( outfile, "a" )
    outfile.write( "input=%i, found=%i, %s\n" % (len(ids), len(found), str(counter)))
    outfile.write( "snps not found: %s\n" % ",".join( set(ids.keys()).difference( found) ) )
    outfile.close()
    
    E.info( "snp import: input=%i, found=%i, %s" % (len(ids), len(found), str(counter)))

    if len(ids) != len(found):
        E.warn( "snps not found: %s" % ",".join( set(ids.keys()).difference( found) ) )

###################################################################
###################################################################
###################################################################
##
###################################################################
@merge( PARAMS["filename_gwas"], "gwas.import" )
def importGWAS( infile, outfile ):
    '''import GWAS intervals. 

    GWAS intervals are defined by their marking SNP. Each 
    snp is looked up in dbsnp and extended by ``gwas_interval_halfwidth``.

    The same SNP can be part of multiple gwas studies.

    The script checks not for overlap between annotations.
    '''

    reader = csv.DictReader( open(infile,"rU") )

    ids = collections.defaultdict( list )
    last_gwas = None

    for row in reader:

        snps, gwas, genes = row["snp"], row["class"], row["genes"]
        if gwas == "": gwas = last_gwas

        snps = snps.split(",")
        for snp in snps:
            s = snp.strip()
            if s == "na" or s == "": continue
            genes=re.sub(" ", "", genes.lower())
            ids[s.strip()].append( (gwas,genes) )

        last_gwas = gwas
            
    tmpf = P.getTempFile()
    tmpf.write( "%s\n" % "\t".join( ("roi_id","snp","class","contig","pos","start","end","genes")) )

    counter = E.Counter()
    
    found = set()
    inf = gzip.open( PARAMS["filename_dbsnp"], "r" )
    nsnps = 0
    halfwidth = PARAMS["gwas_interval_halfwidth"]
    for line in inf:
        data = line[:-1].split("\t")
        contig, snp, pos = data[1], data[4], data[2]
        if snp not in ids: continue
        if contig == None or contig == "": 
            counter.unmapped += 1
            continue
            
        counter.accepted += 1
        pos = int(pos)
        for gwas, genes in ids[snp]:
            nsnps += 1
            tmpf.write( "%s\n" % "\t".join( map(str, 
                                                (nsnps, 
                                                 snp, 
                                                 gwas,
                                                 contig, 
                                                 pos, 
                                                 pos - halfwidth,
                                                 pos + halfwidth,
                                                 genes
                                                 ) )  ))

            
        found.add( snp )

    inf.close()
    tmpf.close()

    table = outfile[:-len(".import")]
    
    tmpfilename = tmpf.name
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
               --index=roi_id \
               --index=class \
               --map=roi_id:str \
               --table=%(table)s \
    < %(tmpfilename)s > %(outfile)s
    '''
    P.run( **locals() )
    
    os.unlink( tmpfilename )

    outfile = open( outfile, "a" )
    outfile.write( "input=%i, found=%i, %s\n" % (len(ids), len(found), str(counter)))
    outfile.write( "snps not found: %s\n" % ",".join( set(ids.keys()).difference( found) ) )
    outfile.close()
    
    E.info( "gwas import: input=%i, found=%i, %s" % (len(ids), len(found), str(counter)))

    if len(ids) != len(found):
        E.warn( "snps not found: %s" % ",".join( set(ids.keys()).difference( found) ) )

###################################################################
###################################################################
###################################################################
##
###################################################################
@merge( PARAMS["filename_gwas"], "gwas_merged.import" )
def importMergedGWAS( infile, outfile ):
    '''import GWAS intervals. 

    GWAS intervals are defined by their marking SNP. Each 
    snp is looked up in dbsnp and extended by ``gwas_interval_halfwidth``.

    The same SNP can be part of multiple gwas studies.

    This method checks for overlap between annotations
    '''

    reader = csv.DictReader( open(infile,"rU") )

    ids = collections.defaultdict( list )
    last_gwas = None

    for row in reader:

        snps, gwas, genes = row["snp"], row["class"], row["genes"]
        if gwas == "": gwas = last_gwas

        snps = snps.split(",")
        for snp in snps:
            s = snp.strip()
            if s == "na" or s == "": continue
            genes=re.sub(" ", "", genes.lower())
            ids[s.strip()].append( (gwas,genes) )

        last_gwas = gwas
            
    counter = E.Counter()
    
    found = set()
    inf = gzip.open( PARAMS["filename_dbsnp"], "r" )
    halfwidth = PARAMS["gwas_interval_halfwidth"]
    # collecat all snp coordinates for each gwas interval
    intervals = collections.defaultdict( list )
    for line in inf:
        data = line[:-1].split("\t")
        contig, snp, pos = data[1], data[4], data[2]
        if snp not in ids: continue
        if contig == None or contig == "": 
            counter.unmapped += 1
            continue
            
        counter.accepted += 1
        pos = int(pos)
        for gwas, genes in ids[snp]:
            intervals[gwas].append( (contig, pos, snp) )
            
        found.add( snp )

    inf.close()
    
    tmpf = P.getTempFile()
    tmpf.write( "%s\n" % "\t".join( ("roi_id","snp","class","contig","pos","start","end","genes")) )


    roi_id = 0
    for gwas, snps in intervals.iteritems():
        
        def iter_overlaps( snps, halfwidth ):

            snps.sort()
            last_contig, last_pos, last_snp = snps[0]
            width = halfwidth * 2
            r = [(last_pos,last_snp)]
            for contig, pos, snp in snps[1:]:
                if last_contig != contig or \
                        pos - last_pos > width:
                    yield last_contig, r
                    r = []
                    last_contig = contig
                    last_pos = pos
                r.append( (pos, snp) )
            yield last_contig, r
                
            
        for contig, overlaps in iter_overlaps( snps, halfwidth ):
            pos = [x[0] for x in overlaps]
            pos.sort()
            start = pos[0] - halfwidth
            end = pos[-1] + halfwidth
            genes = []
            xsnps = [x[1] for x in overlaps]
            for snp in xsnps:
                for gwas, ggenes in ids[snp]:
                    if ggenes: 
                        genes.append( ggenes )
            genes = sorted(set(genes))

            roi_id += 1
            tmpf.write( "%s\n" % "\t".join( map(str, 
                                                (roi_id, 
                                                 ",".join(xsnps), 
                                                 gwas,
                                                 contig, 
                                                 ",".join(map(str,pos)), 
                                                 start,
                                                 end,
                                                 ",".join(genes),
                                                 ) )  ))



    tmpf.close()

    table = outfile[:-len(".import")]
    
    tmpfilename = tmpf.name
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
               --index=roi_id \
               --index=class \
               --map=roi_id:str \
               --table=%(table)s \
    < %(tmpfilename)s > %(outfile)s
    '''
    P.run( **locals() )
    
    os.unlink( tmpfilename )

    outfile = open( outfile, "a" )
    outfile.write( "input=%i, found=%i, %s\n" % (len(ids), len(found), str(counter)))
    outfile.write( "snps not found: %s\n" % ",".join( set(ids.keys()).difference( found) ) )
    outfile.close()
    
    E.info( "gwas import: input=%i, found=%i, %s" % (len(ids), len(found), str(counter)))

    if len(ids) != len(found):
        E.warn( "snps not found: %s" % ",".join( set(ids.keys()).difference( found) ) )

###################################################################
###################################################################
###################################################################
##
###################################################################
@transform( (importRegionsOfInterest, 
             importGWAS, 
             importMergedGWAS ),
            suffix( ".import"),
            "_genes.import")
def importRegionsOfInterestGenes( infile, outfile ):
    '''import association between regions of interest and genes implicated in disease.'''

    intablename = infile[:-len(".import")]
    tablename = outfile[:-len(".import")]

    E.info( "tablenames: %s, %s" % (intablename, tablename))

    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    statement = "SELECT lower(gene_name), gene_id FROM gene_info"
    cc = dbhandle.cursor()
    map_gene2id = dict( cc.execute( statement ).fetchall() )

    statement = "SELECT roi_id, lower(genes) FROM %(intablename)s" % locals()

    tmpf = P.getTempFile()
    tmpf.write("roi_id\tgene_id\tgene_name\n" )
    
    cc = dbhandle.cursor()
    counter = E.Counter()
    notfound = set()

    for roi_id, g in cc.execute( statement ):
        counter.input += 1

        if g == None: 
            counter.nogenes += 1
            continue

        genes = [x.lower() for x in g.split( "," )]

        for gene_name in genes:
            if gene_name in map_gene2id:
                counter.found += 1
                tmpf.write( "%s\n" %\
                            "\t".join( map(str,
                                           ( roi_id, map_gene2id[gene_name], gene_name) ) ) )
            else:
                counter.missed += 1
                notfound.add( gene_name )

    E.info( "importRegionsOfInterestGenes: %s" % (str(counter)))

    outf = open( outfile, "a" )
    outf.write( "%s\n" % (str(counter)))
    outf.write( "gene names not found: %s\n" % ",".join(notfound) )
    outf.close()

    tmpf.close()

    tmpfilename = tmpf.name
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
               --allow-empty
               --index=roi_id 
               --map=roi_id:str 
               --table=%(tablename)s 
    < %(tmpfilename)s 
    >> %(outfile)s
    '''

    P.run()
    
    os.unlink( tmpfilename )
        
###################################################################
###################################################################
###################################################################
##
###################################################################
@merge( PARAMS["filename_snps_of_interest"], "snps_of_interest.import" )
def importSNPsOfInterest( infile, outfile ):
    '''import regions of interest.'''
    
    table = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
               --dialect=excel \
               --index=snp \
               --map=snp:str \
               --map=pos:int \
               --table=%(table)s \
    < %(infile)s 
    > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
## export regions of interest
############################################################
@transform( (importRegionsOfInterest, importGWAS, importMergedGWAS), suffix(".import"), ".bed")
def exportRegionsOfInterest( infile, outfile ):

    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
        
    outs = open(outfile, "w")
    tablename = outfile[:-len(".bed")]

    tracks = [x[0] for x in cc.execute( "SELECT DISTINCT class FROM %(tablename)s" % locals() ).fetchall()]
    
    for track in tracks:
        
        outs.write( "track name=%s\n" % track )

        statement = "SELECT contig, max(0,start), end, roi_id FROM %(tablename)s WHERE class = '%(track)s' ORDER by contig, start" % locals()
        cc.execute( statement )
        for contig, start, end, roi_id in cc:
            outs.write("%s\t%i\t%i\t%s\n" % (contig, start, end, roi_id) )
    outs.close()

############################################################
############################################################
############################################################
## import reference intervals
############################################################
@files( [ ("%s.bed" % x, "%s_intervals.import" % x) for x in TRACKS_REFERENCE ] )
def importReferenceIntervals( infile, outfile ):
    '''import reference intervals'''

    ## add the bed intervals back to table
    track = infile[:-len(".bed")]

    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    tmpfilename = tmpfile.name

    headers = ("AvgVal","DisttoStart","GeneList","Length",
               "PeakCenter","PeakVal","Position",
               "interval_id","nCpGs","nGenes","nPeaks",
               "nProbes","nPromoters", 
               "contig","start","end" )

    avgval,contig,disttostart,end,genelist,length,peakcenter,peakval,position,start,interval_id,ncpgs,ngenes,npeaks,nprobes,npromoters = \
        0,"",0,0,"",0,0,0,0,0,0,0,0,0,0,0,

    tmpfile.write( "\t".join( headers ) + "\n" )

    for line in open( infile, "r" ):
        contig, start, end, interval_id = line[:-1].split("\t")[:4]
        start, end = int(start), int(end)
        length = end - start
        npeaks, peakcenter = 1, start + int(length / 2)
        avgval, peakval, nprobes = 1, 1, 1

        tmpfile.write( "\t".join( map( str, (avgval,disttostart,genelist,length,peakcenter,peakval,position,interval_id,ncpgs,ngenes,npeaks,nprobes,npromoters, contig,start,end) )) + "\n" )

    tmpfile.close()

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=interval_id \
              --table=%(track)s_intervals \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run()

    os.unlink( tmpfile.name )

###################################################################
###################################################################
###################################################################
## Version 1: import from given bed files
###################################################################

if PARAMS["method"] == "intervals":
    ###################################################################
    ###################################################################
    ###################################################################
    ##
    ###################################################################
    @files( PARAMS["filename_intervals"], [ "%s.intervals" % x for x in TRACKS_RAW] )
    def createIntervals( infile, outfiles ):
        '''create intervals files from the input csv file.
        '''
        tmpdir = tempfile.mkdtemp()
        reader = csv.reader( open(infile,"rU") )

        # convert to well formatted file for stats
        output_files_info = IOTools.FilePool( output_pattern = "%s.intervals" )
        output_files_bed = IOTools.FilePool( output_pattern = os.path.join( tmpdir, "%s.bed" ) )

        for row in reader:
            # process header
            row = [ re.sub( "\s", "", x) for x in row ]
            row = [ re.sub( "^#", "n", x) for x in row ]

            if row[0] == 'Sample':
                row = [ re.sub( "^#", "n", x) for x in row ]
                row = [ re.sub( "\s", "", x) for x in row ]
                row[1] = "interval_id"
                # delete coordinates
                del row[2:5]
                output_files_info.setHeader( "\t".join(row[1:])+ "\n" )
            else:
                # convert to 0-based open/closed coordinates
                # start
                row[3] = str(int(row[3]) - 1)
                # PeakCenter
                row[7] = str(int(row[7]) - 1)
                # add chr
                row[2] = "chr%s" % row[2]

                run_id, cell, condition, replicate = row[0].split("_")
                id = "run%s%s%s" % (cell, condition, replicate)

                # output coordinates
                output_files_bed.write( id, "\t".join( row[2:5] + row[1:2] ) + "\n" )

                # delete coordinates
                del row[2:5]
                output_files_info.write( id, "\t".join( row[1:])+ "\n" )

        output_files_bed.close()
        output_files_info.close()

        if PARAMS["filename_chain"]:
            chain = PARAMS["filename_chain"]
            statement = '''liftOver %(filename)s %(chain)s %(basename)s %(basename)s.unmapped'''

        else:
            statement = '''cp %(filename)s %(basename)s'''

        for filename in output_files_bed.keys():
            basename = os.path.basename( filename )
            P.run( **dict( locals().items() + PARAMS.items() ) )

        shutil.rmtree( tmpdir )

    ############################################################
    ############################################################
    ############################################################
    ##
    ############################################################
    @transform( "*.intervals", suffix(".intervals"), "_intervals.import" )
    def importIntervals( infile, outfile ):

        ## add the bed intervals back to table
        track = infile[:-len(".intervals")]

        tmpfile = tempfile.NamedTemporaryFile(delete=False)
        tmpfilename = tmpfile.name
        map_id2coords = {}
        for line in open( track + ".bed", "r" ):
            contig, start, end, id = line[:-1].split("\t")[:4]
            map_id2coords[id] = (contig, start, end )

        for line in open( infile ):
            if line.startswith( "interval_id"):
                tmpfile.write( "%s\t%s\t%s\t%s\n" % ((line[:-1],"contig","start","end") ) )
            else:
                id = line[:-1].split( "\t")[0]
                try:
                    tmpfile.write( "%s\t%s\t%s\t%s\n" % ((line[:-1],) +  map_id2coords[id]))
                except KeyError:
                    P.warn( "%s -> interval %s - omitted due to mapping problems" % (infile, id) )

        tmpfile.close()

        statement = '''
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
                  --index=interval_id \
                  --table=%(track)s_intervals \
        < %(tmpfilename)s > %(outfile)s
        '''

        P.run()

        os.unlink( tmpfile.name )

    # @transform( importIntervals, suffix("_intervals.import"), ".bed" )
    @transform( "*_intervals.import", suffix("_intervals.import"), ".bed" )
    def exportIntervalsAsBed( infile, outfile ):
        _exportIntervalsAsBed( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @files_re( PARAMS["filename_reads"],
               "(\S+).export.txt.gz", 
               r"\1.bam" )
    def buildBAM( infile, outfile ):
        '''build bam formatted files. The files will be sorted
        and indexed.

        The resultant BAM files contain

        1. all the uniquely mapped reads
        2. all the unmapped reads that have passed the quality filter.
        '''

        tmpfilename = P.getTempFilename()

        statement = '''
        /net/cpp-group/src/samtools/b64/samtools-0.1.6/misc/export2sam.pl \
        <(gunzip < %(infile)s) |\
        sed "s/.fa//" |\
        samtools import <( cut -f 1,4 %(genome)s.idx ) - %(tmpfilename)s
        '''

        P.run( **dict( locals().items() + PARAMS.items() ) )

        prefix = outfile[:-len(".bam")]
        pysam.sort( tmpfilename, prefix ) 
        pysam.index( outfile )

        os.unlink( tmpfilename )

    ############################################################
    ############################################################
    ############################################################
    @files_re( buildBAM, 
               "(.*).bam",
               r"\1.readstats" )
    def buildBAMStats( infile, outfile ):
        PipelineChipseq.buildBAMStats( infile, outfile )

    @transform( buildBAMStats,
                regex(r"(run.*).readstats"),
                inputs( (r"\1.bam", r"\1.readstats") ),
                r"\1.norm.bam" )
    def normalizeBAMPerReplicate( infile, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.
        '''
        track = infile[:-len(".bam")]
        PipelineChipseq.buildNormalizedBAM( ((infile,track + ".readstats"),), outfile )

    @follows( buildBAMStats )
    @files( [ ( [ ("%s%s.bam" % (x,y), "%s%s.readstats" % (x,y))  for y in REPLICATES],
                "%s.norm.bam" % x) for x in TRACKS_EXPERIMENTS ] )
    def normalizeBAMPerExperiment( infiles, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.


        Duplicated reads are removed at the same time.

        Merge reads from several replicates.
        '''
        PipelineChipseq.buildNormalizedBAM( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( exportIntervalsAsBed )
    @files( [ ( ["run%s%s.bed" % (y,x) for y in CELLLINES],
                "run%s.bed" % x, 
                x) for x in CONDITIONS ] )
    def combineConditions( infile, outfile, track ):
        '''combine conditions between cell lines. 

        The conditions are merged via intersection.
        '''

        conditions = [ "run%s%s.bed" % (x,track) for x in CELLLINES ]
        PipelineChipseq.intersectBedFiles( conditions, outfile )

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
        PipelineChipseq.subtractBedFiles( infile, subtract, outfile )

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
        PipelineChipseq.intersectBedFiles( replicates, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( (combineReplicates,
                 combineConditions, 
                 combineUnstim, 
                 normalizeBAMPerExperiment, 
                 normalizeBAMPerReplicate),
                suffix(".bed"),
                "_bed.import" )
    def importCombinedIntervals( infiles, outfile ):
        PipelineChipseq.importCombinedIntervals( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    ## master target for this section
    ############################################################
    @follows( importCombinedIntervals )
    def buildIntervals():
        pass

###################################################################
###################################################################
###################################################################
## Version 2: peak calling to get intervals
###################################################################
elif PARAMS["method"] == "bowtie-macs":

    ############################################################
    ############################################################
    ############################################################
    @follows( indexGenome )
    @transform( PARAMS["filename_reads"], 
                regex("(\S+).export.txt.gz"), 
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
        bowtie --sam %(intervals_bowtie_options)s %(intervals_bowtie_index)s - 2>%(outfile)s.log |\
        samtools import %(genome)s - %(tmpfilename)s >& %(outfile)s.log;
        samtools sort %(tmpfilename)s %(prefix)s;
        samtools index %(outfile)s;
        rm -f %(tmpfilename)s
        '''

        P.run( **dict( locals().items() + PARAMS.items() ) )

        if os.path.exists( tmpfilename ):
            os.unlink( tmpfilename )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAM,
                suffix(".bam"),
                ".readstats" )
    def buildBAMStats( infile, outfile ):
        PipelineChipseq.buildBAMStats( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAMStats,
                regex(r"(run.*).readstats"),
                inputs( (r"\1.bam", r"\1.readstats") ),
                r"\1.norm.bam" )
    def normalizeBAMPerReplicate( infiles, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.
        '''
        PipelineChipseq.buildNormalizedBAM( (infiles,), outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAMStats,
                regex(r"(control.*).readstats"),
                inputs( (r"\1.bam", r"\1.readstats") ),
                r"\1.norm.bam" )
    def normalizeBAMControls( infiles, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.
        '''
        PipelineChipseq.buildNormalizedBAM( (infiles,), outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( buildBAMStats )
    @files( [ ( [ ("%s%s.bam" % (x,y), "%s%s.readstats" % (x,y))  for y in REPLICATES],
                "%s.norm.bam" % x) for x in TRACKS_EXPERIMENTS ] )
    def normalizeBAMPerExperiment( infiles, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.

        Merge reads from several replicates.
        '''
        PipelineChipseq.buildNormalizedBAM( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( normalizeBAMControls )
    @transform( (normalizeBAMPerReplicate, normalizeBAMPerExperiment ),
               regex("run(.*).norm.bam"),
               r"run\1.macs" )
    def runMACS( infile, outfile ):
        '''run MACS for peak detection.'''
        PipelineChipseq.runMACS( infile, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runMACS,
                suffix(".macs"),
                "_macs.import" )
    def importMACS( infile, outfile ):
        PipelineChipseq.importMACS( infile, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @merge( runMACS, "macs.summary" )
    def summarizeMACS( infiles, outfile ):
        '''run MACS for peak detection.'''

        PipelineChipseq.summarizeMACS( infiles, outfile )


    ############################################################
    ############################################################
    ############################################################
    @transform( summarizeMACS,
                suffix(".summary"),
                "_summary.import" )
    def importMACSSummary( infile, outfile ):
        '''import macs summary.'''
        PipelineChipseq.importMACSSummary( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( importMACS, suffix("_macs.import"), ".bed" )
    def exportIntervalsAsBed( infile, outfile ):
        PipelineChipseq.exportIntervalsAsBed( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( exportIntervalsAsBed )
    @files( [ ( ["run%s%s.bed" % (y,x) for y in CELLLINES],
                "run%s.bed" % x, 
                x) for x in CONDITIONS ] )
    def combineConditions( infile, outfile, track ):
        '''combine conditions between cell lines. 

        The conditions are merged via intersection.
        '''

        conditions = [ "run%s%s.bed" % (x,track) for x in CELLLINES ]
        PipelineChipseq.intersectBedFiles( conditions, outfile )

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
        PipelineChipseq.subtractBedFiles( infile, subtract, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( (combineConditions, 
                 combineUnstim, 
                 normalizeBAMPerExperiment, 
                 normalizeBAMPerReplicate),
                suffix(".bed"),
                "_bed.import" )
    def importCombinedIntervals( infiles, outfile ):
        PipelineChipseq.importCombinedIntervals( infiles, outfile )

    @follows( exportIntervalsAsBed )
    def combineReplicates( ): pass

    ############################################################
    ############################################################
    ############################################################
    ## master target for this section
    ############################################################
    @follows( importCombinedIntervals, summarizeMACS )
    def buildIntervals():
        pass

###################################################################
###################################################################
###################################################################
## Version 3: peak calling with MACS from ELAND mapped reads
###################################################################
elif PARAMS["method"] == "eland":

    ############################################################
    ############################################################
    ############################################################
    @follows( indexGenome )
    @transform( PARAMS["filename_reads"], 
                regex("(\S+).export.txt.gz"), 
                r"\1.bam" )
    def buildBAM( infile, outfile ):
        '''build bam formatted files. The files will be sorted
        and indexed.

        The resultant BAM files contain

        1. all the uniquely mapped reads
        2. all the unmapped reads that have passed the quality filter.
        '''

        tmpfilename = P.getTempFilename()

        statement = '''
        /net/cpp-group/src/samtools/b64/samtools-0.1.6/misc/export2sam.pl \
        <(gunzip < %(infile)s) |\
        sed "s/.fa//" |\
        samtools import <( cut -f 1,4 %(genome)s.idx ) - %(tmpfilename)s
        '''

        P.run( **dict( locals().items() + PARAMS.items() ) )

        prefix = outfile[:-len(".bam")]
        pysam.sort( tmpfilename, prefix ) 
        pysam.index( outfile )

        os.unlink( tmpfilename )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAM,
                suffix(".bam"),
                ".readstats" )
    def buildBAMStats( infile, outfile ):
        PipelineChipseq.buildBAMStats( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAMStats,
                regex(r"(run.*).readstats"),
                inputs( (r"\1.bam", r"\1.readstats") ),
                r"\1.norm.bam" )
    def normalizeBAMPerReplicate( infile, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.
        '''
        track = infile[:-len(".bam")]
        PipelineChipseq.buildNormalizedBAM( ((infile,track + ".readstats"),), outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAMStats,
                regex(r"(control.*).readstats"),
                inputs( (r"\1.bam", r"\1.readstats") ),
                r"\1.norm.bam" )
    def normalizeBAMControls( infile, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.
        '''
        track = infile[:-len(".bam")]
        PipelineChipseq.buildNormalizedBAM( ((infile,track + ".readstats"),), outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( buildBAMStats )
    @files( [ ( [ ("%s%s.bam" % (x,y), "%s%s.readstats" % (x,y))  for y in REPLICATES],
                "%s.norm.bam" % x) for x in TRACKS_EXPERIMENTS ] )
    def normalizeBAMPerExperiment( infiles, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.

        Merge reads from several replicates.
        '''
        PipelineChipseq.buildNormalizedBAM( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( normalizeBAMControls )
    @transform( (normalizeBAMPerReplicate, normalizeBAMPerExperiment ),
               regex("run(.*).norm.bam"),
               r"run\1.macs" )
    def runMACS( infile, outfile ):
        '''run MACS for peak detection.'''
        PipelineChipseq.runMACS( infile, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runMACS,
                suffix(".macs"),
                "_macs.import" )
    def importMACS( infile, outfile ):
        PipelineChipseq.importMACS( infile, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @merge( runMACS, "macs.summary" )
    def summarizeMACS( infiles, outfile ):
        '''run MACS for peak detection.'''

        PipelineChipseq.summarizeMACS( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( summarizeMACS,
                suffix(".summary"),
                "_summary.import" )
    def importMACSSummary( infile, outfile ):
        '''import macs summary.'''
        PipelineChipseq.importMACSSummary( infile, outfile )


    ############################################################
    ############################################################
    ############################################################
    @transform( importMACS, suffix("_macs.import"), ".bed" )
    def exportIntervalsAsBed( infile, outfile ):
        PipelineChipseq.exportIntervalsAsBed( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( exportIntervalsAsBed )
    @files( [ ( ["run%s%s.bed" % (y,x) for y in CELLLINES],
                "run%s.bed" % x, 
                x) for x in CONDITIONS ] )
    def combineConditions( infile, outfile, track ):
        '''combine conditions between cell lines. 

        The conditions are merged via intersection.
        '''

        conditions = [ "run%s%s.bed" % (x,track) for x in CELLLINES ]
        PipelineChipseq.intersectBedFiles( conditions, outfile )

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
        PipelineChipseq.subtractBedFiles( infile, subtract, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( (combineConditions, 
                 combineUnstim, 
                 normalizeBAMPerExperiment, 
                 normalizeBAMPerReplicate),
                suffix(".bed"),
                "_bed.import" )
    def importCombinedIntervals( infiles, outfile ):
        PipelineChipseq.importCombinedIntervals( infiles, outfile )

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
## Version 4: peak calling with MACS from bowtie, combine by intersection
###################################################################
elif PARAMS["method"] == "bowtie-macs-replicate":

    ############################################################
    ############################################################
    ############################################################
    @follows( indexGenome )
    @transform( PARAMS["filename_reads"], 
                regex("(\S+).export.txt.gz"), 
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
        bowtie --sam %(intevals_bowtie_options)s %(intervals_bowtie_index)s - 2>%(outfile)s.log |\
        samtools import %(genome)s - %(tmpfilename)s >& %(outfile)s.log;
        samtools sort %(tmpfilename)s %(prefix)s;
        samtools index %(outfile)s;
        rm -f %(tmpfilename)s
        '''

        P.run()

        if os.path.exists( tmpfilename ):
            os.unlink( tmpfilename )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAM,
                suffix(".bam"),
                ".readstats" )
    def buildBAMStats( infile, outfile ):
        PipelineChipseq.buildBAMStats( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAMStats,
                regex(r"(run.*).readstats"),
                inputs( (r"\1.bam", r"\1.readstats") ),
                r"\1.norm.bam" )
    def normalizeBAMPerReplicate( infile, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.
        '''
        track = infile[:-len(".bam")]
        PipelineChipseq.buildNormalizedBAM( ((infile,track + ".readstats"),), outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( buildBAMStats )
    @files( [ ( [ ("%s%s.bam" % (x,y), "%s%s.readstats" % (x,y))  for y in REPLICATES],
                "%s.norm.bam" % x) for x in TRACKS_EXPERIMENTS ] )
    def normalizeBAMPerExperiment( infiles, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.


        Duplicated reads are removed at the same time.

        Merge reads from several replicates.
        '''
        PipelineChipseq.buildNormalizedBAM( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAMStats,
                regex(r"(control.*).readstats"),
                inputs( (r"\1.bam", r"\1.readstats") ),
                r"\1.norm.bam" )
    def normalizeBAMControls( infile, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.
        '''
        track = infile[:-len(".bam")]
        PipelineChipseq.buildNormalizedBAM( ((infile,track + ".readstats"),), outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( normalizeBAMControls )
    @transform( normalizeBAMPerReplicate,
                regex("run(.*).norm.bam"),
                r"run\1.macs" )
    def runMACS( infile, outfile ):
        '''run MACS for peak detection.'''
        PipelineChipseq.runMACS( infile, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runMACS,
                suffix(".macs"),
                "_macs.import" )
    def importMACS( infile, outfile ):
        PipelineChipseq.importMACS( infile, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @merge( runMACS, "macs.summary" )
    def summarizeMACS( infiles, outfile ):
        '''run MACS for peak detection.'''

        PipelineChipseq.summarizeMACS( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( summarizeMACS,
                suffix(".summary"),
                "_summary.import" )
    def importMACSSummary( infile, outfile ):
        '''import macs summary.'''
        PipelineChipseq.importMACSSummary( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( importMACS, suffix("_macs.import"), ".bed" )
    def exportIntervalsAsBed( infile, outfile ):
        PipelineChipseq.exportIntervalsAsBed( infile, outfile )

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
        PipelineChipseq.intersectBedFiles( replicates, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @follows( exportIntervalsAsBed, combineReplicates )
    @files( [ ( ["run%s%s.bed" % (y,x) for y in CELLLINES],
                "run%s.bed" % x, 
                x) for x in CONDITIONS ] )
    def combineConditions( infile, outfile, track ):
        '''combine conditions between cell lines. 

        The conditions are merged via intersection.
        '''

        conditions = [ "run%s%s.bed" % (x,track) for x in CELLLINES ]
        PipelineChipseq.intersectBedFiles( conditions, outfile )

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
        PipelineChipseq.subtractBedFiles( infile, subtract, outfile )

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
        PipelineChipseq.importCombinedIntervals( infiles, outfile )

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
## Version 4: peak calling with MACS from bowtie, combine by intersection of peaks, not intervals.
###################################################################
elif PARAMS["method"] == "bowtie-macs-peakreplicate":

    ############################################################
    ############################################################
    ############################################################
    @follows( indexGenome )
    @transform( PARAMS["filename_reads"], 
                regex("(\S+).export.txt.gz"), 
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

        P.run( **dict( locals().items() + PARAMS.items() ) )

        if os.path.exists( tmpfilename ):
            os.unlink( tmpfilename )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAM,
                suffix(".bam"),
                ".reastats" )
    def buildBAMStats( infile, outfile ):
        PipelineChipseq.buildBAMStats( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAMStats,
                regex(r"(run.*).readstats"),
                inputs( (r"\1.bam", r"\1.readstats") ),
                r"\1.norm.bam" )
    def normalizeBAMPerReplicate( infile, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.
        '''
        track = infile[:-len(".bam")]
        PipelineChipseq.buildNormalizedBAM( ((infile,track + ".readstats"),), outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( buildBAMStats )
    @files( [ ( [ ("%s%s.bam" % (x,y), "%s%s.readstats" % (x,y))  for y in REPLICATES],
                "%s.norm.bam" % x) for x in TRACKS_EXPERIMENTS ] )
    def normalizeBAMPerExperiment( infiles, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.


        Duplicated reads are removed at the same time.

        Merge reads from several replicates.
        '''
        PipelineChipseq.buildNormalizedBAM( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAMStats,
                regex(r"(control.*).readstats"),
                inputs( (r"\1.bam", r"\1.readstats") ),
                r"\1.norm.bam" )
    def normalizeBAMControls( infile, outfile ):
        '''build a normalized BAM file such that all
        files have approximately the same number of 
        reads.

        Duplicated reads are removed at the same time.
        '''
        track = infile[:-len(".bam")]
        PipelineChipseq.buildNormalizedBAM( ((infile,track + ".readstats"),), outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( normalizeBAMControls )
    @transform( normalizeBAMPerReplicate,
                regex("run(.*).norm.bam"),
                r"run\1.macs" )
    def runMACS( infile, outfile ):
        '''run MACS for peak detection.'''
        PipelineChipseq.runMACS( infile, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runMACS,
                suffix(".macs"),
                "_macs.import" )
    def importMACS( infile, outfile ):
        PipelineChipseq.importMACS( infile, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @merge( runMACS, "macs.summary" )
    def summarizeMACS( infiles, outfile ):
        '''run MACS for peak detection.'''
        PipelineChipseq.summarizeMACS( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( summarizeMACS,
                suffix(".summary"),
                "_summary.import" )
    def importMACSSummary( infile, outfile ):
        '''import macs summary.'''
        PipelineChipseq.importMACSSummary( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( importMACS, suffix("_macs.import"), ".bed" )
    def exportIntervalsAsBed( infile, outfile ):
        PipelineChipseq.exportPeaksAsBed( infile, outfile )

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
        PipelineChipseq.intersectBedFiles( replicates, outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( exportIntervalsAsBed, combineReplicates )
    @files( [ ( ["run%s%s.bed" % (y,x) for y in CELLLINES],
                "run%s.bed" % x, 
                x) for x in CONDITIONS ] )
    def combineConditions( infile, outfile, track ):
        '''combine conditions between cell lines. 

        The conditions are merged via intersection.
        '''

        conditions = [ "run%s%s.bed" % (x,track) for x in CELLLINES ]
        PipelineChipseq.intersectBedFiles( conditions, outfile )

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
        PipelineChipseq.subtractBedFiles( infile, subtract, outfile )

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
        PipelineChipseq.importCombinedIntervals( infiles, outfile )

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
## Version5: peak calling from bed intervals.
###################################################################
elif PARAMS["method"] == "bed-macs-replicate":

    BAM_SUFFIX = ".bam"
    ############################################################
    ############################################################
    ############################################################
    @transform( "run*.bed.gz",
               suffix(".bed.gz"),
               ".bam" )
    def buildBAM( infile, outfile ):
        '''build bam formatted files from bed files.
        The resultant files are sorted and indexed.
        The BAM files are filled with dummy sequences
        and quality scores.
        '''

        tmpfilename1 = P.getTempFilename()
        tmpfilename2 = P.getTempFilename()

        outf = open( tmpfilename1, "w" )
        fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )
        
        # convert to SAM
        for x,b in enumerate( Bed.iterator( gzip.open( infile ) )):
            l = b.end - b.start
            strand = b.fields[2]
            if strand == "+": flag = "0"
            elif strand == "-": flag = "16"
            else: raise ValueError("unknown strand in %s" % str(b))
            outf.write( "\t".join( 
                ( "read%i" % x,
                  flag,
                  b.contig,
                  str(b.start+1),
                  "10",
                  "%iM" % l,
                  "*", "0", "0",
                  "A" * l,
                  "A" * l,
                  "MD:Z:%i" % l )) + "\n" )
        
        outf.close()
            
        statement = '''
        samtools import <( cut -f 1,4 %(genome)s.idx ) %(tmpfilename1)s %(tmpfilename2)s >& %(outfile)s.log
        '''
        P.run( **dict( locals().items() + PARAMS.items() ) )

        prefix = outfile[:-len(".bam")]
        pysam.sort( tmpfilename2, prefix ) 
        pysam.index( outfile )

        os.unlink( tmpfilename1 )
        os.unlink( tmpfilename2 )

    ############################################################
    ############################################################
    ############################################################
    @transform( buildBAM,
                suffix(".bam"),
                ".macs" )
    def runMACS( infile, outfile ):
        '''run MACS for peak detection.'''
        PipelineChipseq.runMACS( infile, outfile )
        
    ############################################################
    ############################################################
    ############################################################
    @transform( runMACS,
                suffix(".macs"),
                "_macs.import" )
    def importMACS( infile, outfile ):
        PipelineChipseq.importMACS( infile, outfile, suffix = ".bam" )
        
    ############################################################
    ############################################################
    ############################################################
    @merge( runMACS, "macs.summary" )
    def summarizeMACS( infiles, outfile ):
        '''run MACS for peak detection.'''

        PipelineChipseq.summarizeMACS( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( summarizeMACS,
                suffix(".summary"),
                "_summary.import" )
    def importMACSSummary( infile, outfile ):
        '''import macs summary.'''
        PipelineChipseq.importMACSSummary( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @transform( importMACS, suffix("_macs.import"), ".bed" )
    def exportIntervalsAsBed( infile, outfile ):
        PipelineChipseq.exportIntervalsAsBed( infile, outfile )

    ############################################################
    ############################################################
    ############################################################
    @follows( exportIntervalsAsBed )
    def combineReplicates( ): pass
    
    @follows( exportIntervalsAsBed, combineReplicates )
    def combineConditions( ): pass

    @follows( combineConditions )
    def combineUnstim(): pass

    @follows( exportIntervalsAsBed )
    def normalizeBAMPerExperiment(): pass
    
    @follows( exportIntervalsAsBed )
    def normalizeBAMPerReplicate(): pass

    ############################################################
    ############################################################
    ############################################################
    @transform( ( combineReplicates,
                  combineConditions, 
                  combineUnstim ),
                suffix(".bed"),
                "_bed.import" )
    def importCombinedIntervals( infiles, outfile ):
        PipelineChipseq.importCombinedIntervals( infiles, outfile )

    ############################################################
    ############################################################
    ############################################################
    ## master target for this section
    ############################################################
    @follows( exportIntervalsAsBed,
              importCombinedIntervals,
              importMACSSummary )
    def buildIntervals(): pass

###################################################################
###################################################################
###################################################################
## general import
###################################################################

############################################################
############################################################
############################################################
@follows( exportIntervalsAsBed )
@files_re( ["%s.bed" % x for x in TRACKS_RAW],
          combine( "(.*).bed" ),
          "merged.bed" )
def makeMerged( infiles, outfile ):
    '''combine all experiments.

    The replicates are combined using a merge.
    '''
    PipelineChipseq.mergeBedFiles( infiles, outfile )

############################################################
############################################################
############################################################
@files(( 
        ( (PARAMS["filename_map_ensembl2refseq"], 
           PARAMS["filename_ensembl_geneset"], 
           PARAMS["filename_refseq_geneset"]), TARGET_TRANSCRIPTS), 
        ),)
def buildTranscripts( infiles, outfile ):
    '''build a file with all transcripts. 

    remove duplicate transcripts (within the refseq track there usually are).
    '''
    map_ensembl2refseq, ensembl, refseq = infiles
    to_cluster = True

    statement = '''
    cat %(ensembl)s | gunzip |
    awk '$3 == "exon"' |\
    sed "s/hg18_refGene/protein_coding/" |\
    python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s --log=%(outfile)s.log |\
    python %(scriptsdir)s/gtf2gtf.py --remove-duplicates=gene --log=%(outfile)s.log |\
    gzip > %(outfile)s'''
    P.run( **dict( locals().items() + PARAMS.items() ) )

    statement = '''
    cat %(refseq)s | gunzip |
    awk '$3 == "exon"' |\
    sed "s/hg18_refGene/protein_coding/" |\
    python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s --log=%(outfile)s.log |\
    python %(scriptsdir)s/gtf2gtf.py --filter=gene --apply=<(cut -f 3 %(map_ensembl2refseq)s | sort | uniq) --invert-filter |\
    python %(scriptsdir)s/gtf2gtf.py --remove-duplicates=gene --log=%(outfile)s.log |\
    gzip > %(outfile)s'''
    # disabled - only work with ENSEMBL
    # P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@files( PARAMS["filename_ensembl_geneset"], TARGET_ANNOTATOR_GENETERRITORIES )
def buildAnnotatorGeneTerritories( infile, outfile ):
    '''build gene territories.'''

    to_cluster=True
    
    radius = PARAMS['annotator_gene_territories_radius']
    statement = '''
    gunzip < %(infile)s |\
        awk '$2 == "protein_coding"' |\
        python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s |\
        python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --with-utr --log=%(outfile)s.log |\
        python %(scriptsdir)s/gtf2gtf.py --filter=longest-gene --log=%(outfile)s.log |\
        %(scriptsdir)s/gff_sort pos |\
        python %(scriptsdir)s/gtf2gff.py \
                                --genome-file=%(genome)s --log=%(outfile)s.log \
                                --radius=%(radius)s --method=territories > %(outfile)s'''

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@files_re( buildTranscripts, "(.*).gtf.gz", r"\1_gtf.import")
def importTranscripts( infile, outfile ):
    '''import the transrcipt set.'''
    table = "transcripts"

    statement = '''
    gunzip < %(infile)s |\
    python %(scriptsdir)s/gtf2tsv.py |\
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=interval_id \
              --index=transcript_id \
              --index=gene_id \
              --table=%(table)s \
    > %(outfile)s'''

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@files( PARAMS["filename_repeats"], TARGET_REPEATS)
def buildRepeats( infile, outfile ):
    '''import ucsc repeat tracks.
    '''

    instream = gzip.open( infile )
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    tmpfilename = tmpfile.name

    gff = GFF.Entry()
    gff.feature = "exon"
    gff.source = "repeat"
    
    for line in instream:
        if line.startswith("#bin"): continue
        (bin,swScore,milliDiv,milliDel,milliIns,\
         genoName,genoStart,genoEnd,genoLeft,\
         strand,repName,repClass,repFamily,repStart,repEnd,repLeft,id) =\
         line[:-1].split("\t")
         
        if repClass not in REPEATS: continue

        gff.strand = strand
        gff.contig = genoName
        gff.start = int(genoStart)-1
        gff.end = int(genoEnd)
        gff.addAttribute( "name", repName )
        gff.addAttribute( "class", repClass )
        gff.addAttribute( "family", repFamily )
        
        tmpfile.write ( str(gff) + "\n" )
         
    tmpfile.close()

    E.info( "created temporary file" )
    statement = """sort -k1,1 -k4,4n < %(tmpfilename)s |\
	python %(scriptsdir)s/gff2gff.py --merge-features=0,10,0,0 --log=%(outfile)s.log \
	> %(outfile)s"""

    P.run( **dict( locals().items() + PARAMS.items() ) )
    
    os.unlink( tmpfile.name )
        
############################################################
############################################################
############################################################
@files( PARAMS["filename_ensembl_geneset"], TARGET_ANNOTATION )
def buildGeneRegions( infile, outfile ):
    '''annotate genomic regions with reference gene set.

    Only considers protein coding genes. In case of overlapping
    genes, only take the longest (in genomic coordinates).
    Genes not on UCSC contigs are removed.
    '''
    PGeneset.buildGeneRegions( infile, outfile, only_protein_coding = True )

############################################################
############################################################
############################################################
@files( PARAMS["filename_ensembl_geneset"], TARGET_GENESET )
def buildGeneSet( infile, outfile ):
    '''build a collection of exons from the protein-coding
    section of the ENSEMBL gene set. The exons include both CDS
    and UTR.

    The set is filtered in the same way as in :meth:`buildGeneRegions`.
    '''
    PGeneset.buildGeneSet( infile, outfile, only_protein_coding = True )

############################################################
############################################################
############################################################
@files( buildGeneSet, "gene_stats.import" )
def importGeneStats( infile, outfile ):
    '''import the transrcipt set.'''
    PGeneset.importGeneStats( infile, outfile )

############################################################
############################################################
############################################################
@files( [ ("%s.idx" % PARAMS["filename_genome"], "genome.gff"), ] )
def buildGenome( infile, outfile ):

    fasta = IndexedFasta.IndexedFasta( PARAMS["filename_genome"] )

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
@follows( buildTranscripts )
@files( ( (PARAMS["filename_ensembl_geneset"], TARGET_PROMOTORS ),
          (TARGET_TRANSCRIPTS, TARGET_TRANSCRIPTS_PROMOTORS), ) )
def buildPromotorRegions( infile, outfile ):
    '''annotate promotor regions from reference gene set.'''
    statement = """
        gunzip < %(infile)s |\
        python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s --log=%(outfile)s.log |\
        python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=%(promotor_size)s \
                              --genome-file=%(genome)s --log=%(outfile)s.log > %(outfile)s
    """

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@follows( buildTranscripts )
@files( ( (PARAMS["filename_ensembl_geneset"], TARGET_TSS), 
          (TARGET_TRANSCRIPTS, TARGET_TRANSCRIPTS_TSS), ) )
def buildTSSRegions( infile, outfile ):
    '''annotate transcription start sites from reference gene set.

    Similar to promotors, except that the witdth is set to 1.
    '''
    statement = """
        gunzip < %(infile)s |\
        python %(scriptsdir)s/gff2gff.py --sanitize=genome --skip-missing --genome-file=%(genome)s --log=%(outfile)s.log |\
        python %(scriptsdir)s/gtf2gff.py --method=promotors --promotor=1 --genome-file=%(genome)s --log=%(outfile)s.log > %(outfile)s
    """

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@follows( buildPromotorRegions, buildGeneRegions, buildTSSRegions )
@files( [ ("%s.gtf" % x, "%s.bed" % x ) for x in ("ensembl", "promotors", "tss" ) ] )
def exportReferenceAsBed( infile, outfile ):
    '''export a reference gtf file as bed for computing overlap.'''

    bed = Bed.Bed()
    outfile = open( outfile, "w" )
    with open(infile, "r") as inf:
        for gff in GTF.iterator( infile ):
            bed.contig, bed.start, bed.end = gff.contig, gff.start, gff.end
            bed.fields = [ gff.gene_id ]
            outfile.write( "%s\n" % str(bed) )
    outfile.close()

############################################################
############################################################
############################################################
@files_re( (buildBAM, normalizeBAMPerExperiment, normalizeBAMPerReplicate),
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
##
############################################################
@files( ((PARAMS["filename_affymetrix"], "affymetrix_import" ),) )
def importAffymetrixAnnotation( infile, outfile):
    '''import affymetrix probe annotations.
    '''
    table = "affymetrix"

    statement = '''
    gunzip < %(infile)s |\
    sed "s/probeset_id/probeset/" |\
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --dialect=excel \
              --index=probeset \
              --index=mrna_id \
              --table=%(table)s \
    > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
##
############################################################
@files( importAffymetrixAnnotation, "probeset2transcript.table")
def buildProbeset2Transcript( infile, outfile):
    '''build map relating a probeset to a transcript_id, gene_id, ...

    1. only take those that can be mapped to transcripts in table transcripts
    2. check consistency of location with transcripts in table transcripts
    3. remove those with cross hybridyzation signals

    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    E.debug( "collecting exons")
    raise NotImplementError("check fastgtf and gzip.")
    with gzip.open(PARAMS["transcripts"]) as inf:
        exons = GTF.readAndIndex( GTF.iterator( inf ) )

    E.debug( "collecting transcripts")
    cc = dbhandle.cursor()
    statement = "SELECT transcript_id, gene_id, source, contig, min(start), max(end), strand FROM transcripts GROUP BY transcript_id, gene_id, source, contig"
    cc.execute( statement )
    map_transcript2gene = {}
    for transcript_id, gene_id, source, contig, start, end, strand in cc.fetchall():
        map_transcript2gene[transcript_id] = (gene_id, source, contig, start, end, strand )

    counter = E.Counter()

    cc = dbhandle.cursor()
    statement = "SELECT probeset, transcript_cluster_id, seqname, start, stop, mrna_assignment, crosshyb_type FROM affymetrix"
    cc.execute( statement )

    outf_notfound = open( outfile + ".notfound", "w")
    outf_notfound.write("probeset\tcause\tcontig\tstart\tend\tmrna_assignment\n" )

    outf = open( outfile, "w" )
    outf.write( "probeset\tcluster_id\tprobe_contig\tprobe_start\tprobe_end\ttranscript_id\tgene_id\tsource\tcontig\tstart\tend\tstrand\n" )

    genes_found, transcripts_found = set(), set()

    for probeset, transcript_cluster_id, contig, start, end, mrna_assignment, crosshyb_type in cc.fetchall():
        start -= 1 
        if mrna_assignment == "---":
            counter.probesets_unassigned += 1
            continue

        data = [ x.strip().split(" // ") for x in mrna_assignment.split("///") ]
        found, mapped, mismapped, nxhyb = 0, 0, 0, 0

        # collect transcripts by exon overlap
        try:
            overlaps = [ x[2] for x in exons.get( contig, start, end ) ]
        except KeyError:
            overlaps = []

        transcripts = set( [x.transcript_id for x in overlaps ] )

        for transcript_id, seqname, score, direct_probes, possible_probes, xhyb in data:

            if xhyb != "0":
                nxhyb += 1
                continue

            # collect transcripts from annotation file
            if transcript_id.startswith("ENST"):
                found += 1
                    
                try:
                    gene_id, other_source, other_contig, other_start, other_end, other_strand = map_transcript2gene[transcript_id]
                except KeyError:
                    continue

                mapped += 1

                if transcript_id not in transcripts:
                    mismapped += 1
                    continue
                
                if other_contig != contig or \
                        min(other_end,end) - max(other_start, start) < 0:
                    mismapped += 1
                    continue

                outf.write("%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\t%s\t%i\t%i\t%s\n" % \
                               (probeset, transcript_cluster_id, contig, start, end,
                                transcript_id, gene_id, other_source,
                                other_contig, other_start, other_end, other_strand))

                genes_found.add(gene_id)
                transcripts_found.add(transcript_id)

        if found == 0: 
            counter.probesets_notfound += 1
            code = "no ENSEMBL"
        elif mapped == 0:
            counter.probesets_unmapped += 1
            code = "not in transcripts"
        elif mapped > 0 and mapped == mismapped:
            counter.probesets_mismapped += 1
            code = "mismapped"
        elif mapped > 0 and mapped == nxhyb:
            counter.probesets_crosshyb += 1
            code = "cross-hybridized"
        elif mapped > 0 and mapped == nxhyb + mismapped:
            counter.probesets_errors += 1
            code = "errors"
        else:
            counter.probesets_found += 1
            code = None

        if code:
            outf_notfound.write("%s\t%s\t%s\t%i\t%i\t%s\n" % (probeset, code, contig, start, end, mrna_assignment) )

    outf_notfound.close()
    outf.close()

    E.info( "mapping results: %s" % str(counter))
    E.info( "genes found: %i(%5.2f%%): transcripts found %i(%5.2f%%)" %\
                (len(genes_found),
                 100.0 * len(genes_found) / len( set( [x[0] for x in map_transcript2gene.values()]) ),
                 len(transcripts_found),
                 100.0 * len(transcripts_found) / len( map_transcript2gene.keys() ) ) )

############################################################
############################################################
############################################################
##
############################################################
@transform( buildProbeset2Transcript, suffix(".table"), "_table.import" )
def importProbeset2Transcript( infile, outfile):

    table = outfile[:-len("_table.import")] 

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=probeset \
              --index=cluster_id \
              --index=transcript_id \
              --index=gene_id \
              --table=%(table)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
## OBSOLETE
@files( ((importProbeset2Transcript, TARGET_PROBESET),) )
def exportProbesetLocations( infile, outfile ):
    '''export probeset locations.
    
    Probeset locations map the range of a transcript onto
    a probeset.
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    # statement = "SELECT contig, COUNT(DISTINCT contig), min(start), max(end), t.strand, e.probeset, COUNT(DISTINCT t.transcript_id) FROM transcripts AS t, expression as e WHERE t.transcript_id = e.mrna_id GROUP BY t.strand, e.probeset" 
    statement = "SELECT contig, 1, start, end, strand, cluster_id, transcript_id, gene_id, COUNT(DISTINCT gene_id) FROM probeset2transcript GROUP BY cluster_id, transcript_id"
    cc.execute( statement )

    gff = GTF.Entry()
    gff.source = "cluster_id"
    gff.feature = "exon"

    outs = open(outfile, "w")
    noutput, nskipped = 0,0
    max_len = int(PARAMS["max_probeset_length"])
    for contig, ncontigs, start, end, strand, cluster_id, transcript_id, gene_id, ngenes in cc:
        if ncontigs > 1:
            E.warn( "cluster_id %s: multiple contigs" % cluster_id)
            nskipped += 1
            continue
        if ngenes > 1:
            E.warn( "cluster_id %s: multiple genes" % cluster_id)
            nskipped += 1
            continue
        if end - start > max_len:
            E.warn( "cluster_id %s: range (%s:%i..%i) larger than cutoff (%i>%i)" % \
                        (cluster_id, contig, start, end, end-start, max_len))
            nskipped += 1
            continue
            
        gff.contig = contig
        gff.start = start 
        gff.end = end
        gff.strand = strand
        gff.gene_id = transcript_id
        gff.transcript_id = cluster_id
        outs.write( str(gff) + "\n" )
        noutput += 1
    outs.close()

    E.info( "%i cluster_ids exported: %i skipped" %(noutput, nskipped) )

############################################################
############################################################
############################################################
## OBSOLETE
@transform( exportProbesetLocations, suffix(".gtf"), "_gtf.import")
def importProbesetLocations( infile, outfile ):
    '''import probeset locations.'''

    tablename = outfile[:-len(".import")]

    statement = '''
    python %(scriptsdir)s/gtf2tsv.py < %(infile)s |\
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=id \
              --index=probeset \
              --table=%(tablename)s \
    > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
##
############################################################
@split( (PARAMS["filename_expression"], PARAMS["filename_expression2"],importProbeset2Transcript),
        ("expression.skipped", "expression.table",  "expression.map", "exp*.data"))
def buildExpressionTracks( infiles, outfiles ):
    '''build expression tracks from FILENAME_EXPRESSION and FILENAME_EXPRESSION2'''

    # start with a clean sheet
    for f in outfiles:
        try:
            os.remove( f )
        except OSError:
            pass


    outf = open( "expression.map", "w" )
    outf.write("tablename\tcellline\ttime\tstimulus\n" )

    for index, infile in enumerate(infiles):
        map_exp2columns = collections.defaultdict( list )
        if infile == PARAMS["filename_expression"]:
            map_exp2columns = { 'GM00855_00_D3' : (3,4,5),
                                'GM00855_36_D3' : (6,7,8),
                                'GM00855_36_Est' : (9,10,11),
                                'GM00855_36_D3Est' : (12,13,14),
                                'GM00861_00_D3' : (15,16,17),
                                'GM00861_36_D3' : (18,19,20),
                                'GM00861_36_Est' : (21,22,23),
                                'GM00861_36_D3Est' : (24,25,26),
                                'All_00_D3' : (3,4,5,15,16,17),
                                'All_36_D3' : (6,7,8,18,19,20),
                                'All_36_Est' : (9,10,11,21,22,23),
                                'All_36_D3Est' : (12,13,14,24,25,26),
                                }
            
            PExpression.buildExpressionTracks( infile, outfiles, map_exp2columns, suffix = "" )
            for key in map_exp2columns.keys():
                cellline, time, stimulus = key.split("_")
                outf.write( "\t".join( (key, cellline, time, stimulus) ) + "\n" )
                
        elif infile == PARAMS["filename_expression2"]:
            x = 1
            stimulus = "D3"
            # pool lymphoblastoid cell lines (lcl) and others separately

            for cell in ('GMO7019', 'GMO7348', 'GM18054',
                         'Ag09309', 'HL60', 'K562', 
                         'HEPG2', 'AG09319', 'GM12878'):
                for timepoint in ('00', '08', '36'):
                    key = "%s_%s_%s" % (cell, timepoint, stimulus)
                    map_exp2columns[ key ] = (x,x+1)
                    if cell.startswith( "GM" ):
                        key2 = "%s_%s_%s" % ("LCL", timepoint, stimulus)
                    else:
                        key2 = "%s_%s_%s" % ("nonLCL", timepoint, stimulus)
                                            
                    map_exp2columns[key2].extend( [x, x+1] )
                    x += 2

            PExpression.buildExpressionTracks( infile, outfiles, map_exp2columns, suffix = ".%0i" % index )
            for key in map_exp2columns.keys():
                cellline, time, stimulus = key.split("_")
                outf.write( "\t".join( (key, cellline, time, stimulus) ) + "\n" )

    outf.close()
        
############################################################
############################################################
############################################################
##
############################################################
@transform( buildExpressionTracks, suffix(".data"), "_data.import" )
def importExpressionTracks( infile, outfile):

    table = outfile[:-len("_data.import")] + "_levels"

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=cluster_id \
              --table=%(table)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
## 
############################################################
@merge( buildExpressionTracks, "expression.correlation")
def buildExpressionCorrelation( infiles, outfile):
    '''compute correlation between replicates in expression sets.'''

    outf = open(outfile, "w")
    outf.write( "track\treplicate1\treplicate2\t%s\n" % "\t".join(Stats.CorrelationTest().getHeaders()))
    offset = 10
    for filename in infiles:
        if not filename.endswith(".data"): continue
        vals = []
        headers = []
        track = filename[:-len(".data")]
        inf = open( filename, "r")
        for line in inf:
            if line.startswith("#"): continue
            data = line[:-1].split("\t")
            if line.startswith("cluster_id"):
                nreplicates = len(data) - offset
                for x in range(nreplicates):
                    headers.append( data[offset+x] )
                    vals.append( [] )
                continue
            for x in range( 0, nreplicates ):
                vals[x].append( float(data[offset+x]) )

        for x in range(nreplicates-1):
            for y in range( x+1, nreplicates):
                outf.write( "\t".join( (track,
                                        headers[x],
                                        headers[y],
                                        str( Stats.doCorrelationTest( vals[x],
                                                                      vals[y]) ) ) ) + "\n" )
    outf.close()

############################################################
############################################################
############################################################
## 
############################################################
@merge( buildExpressionTracks, "expression_full.correlation")
def buildExpressionFullCorrelation( infiles, outfile):
    '''compute correlation between all replicates in all expression sets.
    '''

    outf = open(outfile, "w")
    outf.write( "track1\treplicate1\ttrack2\ttreplicate2\t%s\n" % "\t".join(Stats.CorrelationTest().getHeaders()))

    nfiles = len(infiles)

    def _readVals( filename ):
        '''read replicate values from file.'''
        vals = []
        headers = []
        offset = 10

        lines = open( filename, "r").readlines()
        data = lines[0][:-1].split("\t")
        nreplicates = len(data) - offset
        for x in range(nreplicates):
            headers.append( data[offset+x] )
            vals.append( [] )
        del lines[0]
        lines.sort()
        for line in lines:
            data = line[:-1].split("\t")
            for x in range( 0, nreplicates ):
                vals[x].append( float(data[offset+x]) )
        return nreplicates, vals, headers
    
    for f1 in range( nfiles - 1):
        filename1 = infiles[f1]
        if not filename1.endswith(".data"): continue
        track1 = filename1[:-len(".data")]

        nreplicates1, vals1, headers1 = _readVals( filename1 )
        
        # do within track correlation
        for x in range(nreplicates1-1):
            for y in range( x+1, nreplicates1):
                outf.write( "\t".join( (track1,
                                        headers1[x],
                                        track1,
                                        headers1[y],
                                        str( Stats.doCorrelationTest( vals1[x],
                                                                      vals1[y]) ) ) ) + "\n" )
        
        for f2 in range(f1+1, nfiles):
            filename2 = infiles[f2]
            if not filename2.endswith(".data"): continue
            track2 = filename2[:-len(".data")]

            nreplicates2, vals2, headers2 = _readVals( filename2 )
            
            for x in range(nreplicates1):
                for y in range(nreplicates2):
                    outf.write( "\t".join( (track1,
                                            headers1[x],
                                            track2,
                                            headers2[y],
                                            str( Stats.doCorrelationTest( vals1[x],
                                                                          vals2[y]) ) ) ) + "\n" )
            
            
    outf.close()

############################################################
############################################################
############################################################
##
############################################################
@transform( (buildExpressionCorrelation, buildExpressionFullCorrelation),
            suffix(".correlation"), "_correlation.import" )
def importExpressionCorrelation( infile, outfile):

    table = outfile[:-len(".import")] 

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=track \
              --table=%(table)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
##
############################################################
@transform( buildExpressionTracks, suffix(".table"), "_table.import" )
def importExpressionProbesets( infile, outfile):

    table = outfile[:-len("_table.import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=cluster_id \
              --index=mrna_id \
              --index=genesymbol \
              --table=%(table)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
##
############################################################
@transform( buildExpressionTracks, suffix(".map"), "_map.import" )
def importExpressionMap( infile, outfile):

    table = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --table=%(table)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
# @follows( importCombinedIntervals, exportIntervalsAsBed, exportProbesetLocations )
@files_re( [ x for x in glob.glob("run*.bed") if not isReplicate( x[:-len("bed")] ) ],
           "(.*).bed", 
           (r"\1.bed", TARGET_PROBESET),
           r"\1.assoc" )
def buildIntervalProbesetAssociations( infiles, outfile ):
    '''build association between probesets and intervals.'''

    to_cluster = True
    infile, probeset = infiles
    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python %(scriptsdir)s/gtf2table.py \
		--counter=neighbours \
		--log=%(outfile)s.log \
		--filename-gff=%(probeset)s \
                --proximal-distance=%(proximal_distance)s \
		--genome-file=%(genome)s |\
        python %(toolsdir)s\table2graph.py \
                --headers=id,probeset \
	> %(outfile)s"""

    P.run()

############################################################
############################################################
############################################################
@files_re( buildIntervalProbesetAssociations,
           "(.*).assoc",
           r"\1_assoc.import")
def importIntervalProbesetAssociations( infile, outfile):
    tablename = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=id \
              --index=probeset \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
# fix: does not work due to exportIntervalsAsBed
# (importCombinedIntervals, exportIntervalsAsBed ),
@transform( [ "%s.bed" % x for x in TRACKS_MOTIFS ],
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
@transform( [ "%s.bed" % x for x in TRACKS_MOTIFS],
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

############################################################
############################################################
############################################################
@merge( (exportReferenceAsBed, exportIntervalsAsBed, 
         ["%s.bed" % x for x in TRACKS_REFERENCE ] ),
        "intervals.overlap" )
def buildOverlap( infiles, outfile ):
    '''compute overlap between intervals.
    '''

    if os.path.exists(outfile): 
        os.rename( outfile, outfile + ".orig" )
        options = "--update=%s.orig" % outfile
    else:
        options = ""

    infiles = " ".join(IOTools.flatten(infiles))

    statement = '''
        python %(scriptsdir)s/diff_bed.py 
               --log=%(outfile)s.log 
               %(options)s 
               %(infiles)s 
        > %(outfile)s
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
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=set1 
              --index=set2 
              --table=%(tablename)s 
    < %(infile)s 
    > %(outfile)s
    '''

    P.run()

############################################################
############################################################
############################################################
@merge( (exportUCSCEncodeTracks, 
         exportReferenceAsBed, 
         ["%s.bed" % x for x in TRACKS_REFERENCE ],
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

    infiles = " ".join(IOTools.flatten(infiles))
    statement = '''
        python %(scriptsdir)s/diff_bed.py 
               --tracks 
               --log=%(outfile)s.log
               %(options)s 
               %(infiles)s 
        > %(outfile)s
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
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
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
    PipelineChipseq.makeIntervalCorrelation( infiles, outfile, "peakval" )

@follows( makeMerged )
@files_re( ["%s.bed" % x for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "avgval.correlation" )
def makeAvgvalCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PipelineChipseq.makeIntervalCorrelation( infiles, outfile, "avgval" )

@follows( makeMerged )
@files_re( ["%s.bed" % x for x in TRACKS_CORRELATION],
           combine("(.*).bed"),
           "length.correlation" )
def makeLengthCorrelation( infiles, outfile ):
    '''compute correlation of interval properties between sets for field peakval.
    '''
    PipelineChipseq.makeIntervalCorrelation( infiles, outfile, "length" )

@transform( ( makePeakvalCorrelation, makeAvgvalCorrelation, makeLengthCorrelation ),
            suffix(".correlation"),
            "_correlation.import")
def importCorrelation( infile, outfile ):

    tablename = outfile[:-len(".import")] 

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
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
    %(motifs_execmeme)s %(tmpname)s -mod oops -dna -revcomp -nmotifs 1 -text -oc %(tmpdir)s > %(outfile)s
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
    maxsize = int(PARAMS["motifs_meme_max_size"])

    E.info( "runMeme: %s: using at most %i sequences for pattern finding" % (track, cutoff) )

    fasta = IndexedFasta.IndexedFasta( PARAMS["genome"] )

    tmpdir = tempfile.mkdtemp( dir = "." )
    tmpfasta =  os.path.join( tmpdir, "in.fa")
    outs = open( tmpfasta , "w" )

    masker = Masker.MaskerDustMasker()

    hw = int(PARAMS["motifs_meme_halfwidth"])
    current_size, nseq = 0, 0
    for peakcenter, id, contig in data[:cutoff]:
        start, end = peakcenter - hw, peakcenter + hw
        id = "%s_%s %s:%i..%i" % (track, str(id), contig, start, end)
        seq = fasta.getSequence( contig, "+", start, end )

        if PARAMS["motifs_meme_masker"] == "repeatmasker":
            # the genome sequence is repeat masked
            masked_seq = seq
        elif PARAMS["motifs_meme_masker"] == "dustmasker":
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
    %(motifs_execmeme)s %(tmpfasta)s -dna -revcomp 
                 -mod %(motifs_meme_model)s 
                 -nmotifs %(motifs_meme_nmotifs)s 
                 -oc %(tmpdir)s 
                 -maxsize %(maxsize)s 
                 %(motifs_meme_options)s 
    > %(outfile)s.log
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
                                       halfwidth = int(PARAMS["motifs_meme_halfwidth"]),
                                       maxsize = int(PARAMS["motifs_meme_max_size"]),
                                       proportion = 0.1 )

    min_sequences = int(nseq / 10.0)
    statement = '''
    %(motifs_execglam2)s -2 -O %(tmpdir)s %(motifs_glam2_options)s -z %(min_sequences)i n %(tmpfasta)s > %(outfile)s.log
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
           %(motifs_exectomtom)s -text -query %(motifs_tomtom_master_motif)s -target %(infile)s > %(outfile)s
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
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
                  --allow-empty \
                  --table=%(tablename)s \
        < %(tmpname)s > %(outfile)s
        '''

        P.run()
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
    def importTomTom(): pass

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
# todo: fix, causes a problem: exportSequences and exportMASTControlSequences,
@follows( makeMotifs, filterMotifs )
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

        of = gzip.open(outfile, "a")
        motif, x = os.path.splitext( motiffile )
        of.write(":: motif = %s ::\n" % motif )
        of.close()

        statement = '''
        cat %(dbfile)s %(controlfile)s 
        | %(motifs_execmast)s %(motiffile)s -stdin -stdout -text -ev %(motifs_mast_evalue)f 
        | gzip 
        >> %(outfile)s
        '''
        P.run()

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
    %(motifs_execbioprospector)s -i %(tmpfasta)s %(motifs_bioprospector_options)s -o %(outfile)s > %(outfile)s.log
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
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
    --allow-empty \
    -b sqlite \
    --index=id \
    --index=motif \
    --index=id,motif \
    --table=%(tablename)s \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

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
##
############################################################
def _runNubiscanMotifSearch( infile, outfile, motif, mode="straight" ):
    '''run NUBISCAN to find matches to nuclear receptors.
    '''

    to_cluster = True
    job_options = "-l mem_free=8000M"

    track = re.sub( ".nubiscan.*", "", outfile )
    tmpfasta = P.getTempFilename( "." )

    if mode == "shuffled":
        nseq = writeSequencesForIntervals( track,
                                           tmpfasta,
                                           full = False,
                                           halfwidth = 100,
                                           shuffled = True )
    elif mode == "shifted":
        nseq = writeSequencesForIntervals( track,
                                           tmpfasta,
                                           full = False,
                                           halfwidth = 100,
                                           offset = -1000000)
    elif mode == "straight":
        nseq = writeSequencesForIntervals( track,
                                           tmpfasta,
                                           full = False,
                                           halfwidth = 100 )
    else:
        raise ValueError("invalid mode `%s`" % mode )


    arrangements = ",".join( PARAMS["motifs_nubiscan_arrangements"] )

    statement = '''python %(scriptsdir)s/run_nubiscan.py \
    --qvalue=%(motifs_nubiscan_qvalue)f \
    --iterations=%(motifs_nubiscan_iterations)i \
    --arrangements=%(arrangements)s \
    --motif=%(motif)s \
    --mask=dust \
    --add-sequence \
    --log=%(outfile)s.log \
    < %(tmpfasta)s > %(outfile)s'''

    P.run()

    os.unlink( tmpfasta )

############################################################
############################################################
############################################################
# todo: fix, causes a problem: exportSequences and exportMASTControlSequences,
@transform( exportMotifSequences,
            suffix(".fasta"),
            ".nubiscan.rxrvdr")
def runNubiscanVDR( infile, outfile ):
    _runNubiscanMotifSearch( infile, outfile, motif = "rxrvdr", mode = "straight" )

@transform( exportMotifSequences,
            suffix(".fasta"),
            ".nubiscan.rxrvdr.shuffled")
def runNubiscanVDRShuffled( infile, outfile ):
    _runNubiscanMotifSearch( infile, outfile, motif = "rxrvdr", mode = "shuffled" )

@transform( exportMotifSequences,
            suffix(".fasta"),
            ".nubiscan.rxrvdr.shifted")
def runNubiscanVDRShifted( infile, outfile ):
    _runNubiscanMotifSearch( infile, outfile, motif = "rxrvdr", mode = "shifted" )

############################################################
############################################################
############################################################
# todo: fix, causes a problem: exportSequences and exportMASTControlSequences,
@transform( exportMotifSequences,
            suffix(".fasta"),
            ".nubiscan.nr")
def runNubiscanNR( infile, outfile ):
    _runNubiscanMotifSearch( infile, outfile, motif = "nr", mode = "straight" )

@transform( exportMotifSequences,
            suffix(".fasta"),
            ".nubiscan.nr.shuffled")
def runNubiscanNRShuffled( infile, outfile ):
    _runNubiscanMotifSearch( infile, outfile, motif = "nr", mode = "shuffled" )

@transform( exportMotifSequences,
            suffix(".fasta"),
            ".nubiscan.nr.shifted")
def runNubiscanNRShifted( infile, outfile ):
    _runNubiscanMotifSearch( infile, outfile, motif = "nr", mode = "shifted" )

############################################################
############################################################
############################################################
# todo: fix, causes a problem: exportSequences and exportMASTControlSequences,
@transform( (runNubiscanVDR, runNubiscanNR),
            suffix( ".nubiscan.rxrvdr" ),
            "_nubiscan.import" )
def importNubiscan( infile, outfile ):

    tablename = outfile[:-len(".import")]
    track = tablename[:-len("_nubiscan")]

    tmpfile = P.getTempFile()
    first = True

    for line in fileinput.input( glob.glob("%s.nubiscan.*" % (track))):
        if fileinput.filename().endswith(".log"):
            fileinput.nextfile()
            continue

        if line.startswith("#"): continue

        if re.match("id\t", line):

            base, motif = os.path.splitext(fileinput.filename())

            if motif == ".shuffled":
                base, motif = os.path.splitext(base)
                motif += "_shuffled"
            elif motif == ".shifted":
                base, motif = os.path.splitext(base)
                motif += "_shifted"

            motif = motif[1:]

            if first:
                tmpfile.write("motif\t%s" % line )
                first = False
            continue

        line = re.sub( "run\S+_", "", line)
        # patch for fixing ids
        data = line.split("\t")
        data[0] = re.sub(" ", "", data[0])
        line = "\t".join(data)

        tmpfile.write("%s\t%s" % (motif, line))
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

    P.run( **dict( locals().items() + PARAMS.items() ) )

    os.unlink( tmpfilename )

############################################################
############################################################
############################################################
@transform( runMAST,
            suffix(".mast.gz"),
            "_mast.import" )
def importMAST( infile, outfile ):
    '''parse mast file and import into database.

    Parse several motif runs and add them to the same
    table.

    Add columns for the control data as well.
    '''

    tablename = outfile[:-len(".import")]
    tmpfile = tempfile.NamedTemporaryFile(delete=False)

    tmpfile.write( MAST.Match().header +\
                   "\tmotif\tcontig" \
                   "\tl_evalue\tl_pvalue\tl_nmatches\tl_length\tl_start\tl_end" \
                   "\tr_evalue\tr_pvalue\tr_nmatches\tr_length\tr_start\tr_end" \
                   "\tmin_evalue\tmin_pvalue\tmax_nmatches" + "\n" )

    lines = gzip.open(infile).readlines()
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
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
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
        cat %(dbfile)s %(controlfile)s | %(motifs_execglam2scan)s -2 -n %(motifs_glam2scan_results)i n %(motiffile)s - >> %(outfile)s
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
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
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
@files( [ ("%s.bed" % x, "%s.annotations" % x ) for x in TRACKS_ANNOTATE] )
def annotateIntervals( infile, outfile ):
    '''classify chipseq intervals according to their location 
    with respect to the gene set.
    '''
    to_cluster = True
    statement = """
    cat < %(infile)s 
        | python %(scriptsdir)s/bed2gff.py --as-gtf 
        | python %(scriptsdir)s/gtf2table.py 
                --counter=position 
                --counter=classifier-chipseq 
                --section=exons 
                --counter=length 
                --log=%(outfile)s.log 
                --filename-gff=%(annotation)s 
                --genome-file=%(genome)s
        > %(outfile)s"""

    P.run()

############################################################
############################################################
############################################################
def straightImport( infile, outfile, options = "" ):
    '''straight import from gtf2table results.

    The table name is given by outfile without the
    ".import" suffix.
    '''

    tablename = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=gene_id \
              %(options)s \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run()

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
@files( [ ("%s.bed" % x, "%s.tracks" % x, "ucsc_encode.bed" ) for x in TRACKS_ANNOTATE] )
def annotateTracks( infile, outfile, trackfile ):
    '''annotate intervals by track overlap
    '''

    statement = """
    python %(scriptsdir)s/bed2table.py 
                --counter=overlap 
                --log=%(outfile)s.log 
                --filename-bed=%(trackfile)s 
                --genome-file=%(genome)s
    < %(infile)s
    | perl -p -e "s/name/gene_id/" 
    > %(outfile)s"""

    P.run()

############################################################
############################################################
############################################################
@transform( annotateTracks, 
            suffix(".tracks"), 
            "_tracks.import" )
def importTracks( infile, outfile ):
    '''import interval annotations: genome architecture
    '''
    straightImport( infile, outfile )

############################################################
############################################################
############################################################
@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, "%s.tss" % x ) for x in TRACKS_ANNOTATE ] )
def annotateTSS( infile, outfile ):

    to_cluster = True
    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
        python %(scriptsdir)s/gtf2table.py \
                --counter=distance-tss \
                --log=%(outfile)s.log \
                --filename-gff=%(transcripts_tss)s \
                --genome-file=%(genome)s
        > %(outfile)s"""

    P.run()

############################################################
############################################################
############################################################
@files( [ ("%s.bed" % x, "%s.repeats" % x ) for x in TRACKS_ANNOTATE ] )
def annotateRepeats( infile, outfile ):
    '''count the overlap between intervals and repeats.'''

    to_cluster = True
    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
        python %(scriptsdir)s/gtf2table.py \
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
@transform( importExpressionTracks, suffix("_data.import"), ".ttest.expdiff" )
def buildExpressionDifferencesTTest( infile, outfile ):
    '''compute stats on expression differences between each set
    and the appropriate unstimulated set.

    A FDR Q-Value is computed as well for each pairwise comparison.
    '''

    nskipped = 0

    track = infile[:-len("_data.import")]

    outs = open( outfile, "w" )
    outs.write( "cluster_id\tmean1\tvar1\tmean2\tvar2\tdf\tpvalue\tqvalue\tdiff\tdiff95lower\tdiff95upper\tfold\tfold95lower\tfold95upper\n" )

    try:
        control, cluster_ids, treatments, controls = PExpression.getExpressionMeasurements( track )
    except P.PipelineError, msg:
        E.warn("skipped: %s" % msg )
        outs.close()
        return

    E.info( "track=%s, control=%s: computing ttest for %i probesets, %i samples, %i control" % (track, control,
                                                                                                len(cluster_ids),
                                                                                                len(treatments),
                                                                                                len(controls)))

    result, nskipped = Expression.WelchsTTest()( cluster_ids, treatments, controls )

    for s in result:
        # not that values are logged (log2), hence the fold-change
        # is exp( mean1 - mean2 )
        fold = math.pow(2,s.mMean1-s.mMean2)
        outs.write ("%s\t%f\t%f\t%f\t%f\t%5.2f\t%e\t%e\t%f\t%f\t%f\t%f\t%f\t%f\n" % \
                    (s.mProbeset, 
                     s.mMean1,
                     s.mSampleVariance1,
                     s.mMean2, 
                     s.mSampleVariance2,
                     s.mDegreesFreedom,
                     s.mPValue, 
                     s.mQValue, 
                     s.mDifference,
                     s.mDifferenceLower,
                     s.mDifferenceUpper,
                     fold,
                     math.pow(2,s.mDifferenceLower),
                     math.pow(2,s.mDifferenceUpper),
                     ) )
    outs.close()

############################################################
############################################################
############################################################
@transform( importExpressionTracks, suffix("_data.import"), ".sam.expdiff" )
def buildExpressionDifferencesSAM( infile, outfile ):
    '''compute stats on expression differences between each set
    and the appropriate unstimulated set using SAM.
    '''

    nskipped = 0

    track = infile[:-len("_data.import")]

    target_path = os.path.join( os.path.abspath( PARAMS["exportdir"] ), "SAM" )
    if not os.path.exists( target_path): 
        try:
            os.makedirs( target_path )
        except OSError: 
            pass

    outs = open( outfile, "w" )
    outs.write( "cluster_id\tmean1\tstd1\tmean2\tstd2\tpvalue\tqvalue\tdiff\tfold\tcalled\n" )

    try:
        control, cluster_ids, treatments, controls = PExpression.getExpressionMeasurements( track )
    except P.PipelineError, msg:
        E.warn("skipped: %s" % msg )
        outs.close()
        return

    E.info( "track=%s, control=%s: computing SAM for %i probesets, %i samples, %i control" % (track, control,
                                                                                              len(cluster_ids),
                                                                                              len(treatments),
                                                                                              len(controls)))

    genes, summary, cutoffs = Expression.SAM()( cluster_ids, treatments, controls,
                                                pattern = os.path.join(target_path, outfile + "%s"),
                                                fdr = float(PARAMS["expression_sam_fdr"]),
                                                ngenes = float(PARAMS["expression_sam_ngenes"]),
                                                ndelta = float(PARAMS["expression_sam_ndelta"]),
                                                npermutations = PARAMS["expression_sam_permutations"],
                                                method = PARAMS["expression_sam_method"],
                                                use_excel_sam = False )

    if summary == None:
        E.warn( "no cutoff when running sam for %s" % infile )

    if cutoffs:
        logs = open( outfile + ".cutoffs", "w" )
        logs.write("\t".join( cutoffs[0]._fields) + "\n" )
        for x in cutoffs:
            logs.write( "\t".join(map(str,x) ) + "\n" )
        logs.close()
    for s in genes:
        outs.write ("%s\t%f\t%f\t%f\t%f\t%e\t%f\t%f\t%f\t%i\n" % \
                        (s.probeset, 
                         s.mean1,
                         s.stddev1,
                         s.mean2,
                         s.stddev2,
                         s.pvalue, 
                         s.qvalue, 
                         s.difference,
                         s.fold,
                         s.called,
                         ) )
    outs.close()

                  
@transform( (buildExpressionDifferencesSAM, buildExpressionDifferencesTTest), 
           suffix(".expdiff"), 
            "_expdiff.import")
def importExpressionDifferences( infile, outfile ):
    '''import expression differences.
    '''

    track, method, control = PExpression.getExpressionMatch( infile )

    if track == control:
        outs = open( outfile, "w" )
        outs.close()
        return

    tablename = "%s_vs_%s_%s" % (track,control,method)
        
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --allow-empty \
              --index=cluster_id \
              --table=%(tablename)s \
         < %(infile)s > %(outfile)s
         '''

    P.run()

@transform( importExpressionDifferences,
            suffix("_expdiff.import"), 
            "_expdiff.export")
def exportDifferentiallyExpressedGenes( infile, outfile ):
    '''import expression differences.
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    track = infile[:-len("_expdiff.import")]
    if track.endswith( ".sam" ):
        field = "qvalue"
        cutoff = float(PARAMS['expression_sam_export_qvalue_cutoff'])
    elif track.endswith( ".ttest" ):
        field = "pvalue"
        cutoff = float(PARAMS['expression_ttest_export_pvalue_cutoff'])

    outs = open(outfile, "w" )

    track, method, control = PExpression.getExpressionMatch( track + ".expdiff" )

    if track == control:
        outs.close()
        return

    tablename = "%s_vs_%s_%s" % (track,control,method)

    statement = '''SELECT d.cluster_id, genesymbol, fold, pvalue, qvalue, ee.description 
                    FROM %(tablename)s as d, expression as ee on ee.cluster_id = d.cluster_id where %(field)s < %(cutoff)f''' % locals()
    cc = dbhandle.cursor()
    cc.execute( statement )

    outs.write( "cluster_id\tgenesymbol\tfold\tpvalue\tqvalue\tdescription\n" )
    for x in cc:
        outs.write( "\t".join( map(str, x) ) + "\n" )
        
    outs.close()

############################################################
############################################################
############################################################
@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, "%s.probetss" % x ) for x in TRACKS_MASTER ] )
def annotateProbesetTSS( infile, outfile ):

    to_cluster = True
    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python %(scriptsdir)s/gtf2table.py \
		--counter=distance-tss \
		--log=%(outfile)s.log \
		--filename-gff=%(probeset)s \
		--genome-file=%(genome)s
	> %(outfile)s"""

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
## count coverage within intervals for each track agains all
## the Unstim tracks
############################################################
@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, "%s.readcounts" % x ) for x in TRACKS_NOT_UNSTIM ] )
def buildIntervalCounts( infile, outfile ):
    '''count read density in bed files.
    '''
    track = outfile[:-len(".readcounts")]

    cellline, condition, replicate = splitTrack( track )

    samfiles_fg, samfiles_bg = [], []

    # save original, because condition is later assigned to
    one = cellline and condition
    two = condition

    if REPLICATES:
        if one:
            for replicate in REPLICATES:
                samfiles_fg.append( FILEPATTERN % locals() + BAM_SUFFIX )
            condition = "Unstim"
            for replicate in REPLICATES:
                samfiles_bg.append( FILEPATTERN % locals() + BAM_SUFFIX )
        elif two:
            for cellline in CELLLINES:
                for replicate in REPLICATES:
                    samfiles_fg.append( FILEPATTERN % locals() + BAM_SUFFIX )
            condition = "Unstim"
            for cellline in CELLLINES:
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
	python %(scriptsdir)s/gtf2table.py \
		--counter=read-coverage \
		--log=%(outfile)s.log \
		--bam-file=%(samfiles_fg)s \
	> %(tmpfile1)s"""

    P.run( **dict( locals().items() + PARAMS.items() ) )

    statement = """
    cat < %(infile)s |\
        python %(scriptsdir)s/bed2gff.py --as-gtf |\
	python %(scriptsdir)s/gtf2table.py \
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

    P.run( **dict( locals().items() + PARAMS.items() ) )
    os.unlink( tmpfile1 )
    os.unlink( tmpfile2 )


@files_re( buildIntervalCounts, "(.*).readcounts", r"\1_readcounts.import" )
def importIntervalCounts( infile, outfile ):
    '''import interval counts.'''
    tablename = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=gene_id \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

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
## export expression calls
############################################################
@merge( importExpressionDifferences, "expression_differences_%s.bed" % PARAMS["version"] )
def viewExpressionDifferences( infiles, outfile ):
    '''export a bed file with genes called as differentially expressed.'''
    dbhandle = sqlite3.connect( PARAMS["database"] )

    outs = open( outfile, "w" )

    for infile in infiles:
        track, method, control = PExpression.getExpressionMatch( infile[:-len(".import")] )
        if track == control: continue
        table = "%(track)s_vs_%(control)s_%(method)s" % locals()
        
        statement = """SELECT DISTINCT g.contig, g.start, g.end, exp.probeset
                          FROM %(table)s AS exp, probeset_gtf AS g
                          WHERE g.gene_id = exp.probeset AND \
                          exp.pvalue < 0.10
                   """ % locals()
        cc = dbhandle.cursor()
        cc.execute(statement)

        outs.write( '''track name="%(track)s_%(method)s" description="%(track)s - %(method)s" visibility=2\n''' % locals() )
        for contig, start, end, probeset in cc:
            if contig == "chrMT": contig = "chrM"
            bed = Bed.Bed()
            bed.contig, bed.start, bed.end = contig, start, end
            bed.fields = [ probeset ]
            outs.write( "%s\n" % str(bed) )

    outs.close()

    dest = os.path.join( PARAMS["ucsc_dir"], outfile ) 
    shutil.copyfile( outfile, dest )
    print "paste the following into the UCSC browser:"
    print '''http://wwwfgu.anat.ox.ac.uk/~andreas/ucsc_tracks/%(outfile)s''' % locals()

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
## get GO assignments
############################################################
@files( [ (None, PARAMS["filename_go"] ), ] )
def createGO( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    
    statement = '''
        python %(scriptsdir)s/GO.py \
                     --filename-dump=%(outfile)s \
                     --host=ensembldb.ensembl.org \
                     --user=anonymous \
                     --database=%(go_database)s \
                     --port=5306 > %(outfile)s.log
        '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
## get GO Slim assignments
############################################################
@files_re( createGO, "(.*).go", r"\1.goslim") 
def createGOSlim( infile, outfile ):
    '''get GO assignments from ENSEMBL'''
    
    statement = '''
        wget http://www.geneontology.org/GO_slims/goslim_goa.obo
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

    statement = '''
        wget http://www.geneontology.org/ontology/gene_ontology.obo
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

    
    statement = '''
        /net/cpp-group/server/lib/perl5/5.10.0/bin/map2slim -outmap go2goslim.map goslim_goa.obo gene_ontology.obo
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

    statement = '''
        python %(scriptsdir)s/GO.py \
                --go2goslim \
                --filename-ontology=gene_ontology.obo \
                --slims=go2goslim.map \
                --log=%(outfile)s.log \
        < %(infile)s > %(outfile)s
        '''

    P.run( **dict( locals().items() + PARAMS.items() ) )
    
############################################################
############################################################
############################################################
##
############################################################
@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, "%s-%s.annodist" % (x,y) ) \
              for x,y in itertools.product( TRACKS_ANNOTATOR, (TARGET_TSS, "intergenic", "intronic", "genic")) ] )
def makeAnnotatorDistance( infile, outfile):
    '''check statistical association between intervals and transcription
    start sites.'''

    track,workspace = outfile[:-len(".annodist")].split("-")

    if workspace == "intronic":
        builder = "gtf-intronic"
        workspace = TARGET_GENESET
    elif workspace == "genic":
        builder = "gtf-genic"
        workspace = TARGET_GENESET
    elif workspace == "intergenic":
        builder = "gtf-intergenic"
        workspace = TARGET_GENESET
    else:
        builder = "gtf-intergenic"

    PAnnotator.makeAnnotatorDistance( infile, outfile, builder, workspace, 
                                      workspace_label = "direction" )


@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, 
           "%s-expression.annodist" % (x) ) \
              for x in TRACKS_MASTER ] )
def makeAnnotatorDistanceExpression( infile, outfile):
    '''check statistical association between intervals and transcription
    start sites.

    Look at different expression sets.
    '''

    track,workspace = outfile[:-len(".annodist")].split("-")
    
    builder = "gtf-intergenic"
    workspace = TARGET_GENESET

    tmpannotations= PAnnotator.buildAnnotatorDistanceAnnotations( "expression" )
                                                    
    PAnnotator.makeAnnotatorDistance( infile, 
                                      outfile, 
                                      builder,
                                      workspace,
                                      workspace_label = "annotation",
                                      annotations = tmpannotations,
                                      )
    
    os.unlink( tmpannotations )

@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, 
           "%s-expressionnopromotor.annodist" % (x) ) \
              for x in TRACKS_MASTER ] )
def makeAnnotatorDistanceWithoutPromotorExpression( infile, outfile):
    '''check statistical association between intervals and transcription
    start sites.

    Look at different expression sets.
    '''

    track,workspace = outfile[:-len(".annodist")].split("-")
    
    builder = "gtf-intergenic"

    tmpfilename = P.getTempFilename( ".")

    # build workspace: add 5kb on either side of
    # a transcription start site.
    statement = '''
    python %(scriptsdir)s/gff2gff.py 
            --genome=genome 
            --extend 
            --add-up-flank=%(promotor_size)i 
            --add-down-flank=%(promotor_size)i 
    < %(tss)s > %(tmpfilename)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

    workspace = tmpfilename

    tmpannotations= PAnnotator.buildAnnotatorDistanceAnnotations( "expression" )
                                                    
    PAnnotator.makeAnnotatorDistance( infile, 
                                      outfile, 
                                      builder,
                                      workspace,
                                      workspace_label = "annotation",
                                      annotations = tmpannotations,
                                      )
    
    os.unlink( tmpannotations )


############################################################
############################################################
############################################################
##
############################################################
@follows( buildGeneRegions, buildAnnotatorGC )
@files( [ ("%s.bed" % x, "%s.annotator_architecture" % x ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorArchitecture( infile, outfile ):
    '''check statistical overlap between intervals and and other genomic features
    defined in the file PARAMS["annotations"].

    Annotator is run with the following parameters:
    1. Segments: the interval track
    2. Annotations:
       1. genomic architecture (PARAMS["annotation"])
       2. promotors (PARAMS["promotors"])
    3. Workspace: the full genome
    '''
    PAnnotator.makeAnnotatorArchitecture( infile, outfile )

@follows( buildGeneRegions, buildAnnotatorGC )
@files( [ (("%s.bed" % x, "gwas.bed" ), "%s.annotator_gwas" % x ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorGWAS( infiles, outfile ):
    '''check statistical overlap between intervals and and other genomic features
    defined in the file PARAMS["gwas"].
    '''
    PAnnotator.makeAnnotatorRegionsOfInterest( infiles, outfile )

@follows( buildGeneRegions, buildAnnotatorGC )
@files( [ (("%s.bed" % x, "gwas.bed" ), "%s.annotator_gwas_with_motif" % x ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorGWASWithMotif( infiles, outfile ):
    '''check statistical overlap between intervals and and other genomic features
    defined in the file PARAMS["gwas"].
    '''
    PAnnotator.makeAnnotatorRegionsOfInterest( infiles, 
                                               outfile,
                                               with_motif = PARAMS['annotator_master_motif'] )

@follows( buildGeneRegions, buildAnnotatorGC )
@files( [ (("%s.bed" % x, "gwas.bed" ), "%s.annotator_gwas_without_motif" % x ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorGWASWithoutMotif( infiles, outfile ):
    '''check statistical overlap between intervals and and other genomic features
    defined in the file PARAMS["gwas"].
    '''
    PAnnotator.makeAnnotatorRegionsOfInterest( infiles, 
                                               outfile,
                                               without_motif = PARAMS['annotator_master_motif'] )

@follows( buildGeneRegions, buildAnnotatorGC )
@files( [ (("%s.bed" % x, "selection.bed" ), "%s.annotator_selection" % x ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorSelection( infiles, outfile ):
    '''check statistical overlap between intervals and and other genomic features
    defined in the file PARAMS["selection"].
    '''
    PAnnotator.makeAnnotatorRegionsOfInterest( infiles, outfile )
    
############################################################
############################################################
############################################################
##
############################################################
@follows( exportUCSCEncodeTracks, buildAnnotatorGC )
@files( [ (("%s.bed" % x, "ucsc_encode.bed"), "%s.annotator_tracks" % x ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorTracks( infiles, outfile ):
    PAnnotator.makeAnnotatorTracks( infiles, outfile )
    
############################################################
############################################################
############################################################
##
############################################################
@follows( exportRegionsOfInterest, buildAnnotatorGC )
@files( [ (("%s.bed" % x, "roi.bed"), "%s.annotator_roi" % x ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorRegionsOfInterest( infiles, outfile ):
    PAnnotator.makeAnnotatorRegionsOfInterest( infiles, outfile )

############################################################
############################################################
############################################################
##
############################################################
@follows( exportUCSCEncodeTracks, buildAnnotatorGC )
@files( [ (("%s.bed" % x, "ucsc_encode.bed"), "%s.annotator_promotors" % x ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorPromotors( infiles, outfile ):
    '''check statistical overlap between intervals and selected ucsc tracks

    Annotator is run with the following parameters:
    1. Segments: the interval track
    2. Annotations:
       1. ucsc encode features
    3. Workspace:
       1. the mappable part of the genome
       2. GC controlled
    '''

    infile, infile_annotations = infiles
    tmpdir = P.getTempDir( dir = os.getcwd() )

    annotations = PAnnotator.buildAnnotatorAnnotations( tmpdir, outfile, bedfiles=(infile_annotations,) )
    
    workspaces, synonyms = PAnnotator.buildAnnotatorWorkSpace( tmpdir, outfile,
                                                    workspaces = ("mappable", "promotors"),
                                                    gc_control = True )

    segments = PAnnotator.buildAnnotatorSegments( tmpdir, infile, outfile )

    PAnnotator.runAnnotator( tmpdir, outfile, annotations, segments, workspaces, synonyms )

    shutil.rmtree( tmpdir )

############################################################
############################################################
############################################################
##
############################################################
@merge( makeAnnotatorTracks,
        "annotator_tracks.import" )
def importAnnotatorTracks( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_tracks", "genome", "all", "all", "annotator" )

############################################################
############################################################
############################################################
##
############################################################
@merge( makeAnnotatorArchitecture,
        "annotator_architecture.import" )
def importAnnotatorArchitecture( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_architecture", "genome", "all", "all", "annotator" )

############################################################
############################################################
############################################################
##
############################################################
@merge( makeAnnotatorRegionsOfInterest,
        "annotator_roi.import" )
def importAnnotatorRegionsOfInterest( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_roi", "genome", "all", "all", "annotator" )

@merge( makeAnnotatorGWAS,
        "annotator_gwas.import" )
def importAnnotatorGWAS( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_gwas", "genome", "all", "all", "annotator" )

@merge( makeAnnotatorSelection,
        "annotator_selection.import" )
def importAnnotatorSelection( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_selection", "genome", "all", "all", "annotator" )

@follows( importAnnotatorRegionsOfInterest,
          importAnnotatorGWAS,
          importAnnotatorSelection )
def annotator_regions(): pass
    
    
############################################################
############################################################
############################################################
##
############################################################
@merge( makeAnnotatorPromotors,
        "annotator_promotors.import" )
def importAnnotatorPromotors( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_promotors", "genome", "all", "all", "annotator" )

############################################################
############################################################
############################################################
## GO annotator analysis
## scan within promotors only
## scan within gene territories only
############################################################
@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, "%s.annotator_promotors_go" % x,
           PARAMS["filename_go"], "promotors" ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorPromotorsGO( infile, outfile, gofile, workspace ):
    PAnnotator.makeAnnotatorGO( infile, outfile, gofile, workspace )

@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, "%s.annotator_promotors_goslim" % x,
           PARAMS["filename_goslim"], "promotors" ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorPromotorsGOSlim( infile, outfile, gofile, workspace ):
    PAnnotator.makeAnnotatorGO( infile, outfile, gofile, workspace )

@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, "%s.annotator_territories_go" % x,
           PARAMS["filename_go"], "gene-territories" ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorTerritoriesGO( infile, outfile, gofile, workspace ):
    PAnnotator.makeAnnotatorGO( infile, outfile, gofile,workspace )

@follows( buildGeneRegions )
@files( [ ("%s.bed" % x, "%s.annotator_territories_goslim" % x,
           PARAMS["filename_goslim"], "gene-territories" ) for x in TRACKS_ANNOTATOR ] )
def makeAnnotatorTerritoriesGOSlim( infile, outfile, gofile, workspace ):
    PAnnotator.makeAnnotatorGO( infile, outfile, gofile, workspace )

############################################################
############################################################
############################################################
##
############################################################
@merge( makeAnnotatorPromotorsGO, "annotator_promotors_go.import" )
def importAnnotatorPromotorsGO( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile,
                                      "annotator_promotors_go", "promotors",
                                      "all", "all",
                                      fdr_method = "annotator-estimate" )

@merge( makeAnnotatorPromotorsGOSlim, "annotator_promotors_goslim.import" )
def importAnnotatorPromotorsGOSlim( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile,
                                      "annotator_promotors_goslim", "promotors",
                                      "all", "all",
                                      fdr_method = "annotator" )

@merge( makeAnnotatorTerritoriesGO, "annotator_territories_go.import" )
def importAnnotatorTerritoriesGO( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile,
                                      "annotator_territories_go", "territories",
                                      "all", "all",
                                      fdr_method = "annotator-estimate" )

@merge( makeAnnotatorTerritoriesGOSlim, "annotator_territories_goslim.import" )
def importAnnotatorTerritoriesGOSlim( infiles, outfile ):
    PAnnotator.genericImportAnnotator( infiles, outfile,
                                      "annotator_territories_goslim", "territories",
                                      "all", "all",
                                      fdr_method = "annotator" )


############################################################
############################################################
############################################################
## optional annotator analysis:
##
## GO annotator analysis for regions of interest
## Note that in contrast to other annotator analyses, the 
## segments used are those for regions of interest on not
## the ChIPSeq intervals.
############################################################
if PARAMS["annotator_master_roi"]:

    @follows( buildGeneRegions )
    @files( [ ("roi.bed", 
               "%s.annotator_roi_goslim" % x ,
               x,
               PARAMS["filename_goslim"], 
               "gene-territories" ) for x in TRACKS_ROI ] )
    def makeAnnotatorROIGOSlim( infile, outfile, roi_class, gofile, workspace ):
        '''do GO annotator analysis on selection intervals.
        '''
        PAnnotator.makeAnnotatorROIGO( roi_class, outfile, gofile, workspace, overlap=None )

    @follows( buildGeneRegions )
    @files( [ ("roi.bed", 
               "%s.annotator_roi_go" % x ,
               x,
               PARAMS["filename_go"], 
               "gene-territories" ) for x in TRACKS_ROI ] )
    def makeAnnotatorROIGO( infile, outfile, roi_class, gofile, workspace ):
        '''do GO annotator analysis on selection intervals.
        '''
        PAnnotator.makeAnnotatorROIGO( roi_class, outfile, gofile, workspace, overlap=None )

    @follows( buildGeneRegions )
    @files( [ ("roi.bed", 
               "%s.annotator_roi_overlap_goslim" % x ,
               x,
               PARAMS["filename_goslim"], 
               "gene-territories" ) for x in TRACKS_ROI ] )
    def makeAnnotatorROIOverlapGOSlim( infile, outfile, roi_class, gofile, workspace ):
        '''do GO annotator analysis on selection intervals.
        '''
        PAnnotator.makeAnnotatorROIGO( roi_class, outfile, gofile, workspace, overlap=PARAMS["annotator_master_roi"] )

    @follows( buildGeneRegions )
    @files( [ ("roi.bed", 
               "%s.annotator_roi_overlap_go" % x ,
               x,
               PARAMS["filename_go"], 
               "gene-territories" ) for x in TRACKS_ROI ] )
    def makeAnnotatorROIOverlapGO( infile, outfile, roi_class, gofile, workspace ):
        '''do GO annotator analysis on selection intervals.
        '''
        PAnnotator.makeAnnotatorROIGO( roi_class, outfile, gofile, workspace, overlap=PARAMS["annotator_master_roi"] )

    @merge( makeAnnotatorROIGO, "annotator_roi_go.import" )
    def importAnnotatorROIGO( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile,
                                          "annotator_roi_go", "territories",
                                          "all", "all",
                                          fdr_method = "annotator-estimate" )

    @merge( makeAnnotatorROIGOSlim, "annotator_roi_goslim.import" )
    def importAnnotatorROIGOSlim( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile,
                                          "annotator_roi_goslim", "territories",
                                          "all", "all",
                                          fdr_method = "annotator" )


    @merge( makeAnnotatorROIOverlapGO, "annotator_roi_overlap_go.import" )
    def importAnnotatorROIOverlapGO( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile,
                                          "annotator_roi_overlap_go", "territories",
                                          "all", "all",
                                          fdr_method = "annotator-estimate" )


    @merge( makeAnnotatorROIOverlapGOSlim, "annotator_roi_overlap_goslim.import" )
    def importAnnotatorROIOverlapGOSlim( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile,
                                          "annotator_roi_overlap_goslim", "territories",
                                          "all", "all",
                                          fdr_method = "annotator" )

    @follows( makeAnnotatorROIGO,            importAnnotatorROIGO,          
              makeAnnotatorROIGOSlim,        importAnnotatorROIGOSlim,      
              makeAnnotatorROIOverlapGO,     importAnnotatorROIOverlapGO,   
              makeAnnotatorROIOverlapGOSlim, importAnnotatorROIOverlapGOSlim, )
    def annotator_roi():
        pass
else:
    @follows(buildGeneSet)
    def annotator_roi():
        pass

############################################################
############################################################
############################################################
## optional annotator section: 
## run annotator analysis on interval sets that overlap/do not
##          overlap motifs.
############################################################
if PARAMS["annotator_master_motif"]:

    @follows( buildGeneRegions, buildAnnotatorGC )
    @files( [ ("%s.bed" % x, "%s.annotator_architecture_with_motif" % x ) for x in TRACKS_ANNOTATOR ] )
    def makeAnnotatorArchitectureWithMotif( infile, outfile ):
        '''check statistical overlap between intervals and and other genomic features
        defined in the file PARAMS["annotations"].

        Annotator is run with the following parameters:
        1. Segments: the interval track, only those with motifs
        2. Annotations:
           1. genomic architecture (PARAMS["annotation"])
           2. promotors (PARAMS["promotors"])
        3. Workspace: the full genome
        '''
        PAnnotator.makeAnnotatorArchitecture( infile, outfile, with_motif = PARAMS['annotator_master_motif'] )

    @follows( buildGeneRegions, buildAnnotatorGC )
    @files( [ ("%s.bed" % x, "%s.annotator_architecture_without_motif" % x ) for x in TRACKS_ANNOTATOR ] )
    def makeAnnotatorArchitectureWithoutMotif( infile, outfile ):
        '''check statistical overlap between intervals and and other genomic features
        defined in the file PARAMS["annotations"].

        Annotator is run with the following parameters:
        1. Segments: the interval track, only those without motifs
        2. Annotations:
           1. genomic architecture (PARAMS["annotation"])
           2. promotors (PARAMS["promotors"])
        3. Workspace: the full genome
        '''
        PAnnotator.makeAnnotatorArchitecture( infile, outfile, without_motif = PARAMS['annotator_master_motif'] )


    @follows( exportUCSCEncodeTracks, buildAnnotatorGC )
    @files( [ (("%s.bed" % x, "ucsc_encode.bed"), "%s.annotator_tracks_with_motif" % x ) for x in TRACKS_ANNOTATOR ] )
    def makeAnnotatorTracksWithMotif( infiles, outfile ):
        '''check statistical overlap between intervals and selected ucsc tracks

        Annotator is run with the following parameters:
        1. Segments: the interval track, but stratified according to
                     whether they have or do not have a motif.
        2. Annotations:
           1. ucsc encode features
           2. disease intervals (regions of interest)
        3. Workspace: the full genome
        '''

        PAnnotator.makeAnnotatorTracks( infiles, outfile, with_motif = PARAMS['annotator_master_motif'] )

    ############################################################
    ############################################################
    ############################################################
    ##
    ############################################################
    @follows( exportUCSCEncodeTracks, buildAnnotatorGC )
    @files( [ (("%s.bed" % x, "ucsc_encode.bed"), "%s.annotator_tracks_without_motif" % x ) for x in TRACKS_ANNOTATOR ] )
    def makeAnnotatorTracksWithoutMotif( infiles, outfile ):
        '''check statistical overlap between intervals and selected ucsc tracks

        Annotator is run with the following parameters:
        1. Segments: the interval track, but stratified according to
                     whether they have or do not have a motif.
        2. Annotations:
           1. ucsc encode features
           2. disease intervals (regions of interest)
        3. Workspace: the full genome
        '''

        PAnnotator.makeAnnotatorTracks( infiles, outfile, without_motif = PARAMS['annotator_master_motif'] )

    @follows( exportRegionsOfInterest, buildAnnotatorGC )
    @files( [ (("%s.bed" % x, "roi.bed"), "%s.annotator_roi_with_motif" % x ) for x in TRACKS_ANNOTATOR ] )
    def makeAnnotatorRegionsOfInterestWithMotif( infiles, outfile ):
        PAnnotator.makeAnnotatorRegionsOfInterest( infiles, outfile, with_motif = PARAMS['annotator_master_motif'] )

    @follows( exportRegionsOfInterest, buildAnnotatorGC )
    @files( [ (("%s.bed" % x, "roi.bed"), "%s.annotator_roi_without_motif" % x ) for x in TRACKS_ANNOTATOR ] )
    def makeAnnotatorRegionsOfInterestWithoutMotif( infiles, outfile ):
        PAnnotator.makeAnnotatorRegionsOfInterest( infiles, outfile, without_motif = PARAMS['annotator_master_motif'] )

    @merge( makeAnnotatorArchitectureWithMotif,
            "annotator_architecture_with_motif.import" )
    def importAnnotatorArchitectureWithMotif( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_architecture_with_motif", "genome", "all", "all", "annotator" )

    @merge( makeAnnotatorArchitectureWithoutMotif,
            "annotator_architecture_without_motif.import" )
    def importAnnotatorArchitectureWithoutMotif( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_architecture_without_motif", "genome", "all", "all", "annotator" )

    @merge( makeAnnotatorTracksWithMotif,
            "annotator_tracks_with_motif.import" )
    def importAnnotatorTracksWithMotif( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_tracks_with_motif", "genome", "all", "all", "annotator" )

    @merge( makeAnnotatorTracksWithoutMotif,
            "annotator_tracks_without_motif.import" )
    def importAnnotatorTracksWithoutMotif( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_tracks_without_motif", "genome", "all", "all", "annotator" )

    @merge( makeAnnotatorRegionsOfInterestWithMotif,
            "annotator_roi_with_motif.import" )
    def importAnnotatorRegionsOfInterestWithMotif( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_roi_with_motif", "genome", "all", "all", "annotator" )

    @merge( makeAnnotatorRegionsOfInterestWithoutMotif,
            "annotator_roi_without_motif.import" )
    def importAnnotatorRegionsOfInterestWithoutMotif( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_roi_without_motif", "genome", "all", "all", "annotator" )

    @follows( makeAnnotatorArchitectureWithMotif, importAnnotatorArchitectureWithMotif,
              makeAnnotatorArchitectureWithoutMotif, importAnnotatorArchitectureWithoutMotif,
              makeAnnotatorTracksWithMotif, importAnnotatorTracksWithMotif,
              makeAnnotatorTracksWithoutMotif, importAnnotatorTracksWithoutMotif,
              makeAnnotatorRegionsOfInterestWithMotif, importAnnotatorRegionsOfInterestWithMotif,
              makeAnnotatorRegionsOfInterestWithoutMotif, importAnnotatorRegionsOfInterestWithoutMotif )
    def annotator_motifs():
        pass

else:
    @follows(buildGeneSet)
    def annotator_motifs():
        pass

############################################################
############################################################
############################################################
## optional annotator section: 
## run annotator analysis on top/bottom of interval sets 
############################################################

if PARAMS["annotator_proportion"]:
    @follows( exportRegionsOfInterest, buildAnnotatorGC )
    @files( [ (("%s.bed" % x, "roi.bed"), "%s.annotator_roi_top" % x ) for x in TRACKS_ANNOTATOR ] )
    def makeAnnotatorRegionsOfInterestTop( infiles, outfile ):
        PAnnotator.makeAnnotatorRegionsOfInterest( infiles, outfile, proportion = PARAMS['annotator_proportion'] )

    @follows( exportRegionsOfInterest, buildAnnotatorGC )
    @files( [ (("%s.bed" % x, "roi.bed"), "%s.annotator_roi_bottom" % x ) for x in TRACKS_ANNOTATOR ] )
    def makeAnnotatorRegionsOfInterestBottom( infiles, outfile ):
        PAnnotator.makeAnnotatorRegionsOfInterest( infiles, outfile, proportion = -PARAMS['annotator_proportion'] )

    @merge( makeAnnotatorRegionsOfInterestTop,
            "annotator_roi_top.import" )
    def importAnnotatorRegionsOfInterestTop( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_roi_top", "genome", "all", "all", "annotator" )

    @merge( makeAnnotatorRegionsOfInterestBottom,
            "annotator_roi_bottom.import" )
    def importAnnotatorRegionsOfInterestBottom( infiles, outfile ):
        PAnnotator.genericImportAnnotator( infiles, outfile, "annotator_roi_bottom", "genome", "all", "all", "annotator" )

    @follows( makeAnnotatorRegionsOfInterestTop,
              importAnnotatorRegionsOfInterestTop,
              makeAnnotatorRegionsOfInterestBottom,
              importAnnotatorRegionsOfInterestBottom )
    def annotator_proportion():
        pass
else:
    @follows(buildGeneSet)
    def annotator_proportion():
        pass

############################################################
############################################################
############################################################
##
############################################################
@follows( buildGeneRegions )
@files( [ ( ( "%s.bed" % x, TARGET_TRANSCRIPTS_TSS), "%s.distance" % x ) for x in TRACKS_ANNOTATE ] )
def annotateTSSIntervalDistance( infiles, outfile ):
    '''annotate the associations between a tss and of a transcript.
    Only the closest distance is recorded in either direction.'''

    infile_bed, infile_tss = infiles

    to_cluster = True

    statement = """
	python %(scriptsdir)s/gtf2table.py \
		--counter=distance \
		--log=%(outfile)s.log \
                --reporter=transcripts \
		--filename-gff=<( python %(scriptsdir)s/bed2gff.py --as-gtf < %(infile_bed)s) \
                --filename-format=gtf \
		--genome-file=%(genome)s \
	< %(infile_tss)s > %(outfile)s"""

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@files_re( annotateTSSIntervalDistance,
           "(.*).distance",
           r"\1_distance.import")
def importTSSIntervalDistance( infile, outfile):
    tablename = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=transcript_id \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
##
############################################################
@follows( buildGeneRegions )
@files( [ ( ( "%s.bed" % x, TARGET_TRANSCRIPTS_TSS), "%s.assoc" % x ) for x in TRACKS_ANNOTATE ] )
def annotateTSSIntervalAssociations( infiles, outfile ):
    '''annotate the associations between tss and of a transcript.
    All distances within *proximal_distance* are recorded.'''

    infile_bed, infile_tss = infiles

    to_cluster = True

    statement = """
	python %(scriptsdir)s/gtf2table.py \
		--counter=neighbours \
		--log=%(outfile)s.log \
                --reporter=transcripts \
                --proximal-distance=%(proximal_distance)i \
		--filename-gff=<( python %(scriptsdir)s/bed2gff.py --as-gtf < %(infile_bed)s) \
		--genome-file=%(genome)s \
	< %(infile_tss)s |\
        python %(toolsdir)s/table2graph.py \
                --headers=transcript_id,id \
    > %(outfile)s"""

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@files_re( annotateTSSIntervalAssociations,
           "(.*).assoc",
           r"\1_assoc.import")
def importTSSIntervalAssociations( infile, outfile):
    tablename = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=id \
              --index=transcript_id \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

############################################################
############################################################
############################################################
@files( PARAMS["filename_t1d"], "t1d.bed" )
def buildDiseaseIntervals( infile, outfile):

    reader = csv.reader( open(infile,"rU") )
    fasta = IndexedFasta.IndexedFasta( FILENAME_GENOME )

    outs = open( outfile, "w" )

    for row in reader:
        if row[0] == "Locus": 
            headers = row
            continue
        loc, chrom, start, end = row[:4]
        if loc == "": continue
        start, end = int(start)-1, int(end)
        
        contig = fasta.getToken( re.sub("[pq]\d[.0-9]*", "", chrom) )

        outs.write("%s\t%i\t%i\t%s\n" % (contig, start, end, loc) )

    outs.close()

###################################################################
###################################################################
###################################################################
## import data from Wang et. al
###################################################################
@files( ((PARAMS["filename_expression_wang"], "expression_wang.import"), ) )
def importExpressionWang( infile, outfile ):
    '''import gene expression data from Wang et al. (2005) 
    
    :pmid:`16002434`
    '''

    headers = (
        ("Probe Set ID", "probeset"),
        ("Accession", "accession"),
        ("Chromosomal Location", "chromosome"),
        ("Gene Title", "genetitle"),
        ("Gene Symbol", "genesymbol"),
        ("Entrez Gene", "entrezgene"),
        ("Score", "score"),
        ("Q-value", "qvalue"),
        ("Gene Symbol (May2005)", "genesymbol1"),
        ("Entrez Gene (May 2005)", "entrezgene1"),
        ("Unigene (Release # 184)", "unigene",),
        ("DR3 cons.","d3_cons" ),
        ("DR3 1 mism.","d3_mismatch",),
        ("ER6 cons.", "er6_cons"), )

    old_headers = set( [x[0] for x in headers] )
    new_headers = [x[1] for x in headers]
    take = []
    outs = P.getTempFile()
    first = True

    reader = csv.reader( open(infile,"rU") )
    for row in reader:
        if first:
            first = False
            for x, old_header in enumerate(row ):
                if old_header in old_headers: take.append( x )
            outs.write("\t".join(new_headers)+ "\n")
        else:
            new_row = []
            for x in take:
                if row[x].strip() != "---":
                    new_row.append(row[x].strip())
                else:
                    new_row.append("")
            outs.write("\t".join( new_row )+ "\n")
            
    outs.close()
    tmpname = outs.name
    tablename = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=probeset \
              --index=genesymbol \
              --table=%(tablename)s \
    < %(tmpname)s > %(outfile)s
    '''
    P.run()

    os.unlink( tmpname )

###################################################################
###################################################################
###################################################################
## import data from Wang et. al
###################################################################
@files( ((PARAMS["filename_expression_kovalenko"], "expression_kovalenko.import"), ) )
def importExpressionKovalenko( infile, outfile ):
    '''import gene expression data from Kovalenko et al. (2010) 

    Converts negative fold change into ratios and % qvalues
    into qvalues.
    '''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    map_name2id = dict( cc.execute( "SELECT UPPER(gene_name), gene_id FROM gene_info" ).fetchall() )

    headers = (
        ("affy_id", "probeset"),
        ("gene_name", "gene_id"),
        ("fold_6", "fold_6"),
        ("qvalue_6", "qvalue_6"),
        ("fold_24", "fold_24"),
        ("qvalue_24", "qvalue_24"),
        ("fold_48", "fold_48"),
        ("qvalue_48", "qvalue_48"),
        )

    old_headers = set( [x[0] for x in headers] )
    new_headers = [x[1] for x in headers]
    take = []
    outs = P.getTempFile()
    first = True

    counts = E.Counter()

    reader = csv.reader( open(infile,"rU") )
    notfound = set()
    for row in reader:
        if first:
            first = False
            for x, old_header in enumerate(row ):
                if old_header in old_headers: take.append( x )
            outs.write("\t".join(new_headers)+ "\n")
        else:
            skip = False
            new_row = []
            for x in take:
                if headers[x][0].startswith("fold"):
                    foldval = float( row[x] )
                    if foldval < 0: 
                        row[x] = str( -1.0 / foldval )
                elif headers[x][0].startswith("qvalue"):
                    qvalue = float( row[x] )
                    row[x] = str(qvalue / 100.0)
                elif headers[x][0].startswith("gene_name"):
                    parts = [ c.upper().strip() for c in row[x].split("//") ]

                    for part in parts:
                        if part in map_name2id:
                            row[x] = map_name2id[part].upper()
                            counts.found += 1
                            break
                    else:
                        counts.notfound += 1
                        notfound.add( row[x] )
                        skip = True
                        break


                if row[x].strip() != "---":
                    new_row.append(row[x].strip())
                else:
                    new_row.append("")
                    
            if skip: continue
            outs.write("\t".join( new_row )+ "\n")

    outs.close()
    tmpname = outs.name
    tablename = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              --index=probeset \
              --index=gene_id \
              --table=%(tablename)s \
    < %(tmpname)s > %(outfile)s
    '''
    P.run()

    os.unlink( tmpname )

    E.info( "%s" % str(counts))
    if notfound:
        E.warn("not found: %s" % str(notfound) )

############################################################
############################################################
############################################################
## run GO analysis on differentially expressed genes
############################################################
@transform( importExpressionDifferences,
            suffix( ".import"),
            ".goslim" )
def runExpressionGOSlim( infile, outfile ):
    '''run GOSlim analysis on expression tracks.'''
    PExpression.runGO( infile, outfile, PARAMS["filename_goslim"] )

@transform( importExpressionDifferences,
            suffix( ".import"),
            ".go" )
def runExpressionGO( infile, outfile ):
    '''run GO analysis on expression tracks.'''
    PExpression.runGO( infile, outfile, PARAMS["filename_go"] )

############################################################
############################################################
############################################################
## import GO analyses
############################################################
@transform( runExpressionGO,
            suffix(".go"),
            ("_go.import") )
def importExpressionGO( infile, outfile ):
    PExpression.importGO( infile, outfile, "go" )

@transform( runExpressionGOSlim,
            suffix(".goslim"),
            ("_goslim.import") )
def importExpressionGOSlim( infile, outfile ):
    PExpression.importGO( infile, outfile, "goslim" )

############################################################
############################################################
############################################################
@transform( buildGeneSet, suffix(".gtf"), "_roi.import" )
def importOverlapRegionsOfInterestEnsembl( infile, outfile ):
    '''compute overlap of intervals with regions of interest.'''

    tmpfilename = P.getTempFilename( "." )
    
    PGeneset.buildOverlapWithEnsembl( infile,  
                                      tmpfilename,
                                      "roi.bed" )

    tablename = outfile[:-len(".import")]

    statement = '''
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
            --allow-empty \
            --index=interval_id \
            --index=roi_id \
            --table=%(tablename)s \
        < %(tmpfilename)s > %(outfile)s
        '''

    P.run()
    
    os.unlink( tmpfilename )
    os.unlink( tmpfilename+".log" )

############################################################
############################################################
############################################################
@transform( buildGeneSet, suffix(".gtf"), "_gwas.import" )
def importGWASEnsembl( infile, outfile ):
    '''compute overlap of intervals with regions of interest.'''

    tmpfilename = P.getTempFilename( "." )
    
    PGeneset.buildOverlapWithEnsembl( infile,  
                                      tmpfilename,
                                      "gwas.bed" )

    tablename = outfile[:-len(".import")]

    statement = '''
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
            --allow-empty \
            --index=interval_id \
            --index=roi_id \
            --table=%(tablename)s \
        < %(tmpfilename)s > %(outfile)s
        '''

    P.run()
    
    os.unlink( tmpfilename )
    os.unlink( tmpfilename+".log" )

############################################################
############################################################
############################################################
@follows( buildGeneRegions, importRegionsOfInterest )
@files( [ ( ("%s.bed" % x, "%s.bed" % y), "%s_%s.import" % (x,y))  
        for x,y in itertools.product( TRACKS_ANNOTATOR, TRACKS_REGIONS) ] )
def importOverlapRegions( infiles, outfile ):
    '''compute overlap of infile with bedfile.
    '''

    infile, bedfile = infiles

    tablename = outfile[:-len(".import")]
    to_cluster = True
    statement = '''
        python %(scriptsdir)s/bed2graph.py \
            --output=name \
            --log=%(outfile)s \
            %(infile)s %(bedfile)s |\
        sed "s/name1/interval_id/; s/name2/roi_id/" |\
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
            --allow-empty \
            --index=interval_id \
            --index=roi_id \
            --table=%(tablename)s \
        > %(outfile)s
        '''
    
    P.run()

# ############################################################
# ############################################################
# ############################################################
# @follows( buildGeneRegions, importRegionsOfInterest )
# @files( [ ("%s.bed" % x, "%s_roi.import" % x ) for x in TRACKS_ANNOTATOR ] )
# def importOverlapRegionsOfInterest( infile, outfile ):
#     '''compute overlap of intervals with regions of interest.'''
#     importOverlapRegions( infile, outfile, "roi.bed" )
    
# ############################################################
# ############################################################
# ############################################################
# @follows( buildGeneRegions, importGWAS )
# @files( [ ("%s.bed" % x, "%s_gwas.import" % x ) for x in TRACKS_ANNOTATOR ] )
# def importOverlapGWAS( infile, outfile ):
#     '''compute overlap of intervals with regions of interest.'''
#     importOverlapRegions( infile, outfile, "gwas.bed" )


###################################################################
###################################################################
###################################################################
## compute the coverage of snps
###################################################################
@follows( importSNPsOfInterest, runMACS )
@transform(  "*.bam", 
             regex( r"(.*).bam"),
             inputs( (r"\1.bam", "snps_of_interest.import",) ), 
             r"\1.snpcoverage" )
def buildSNPCoverage( infiles, outfile ):

    infile, other = infiles
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()

    tmpfilename = P.getTempFilename()

    outf = open(outfile, "w")
    fields = pysam.Pileup.PileupSubstitution._fields
    outf.write( "snp\tnreads\t%s\n" % "\t".join( fields ))

    samfile = pysam.Samfile( infile, "rb" )
    
    for snp, contig, pos, reference in cc.execute( "SELECT snp, contig, pos, reference FROM snps_of_interest"):

        # count reads in neighbourhood
        start = pos - PARAMS["snp_halfwidth"]
        end = pos + PARAMS["snp_halfwidth"] + 1
        nreads = len( list(samfile.fetch( contig, start, end ) ))

        # collect pileup stats
        start = pos
        end = pos + 1
        statement = '''
        samtools pileup -f genome.fa -c -S <(samtools view -h %(infile)s %(contig)s:%(start)i-%(end)i ) 2>/dev/null > %(tmpfilename)s
        '''
        P.run( **dict( locals().items() + PARAMS.items() ) )

        for s in pysam.Pileup.iterate(open( tmpfilename, "r")):
            if s.position == pos:
                if s.reference_base != reference:
                    E.warn("mismatch for snp %s: expected %s, but got %s" % (snp, reference, s.reference_base))
                outf.write( "%s\t%i\t%s\n" % (snp, 
                                              nreads,
                                              "\t".join( map(str,s)) ) )
                break
        else:
            outf.write( "%s\t%i\t%s\n" % (snp, 
                                          nreads,
                                          "\t".join( [""] * len(fields)) ) )

    os.unlink( tmpfilename )

############################################################
############################################################
############################################################
@transform( buildSNPCoverage, suffix( ".snpcoverage" ), "_snpcoverage.import")
def importSNPCoverage( infile, outfile ):

    tablename = outfile[:-len(".import")]
    
    statement = '''
    sed "s/chromosome/contig/" < %(infile)s |\
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \ 
              --allow-empty \
              --index=snp \
              --index=genesymbol \
              --table=%(tablename)s \
    > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
##
############################################################

@follows( indexGenome, buildGenome,
          buildGeneRegions, buildGeneSet, buildPromotorRegions, buildTSSRegions,
          buildRepeats, buildTranscripts, buildTSSRegions,
          importTranscripts,
          importRegionsOfInterest, exportRegionsOfInterest,
          importSNPsOfInterest,
          exportReferenceAsBed,
          importUCSCEncodeTracks, exportUCSCEncodeTracks,
          importTranscriptInformation,
          importGeneInformation, importGeneStats,
          importGWAS,
          buildAnnotatorGeneTerritories,
          importRegionsOfInterestGenes )
def prepare():
    pass

@follows( buildIntervals, importReferenceIntervals )
def intervals():
    pass

@follows( buildBigwig, viewIntervals, viewBigwig, viewExpressionDifferences )
def export():
    pass

@follows( importAffymetrixAnnotation, buildProbeset2Transcript, importProbeset2Transcript,
          buildExpressionTracks, importExpressionTracks,
          importExpressionProbesets, importExpressionMap,
          buildExpressionCorrelation, importExpressionCorrelation,
          importExpressionDifferences, exportDifferentiallyExpressedGenes,
          runExpressionGOSlim, runExpressionGO,
          importExpressionGOSlim, importExpressionGO,
          )
def expression():
    '''run expression data related analyses.'''
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

@follows( runNubiscanNR, runNubiscanVDR,
          runNubiscanNRShuffled, runNubiscanVDRShuffled,
          runNubiscanNRShifted, runNubiscanVDRShifted,
          importNubiscan)
def run_nubiscan():
    pass

@follows( importCorrelation, 
          importOverlap,
          importUCSCOverlap,
          reproducibility)
def correlation():
    pass

@follows( annotateIntervals, importAnnotations, 
          annotateTSS, importTSS, 
          annotateRepeats, importRepeats,
          annotateTSSIntervalAssociations, importTSSIntervalAssociations,
          annotateTSSIntervalDistance, importTSSIntervalDistance, 
          annotateTracks, importTracks,
          buildIntervalCounts, importIntervalCounts,
          buildSNPCoverage, importSNPCoverage,
          importOverlapRegions,
          importOverlapRegionsOfInterestEnsembl,)
def annotate():
    pass

@follows( buildAnnotatorGC,
          makeAnnotatorDistance,
          makeAnnotatorArchitecture, importAnnotatorArchitecture,
          makeAnnotatorTracks, importAnnotatorTracks,
          makeAnnotatorPromotorsGO, importAnnotatorPromotorsGO, 
          makeAnnotatorPromotorsGOSlim, importAnnotatorPromotorsGOSlim,
          makeAnnotatorTerritoriesGO, importAnnotatorTerritoriesGO,
          makeAnnotatorTerritoriesGOSlim, importAnnotatorTerritoriesGOSlim,
          annotator_regions,
          annotator_proportion,
          annotator_motifs,
          annotator_roi,
          )
def annotator():
    pass

@follows( prepare, 
          intervals, 
          build_motifs, 
          find_motifs,
          expression,
          correlation, 
          annotate, 
          annotator 
          )
def full():
    pass

@files( ((".", "clean.log"),))
def clean( infile, outfile):
    '''remove all files not necessary for document creation.'''
    
    patterns = ("*.controlfasta",
                "*.mast.gz",
                "*.norm.bam",
                "*.norm.bam.bai",
                )
    
    cleaned = P.clean( patterns, dry_run = False )

    outf = open( outfile, "a" )

    outf.write( "filename\tsize\taccess\tmodified\tcreated\n" )
    
    for filename, statinfo in cleaned:
        outf.write( "\t".join( map(str, (filename,
                                         statinfo.st_size,
                                         statinfo.st_atime,
                                         statinfo.st_mtime,
                                         statinfo.st_ctime)) ) + "\n" )
    
    outf.close()


if __name__== "__main__":
    #P.checkFiles( ("genome.fasta", "genome.idx" ) )
    #P.checkExecutables( ("liftOver",) )
    sys.exit( P.main(sys.argv) )


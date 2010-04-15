################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_hvc.py 2861 2010-02-23 17:36:32Z andreas $
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
:Release: $Id: pipeline_hvc.py 2861 2010-02-23 17:36:32Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

snp annotation pipeline.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

"""
from ruffus import *
import sys, glob, subprocess, re, csv, os
import optparse, shutil, collections
import sqlite3
import numpy

import Experiment as E
import Pipeline as P
import IOTools

## global options
PARAMS = P.getParameters()

import GTF

TRACKS = ("Down","NonDiff","Up","UpAndDown")

@files( ((PARAMS["filename_segments"], "segments.import"), ))
def importSegments( infile, outfile ):
    tablename = outfile[:-len(".import")]

    statement = '''
    python %(scriptsdir)s/gtf2tab.py < %(infile)s |\
    csv2db.py %(csv2db_options)s \
              --index=gene_id \
              --table=%(tablename)s \
    > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

@files( ((PARAMS["filename_fasta"], "identifiers.import" ), ))
def importIdentifiers( infile, outfile):

    tablename = outfile[:-len(".import")]
    
    statement = '''
    grep ">" %(infile)s |\
    sed "s/>//; s/ .*//;" |\
    awk 'BEGIN {printf("gene_id\\n");} {print;}' |\
    csv2db.py %(csv2db_options)s \
              --index=gene_id \
              --table=%(tablename)s \
    > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

@merge( glob.glob( "../data/HVC*.csv" ), "assignments.import" )
def importAssignments( infiles, outfile):

    genes = collections.defaultdict( set )
    tags = []
    for infile in infiles:
        tag = re.sub(".*HVC", "", infile[:-len(".csv")])
        reader = csv.DictReader( open(infile,"rU") )
        for row in reader:
            try:
                gene_id = row["ESTIMA_ID"]
            except KeyError:
                print "Parsing error in line %s" % str(row)
                raise
            if gene_id == "": continue
            # remove ".A", ".B" and ".M" suffixes
            if gene_id[-2] == ".":
                gene_id = gene_id[:-2]

            genes[gene_id].add( tag )
            
        tags.append(tag)

    outf = P.getTempFile()
    outf.write( "gene_id\t%s\n" % "\t".join( tags) )
    for gene, present in genes.iteritems():
        x = []
        for t in tags:
            if t in present: x.append( "1" )
            else: x.append( "0" )
        outf.write( "%s\t%s\n" % (gene, "\t".join(x) ))

    outf.close()

    tablename = outfile[:-len(".import")]
    tmpfilename = outf.name
    
    statement = '''
    csv2db.py %(csv2db_options)s \
              --index=gene_id \
              --table=%(tablename)s \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

    os.unlink( tmpfilename )

@merge( (importSegments, importIdentifiers, importAssignments), "counts.table" )
def outputCounts( infiles, outfile ):
    '''number of hvc markers imported.'''
    
    dbhandle = sqlite3.connect( PARAMS["database"] )

    outf = open( outfile, "w")
    outf.write("track\timport\tmapped\t\twith_id\t\n" )
    
    for track in TRACKS:
        cc = dbhandle.cursor()
        statement = """SELECT COUNT(DISTINCT a.gene_id) FROM assignments AS a WHERE a.%(track)s""" % locals()
        result1 = cc.execute( statement ).fetchone()[0]
        
        statement = """SELECT COUNT(DISTINCT a.gene_id) FROM assignments AS a, segments as s WHERE s.gene_id = a.gene_id AND a.%(track)s""" % locals()
        result2 = cc.execute( statement ).fetchone()[0]

        statement = """SELECT COUNT(DISTINCT a.gene_id) FROM assignments AS a, identifiers as s WHERE s.gene_id = a.gene_id AND a.%(track)s""" % locals()
        result3 = cc.execute( statement ).fetchone()[0]

        outf.write( "\t".join( map( str, (track, result1,
                                          result2,
                                          "%5.2f" % (100.0 * result2/result1),
                                          result3,
                                          "%5.2f" % (100.0 * result3/result1)))) + "\n" )

    outf.close()
    
@merge( (importSegments, importAssignments), "tracks.gtf" )
def exportLocations( infiles, outfile ):

    outf = open( outfile, "w" )

    dbhandle = sqlite3.connect( PARAMS["database"] )

    for track in TRACKS:
        outf.write( 'track name=%(track)s description="%(track)s"\n' % locals() )
    
        cc = dbhandle.cursor()
        # ignores the attributes
        statement = """SELECT DISTINCT contig, end, feature, frame, s.gene_id, score, source, start, strand, transcript_id \
                         FROM segments AS s, assignments AS a WHERE a.gene_id = s.gene_id and a.%(track)s""" % locals()
        cc.execute( statement )

        for row in cc:
            gtf = GTF.Entry()
            gtf.contig, gtf.end, gtf.feature, gtf.frame, gtf.gene_id, gtf.score, gtf.source, gtf.start, gtf.strand, gtf.transcript_id =\
                         row
            outf.write( str(gtf) + "\n" )

    outf.close()    

@follows( importSegments, importAssignments )
@files( [(importSegments, "%s.gtf" % x) for x in TRACKS] )
def exportSegments( infiles, outfile ):

    track = outfile[:-len(".gtf")]
    
    outf = open( outfile, "w" )

    dbhandle = sqlite3.connect( PARAMS["database"] )

    cc = dbhandle.cursor()
    # ignores the attributes
    statement = """SELECT DISTINCT contig, end, feature, frame, s.gene_id, score, source, start, strand, transcript_id \
    FROM segments AS s, assignments AS a WHERE a.gene_id = s.gene_id and a.%(track)s""" % locals()

    cc.execute( statement )

    for row in cc:
        gtf = GTF.Entry()
        gtf.contig, gtf.end, gtf.feature, gtf.frame, gtf.gene_id, gtf.score, gtf.source, gtf.start, gtf.strand, gtf.transcript_id =\
                     row
        outf.write( str(gtf) + "\n" )

    outf.close()    

@files( PARAMS["genome"]+".fasta", "windows.bed") 
def buildWindows( infile, outfile ):
    
    statement = '''python %(scriptsdir)s/fasta2bed.py \
    --method=fixed-width-windows \
    --window-size=%(window_size)s \
    --log=%(outfile)s.log \
    < %(genome)s.fasta > %(outfile)s'''

    P.run( **dict( locals().items() + PARAMS.items() ) )    

@follows( buildWindows)
@transform( exportSegments, suffix(".gtf"), ".overlap" )
def computeOverlap( infile, outfile ):
    '''compute overlap between markers and the windows used for clustering.
    '''
    
    statement = '''python %(scriptsdir)s/gff2table.py 
    --filename-windows=<(python %(scriptsdir)s/bed2gff.py < windows.bed) 
    --decorator=counts
    --filename-data=%(infile)s \
    --skip-empty \
    --is-gtf \
    --log=%(outfile)s.log \
    < %(genome)s.fasta > %(outfile)s'''

    P.run( **dict( locals().items() + PARAMS.items() ) )    

@transform( computeOverlap, suffix(".overlap"), "_overlap.import" )
def importOverlap( infile, outfile ):

    tablename = outfile[:-len(".import")]
    
    statement = '''
    csv2db.py %(csv2db_options)s \
              --table=%(tablename)s \
    < %(infile)s > %(outfile)s
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

@merge( (importOverlap), "overlaps.table" )
def outputOverlapHistogram( infiles, outfile ):
    '''number of hvc markers imported.'''
    
    dbhandle = sqlite3.connect( PARAMS["database"] )

    max_counts = 100
    outf = open( outfile, "w")

    data = numpy.zeros( ( len(TRACKS), max_counts ) )
    observed_max = 0
    
    for x, track in enumerate(TRACKS):
        cc = dbhandle.cursor()
        statement = """SELECT ngenes, COUNT(*) FROM %(track)s_overlap GROUP by ngenes""" % locals()
        cc.execute( statement )
        for ngenes, c in cc:
            observed_max = max( ngenes, observed_max )
            data[x,ngenes] = c

    outf.write("counts\t%s\n" % "\t".join(TRACKS))
    for y in range(1, observed_max + 1 ):
        outf.write( "%i\t%s\n" % (y, "\t".join( map(str, data[:,y] ))) )

    outf.close()

@transform( importOverlap, suffix("_overlap.import"), "_clusters.bed" )
def exportClusters( infile, outfile ):
    '''number of hvc markers imported.'''

    tablename = infile[:-len(".import")]
    dbhandle = sqlite3.connect( PARAMS["database"] )

    outf = open( outfile, "w")
    cc = dbhandle.cursor()
    statement = """SELECT contig, start, end, ngenes FROM %(tablename)s WHERE ngenes >= 2""" % locals()
    cc.execute( statement )
    id = 0
    for contig, start, end, ngenes in cc:
        id += 1
        outf.write("%s\t%i\t%i\t%i\t%i\n" % (contig, start, end, id, ngenes ) )

    outf.close()

@merge( (PARAMS["filename_segments"], PARAMS["filename_ensembl"]),
        "ensembl.diff" )
def diffEnsembl( infiles, outfile ):
    '''compute overlap between markers and ENSEMBL genes.
    '''

    to_cluster = True
    
    infile1, infile2 = infiles
    
    statement = '''
    python %(scriptsdir)s/diff_gtf.py \
        -p --write-equivalent --ignore-strand --output-pattern=%(outfile)s.%%s
    %(infile1)s %(infile2)s > %(outfile)s
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

@follows( buildWindows, diffEnsembl )
@transform( exportSegments, suffix(".gtf"), ".coding_overlap" )
def computeOverlapCoding( infile, outfile ):
    '''compute overlap between coding markers and windows.

    This is done by setting the gene_id and transcript_id of markers to the ENSEMBL gene id
    and transcript_id that it overlaps with. Markers not overlapping an ENSEMBL gene id
    are removed.
    '''
    
    to_cluster = True
    tmpfilename = P.getTempFilename( dir = "." )
    
    statement = '''python %(scriptsdir)s/gtf2gtf.py
    --rename=gene \
    --apply=ensembl.diff.genes_ovl \
    < %(infile)s > %(tmpfilename)s
    '''
    
    P.run( **dict( locals().items() + PARAMS.items() ) )

    statement = '''python %(scriptsdir)s/gff2table.py 
    --filename-windows=<(python %(scriptsdir)s/bed2gff.py < windows.bed) 
    --decorator=counts
    --filename-data=%(tmpfilename)s \
    --skip-empty \
    --is-gtf \
    --log=%(outfile)s.log \
    < %(genome)s.fasta > %(outfile)s'''

    P.run( **dict( locals().items() + PARAMS.items() ) )

    os.unlink( tmpfilename )

@follows( buildWindows, diffEnsembl )
@transform( exportSegments, suffix(".gtf"), ".go_overlap" )
def computeOverlapGO( infile, outfile ):
    '''compute overlap between codingmarkers and windows.
    Only markers of certain GO categories are counted.

    This is done by setting the gene_id and transcript_id of markers of the
    ENSEMBEL gene that it overlaps with. This list is filtered first to
    keep only those ids with valid GO associations
    
    '''
    
    to_cluster = False

    filter_goid = set(IOTools.readList( open( PARAMS["filename_gofilter"] ) ))
    filter_genes = set()

    E.info( "number of goids: %i" % len(filter_goid))
    
    for l in open( PARAMS["filename_go"]):
        f, id, goid, desc, evd = l[:-1].split("\t")[:5]
        if goid in filter_goid:
            filter_genes.add( id )

    tmpfile1 = P.getTempFile( dir = "." )

    for line in open("ensembl.diff.genes_ovl" ):

        a,b = line[:-1].split( "\t" )
        if b not in filter_genes: continue
        tmpfile1.write(line)
        
    E.info( "number of genes taken: %i" % len(filter_genes))
    
    tmpfile1.close()
    tmpfilename1 = tmpfile1.name

    tmpfilename = P.getTempFilename( dir = "." )
    
    statement = '''python %(scriptsdir)s/gtf2gtf.py
    --rename=gene \
    --apply=%(tmpfilename1)s \
    < %(infile)s > %(tmpfilename)s
    '''
    
    P.run( **dict( locals().items() + PARAMS.items() ) )
    
    statement = '''python %(scriptsdir)s/gff2table.py 
    --filename-windows=<(python %(scriptsdir)s/bed2gff.py < windows.bed) 
    --decorator=counts
    --filename-data=%(tmpfilename)s \
    --skip-empty \
    --is-gtf \
    --log=%(outfile)s.log \
    < %(genome)s.fasta > %(outfile)s'''

    P.run( **dict( locals().items() + PARAMS.items() ) )

    os.unlink( tmpfilename )

@transform( computeOverlapCoding, suffix(".coding_overlap"), ".clusters_level1" )
def filterClustersLevel1( infile, outfile ):
    '''filter clusters.

    Only keep clusters that:
    1. have at least 2 ESTs
    2. have only one protein coding EST
    '''

    track = infile[:-len(".coding_overlap")]
    infile_clusters = track + ".overlap"

    inf_all, inf_filter = open(infile_clusters, "r"), open( infile, "r")
    
    def _parse(l):
        contig, start, end, ngenes = l[:-1].split("\t")[:4]
        return (contig, int(start), int(end)), int(ngenes)
    
    this = None
    # skip headers
    other_line = inf_filter.readline()
    line = inf_all.readline()

    outf = open(outfile, "w")
    outf.write( line[:-1] + "\tnprotein_coding\n" )
    
    while 1:
        other_line = inf_filter.readline()
        if not other_line: break
        if other_line.startswith("#"): continue
        
        other, other_ngenes = _parse(other_line)
        
        while this != other:
            line = inf_all.readline()
            if not line: break
            if line.startswith("#"): continue
            this, ngenes = _parse( line )

        if not line: break
        
        if ngenes >= 2 and other_ngenes == 1:
            outf.write( line[:-1] + "\t%i\n" % other_ngenes)

    outf.close()

@transform( computeOverlapGO, suffix(".go_overlap"), ".clusters_level2" )
def filterClustersLevel2( infile, outfile ):
    '''filter clusters.

    Only keep clusters that:
    1. have at least one gene with GO annotation
    '''

    track = infile[:-len(".go_overlap")]
    infile_clusters = track + ".clusters_level1"

    inf_all, inf_filter = open(infile_clusters, "r"), open( infile, "r")
    
    def _parse(l):
        contig, start, end, ngenes = l[:-1].split("\t")[:4]
        return (contig, int(start), int(end)), int(ngenes)
    
    this = None
    # skip headers
    other_line = inf_filter.readline()
    line = inf_all.readline()

    outf = open(outfile, "w")
    outf.write( line[:-1] + "\twith_go\n" )
    
    while 1:
        other_line = inf_filter.readline()
        if not other_line: break
        if other_line.startswith("#"): continue
        
        other, other_ngenes = _parse(other_line)
        
        while this != other:
            line = inf_all.readline()
            if not line: break
            if line.startswith("#"): continue
            this, ngenes = _parse( line )

        if not line: break
        
        if other_ngenes == 1:
            outf.write( line[:-1] + "\t%i\n" % other_ngenes)

    outf.close()

@follows( importSegments, outputCounts,
          exportLocations, buildWindows,
          computeOverlap, importOverlap, outputOverlapHistogram,
          exportClusters, diffEnsembl )
def full():
    pass

if __name__== "__main__":
    P.checkFiles( ("genome.fasta", "genome.idx" ) )
    sys.exit( P.main(sys.argv) )

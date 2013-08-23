################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_kamilah.py 2869 2010-03-03 10:20:13Z andreas $
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
:Release: $Id: pipeline_kamilah.py 2869 2010-03-03 10:20:13Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Kamilah pipeline

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
import csv
import gzip
from ruffus import *
import sqlite3

import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGAT.GFF as GFF
import CGAT.GTF as GTF
import CGAT.Database as Database

import PipelineGeneset as PGeneset
import PipelineAnnotator as PAnnotator

if os.path.exists("conf.py"):
    execfile("conf.py")
else:
    EXPERIMENTAL_TRACKS = []

PARAMS = P.getParameters(
    defaults = { 'ancestral_repeats_filename': None } )

PARAMS.update( {
    "annotation" : "regions.gff",
    "genes": "genes.gtf.gz",
    "merged": "merged.gtf.gz",
    "transcripts": "transcripts.gtf.gz" } )

def getSourceTrack( track, all_tracks ):
    '''if track is in derived tracks, get the source track.

    returns None if track is not a derived track
    '''
    if len(all_tracks) == 0: return None

    # get rid of any extensions
    all_tracks = [ re.sub("\..*$", "", x) for x in all_tracks ]
    track = re.sub("\..*$", "", track )

    if len(all_tracks) == 1:
        if len(os.path.commonprefix( (track, all_tracks[0]))) > 0:
            return all_tracks[0]
        else:
            return None

    # get all tracks with a common prefix of length 3 or more
    prefixes = [ t for t in all_tracks if len(os.path.commonprefix( (track, t) ) ) > 3 ]

    prefixes.sort( key = lambda x: len(x) )
    # return the shortest
    return prefixes[0]

def getRelatedTracks( track, all_tracks ):
    '''return tracks in ``all_tracks`` that are related to ``track`` (including itself)

    related tracks are build by merging or slicing another track.
    '''

    source = getSourceTrack( track, all_tracks )

    if not source: source = track

    related = set([x for x in all_tracks if x.startswith( source ) ])

    if track not in related: related.add( track )
    
    for x in related:
        if x in EXPERIMENTAL_TRACKS:
            related.add( PARAMS["merged"] )
            break

    return list(related)

###################################################################
###################################################################
@files( [ (x, "%sLinc.import" % x[:-len(".gtf.gz")]) for x in EXPERIMENTAL_TRACKS] )
def importLincRNA( infile, outfile ):
    '''build a linc RNA set.

        * no coding potential
        * unknown and intergenic transcripts
        * no overlap with ``linc_exclude`` (usually: human refseq)
        * at least ``linc_min_length`` bp in length
        * at least ``linc_min_reads`` reads in transcript

    '''

    table = outfile[:-len(".import")]
    track = table[:-len("Linc")]

    dbhandle = sqlite3.connect( PARAMS["database"] )

    Database.executewait( dbhandle, '''DROP TABLE IF EXISTS %(table)s''' % locals())
    Database.executewait( dbhandle, '''CREATE TABLE %(table)s (gene_id TEXT)''' % locals())
    Database.executewait( dbhandle, '''CREATE INDEX %(table)s_index1 ON %(table)s (gene_id)''' % locals())

    joins, wheres = [], ["1"]

    if PARAMS["linc_min_reads"] > 0:
        joins.append( ", %(track)s_coverage as cov" % locals() )
        wheres.append( "cov.gene_id = m.gene_id2 AND cov.nmatches >= %(i)" % PARAMS["linc_min_reads"] )

    if PARAMS["linc_exclude"] > 0:
        joins.append( "LEFT JOIN %s_vs_%s_ovl as ovl on ovl.gene_id2 = a.gene_id" %\
                      (PARAMS["linc_exclude"], track ) )
        wheres.append( "ovl.gene_id1 IS NULL" )

    wheres = " AND ".join( wheres )
    joins = " ".join( joins )

    statement = '''INSERT INTO %(table)s 
                   SELECT DISTINCT(a.gene_id) FROM 
                          %(track)s_annotation as a
                          %(joins)s
                          LEFT JOIN %(track)s_coding AS c on c.gene_id = a.gene_id 
                    WHERE is_unknown 
                          AND is_intergenic 
                          AND exons_sum >= %(linc_min_length)i
                          AND (c.is_coding IS NULL or not c.is_coding) 
                          AND %(wheres)s
                     ''' % dict( PARAMS.items() + locals().items() )

    E.debug( "statement to build lincRNA: %s" % statement)

    Database.executewait( dbhandle, statement % locals())

    dbhandle.commit()
    
    cc = dbhandle.cursor()
    result = cc.execute("SELECT COUNT(*) FROM %(table)s" % locals() ).fetchall()[0][0]
    E.info( "build lincRNA set for %s: %i entries" % ( track, result ))

    outgtf = "%s.gtf.gz" % table

    E.info( "creating gtf file `%s`" % outgtf )
    # output gtf file
    statement = '''%(cmd-sql)s %(database)s "SELECT g.* FROM %(track)s_gtf as g, %(table)s AS t
                          WHERE t.gene_id = g.gene_id" 
                | python %(scriptsdir)s/gtf2tsv.py --invert --log=%(outfile)s
                | gzip
                > %(outgtf)s'''
    
    P.run()

###################################################################
###################################################################
@files( PARAMS["ancestral_repeats_filename"], "repeats.gtf.gz" )
def buildRepeatTrack( infile, outfile ):
    '''build a repeat track as negative control.'''

    nrepeats = 0
    for gff in GFF.iterator( gzip.open(infile, "r" ) ): nrepeats+=1
    sample = set( random.sample( xrange( nrepeats), PARAMS["ancestral_repeats_samplesize"]) )

    outf = gzip.open( outfile, "w" )
    gtf = GTF.Entry()
    for x,gff in enumerate( GFF.iterator( gzip.open(infile, "r" ) ) ):
        if not x in sample: continue
        gtf.fromGFF( gff, "%08i" % x, "%08i" % x )
        outf.write( "%s\n" % str(gtf) )
    outf.close()

    E.debug( "created sample of %i repeats out of %i in %s" % (len(sample), nrepeats, outfile))

###################################################################
###################################################################
@files( PARAMS["ancestral_repeats_filename"], 
        "introns.gtf.gz" )
def buildIntronTrack( infile, outfile ):
    '''build an intron track as negative control.

    The ``intron`` track is build from ancestral repeats. Here,
    repeats are taken as exons, and artifical genes are constructed
    from 2 to 5 ancestral repeats at least 100 bases apart, but
    within 10 kb of each other.

    '''
    to_cluster = True

    statement = '''gunzip 
        < %(infile)s 
	| %(scriptsdir)s/gff_sort pos
	| python %(scriptsdir)s/gff2gff.py --join=100,10000,2,5 --log=%(outfile)s.log 
	| python %(scriptsdir)s/gtf2gtf.py 
                --filter=gene --sample-size=%(ancestral_repeats_samplesize)i --log=%(outfile)s.log 
        | gzip
        > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
@files( [ ("%s.fasta" % PARAMS["genome"], "%s.table" % PARAMS["genome"]), ] )
def buildGenomeInformation( infile, outfile ):
    '''compute genome composition information.'''

    to_cluster = True

    statement = '''
    python %(scriptsdir)s/analyze_codonbias_shannon.py 
        --sections=length,na
    < %(infile)s
    > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
@transform( buildGenomeInformation, suffix(".table"), "_table.import" )
def importGenomeInformation( infile, outfile ):
    '''import genome information.'''
    table = outfile[:-len(".import")]

    statement = '''
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --table=%(table)s 
        < %(infile)s
        > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
@files( ( ( (PARAMS["repeats_filename"], PARAMS["genome"] + ".idx" ),
            "repeats_table.import" ), ) )
def importRepeatInformation( infiles, outfile ):
    '''import genome information.'''
    
    to_cluster = True

    table = outfile[:-len(".import")]

    repeatsfile, indexfile = infiles

    tmpfilename = P.getTempFilename( "." )

    statement = '''awk '{printf("%%s\\t0\\t%%i\\n", $1, $4)}' < %(indexfile)s > %(tmpfilename)s'''
    P.run()

    statement = '''
        gunzip < %(repeatsfile)s 
        | python %(scriptsdir)s/gff2bed.py -v 0 
        | coverageBed -a stdin -b %(tmpfilename)s
        | awk 'BEGIN { printf("contig\\tstart\\tend\\tnover_entries\\tnover_bases\\tlength\\tpover\\n" );} {print;}'
        |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --table=%(table)s 
        > %(outfile)s
    '''
    P.run()

    os.unlink( tmpfilename )

###################################################################
###################################################################
@files( [ (PARAMS["%s_merge" % x], "%s.gtf.gz" % x) for x in P.asList(PARAMS["merge"])] +\
            [ (EXPERIMENTAL_TRACKS, PARAMS["merged"] ) ] )
def buildMergedTracks( infiles, outfile ):
    '''merge tracks.'''

    infiles = " ".join(infiles)
    statement = '''
	zcat %(infiles)s 
        | python %(scriptsdir)s/gff2psl.py 
                 --log=%(outfile)s.log 
                 --is-gtf 
                 --allow-duplicates 
	| python %(scriptsdir)s/psl2psl.py 
                 --log=%(outfile)s.log 
                 --method=rename-query 
                 --unique 
                 --output-filename-map=%(outfile)s.id2new 
	| sort -k 14,14 -k16,16n 
        | %(cmd-farm)s
		--split-at-column=14 
		--output-header 
		--renumber="M%%06i" 
		--renumber-column=":id" 
		--log=%(outfile)s.log 
		--subdirs 
		--max-files=60 
	"python %(scriptsdir)s/psl2assembly.py 
		--method=locus 
		--threshold-merge-distance=0 
		--threshold-merge-overlap=1 
		--staggered=all 
		--log=%(outfile)s.log 
		--genome=genome 
		--output-filename-pattern=%%DIR%%%(outfile)s.%%s" 
	> %(outfile)s.new2locus
    '''
    P.run()

    statement = '''
	python %(scriptsdir)s/psl2gff.py 
		--log=%(outfile)s.log 
		--as-gtf
	< %(outfile)s.locus.psl 
        | gzip
        > %(outfile)s
    '''
    P.run()

###################################################################
###################################################################
@files( None, PARAMS["geneset_gtf"] )
def exportEnsembl( infile, outfile ):
    '''export gtf file with ensembl transcripts.
    '''
    tmpfile = P.getTempFilename()

    statement = '''
         perl %(scriptsdir)s/ensembl2gtf.pl 
            -dbname %(mysql_database_ensembl)s
            -host %(mysql_host)s
            -user %(mysql_user)s
            -dbpass %(mysql_pass)s
            -dnadbname %(mysql_database_ensembl)s
            -dnahost %(mysql_host)s
            -dnauser %(mysql_user)s
            -dnapass %(mysql_pass)s
            -gtffile %(tmpfile)s
            -schema '%(ensembl_schema)s' 
            -coordsystem %(ensembl_coordsystem)s
            -genetypes %(ensembl_genetypes)s > %(outfile)s.log'''

    P.run()

    statement = 'gzip < %(tmpfile)s > %(outfile)s'
    P.run()

    os.unlink( tmpfile )

###################################################################
###################################################################
###################################################################
## gene set section
###################################################################
############################################################
############################################################
############################################################
@files( exportEnsembl, PARAMS['annotation'] )
def buildGeneRegions( infile, outfile ):
    '''annotate genomic regions with reference gene set.

    Only considers protein coding genes. In case of overlapping
    genes, only take the longest (in genomic coordinates).
    Genes not on UCSC contigs are removed.
    '''
    PGeneset.buildGeneRegions( infile, outfile )

############################################################
############################################################
############################################################
@follows( buildGeneRegions )
@files( exportEnsembl, PARAMS['genes'] )
def buildGenes( infile, outfile ):
    '''build a collection of exons from the protein-coding
    section of the ENSEMBL gene set. The exons include both CDS
    and UTR.

    The set is filtered in the same way as in :meth:`buildGeneRegions`.
    '''
    PGeneset.buildProteinCodingGenes( infile, outfile )

############################################################
############################################################
############################################################
@files( exportEnsembl, "gene_info.import" )
def importGeneInformation( infile, outfile ):
    '''import the transcript set.'''
    PGeneset.importGeneInformation( infile, 
                                    outfile, 
                                    only_proteincoding = PARAMS["geneset_only_proteincoding"] )

############################################################
############################################################
############################################################
@files( buildGenes, "gene_stats.import" )
def importGeneStats( infile, outfile ):
    '''import the transcript set.'''

    PGeneset.importGeneStats( infile, outfile )

############################################################
############################################################
############################################################
@files( exportEnsembl, PARAMS["transcripts"] )
def buildTranscripts( infile, outfile ):
    '''build a collection of transcripts from the protein-coding
    section of the ENSEMBL gene set.

    Only CDS are used.
    '''
    PGeneset.buildProteinCodingTranscripts( infile, outfile )

############################################################
############################################################
############################################################
@transform( buildTranscripts, suffix(".gtf.gz"), "_gtf.import" )
def importTranscripts( infile, outfile ):
    '''import the transcript set.'''
    PGeneset.importTranscripts( infile, outfile )

############################################################
############################################################
############################################################
@files( buildTranscripts, "transcript_stats.import" )
def importTranscriptStats( infile, outfile ):
    '''import the transcript set.'''

    PGeneset.importTranscriptStats( infile, outfile )

############################################################
############################################################
############################################################
@files( exportEnsembl, "transcript_info.import" )
def importTranscriptInformation( infile, outfile ):
    '''import the transcript set.'''
    PGeneset.importTranscriptInformation( infile, outfile,
                                          only_proteincoding = PARAMS["geneset_only_proteincoding"] )

############################################################
############################################################
############################################################
@files( PARAMS["geneset_pep"], "protein_stats.import" )
def importProteinStats( infile, outfile ):
    '''import the transcript set.'''

    PGeneset.importProteinStats( infile, outfile )


############################################################
############################################################
############################################################
@files( (PARAMS["rates_alignment"], 
         PARAMS["genes"]), 
        "alignment_filtered.psl.gz" )
def buildFilteredAlignment( infiles, outfile ):
    '''build a genomic alignment without exons
    from the reference gene set.

    This alignment is used for rate computation within introns.
    '''
    infile_alignment, infile_genes = infiles
    to_cluster = True
    statement = '''gunzip 
	< %(infile_alignment)s 
	| python %(scriptsdir)s/psl2psl.py 
		--method=filter-remove 
		--filter-query=<(gunzip < %(infile_genes)s)
		--log=%(outfile)s.log 
        > %(outfile)s '''

    P.run()

@transform(buildFilteredAlignment, 
           suffix(".psl.gz"), 
           ".stats" )
def buildAlignmentStats( infile, outfile ):

    to_cluster = True
    
    statement = '''
	python %(scriptsdir)s/psl2stats.py 
	    --log=%(outfile)s.log 
	< %(infile)s
        | gzip
        > %(outfile)s'''
    P.run()

###################################################################
###################################################################
###################################################################
@transform( TRACKS, suffix(".gtf.gz"), "_gtf.import" )
def importGTF( infile, outfile ):
    '''import gtf files.'''
    
    table = outfile[:-len(".import")]

    statement = '''gunzip
        < %(infile)s
        | python %(scriptsdir)s/gtf2tsv.py 
        |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=gene_id 
              --map=gene_id:str 
              --index=transcript_id 
              --map=transcript_id:str 
              --table=%(table)s 
        > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
@transform( TRACKS,
            regex("(.*).gtf.gz"), 
            inputs( (r"\1.gtf.gz", buildGeneRegions) ),
            r"\1.annotation")
def makeAnnotations( infiles, outfile ):
    '''annotate transcripts by location (intergenic, intronic, ...)'''
    
    infile, annotation = infiles

    statement = '''gunzip 
    < %(infile)s 
    | %(scriptsdir)s/gff_sort gene 
    | %(cmd-farm)s --split-at-column=1 --output-header --log=%(outfile)s.log --max-files=60 
	"python %(scriptsdir)s/gtf2table.py 
		--counter=position 
		--counter=classifier 
		--section=exons 
		--section=introns 
		--counter=length 
		--counter=splice 
		--counter=composition-na 
		--counter=splice-comparison 
		--log=%(outfile)s.log 
                --filename-format=gff
		--filename-gff=%(annotation)s 
		--genome-file=%(genome)s" 
    > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
@transform( makeAnnotations, suffix(".annotation"), "_annotation.import")
def importAnnotations( infile, outfile ):
    '''import annotations'''
    
    table = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=gene_id 
              --map=gene_id:str 
              --table=%(table)s 
    < %(infile)s
    > %(outfile)s'''

    P.run()

########################################################
########################################################
########################################################
@transform( TRACKS,
            regex("(.*).gtf.gz"), 
            inputs( (r"\1.gtf.gz", buildGeneRegions)),
            r"\1.overrun")
def makeOverrun( infiles, outfile ):
    '''compute intron overrun.'''

    infile, annotation = infiles

    statement = '''gunzip 
    < %(infile)s 
    | %(scriptsdir)s/gff_sort gene 
    | %(cmd-farm)s --split-at-column=1 --output-header --log=%(outfile)s.log --max-files=60 
	"python %(scriptsdir)s/gtf2table.py 
                --filename-format=gff
		--counter=overrun 
		--log=%(outfile)s.log 
                --filename-format=gff 
		--filename-gff=%(annotation)s" 
    > %(outfile)s 
    '''
    P.run()

############################################################
############################################################
############################################################
@transform( makeOverrun, suffix(".overrun"), "_overrun.import")
def importOverrun( infile, outfile ):
    '''import annotations'''
    
    table = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=gene_id 
              --map=gene_id:str 
              --table=%(table)s 
    < %(infile)s
    > %(outfile)s'''

    P.run()

########################################################
@transform( TRACKS, 
            regex("(.*).gtf.gz"), 
            inputs( (r"\1.gtf.gz", buildGenes) ),
            r"\1.distances")
def makeDistances( infiles, outfile ):
    '''compute intron overrun.'''

    infile, annotation = infiles

    statement = '''gunzip
    < %(infile)s 
    | %(scriptsdir)s/gff_sort gene 
    | %(cmd-farm)s --split-at-column=1 --output-header --log=%(outfile)s.log --max-files=60 
	"python %(scriptsdir)s/gtf2table.py 
		--counter=distance-genes 
		--log=%(outfile)s.log 
		--filename-gff=<( gunzip < %(annotation)s ) " 
    > %(outfile)s 
    '''
    P.run()

############################################################
@transform(  makeDistances, suffix(".distances"), "_distances.import")
def importDistances( infile, outfile ):
    '''import annotations'''
    
    table = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=gene_id 
              --index=closest_id
              --map=gene_id:str 
              --table=%(table)s 
    < %(infile)s 
    > %(outfile)s'''

    P.run()


########################################################
@transform( TRACKS, 
            regex("(.*).gtf.gz"), 
            r"\1.segments")
def makeSegments( infile, outfile ):
    '''compute intron overrun.'''

    to_cluster = True

    statement = '''gunzip < %(infile)s 
    | %(scriptsdir)s/gff_sort pos 
    | python %(scriptsdir)s/gff2histogram.py 
		--method=values 
		--output-filename-pattern="%(outfile)s.%%s"
		--force 
		--log=%(outfile)s.log 
    > %(outfile)s 
    '''
    P.run()

    statement = '''gunzip 
    < %(infile)s 
    | %(scriptsdir)s/gff_sort gene-pos 
    | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts 
    | %(scriptsdir)s/gff_sort pos 
    | python %(scriptsdir)s/gff2histogram.py 
		--method=values 
		--force 
		--output-filename-pattern="%(outfile)s_genes.%%s" 
		--log=%(outfile)s.log
    >> %(outfile)s'''
    P.run()

############################################################
@transform(  makeSegments, suffix(".segments"), "_segments.import")
def importSegments( infile, outfile ):
    '''import segments'''
    
    table = outfile[:-len(".import")]
    
    for x in (".distances", ".sizes", ".overlaps", 
              "_genes.distances", "_genes.sizes", "_genes.overlaps" ):
        y = re.sub("\.", "_", x)
        statement = '''
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=gene_id 
              --map=gene_id:str 
              --table=%(table)s%(y)s 
        < %(infile)s%(x)s
        >> %(outfile)s'''

        P.run()

########################################################
@transform( TRACKS, 
            regex("(.*).gtf.gz"), 
            inputs( (r"\1.gtf.gz", PARAMS["ancestral_repeats_filename"]) ),
            r"\1.repeats")
def makeRepeats( infiles, outfile ):
    '''compute overlap with repeats
    '''

    infile, annotation = infiles

    statement = '''gunzip
    < %(infile)s 
    | %(scriptsdir)s/gff_sort gene 
    | %(cmd-farm)s --split-at-column=1 --output-header --log=%(outfile)s.log --max-files=60 
	"python %(scriptsdir)s/gtf2table.py 
                --filename-format=gff
		--counter=overlap 
		--log=%(outfile)s.log 
		--filename-gff=<( gunzip < %(annotation)s )" 
    > %(outfile)s 
    '''
    P.run()

############################################################
@transform(  makeRepeats, suffix(".repeats"), "_repeats.import")
def importRepeats( infile, outfile ):
    '''import repeat overlap'''
    
    table = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=gene_id 
              --map=gene_id:str 
              --table=%(table)s 
    < %(infile)s 
    > %(outfile)s'''

    P.run()


############################################################
@files( [ ( 
            (x, "%s" % y),
            "%s_vs_%s.diff" % (x[:-len(".gtf.gz")], y[:-len(".gtf.gz")])) \
              for x,y in itertools.product( OVERLAP_TRACKS, TRACKS) ] )
def makeDifference( infiles, outfile ):
    '''compute the difference between two sets

    This target computes three differences
    1. exon overlap %.diff.{diff,overlap,total,genes_ovl,genes_total,genes_uniq1,genes_uniq2}
    An intronic transcript will no overlap with an exonic transrcipt.
    2. overlap statistics on each gene in set 1: %.diff
    3. intron/exon overlap %_genes.diff.{diff,overlap,total,genes_ovl,genes_total,genes_uniq1,genes_uniq2}
    An intronic transcript will overlap with an exonic transcript.
    '''
    
    to_cluster = True

    first, last = [x[:-len(".gtf.gz")] for x in infiles]

    statement = '''
         python %(scriptsdir)s/diff_gtf.py 
		--log=%(outfile)s.log 
		-p --write-equivalent --ignore-strand --output-pattern="%(outfile)s.%%s" 
		<(gunzip < %(first)s.gtf.gz)
		<(gunzip < %(last)s.gtf.gz)
	> %(outfile)s.log'''
    P.run()

    statement = '''python %(scriptsdir)s/diff_gtf.py 
		--log=%(outfile)s.log 
		-p --write-equivalent --ignore-strand --output-pattern="%(outfile)s_genes.%%s" 
		<( gunzip < %(first)s.gtf.gz | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --log=%(outfile)s.log ) 
		<( gunzip < %(last)s.gtf.gz | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts --log=%(outfile)s.log )
	> %(outfile)s.log'''
    P.run()
    
    # note that the option --counter=overlap-transcripts
    # has been set to --counter=overlap
    statement = '''gunzip
                < %(first)s.gtf.gz 
                | python %(scriptsdir)s/gtf2table.py 
		--log=%(outfile)s.log 
		--counter=length 
		--counter=overlap
		--counter=coverage 
		--filename-gff=<( gunzip < %(last)s.gtf.gz)
                > %(outfile)s'''
    
    P.run()

############################################################
@transform(  makeDifference, suffix(".diff"), "_diff.import")
def importDifference( infile, outfile ):
    '''import overlap'''

    table = outfile[:-len("_diff.import")]

    statement = '''
    	grep -v "\\bna\\b" 
        < %(infile)s 
	| python %(toolsdir)s/csv_cut.py --large --remove cov_values 
        |python %(scriptsdir)s/csv2db.py --allow-empty 
               %(csv2db_options)s 
              --index=gene_id 
              --map=gene_id:str
              --table=%(table)s    
        > %(outfile)s
    '''
    P.run()

    i = infile + ".genes_ovl"
    if os.path.exists( i ):
        statement = '''
       python %(scriptsdir)s/csv2db.py --allow-empty 
               %(csv2db_options)s 
               --map=gene_id1:str 
               --map=gene_id2:str 
               --index=gene_id1 
               --index=gene_id2
              --table=%(table)s_ovl
        < %(i)s
        > %(outfile)s'''
        P.run()
    else:
        E.warn( "file %s does not exist: not imported" % i )

    i = infile + "_genes.genes_ovl"
    if os.path.exists( i ):
        statement = '''
       python %(scriptsdir)s/csv2db.py --allow-empty 
               %(csv2db_options)s 
               --map=gene_id1:str 
               --map=gene_id2:str 
               --index=gene_id1 
               --index=gene_id2
              --table=%(table)s_geneovl
        < %(i)s
        > %(outfile)s'''
        P.run()
    else:
        E.warn( "file %s does not exist: not imported" % i )

############################################################
def _makeOverlap( infiles, outfile, subset = "all" ):
    '''compute overlaps between sets.'''

    tracks = [x[:-len(".gtf.gz")] for x in infiles ]
    to_cluster = True
    
    if os.path.exists( outfile ):
        shutil.move( outfile, outfile + ".old" )
        extra_options = "--update=%(outfile)s.old" % locals()
    else:
        extra_options = ""

    if subset == "all": where = "1" 
    else: where = "is_%s" % subset

    tracks = " ".join(tracks)

    statement = '''
	for d in %(tracks)s; do 
                gunzip < ${d}.gtf.gz
		| python %(scriptsdir)s/gtf2gtf.py 
		    --apply=<( %(cmd-sql)s %(database)s "SELECT gene_id FROM ${d}_annotation WHERE %(where)s" )
		    --log=%(outfile)s.log 
		    --filter=gene 
		> %(outfile)s_tmp_${d}.xgtf ;
	done;
	python %(scriptsdir)s/diff_gtfs.py
		%(extra_options)s 
		--pattern-id='%(outfile)s_tmp_(.*).xgtf'
		--log=%(outfile)s.log 
                %(outfile)s_tmp_*.xgtf 
        > %(outfile)s
    '''
    P.run()
    
    statement = "rm -f %(outfile)s_tmp*" 
    P.run()

@merge( TRACKS, "overlap_all.table")
def makeOverlapAll( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "all" )

@merge( TRACKS, "overlap_unknown.table")
def makeOverlapUnknown( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "unknown" )

@merge( TRACKS, "overlap_ambiguous.table")
def makeOverlapAmbiguous( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "ambiguous" )

@merge( TRACKS, "overlap_known.table")
def makeOverlapKnown( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "known" )

@merge( TRACKS, "overlap_pc.table")
def makeOverlapPC( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "pc" )

@transform( (makeOverlapAll,
             makeOverlapUnknown,
             makeOverlapAmbiguous,
             makeOverlapKnown,
             makeOverlapPC ),
            suffix(".table"),
            "_table.import" )
def importOverlap( infile, outfile ):
    '''import results of overlap computation.'''
    
    tablename = outfile[:-len("_table.import")]
    statement = '''
	grep -v "\\bna\\b" 
        < %(infile)s 
        |python %(scriptsdir)s/csv2db.py %(csv2db_options)s
             --map set1:str 
             --map set2:str 
             --index=set1 
             --index=set2 
             --table=%(tablename)s
        > %(outfile)s
    '''
    P.run()

@transform( TRACKS, 
            suffix(".gtf.gz"), 
            "_merged.overlap")
def makeMergedOverlap( infile, outfile ):
    '''compute overlap of a track with the merged set.'''
    
    pass

########################################################
@transform( TRACKS,
            suffix(".gtf.gz"), 
            "_norepeats.fasta")
def exportRepeatMaskedSequences( infile, outfile ):
    '''output sequences masking repeats.'''

    to_cluster = True
    statement = '''gunzip < %(infile)s 
        | %(scriptsdir)s/gff_sort gene 
	| python %(scriptsdir)s/gff2fasta.py 
		--is-gtf 
		--genome-file=%(genome)s
		--remove-masked-regions 
		--filename-masks=<(gunzip < %(repeats_filename)s )
		--min-length=%(codingpotential_min_length)i 
		--max-length=%(codingpotential_max_length)i 
		--log=%(outfile)s.log 
	> %(outfile)s
    '''
    P.run()

########################################################
@transform( exportRepeatMaskedSequences,
            suffix("_norepeats.fasta"),
            ".coding")
def makeCodingPotential( infile, outfile ):
    '''run CPC to predict coding potential.'''

    statement = '''
	cpc.sh %(infile)s 
               %(outfile)s.forward.table 
               %(outfile)s.tmp.dir 
               %(outfile)s.forward.evidence 
               %(codingpotential_database)s > %(outfile)s.log'''
    P.run()

    tmpfilename = P.getTempFilename( "." )
    statement = '''python %(toolsdir)s/fasta2fasta.py 
                          --method=reverse-complement -v 0
                  < %(infile)s
                  > %(tmpfilename)s'''
    P.run()
    
    statement = '''
	cpc.sh %(tmpfilename)s 
               %(outfile)s.reverse.table 
               %(outfile)s.tmp.dir 
               %(outfile)s.reverse.evidence 
               %(codingpotential_database)s >> %(outfile)s.log'''
    P.run()
    
    outf = open(outfile, "w")
    outf.write( "gene_id\tlength\tf_iscoding\tf_value\tf_orfstart\tf_orfend\tf_orfval1\tf_orfval2\tf_orf\tr_iscoding\tr_value\tr_orfstart\tr_orfend\tr_orfval1\tr_orfval2\tr_orf\n")
    outf.close()

    to_cluster = True
    
    statement = '''
	python %(toolsdir)s/combine_tables.py -v 0 
            %(outfile)s.forward.table 
            %(outfile)s.forward.evidence.orf 
            %(outfile)s.reverse.table 
            %(outfile)s.reverse.evidence.orf |\
	cut -f 1,2,3,4,6,7,8,9,10,12,13,15- 
        >> %(outfile)s
    '''
    P.run()

    # save space by compressing the result of the homology searches
    E.info( "compressing CPC output" )
    statement ='''rm -f %(outfile)s.*homo.gz; gzip %(outfile)s.*homo'''
    P.run()

    os.unlink( tmpfilename )

########################################################
@transform(  makeCodingPotential, suffix(".coding"), "_coding.import")
def importCodingPotential( infile, outfile ):
    '''import annotations'''
    
    table = outfile[:-len(".import")]

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --allow-empty
              --index=gene_id 
              --map=gene_id:str 
              --table=%(table)s 
    < %(infile)s 
    > %(outfile)s'''

    P.run()

    # set the is_coding flag
    dbhandle = sqlite3.connect( PARAMS["database"] )
    Database.executewait( dbhandle, '''ALTER TABLE %(table)s ADD COLUMN is_coding INTEGER''' % locals())
    Database.executewait( dbhandle, '''UPDATE %(table)s SET is_coding = (f_iscoding == 'coding') OR (r_iscoding == 'coding') ''' % locals())
    dbhandle.commit()

########################################################
@transform( TRACKS, 
            suffix(".gtf.gz"), 
            "_introns.gtf.gz")
def exportIntrons( infile, outfile ):
    '''convert the exons to introns.

    10 bp are truncated on either end of an intron.

    Introns are filtered by a minimum length (100).
    '''

    to_cluster = True
    
    statement = '''gunzip
        < %(infile)s 
	| %(scriptsdir)s/gff_sort gene 
	| python %(scriptsdir)s/gtf2gtf.py 
               --exons2introns 
               --intron-min-length=100 
               --intron-border=10 
               --log=%(outfile)s.log
	| python %(scriptsdir)s/gtf2gtf.py 
              --set-transcript-to-gene 
              --log=%(outfile)s.log 
	| perl -p -e 's/intron/exon/'
        | gzip
        > %(outfile)s
    '''
    P.run()


########################################################
@transform( TRACKS + [exportIntrons,], 
            suffix(".gtf.gz"), 
            ".fasta")
def exportSequences( infile, outfile ):
    '''collect sequences from a gtf file.'''
    
    prefix = outfile[:-len(".fasta")]

    to_cluster = True
    statement = '''gunzip 
        < %(infile)s
	| %(scriptsdir)s/gff_sort gene  
	| python %(scriptsdir)s/gff2fasta.py 
		--is-gtf 
		--genome-file=genome 
		--log=%(outfile)s.log 
	| python %(toolsdir)s/index_fasta.py --force %(prefix)s - 
        > %(outfile)s.log'''

    P.run()

############################################################
@transform( (TRACKS, exportIntrons, exportSequences), 
            regex("(.*).fasta"), 
            inputs( (r"\1.gtf.gz", r"\1.fasta", PARAMS["rates_alignment"]) ),
            r"\1.rates.gz")
def makeRates( infiles, outfile ):
    '''compute nucleotide substitution rates for transcripts from a gtf file - 
    this applies only for transcripts mapped onto the reference genome.

    Sequences from the transcripts are mapped onto the rate genome.
    
    Softmasked sequence will be ignored unless track is in CONTROL_TRACKS.

    The longest contiguous block is selected ignoring matches to other parts of the genome.
    '''


    infile_gtf, infile_sequences, alignment = infiles

    track = infile_gtf[:len(".gtf.gz")]
    
    if infile_gtf in CONTROL_TRACKS:
        # when aligning repeats, do not mask lower case characters
        mask = ""
    else:
        mask = "--mask-lowercase"

    statement = '''gunzip 
        < %(infile_gtf)s 
	| %(scriptsdir)s/gff_sort gene 
	| python %(scriptsdir)s/gff2psl.py 
               --is-gtf 
               --genome-file=%(genome)s 
               --log=%(outfile)s.log 
	| pslMap stdin <(gunzip < %(alignment)s ) stdout 
	| sort -k10,10 -k14,14 -k9,9 -k12,12n 
	| python %(scriptsdir)s/psl2psl.py 
                --method=merge 
                --log=%(outfile)s.log 
	| python %(scriptsdir)s/psl2psl.py 
                --method=select-query 
                --select=most-nmatches 
                --log=%(outfile)s.log 
	| python %(scriptsdir)s/psl2psl.py 
                --method=add-sequence 
                --filename-target=%(rates_genome)s
                --filename-queries=%(infile_sequences)s
                --log=%(outfile)s.log 
	| %(cmd-farm)s 
                --split-at-lines=10000 
                --output-header 
                --log=%(outfile)s.log 
	   "python %(scriptsdir)s/psl2table.py 
                %(mask)s
                --method=counts 
                --method=baseml 
                --baseml-model=REV" 
	| gzip 
        > %(outfile)s
    '''
    
    P.run()

############################################################
@transform(  makeRates, suffix(".rates.gz"), "_rates.import")
def importRates( infile, outfile ):
    '''import rates.

    Select the longest stretch for each transcript.
    '''

    track = outfile[:-len(".import")]
    statement = '''
    gunzip 
    < %(infile)s 
    | python %(toolsdir)s/csv_cut.py 
          --large --remove qStarts tStarts blockSizes qSequence tSequence --log=%(outfile)s 
    | csort -k:qName: -k:aligned:rn 
    | perl -p -e "s/qName/gene_id/" 
    | awk '{if (l==$10) {next;} l = $10; print; }' 
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s --map gene_id:str --table=%(track)s --index=gene_id --allow-empty
    > %(outfile)s
    '''

    P.run()

########################################################
@transform( TRACKS,
            regex("(.*).gtf.gz"), 
            r"\1.repeats_rates.gz")
def makeRepeatsRates( infile, outfile ):
    '''collect ancestral repeats around transcripts and
    estimate neutral rate.'''

    to_cluster = True

    statement = '''gunzip 
    < %(infile)s 
    | python %(scriptsdir)s/gtf2table.py 
        --proximal-distance=%(ancestral_repeats_max_distance)i 
        --section=exons 
        --counter=proximity-exclusive
        --filename-format=gff 
        --filename-gff=<(gunzip < %(ancestral_repeats_filename)s )
    | gzip 
    > %(outfile)s'''
    
    P.run()


############################################################
@transform(  makeRepeatsRates, suffix(".repeats_rates.gz"), "_repeats_rates.import")
def importRepeatsRates( infile, outfile ):
    '''import repeat overlap'''
    
    table = outfile[:-len(".import")]

    statement = '''gunzip 
    < %(infile)s 
    | awk '$4 > 0'
    | python %(toolsdir)s/csv_cut.py --remove exons_lengths exons_values
    |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              --index=gene_id 
              --map=gene_id:str 
              --table=%(table)s 
              --allow-empty
    > %(outfile)s'''

    P.run()

############################################################
@transform( TRACKS, suffix(".gtf.gz"), "_evol.import") 
def importSummary( infile, outfile ):
    '''import several rates into a single convenience table.
    '''

    stmt_select = []
    stmt_from = []
    stmt_where = ["1"]

    track = infile[:-len(".gtf.gz")]

    tablename = "%s_evol" % track

    if os.path.exists( "%s_rates.import" % track ):
        stmt_select.append( "a.distance AS ks, a.aligned AS aligned" )
        stmt_from.append('''LEFT JOIN %(track)s_rates AS a
                     ON r.gene_id = a.gene_id AND 
                     a.aligned >= %(rates_min_aligned)i AND 
                     a.distance <= %(rates_max_rate)f''' )

    if os.path.exists( "%s_coverage.import" % track ):
        stmt_select.append("cov.nmatches AS nreads, cov.mean AS meancoverage" )
        stmt_from.append("LEFT JOIN %(track)s_coverage AS cov ON r.gene_id = cov.gene_id" )

    if os.path.exists( "%s_repeats_gc.import" % track ):
        stmt_select.append("ar_gc.exons_mean AS repeats_gc" )
        stmt_from.append("LEFT JOIN %(track)s_repeats_gc AS ar_gc ON r.gene_id = ar_gc.gene_id" )

    if os.path.exists( "%s_repeats_rates.import" % track ):
        stmt_select.append("ar.exons_length AS ar_aligned, ar.exons_median AS ka, a.distance/ar.exons_median AS kska" )
        stmt_from.append('''LEFT JOIN %(track)s_repeats_rates AS ar 
                     ON r.gene_id = ar.gene_id AND 
                     ar.exons_nval >= %(rates_min_repeats)i''' )

    if os.path.exists( "%s_introns_rates.import" % track ):
        stmt_select.append("ir.aligned AS ir_aligned, ir.distance AS ki, a.distance/ir.distance AS kski" )
        stmt_from.append('''LEFT JOIN %(track)s_introns_rates AS ir 
                            ON r.gene_id = ir.gene_id AND 
                            ir.aligned >= %(rates_min_aligned)i''' )

    x = locals()
    x.update( PARAMS )
    stmt_select = ", ".join( stmt_select ) % x 
    stmt_from = " ".join( stmt_from ) % x
    stmt_where = " AND ".join( stmt_where ) % x

    dbhandle = sqlite3.connect( PARAMS["database"] )

    Database.executewait( dbhandle, "DROP TABLE IF EXISTS %(tablename)s " % locals() )

    statement = '''
    CREATE TABLE %(tablename)s AS
    SELECT
         CAST(r.gene_id AS TEXT) AS gene_id, 
				r.exons_sum as length, 
				r.exons_pGC as pgc,
                                %(stmt_select)s
	FROM 
          %(track)s_annotation AS r 
	  %(stmt_from)s
        WHERE %(stmt_where)s
    ''' % locals()

    Database.executewait( dbhandle, statement)
    dbhandle.commit()
    P.touch(outfile)

############################################################
@files( PARAMS["genome"] + ".fasta", "genome_gc.bed" )
def buildGenomeGCSegmentation( infile, outfile ):
    PAnnotator.buildGenomeGCSegmentation( infile, outfile )

############################################################
@files( buildGenomeGCSegmentation, "genome_gc_regions.bed" )
def buildAnnotatorGC( infile, outfile ):
    '''compute G+C regions.'''
    PAnnotator.buildAnnotatorGC( infile, outfile )

############################################################
@files( buildAnnotatorGC, PARAMS["annotator_gc_workspace"] )
def buildAnnotatorGCWorkspace( infile, outfile ):
    '''compute G+C regions.'''
    PAnnotator.buildAnnotatorGCWorkspace( infile, outfile )

############################################################
@files( [( TRACKS, "%s.annotations" % x, x) for x in ("known", "unknown", "all", "intronic", "intergenic" ) ] )
def buildAnnotatorGeneSetAnnotations( infiles, outfile, slice ):
    PAnnotator.buildGeneSetAnnotations( infiles, outfile, slice )

############################################################
def makeAnnotatorArchitecture( TRACKS ):
    '''compute annotator overlap with architecture.'''
    pass

############################################################
@files( [ (track, "%s.%s.sets.annotator" % (track[:-len(".gtf.gz")],slice), slice) \
             for track, slice in list( itertools.product( EXPERIMENTAL_TRACKS + DERIVED_TRACKS, 
                                                          ("known", "unknown", "all", "intronic", "intergenic" ))) ] )
def makeAnnotatorGeneSets( infile, outfile, slice ):
    '''compute annotator overlap between sets.
    '''
    
    workspaces = ("genomic", "alignable", slice )

    track = infile[:-len(".gtf.gz")]

    infiles = ANNOTATOR_TRACKS

    related = getRelatedTracks( infile, infiles )

    if related:
        E.info("removing related tracks %s from %s" % \
                   ( related, infile ) )
        related = set(related)
        infiles = [x for x in TRACKS if x not in related ]
        
    tmpdir = tempfile.mkdtemp( dir = os.getcwd() )

    annotations = os.path.join( tmpdir, "annotations")
    PAnnotator.buildGeneSetAnnotations( infiles,
                                        annotations,
                                        slice )

    segments = PAnnotator.buildAnnotatorSlicedSegments( tmpdir, 
                                                        outfile, 
                                                        track, 
                                                        slice )

    if not segments:
        E.warn( "no segments for %s - no annotator results" % outfile )
        shutil.rmtree( tmpdir )
        P.touch( outfile )
        return

    workspaces, synonyms = PAnnotator.buildAnnotatorWorkSpace( tmpdir, 
                                                               outfile,
                                                               workspaces = workspaces,
                                                               gc_control = True )
    
    PAnnotator.runAnnotator( tmpdir, 
                             outfile, 
                             annotations, 
                             segments, 
                             workspaces, 
                             synonyms )

    shutil.rmtree( tmpdir )

@merge( makeAnnotatorGeneSets, "annotator_genesets.import" )
def importAnnotatorGeneSets( infiles, outfile ):
    '''import genesets.'''

    PAnnotator.importAnnotator( infiles, 
                                outfile, 
                                regex_id = "(.*).sets.annotator",
                                table = "annotator_sets", 
                                fdr_method = "annotator",
                                with_slice = True )


############################################################
@follows( buildRepeatTrack,
          buildIntronTrack )
def setup():
    pass

@follows( importTranscripts,
          importTranscriptInformation,
          importGeneStats,
          importGeneInformation,
          importGenomeInformation,
          importRepeatInformation,
          buildGenes,
          buildGeneRegions,
          buildMergedTracks,
          buildFilteredAlignment,
          importGTF)
def prepare():
    pass

@follows( importDifference,
          importOverlap
          )
def difference(): pass

@follows( importAnnotations,
          importOverrun,
          importDistances,
          importRepeats,
          importSegments,
          importCodingPotential,
          )
def annotate():
    pass

@follows( importRates, importRepeatsRates )
def rates():
    pass

@follows( buildAnnotatorGCWorkspace, 
          importAnnotatorGeneSets)
def annotator():
    pass

@follows( annotate, rates, difference )
def full(): pass

@follows( annotate, rates, importSummary )
def summary(): pass

if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

    

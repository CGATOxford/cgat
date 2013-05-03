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
===============================
Transcriptome analysis pipeline
===============================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

The transcriptome analysis pipeline examines various gene/transcript sets
and attempts to interprete and compare them. Among the quantities that it 
computes are:

   * coding potential (using CPC)
   * overlap with repeats
   * substitution rates
   * conservation 
      * with respect to ancestral repeat rates)
      * TODO: conservation (with respect to UCSC tracks)
   * annotation with respect to a reference gene set, for example
      * location (intronic, exonic, intergenic, ...)
      * distance to closest gene
      * intron overrun


Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

The sphinxreport report requires a :file:`conf.py` and :file:`sphinxreport.ini` file 
(see :ref:`PipelineDocumenation`). To start with, use the files supplied with the
:ref:`Example` data.

Input
-----

The pipeline works from gene sets in :term:`gtf` format within the :term:`working directory`.

Requirements
------------

The pipeline requires the results from :doc:`pipeline_annotations` and :doc:`pipeline_ancestral_repeats`.

On top of the default CGAT setup, the pipeline requires the following software to be in the 
path:

+--------------------+-------------------+------------------------------------------------+
|*Program*           |*Version*          |*Purpose*                                       |
+--------------------+-------------------+------------------------------------------------+
|CPC                 |>=0.9-r2           |prediction of coding potential                  |
+--------------------+-------------------+------------------------------------------------+

Example
=======

Example data is available at http://www.cgat.org/~andreas/sample_data/pipeline_transcriptome.tgz.
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_transcriptome.tgz
   tar -xvzf pipeline_transcriptome.tgz
   cd pipeline_transcriptome
   python <srcdir>/pipeline_transcriptome.py make full

.. note:: 
   For the pipeline to run, install the :doc:`pipeline_annotations` and :doc:`pipeline_ancestral_repeats`
   as well.

Glossary
========

.. glossary::

Code
====

"""
import sys, tempfile, optparse, shutil, itertools, csv, math, random, re, glob, os, shutil, collections
import csv, gzip
from ruffus import *
import sqlite3

import Experiment as E
import Pipeline as P
import IndexedFasta, IndexedGenome, FastaIterator, Genomics
import IOTools, Database, GFF, GTF
import Database

import PipelineGeneset as PGeneset
import PipelineAnnotator as PAnnotator

# load options from the config file
import Pipeline as P
P.getParameters( 
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini" ] )

PARAMS = P.PARAMS

USECLUSTER = True

## link up with annotations
PARAMS_ANNOTATIONS = P.peekParameters( PARAMS["annotations_dir"],
                                       "pipeline_annotations.py" )

## link up with ancestral repeats
PARAMS_ANCESTRAL_REPEATS = P.peekParameters( PARAMS["ancestral_repeats_dir"],
                                            "pipeline_ancestral_repeats.py" )

###################################################################
###################################################################
## Helper functions mapping tracks to conditions, etc
###################################################################
import PipelineTracks

# collect sra nd fastq.gz tracks
TRACKS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    glob.glob( "*.gtf.gz" ), "(\S+).gtf.gz", exclude=("repeats.gtf.gz", "introns.gtf.gz", "merged.gtf.gz") )

TRACKS_CONTROL = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    ( "repeats.gtf.gz", "introns.gtf.gz"), "(\S+).gtf.gz" )

TRACKS_META = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    ( "merged.gtf.gz",), "(\S+).gtf.gz" )

TRACKS_GENESETS = PipelineTracks.Tracks( PipelineTracks.Sample ).loadFromDirectory( 
    ( "genes.gtf.gz",), "(\S+).gtf.gz" )

# collection of all tracks including controls
TRACKS_WITH_CONTROLS = TRACKS + TRACKS_CONTROL

# # TRACKS to work on
# REFERENCE_TRACKS= [ "guttmanExons.gtf.gz", 
#                     "guttmanK4K36.gtf.gz",
#                     "khalilExons.gtf.gz",
#                     "khalilK4K36.gtf.gz",
#                     "riken.gtf.gz",
#                     "rikenFiltered.gtf.gz",
#                     "lincHinv.gtf.gz",
#                     "lincOrang.gtf.gz",
#                     "humanRefseq.gtf.gz" ]
# GENESET_TRACKS = [ "genes.gtf.gz" ]
# EXPERIMENTAL_TRACKS = ["kamilah.gtf.gz" ]
# DERIVED_TRACKS = [ x for x in ("kamilahLinc.gtf.gz", ) if os.path.exists(x) ]
# OVERLAP_TRACKS = GENESET_TRACKS + META_TRACKS + ["humanRefseq.gtf.gz",  "lincPrimates.gtf.gz" ]
# TRACKS = CONTROL_TRACKS + GENESET_TRACKS + EXPERIMENTAL_TRACKS + DERIVED_TRACKS + META_TRACKS + REFERENCE_TRACKS
# ANNOTATOR_TRACKS = TRACKS + DERIVED_TRACKS

TRACKS_OVERLAP = TRACKS_META + TRACKS_GENESETS

###################################################################
###################################################################
###################################################################
##
###################################################################
if os.path.exists("pipeline_conf.py"): 
    L.info( "reading additional configuration from pipeline_conf.py" )
    execfile("pipeline_conf.py")

###################################################################
###################################################################
###################################################################
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
@files( ( (os.path.join( PARAMS["ancestral_repeats_dir"],
                         PARAMS_ANCESTRAL_REPEATS["interface_rates_query_gff"]), 
           "repeats.gtf.gz"),) )
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
@files( ( (os.path.join( PARAMS["ancestral_repeats_dir"],
                         PARAMS_ANCESTRAL_REPEATS["interface_rates_query_gff"]), 
           "introns.gtf.gz"),) )
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
@files( ( ( (PARAMS["repeats_filename"], PARAMS["genome"] + ".idx" ),
            "repeats_table.load" ), ) )
def loadRepeatInformation( infiles, outfile ):
    '''load genome information.'''
    
    to_cluster = True

    table = outfile[:-len(".load")]

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
@merge( TRACKS.getTracks( "%s.gtf.gz"), "merged.gtf.gz" )
def buildMergedTracks( infiles, outfile ):
    '''merge gtf transcripts across all data sets.'''

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
		--genome=%(genome_dir)s/%(genome)s
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


############################################################
############################################################
############################################################
@files( ( os.path.join( PARAMS["ancestral_repeats_dir"],
                        PARAMS_ANCESTRAL_REPEATS["interface_alignment_psl"]),
          os.path.join( PARAMS["annotations_dir"],
                       PARAMS_ANNOTATIONS["interface_geneset_flat_gtf"])), 
        "alignment_filtered.psl.gz" )
def buildFilteredAlignment( infiles, outfile ):
    '''build a genomic alignment without exons from the reference gene set.

    This alignment is used for rate computation within introns.
    '''
    infile_alignment, infile_genes = infiles
    to_cluster = USECLUSTER

    statement = '''gunzip 
	< %(infile_alignment)s 
	| python %(scriptsdir)s/psl2psl.py 
		--method=filter-remove 
		--filter-query=<(gunzip < %(infile_genes)s)
		--log=%(outfile)s.log 
        | gzip
        > %(outfile)s '''

    P.run()

############################################################
############################################################
############################################################
@transform(buildFilteredAlignment, 
           suffix(".psl.gz"), 
           ".stats" )
def buildAlignmentStats( infile, outfile ):
    '''build alignment statistics.'''

    to_cluster = USECLUSTER
    
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
@transform( TRACKS.getTracks( "%s.gtf.gz" ),
            suffix(".gtf.gz"), 
            "_gtf.load" )
def loadGTF( infile, outfile ):
    '''load gtf files.'''
    
    table = P.toTable( outfile )

    to_cluster = USECLUSTER

    statement = '''gunzip
        < %(infile)s
        | python %(scriptsdir)s/gtf2tab.py 
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
@transform( TRACKS_WITH_CONTROLS.getTracks( "%s.gtf.gz" ),
            suffix(".gtf.gz"), 
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                      PARAMS_ANNOTATIONS["interface_annotation_gff"])), 
            ".annotation.gz")
def buildAnnotations( infiles, outfile ):
    '''annotate transcripts by location (intergenic, intronic, ...)'''
    
    infile, annotation = infiles

    statement = '''gunzip 
    < %(infile)s 
    | python %(scriptsdir)s/gtf2gtf.py --sort=gene
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
		--genome-file=%(genome_dir)s/%(genome)s"
    | gzip
    > %(outfile)s
    '''
    P.run()

############################################################
############################################################
############################################################
@transform( buildAnnotations, suffix(".annotation.gz"), "_annotation.load")
def loadAnnotations( infile, outfile ):
    '''load annotations'''
    P.load( infile, outfile, "--index=gene_id --map=gene_id:str")
    
########################################################
########################################################
########################################################
@transform( TRACKS_WITH_CONTROLS.getTracks( "%s.gtf.gz" ),
            suffix(".gtf.gz"), 
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_annotation_gff"])), 
            ".overrun.gz")
def makeOverrun( infiles, outfile ):
    '''compute intron overrun.'''

    infile, annotation = infiles

    statement = '''gunzip 
    < %(infile)s 
    | python %(scriptsdir)s/gtf2gtf.py --sort=gene
    | %(cmd-farm)s --split-at-column=1 --output-header --log=%(outfile)s.log --max-files=60 
	"python %(scriptsdir)s/gtf2table.py 
                --filename-format=gff
		--counter=overrun 
		--log=%(outfile)s.log 
                --filename-format=gff 
		--filename-gff=%(annotation)s" 
    | gzip
    > %(outfile)s 
    '''
    P.run()

############################################################
############################################################
############################################################
@transform( makeOverrun, suffix(".overrun.gz"), "_overrun.load")
def loadOverrun( infile, outfile ):
    '''load annotations'''
    P.load( infile, outfile, "--index=gene_id --map=gene_id:str")    

########################################################
########################################################
########################################################
@transform( TRACKS_WITH_CONTROLS.getTracks( "%s.gtf.gz" ),
            suffix(".gtf.gz"), 
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_geneset_flat_gtf"])), 
            ".distances")
def makeDistances( infiles, outfile ):
    '''compute intron overrun.'''

    infile, annotation = infiles

    statement = '''gunzip
    < %(infile)s 
    | python %(scriptsdir)s/gtf2gtf.py --sort=gene
    | %(cmd-farm)s --split-at-column=1 --output-header --log=%(outfile)s.log --max-files=60 
	"python %(scriptsdir)s/gtf2table.py 
		--counter=distance-genes 
		--log=%(outfile)s.log 
		--filename-gff=<( gunzip < %(annotation)s ) " 
    > %(outfile)s 
    '''
    P.run()

########################################################
########################################################
############################################################
@transform(  makeDistances, suffix(".distances"), "_distances.load")
def loadDistances( infile, outfile ):
    '''load annotations'''
    P.load( infile, outfile, "--index=gene_id --map=gene_id:str --index=closest_id --map=closest_id:str")    
    table = outfile[:-len(".load")]

########################################################
########################################################
########################################################
@transform( TRACKS_WITH_CONTROLS.getTracks( "%s.gtf.gz" ),
            suffix(".gtf.gz"), 
            ".segments")
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
    | python %(scriptsdir)s/gtf2gtf.py --sort=position+gene
    | python %(scriptsdir)s/gtf2gtf.py --merge-transcripts 
    | python %(scriptsdir)s/gtf2gtf.py --sort=gene
    | python %(scriptsdir)s/gff2histogram.py 
		--method=values 
		--force 
		--output-filename-pattern="%(outfile)s_genes.%%s" 
		--log=%(outfile)s.log
    >> %(outfile)s'''
    P.run()

############################################################
@transform(  makeSegments, suffix(".segments"), "_segments.load")
def loadSegments( infile, outfile ):
    '''load segments'''
    
    table = outfile[:-len(".load")]
    
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
@transform( TRACKS_WITH_CONTROLS.getTracks( "%s.gtf.gz" ),
            suffix(".gtf.gz"), 
            add_inputs( os.path.join( PARAMS["annotations_dir"],
                                      PARAMS_ANNOTATIONS["interface_repeats_gff"])), 
            ".repeats.gz")
def makeRepeats( infiles, outfile ):
    '''compute overlap with repeats
    '''

    infile, annotation = infiles

    statement = '''gunzip
    < %(infile)s 
    | python %(scriptsdir)s/gtf2gtf.py --sort=gene
    | %(cmd-farm)s --split-at-column=1 --output-header --log=%(outfile)s.log --max-files=60 
	"python %(scriptsdir)s/gtf2table.py 
                --filename-format=gff
		--counter=overlap 
		--log=%(outfile)s.log 
		--filename-gff=<( gunzip < %(annotation)s )" 
    | gzip
    > %(outfile)s 
    '''
    P.run()

############################################################
@transform(  makeRepeats, suffix(".repeats.gz"), "_repeats.load")
def loadRepeats( infile, outfile ):
    '''load repeat overlap'''
    P.load( infile, outfile, "--index=gene_id --map=gene_id:str")

############################################################
@files( [ ( 
            ("%s.gtf.gz" % x, "%s.gtf.gz" % y),
            "%s_vs_%s.diff" % (x,y)) for x,y in itertools.product( TRACKS_OVERLAP, TRACKS) ] )
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
@transform(  makeDifference, suffix(".diff"), "_diff.load")
def loadDifference( infile, outfile ):
    '''load overlap'''

    table = outfile[:-len("_diff.load")]

    statement = '''
        cat < %(infile)s 
	| python %(toolsdir)s/csv_cut.py --large --remove cov_values 
    	| grep -v "\\bna\\b" 
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
        E.warn( "file %s does not exist: not loaded" % i )

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
        E.warn( "file %s does not exist: not loaded" % i )

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

@merge( TRACKS.getTracks( "%s.gtf.gz" ), "overlap_all.table")
def makeOverlapAll( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "all" )

@merge( TRACKS.getTracks( "%s.gtf.gz" ), "overlap_unknown.table")
def makeOverlapUnknown( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "unknown" )

@merge( TRACKS.getTracks( "%s.gtf.gz" ), "overlap_ambiguous.table")
def makeOverlapAmbiguous( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "ambiguous" )

@merge( TRACKS.getTracks( "%s.gtf.gz" ), "overlap_known.table")
def makeOverlapKnown( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "known" )

@merge( TRACKS.getTracks( "%s.gtf.gz" ), "overlap_pc.table")
def makeOverlapPC( infiles, outfile ):
    '''compute overlaps between sets.'''
    _makeOverlap( infiles, outfile, subset = "pc" )

@transform( (makeOverlapAll,
             makeOverlapUnknown,
             makeOverlapAmbiguous,
             makeOverlapKnown,
             makeOverlapPC ),
            suffix(".table"),
            "_table.load" )
def loadOverlap( infile, outfile ):
    '''load results of overlap computation.'''
    
    tablename = outfile[:-len("_table.load")]
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

@transform( TRACKS.getTracks( "%s.gtf.gz" ), 
            suffix(".gtf.gz"), 
            "_merged.overlap")
def makeMergedOverlap( infile, outfile ):
    '''compute overlap of a track with the merged set.'''
    
    pass

########################################################
@transform( TRACKS.getTracks( "%s.gtf.gz" ),
            suffix(".gtf.gz"), 
            "_norepeats.fasta")
def buildRepeatMaskedSequences( infile, outfile ):
    '''output sequences masking repeats.

    .. todo: assess which repeats to use (+low complexity?)
    '''

    to_cluster = True

    repeats = os.path.join( PARAMS["annotations_dir"],
                            PARAMS_ANNOTATIONS["interface_repeats_gff"])

    statement = '''gunzip < %(infile)s 
        | python %(scriptsdir)s/gtf2gtf.py --sort=gene
	| python %(scriptsdir)s/gff2fasta.py 
		--is-gtf 
		--genome-file=%(genome_dir)s/%(genome)s
		--remove-masked-regions 
		--filename-masks=<(gunzip < %(repeats)s )
		--min-length=%(codingpotential_min_length)i 
		--max-length=%(codingpotential_max_length)i 
		--log=%(outfile)s.log 
	> %(outfile)s
    '''
    P.run()

########################################################
########################################################
########################################################
@transform( buildRepeatMaskedSequences,
            suffix("_norepeats.fasta"),
            ".blast.gz")
def runAgainstProteinDatabase( infile, outfile ):
    '''run blastx.
    
    search on both strands. Note that the CPC default is: only forward strand
    '''

    statement = '''
    cat %(infile)s | 
       %(cmd-farm)s 
                --split-at-regex="^>(\S+)"
                --chunksize=100
                --output-header 
                --log=%(outfile)s.log 
    blastx -strand both -evalue 1e-10 -ungapped -threshold 14 -db /ifs/apps/bio/cpc-0.9-r2/bin/../share/prot_db
    | gzip 
    > %(outfile)s
    '''
    
    P.run()

########################################################
########################################################
########################################################
@transform( buildRepeatMaskedSequences,
            suffix("_norepeats.fasta"),
            ".frame.gz")
def runFrameFinder( infile, outfile ):
    '''run FrameFinder

    search on both strands (-r TRUE). Note that CPC default is: only forward strand.

    '''
    cpc_dir = "/ifs/apps/bio/cpc-0.9-r2"
    statement = '''
    cat %(infile)s |
    %(cpc_dir)s/libs/estate/bin/framefinder
    -r TRUE -w %(cpc_dir)s/data/framefinder.model /dev/stdin
    | gzip
     > %(outfile)s
    '''
    
    P.run()

########################################################
########################################################
########################################################
@follows( runFrameFinder )
@transform( runAgainstProteinDatabase,
            suffix(".blast.gz"),
            ".coding.gz")
def buildCodingPotential( infile, outfile ):
    '''run CPC analysis as in the cpc script.

    This module runs framefinder and blastx on both strands.
    It seems to work, but I have not thoroughly tested it.
    I expect that the false positive rate increases (i.e.,
    predicting non-coding as coding) in cases where the best 
    framefinder match and the best blast match are on opposite
    strands. In the original CPC, these would be separated.
    '''

    try:
        cpc_dir = os.environ["CPC_HOME"]
    except KeyError:
        raise ValueError("CPC_HOME environment variable is not set. ")

    tmpdir = P.getTempDir( ".")
    track = P.snip( outfile, ".coding.gz" )

    # extract features for frame finder
    # replaces extract_framefinder_feats.pl to parse both strands
    with open( os.path.join(tmpdir, "ff.feat"), "w") as outf:
        outf.write( "\t".join(("QueryID", "CDSLength", "Score", "Used", "Strict")) + "\n")
        for line in IOTools.openFile( "%s.frame.gz" % track ):
            if line.startswith(">"):
                try:
                    ( id, start, end, score, used, mode, tpe) = \
                        re.match(
                        ">(\S+).*framefinder \((\d+),(\d+)\) score=(\S+) used=(\S+)% \{(\S+),(\w+)\}", line ).groups()
                except AttributeError:
                    raise ValueError( "parsing error in line %s" % line )
                length = int(end) - int(start) + 1
                strict = int(tpe == "strict")
                outf.write( "\t".join( (id, str(length), used, str(strict )) )+ "\n")

    to_cluster = USECLUSTER

    # extract features and prepare svm data
    s = []
            
    s.append( '''
    zcat %(infile)s
    | perl %(cpc_dir)s/libs/blast2table.pl 
    | tee %(tmpdir)s/blastx.table
    | perl %(cpc_dir)s/bin/extract_blastx_features.pl
    > %(tmpdir)s/blastx.feat1;
    ''' )

    s.append( '''
    cat %(track)s_norepeats.fasta 
    | perl %(cpc_dir)s/bin/add_missing_entries.pl
       %(tmpdir)s/blastx.feat1 
    > %(tmpdir)s/blastx.feat;
    ''')

    # step 2 - prepare data
    s.append( '''
    perl %(cpc_dir)s/bin/feat2libsvm.pl -c 2,4,6 NA NA %(tmpdir)s/blastx.feat 
    > %(tmpdir)s/blastx.lsv;
    ''' )
    
    s.append( '''
    perl %(cpc_dir)s/bin/feat2libsvm.pl -c 2,3,4,5 NA NA %(tmpdir)s/ff.feat 
    > %(tmpdir)s/ff.lsv;
    ''' )

    s.append( '''
    perl -w %(cpc_dir)s/bin/lsv_cbind.pl %(tmpdir)s/blastx.lsv %(tmpdir)s/ff.lsv 
    > %(tmpdir)s/test.lsv;
    ''' )

    s.append( '''
    %(cpc_dir)s/libs/libsvm/libsvm-2.81/svm-scale 
               -r %(cpc_dir)s/data/libsvm.range  
               %(tmpdir)s/test.lsv 
    > %(tmpdir)s/test.lsv.scaled;
    ''' )
    
    # step 3: prediction
    m_libsvm_model0=os.path.join( cpc_dir, "data/libsvm.model0") # standard
    m_libsvm_model=os.path.join( cpc_dir, "data/libsvm.model") # Prob
    m_libsvm_model2=os.path.join( cpc_dir, "data/libsvm.model2" )	# Prob + weighted version
    m_libsvm_range=os.path.join( cpc_dir, "data/libsvm.range" )

    s.append( '''
               %(cpc_dir)s/libs/libsvm/libsvm-2.81/svm-predict2
               %(tmpdir)s/test.lsv.scaled 
               %(m_libsvm_model0)s 
               %(tmpdir)s/test.svm0.predict 
    > %(tmpdir)s/test.svm0.stdout 2> %(tmpdir)s/test.svm0.stderr;
    ''' )

    s.append( '''
    printf "gene_id\\tlength\\tresult\\tvalue\\n" 
    | gzip > %(outfile)s;
    cat %(tmpdir)s/test.svm0.predict  
    | perl -w %(cpc_dir)s/bin/predict.pl %(track)s_norepeats.fasta 
    | gzip >> %(outfile)s;
    ''' )

    # generate reports
    s.append( '''cat %(tmpdir)s/blastx.feat
    | perl -w %(cpc_dir)s/bin/generate_plot_features.pl %(tmpdir)s/blastx.table <( zcat %(track)s.frame.gz) 
    | perl -w %(cpc_dir)s/bin/split_plot_features_by_type.pl %(outfile)s.homology %(outfile)s.orf;
    gzip %(outfile)s.orf %(outfile)s.homology;
    ''' )

    # now run it all
    statement = " checkpoint; ".join( s )
    P.run()

    # clean up
    shutil.rmtree( tmpdir )

########################################################
@transform(  buildCodingPotential, suffix(".coding.gz"), "_coding.load")
def loadCodingPotential( infile, outfile ):
    '''load annotations'''
    
    table = P.toTable( outfile )

    statement = '''
    gunzip < %(infile)s 
    | python %(scriptsdir)s/csv2db.py 
              %(csv2db_options)s 
              --allow-empty
              --index=gene_id 
              --map=gene_id:str 
              --table=%(table)s 
    > %(outfile)s'''

    P.run()

    # set the is_coding flag
    dbhandle = sqlite3.connect( PARAMS["database"] )
    Database.executewait( dbhandle, '''ALTER TABLE %(table)s ADD COLUMN is_coding INTEGER''' % locals())
    Database.executewait( dbhandle, '''UPDATE %(table)s SET is_coding = (result == 'coding')''' % locals())
    dbhandle.commit()

########################################################
@transform( TRACKS_WITH_CONTROLS.getTracks( "%s.gtf.gz" ),
            suffix(".gtf.gz"), 
            ".fasta")
def exportSequences( infile, outfile ):
    '''collect sequences from a gtf file.'''
    
    prefix = outfile[:-len(".fasta")]

    to_cluster = True
    statement = '''gunzip 
        < %(infile)s
        | python %(scriptsdir)s/gtf2gtf.py --sort=gene
	| python %(scriptsdir)s/gff2fasta.py 
		--is-gtf 
		--genome-file=%(genome_dir)s/%(genome)s
		--log=%(outfile)s.log 
	| python %(toolsdir)s/index_fasta.py --force %(prefix)s - 
        > %(outfile)s.log'''

    P.run()

############################################################
@transform( exportSequences, 
            regex("(.*).fasta"), 
            add_inputs( r"\1.gtf.gz", 
                        os.path.join( PARAMS["ancestral_repeats_dir"],
                                      PARAMS_ANCESTRAL_REPEATS["interface_alignment_psl"] )),
            r"\1.rates.gz")
def makeRates( infiles, outfile ):
    '''compute nucleotide substitution rates for transcripts from a gtf file - 
    this applies only for transcripts mapped onto the reference genome.

    Sequences from the transcripts are mapped onto the rate genome.
    
    Softmasked sequence will be ignored unless track is in CONTROL_TRACKS.

    The longest contiguous block is selected ignoring matches to other parts of the genome.
    '''


    infile_sequences, infile_gtf, alignment = infiles

    track = P.snip( infile_sequences, ".fasta" )
    
    if track in TRACKS_CONTROL:
        # when aligning repeats, do not mask lower case characters
        mask = ""
    else:
        mask = "--mask-lowercase"

    # locate target genome from ancestral repeats ini file
    target_genome = os.path.join( PARAMS["genome_dir"],
                                  PARAMS_ANCESTRAL_REPEATS["target"] )

    statement = '''gunzip 
        < %(infile_gtf)s 
        | python %(scriptsdir)s/gtf2gtf.py --sort=gene
	| python %(scriptsdir)s/gff2psl.py 
               --is-gtf 
               --genome-file=%(genome_dir)s/%(genome)s 
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
                --filename-target=%(target_genome)s
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
@transform(  makeRates, suffix(".rates.gz"), "_rates.load")
def loadRates( infile, outfile ):
    '''load rates.

    Select the longest stretch for each transcript.
    '''

    track = outfile[:-len(".load")]
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
@transform( TRACKS_WITH_CONTROLS.getTracks( "%s.gtf.gz" ),
            regex("(.*).gtf.gz"), 
            r"\1.repeats_rates.gz")
def makeRepeatsRates( infile, outfile ):
    '''collect ancestral repeats around transcripts and
    estimate neutral rate.'''

    to_cluster = True

    ancestral_repeats_filename = os.path.join( PARAMS["ancestral_repeats_dir"],
                                               PARAMS_ANCESTRAL_REPEATS["interface_rates_query_gff"])

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
@transform(  makeRepeatsRates, suffix(".repeats_rates.gz"), "_repeats_rates.load")
def loadRepeatsRates( infile, outfile ):
    '''load repeat overlap'''
    
    table = outfile[:-len(".load")]

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
@transform( TRACKS_WITH_CONTROLS.getTracks( "%s.gtf.gz" ),
            suffix(".gtf.gz"), 
            "_evol.load") 
def loadSummary( infile, outfile ):
    '''load several rates into a single convenience table.
    '''

    stmt_select = []
    stmt_from = []
    stmt_where = ["1"]

    track = infile[:-len(".gtf.gz")]

    tablename = "%s_evol" % track

    if os.path.exists( "%s_rates.load" % track ):
        stmt_select.append( "a.distance AS ks, a.aligned AS aligned" )
        stmt_from.append('''LEFT JOIN %(track)s_rates AS a
                     ON r.gene_id = a.gene_id AND 
                     a.aligned >= %(rates_min_aligned)i AND 
                     a.distance <= %(rates_max_rate)f''' )

    if os.path.exists( "%s_coverage.load" % track ):
        stmt_select.append("cov.nmatches AS nreads, cov.mean AS meancoverage" )
        stmt_from.append("LEFT JOIN %(track)s_coverage AS cov ON r.gene_id = cov.gene_id" )

    if os.path.exists( "%s_repeats_gc.load" % track ):
        stmt_select.append("ar_gc.exons_mean AS repeats_gc" )
        stmt_from.append("LEFT JOIN %(track)s_repeats_gc AS ar_gc ON r.gene_id = ar_gc.gene_id" )

    if os.path.exists( "%s_repeats_rates.load" % track ):
        stmt_select.append("ar.exons_length AS ar_aligned, ar.exons_median AS ka, a.distance/ar.exons_median AS kska" )
        stmt_from.append('''LEFT JOIN %(track)s_repeats_rates AS ar 
                     ON r.gene_id = ar.gene_id AND 
                     ar.exons_nval >= %(rates_min_repeats)i''' )

    if os.path.exists( "%s_introns_rates.load" % track ):
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

# ############################################################
# @files( PARAMS["genome"] + ".fasta", "genome_gc.bed" )
# def buildGenomeGCSegmentation( infile, outfile ):
#     PAnnotator.buildGenomeGCSegmentation( infile, outfile )

# ############################################################
# @files( buildGenomeGCSegmentation, "genome_gc_regions.bed" )
# def buildAnnotatorGC( infile, outfile ):
#     '''compute G+C regions.'''
#     PAnnotator.buildAnnotatorGC( infile, outfile )

# ############################################################
# @files( buildAnnotatorGC, PARAMS["annotator_gc_workspace"] )
# def buildAnnotatorGCWorkspace( infile, outfile ):
#     '''compute G+C regions.'''
#     PAnnotator.buildAnnotatorGCWorkspace( infile, outfile )

# ############################################################
# @files( [( TRACKS, "%s.annotations" % x, x) for x in ("known", "unknown", "all", "intronic", "intergenic" ) ] )
# def buildAnnotatorGeneSetAnnotations( infiles, outfile, slice ):
#     PAnnotator.buildGeneSetAnnotations( infiles, outfile, slice )

# ############################################################
# def makeAnnotatorArchitecture( TRACKS ):
#     '''compute annotator overlap with architecture.'''
#     pass

# ############################################################
# @files( [ (track, "%s.%s.sets.annotator" % (track[:-len(".gtf.gz")],slice), slice) \
#              for track, slice in list( itertools.product( EXPERIMENTAL_TRACKS + DERIVED_TRACKS, 
#                                                           ("known", "unknown", "all", "intronic", "intergenic" ))) ] )
# def makeAnnotatorGeneSets( infile, outfile, slice ):
#     '''compute annotator overlap between sets.
#     '''
    
#     workspaces = ("genomic", "alignable", slice )

#     track = infile[:-len(".gtf.gz")]

#     infiles = ANNOTATOR_TRACKS

#     related = getRelatedTracks( infile, infiles )

#     if related:
#         E.info("removing related tracks %s from %s" % \
#                    ( related, infile ) )
#         related = set(related)
#         infiles = [x for x in TRACKS if x not in related ]
        
#     tmpdir = tempfile.mkdtemp( dir = os.getcwd() )

#     annotations = os.path.join( tmpdir, "annotations")
#     PAnnotator.buildGeneSetAnnotations( infiles,
#                                         annotations,
#                                         slice )

#     segments = PAnnotator.buildAnnotatorSlicedSegments( tmpdir, 
#                                                         outfile, 
#                                                         track, 
#                                                         slice )

#     if not segments:
#         E.warn( "no segments for %s - no annotator results" % outfile )
#         shutil.rmtree( tmpdir )
#         P.touch( outfile )
#         return

#     workspaces, synonyms = PAnnotator.buildAnnotatorWorkSpace( tmpdir, 
#                                                                outfile,
#                                                                workspaces = workspaces,
#                                                                gc_control = True )
    
#     PAnnotator.runAnnotator( tmpdir, 
#                              outfile, 
#                              annotations, 
#                              segments, 
#                              workspaces, 
#                              synonyms )

#     shutil.rmtree( tmpdir )

# @merge( makeAnnotatorGeneSets, "annotator_genesets.load" )
# def loadAnnotatorGeneSets( infiles, outfile ):
#     '''load genesets.'''

#     PAnnotator.loadAnnotator( infiles, 
#                                 outfile, 
#                                 regex_id = "(.*).sets.annotator",
#                                 table = "annotator_sets", 
#                                 fdr_method = "annotator",
#                                 with_slice = True )


############################################################
@follows( buildRepeatTrack,
          buildIntronTrack )
def setup():
    pass

@follows( buildMergedTracks,
          buildFilteredAlignment,
          loadGTF)
def prepare():
    pass

@follows( loadDifference,
          loadOverlap,
          #buildAnnotatorGCWorkspace, 
          #loadAnnotatorGeneSets,
          )
def difference(): pass

@follows( loadAnnotations,
          loadOverrun,
          loadDistances,
          loadRepeats,
          loadSegments,
          loadCodingPotential,
          )
def annotation():
    pass

@follows( loadRates, loadRepeatsRates )
def rates():
    pass

@follows( annotation, rates, difference )
def full(): pass

@follows( annotation, rates, loadSummary )
def summary(): pass

@follows( mkdir( "report" ) )
def build_report():
    '''build report from scratch.'''

    E.info( "starting documentation build process from scratch" )
    P.run_report( clean = True )


if __name__== "__main__":
    sys.exit( P.main(sys.argv) )

    

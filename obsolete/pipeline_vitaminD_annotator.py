################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
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
'''
pipeline_vitaminD_annotator.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python pipeline_vitaminD_annotator.py --help

Type::

   python pipeline_vitaminD_annotator.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import os
import tempfile
import collections
import shutil

import CGAT.Experiment as E
import CGAT.Pipeline as P
import sqlite3
import gff2annotator

PARAMS = P.getParameters()

############################################################
############################################################
############################################################
## Annotator utility functions
############################################################
def buildAnnotatorWorkSpace( tmpdir, outfile, workspaces=("genomic",), gc_control = False ):
    '''write genomic workspace.'''

    to_cluster = True
    job_options = "-l mem_free=4000M"
    
    tmpworkspaces = []

    if gc_control:
        tmpworkspace = os.path.join( tmpdir, "workspace_gc" )
        tmpsynonyms = os.path.join( tmpdir, "synonyms" )
        tmpworkspaces.append( tmpworkspace )

        statement = '''
        awk '{ printf("%%s\\t%%s\\t%%s\\t%%s.%%s\\n", $1,$2,$3,$1,$4)}' < %(annotator_gc)s |\
        python %(scriptsdir)s/bed2gff.py |\
 	python %(scriptsdir)s/gff2annotator.py \
		--section=workspace \
                --output-filename-synonyms=%(tmpsynonyms)s \
                --max-length=0 \
		--log=%(outfile)s.log \
        > %(tmpworkspace)s''' 
        
        P.run()
    else:
        tmpsynonyms = None

    for workspace in workspaces:
        tmpworkspace = os.path.join( tmpdir, "workspace_%s" % workspace )

        if workspace == "genomic":
            statement = '''
            python %(scriptsdir)s/gff2annotator.py \
                    --section=workspace \
                    --max-length=0 \
                    --log=%(outfile)s.log \
            < %(genome)s.gff > %(tmpworkspace)s
            '''
        elif workspace == "promotors":
            statement = '''
            python %(scriptsdir)s/gff2annotator.py \
                    --section=workspace \
                    --max-length=0 \
                    --log=%(outfile)s.log \
            < %(promotors)s > %(tmpworkspace)s
            '''
        elif workspace == "promotors-go":
            # promotors with GO categories
            statement = '''cat < %(promotors)s |\
            python %(scriptsdir)s/gtf2gtf.py \
            --filter=gene \
            --apply=<(cut -f 2 < %(gofile)s | sort | uniq ) |\
            python %(scriptsdir)s/gff2annotator.py \
            --section=workspace \
            --max-length=0 \
            --log=%(outfile)s.log \
            > %(tmpworkspace)s
            '''
        elif workspace == "gene-territories":
            statement = '''
            python %(scriptsdir)s/gff2annotator.py \
                    --section=workspace \
                    --max-length=0 \
                    --log=%(outfile)s.log \
            < %(annotator_geneterritories)s > %(tmpworkspace)s
            '''
        elif workspace == "mappable":
            statement = '''
            python %(scriptsdir)s/bed2gff.py < %(annotator_mappability)s |\
            python %(scriptsdir)s/gff2annotator.py \
                    --section=workspace \
                    --max-length=0 \
                    --log=%(outfile)s.log \
            > %(tmpworkspace)s
            '''
        else:
            raise P.PipelineError("unknown workspace '%s'" % workspace )

        P.run( **dict( locals().items() + PARAMS.items() ) )
        tmpworkspaces.append( tmpworkspace )

    return tmpworkspaces, tmpsynonyms

############################################################
############################################################
############################################################
##
############################################################
def buildAnnotatorAnnotations( tmpdir, outfile,
                               annotations=None,
                               bedfiles = None,
                               gofile = None):
    '''write annotations.'''

    tmpannotations = os.path.join( tmpdir, "annotations" )
    to_cluster = True
    job_options = "-l mem_free=4000M"
    
    if annotations == "architecture":
        statement = '''
         cat %(promotors)s %(annotation)s |\
         python %(scriptsdir)s/gff2annotator.py \
         	--section=annotations-gff \
         	--log=%(outfile)s.log \
         > %(tmpannotations)s
        '''
    elif annotations=="go":
        statement = '''
        cat %(annotator_geneterritories)s |\
        python %(scriptsdir)s/gff2annotator.py \
        --section=annotations-go \
        --input-filename-map=<(cut -f 2,4 < %(gofile)s) \
        --log=%(outfile)s.log \
        > %(tmpannotations)s
        '''
    elif bedfiles:
        bedfiles = " ".join(bedfiles)
        statement = '''
        cat %(bedfiles)s |\
        python %(scriptsdir)s/bed2annotator.py \
        --max-length=0 \
        --merge \
        --section=annotations \
        --log=%(outfile)s.log \
        > %(tmpannotations)s
        '''
    else:
        raise P.PipelineError("unknown annotations '%s'" % workspace )
    
    P.run()

    return tmpannotations

############################################################
############################################################
############################################################
##
############################################################
def buildAnnotatorSegments( tmpdir, infile, outfile ):
    '''convert segments in bed format to annotator format
    from infile to outfile.
    '''

    tmpsegments = os.path.join( tmpdir, "segments" )
    to_cluster = True
    
    statement = '''
        python %(scriptsdir)s/bed2gff.py < %(infile)s |\
	python %(scriptsdir)s/gff2annotator.py --log=%(outfile)s.log --section=segments > %(tmpsegments)s \
    '''

    P.run( **dict( locals().items() + PARAMS.items() ) )

    return tmpsegments

############################################################
############################################################
############################################################
##
############################################################
def buildAnnotatorSegmentsFromDatabase( tmpdir, track,
                                        outfile,
                                        with_motif = None,
                                        without_motif = None,
                                        proportion = None):
    '''output segments annotator format

    with_motif
       only output segments not matching motif
    without_motif
       only output segments not matching motif
    proportion
       only output top x percent of segments (peakval). If the
       proportion is negative, the bottom x percent are output.
       For example, proportion=-33 outputs the third of intervals
       with the smallest peakval.
    '''

    tmpsegments = os.path.join( tmpdir, "segments" )
    to_cluster = True

    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    if with_motif:
        statement = '''
        SELECT i.contig, i.start, i.end, i.interval_id, i.peakval
        FROM %(track)s_intervals AS i, %(track)s_mast AS m
        WHERE i.interval_id = m.id
             AND m.motif = '%(with_motif)s' 
             AND m.nmatches > 0
        ORDER by i.contig, i.start''' % locals()
    elif without_motif:
        statement = '''
        SELECT i.contig, i.start, i.end, i.interval_id, i.peakval
        FROM %(track)s_intervals AS i, %(track)s_mast AS m
        WHERE i.interval_id = m.id
             AND m.motif = '%(without_motif)s' 
             AND m.nmatches = 0
        ORDER by i.contig, i.start''' % locals()
    elif proportion != None:
        statement = '''SELECT COUNT(*) FROM %(track)s_intervals AS i'''
        cc = dbhandle.cursor()
        cc.execute( statement % locals() )
        total = cc.fetchall()[0][0]
        cutoff = int(float(total) * (abs(proportion) / 100.0))
        if proportion > 0:
            statement = '''
            SELECT i.contig, i.start, i.end, i.interval_id, i.peakval
            FROM %(track)s_intervals AS i
            ORDER BY i.peakval DESC LIMIT %(cutoff)i
            ''' % locals()
        else:
            statement = '''
            SELECT i.contig, i.start, i.end, i.interval_id, i.peakval
            FROM %(track)s_intervals AS i
            ORDER BY i.peakval ASC LIMIT %(cutoff)i
            ''' % locals()
            
    cc = dbhandle.cursor()
    try:
        cc.execute( statement )
    except sqlite3.OperationalError, msg:
        E.warn( "error in sql statement: %s" % msg)
        return None

    contigs = collections.defaultdict( list )
    for result in cc:
        contig, start, end, interval_id,peakval = result
        contigs[contig].append( (start,end) )

    outs = open(tmpsegments, "w" )
    gff2annotator.outputSegments( outs, contigs,
                                  section = "segments" )
    outs.close()
    
    return tmpsegments

############################################################
############################################################
############################################################
##
############################################################
def runAnnotator( tmpdir, outfile, 
                  tmpannotations, 
                  tmpsegments, 
                  tmpworkspaces, 
                  tmpsynonyms,
                  options = ""):
    '''run annotator.'''

    if tmpsegments == None:
        E.warn( "no segments - annotator not run for %s" % outfile )
        return

    if tmpannotations == None:
        E.warn( "no annotations - annotator not run for %s" % outfile )
        return
 
    to_cluster = True
    job_queue = "medium_jobs.q"
    job_options = "-l mem_free=8000M"

    workspace_options = ""
    for x,workspace in enumerate( tmpworkspaces ):
        if x == 0:
            workspace_options += " -workspace %s" % workspace
        else:
            workspace_options += " -workspace%i %s" % (x+1, workspace)

    if tmpsynonyms:
        workspace_options += " -synonyms %s" % tmpsynonyms
          

    statement = '''
    java -Xmx8000M -cp %(annotator_dir)s/commons-cli-1.0.jar:%(annotator_dir)s/Annotator.jar app.Annotator \
    -verbose 4 -iterations %(annotator_iterations)s \
    -annotation %(tmpannotations)s \
    -segments %(tmpsegments)s \
    %(workspace_options)s \
    %(options)s \
    > %(outfile)s '''
    
    P.run()

############################################################
############################################################
############################################################
## import annotator GO results
############################################################
def genericImportAnnotator( infiles, outfile, table, workspace, slice, subset, fdr_method ):
    '''generic import of annotator results.

    Assumes that the suffix of all infiles is the same.
    '''

    infile = " ".join(infiles)
    x, suffix = os.path.splitext( infiles[0] )

    tmpfilename = P.getTempFilename()

    statement = '''
	python %(scriptsdir)s/annotator.py \
		--method=fdr-table \
		--fdr-method=%(fdr_method)s \
		--log=%(outfile)s.log \
                --regex-id="(.*)%(suffix)s" \
                %(infile)s > %(tmpfilename)s
        '''
    P.run()

    tmpfile = P.getTempFile()
    
    for line in open( tmpfilename, "r" ):
        if line.startswith("id"):
            line = "subset\tworkspace\tslice\t" + re.sub("^id", "track", line)
        else:
            line = "%s\t%s\t%s\t%s" % (subset, workspace, slice, line)
        tmpfile.write(line)
    tmpfile.close()
    tmpfilename2 = tmpfile.name
        
    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
            --table=%(table)s 
    < %(tmpfilename2)s > %(outfile)s'''

    P.run()
    os.unlink( tmpfilename )
    os.unlink( tmpfilename2 )

############################################################
############################################################
############################################################
##
############################################################
def makeAnnotatorGO( infile, outfile, gofile, workspace ):
    '''check statistical overlap between intervals and genomic
    segements having GO assignments.

    worspace should be ``promotors`` or ``gene-territories``
    '''

    to_cluster = True
    # require 4Gb of free memory
    job_options = "-l mem_free=4000M"

    tmpdir = tempfile.mkdtemp( dir = os.getcwd() )

    annotations = buildAnnotatorAnnotations( tmpdir, 
                                             outfile, 
                                             annotations="go", 
                                             gofile = gofile )

    # take only those promotors with GO categories for workspace
    workspaces, synonyms = buildAnnotatorWorkSpace(
        tmpdir, outfile,
        workspaces = ("mappable", workspace),
        gc_control=True )
    
    segments = buildAnnotatorSegments( tmpdir, infile, outfile )

    runAnnotator( tmpdir, outfile, annotations, segments, workspaces, synonyms )

############################################################
############################################################
############################################################
##
############################################################
def buildAnnotatorSegmentsROI( tmpdir, roi_class, outfile, overlap = None ):
    '''convert segments in bed format to annotator format
    from infile to outfile.
    '''

    tmpsegments = os.path.join( tmpdir, "segments" )
    to_cluster = True

    dbhandle = sqlite3.connect( PARAMS["database"] )

    if overlap:
            statement = '''
            SELECT roi.contig, roi.start, roi.end
            FROM regions_of_interest AS roi,
                 %(overlap)s_intervals AS i
            WHERE roi.class='%(roi_class)s' AND 
                  i.contig = roi.contig AND
                  min(roi.end, i.end) - max(roi.start, i.start) > 0
        '''
    else:
        statement = '''
            SELECT roi.contig, roi.start, roi.end
            FROM regions_of_interest AS roi
            WHERE class='%(roi_class)s'
        '''

    cc = dbhandle.cursor()
    cc.execute( statement % locals() )

    noutput = 0
    contigs = collections.defaultdict( list )
    for result in cc:
        contig, start, end = result
        contigs[contig].append( (start,end) )
        noutput += 1

    E.info("segments for roi_class `%s` and overlap `%s`: %i" % (roi_class, overlap, noutput))
        
    outs = open(tmpsegments, "w" )
    gff2annotator.outputSegments( outs, contigs,
                                  section = "segments" )
    outs.close()

    if noutput == 0:
        return None
    else:
        return tmpsegments

############################################################
############################################################
############################################################
##
############################################################
def makeAnnotatorROIGO( roi_class, outfile, gofile, workspace, overlap = None ):
    '''check statistical overlap between intervals and genomic
    segements having GO assignments.

    worspace should be ``promotors`` or ``gene-territories``
    '''

    to_cluster = True
    # require 4Gb of free memory
    job_options = "-l mem_free=4000M"

    tmpdir = tempfile.mkdtemp( dir = os.getcwd() )

    segments = buildAnnotatorSegmentsROI( tmpdir, 
                                          roi_class,
                                          outfile, 
                                          overlap = overlap )

    if segments == None: 
        E.info("no segments for roi_class `%s` and overlap `%s` - no computation." %(roi_class,
                                                                                     overlap ))
        return
    annotations = buildAnnotatorAnnotations( tmpdir, 
                                             outfile, 
                                             annotations="go", 
                                             gofile = gofile )

    # take only those promotors with GO categories for workspace
    workspaces, synonyms = buildAnnotatorWorkSpace(
        tmpdir, outfile,
        workspaces = (workspace,),
        gc_control=True )
    
    # these are large segments, so increase bucket size
    runAnnotator( tmpdir, outfile, annotations, segments, workspaces, synonyms,
                  "-bucketsize 100" )

############################################################
############################################################
############################################################
##
############################################################
def makeAnnotatorArchitecture( infile, outfile, **kwargs ):
    '''check statistical overlap between intervals and and other genomic features
    defined in the file PARAMS["annotations"].

    Annotator is run with the following parameters:
    1. Segments: the interval track
    2. Annotations:
       1. genomic architecture (PARAMS["annotation"])
       2. promotors (PARAMS["promotors"])
    3. Workspace: the full genome
    '''
    
    tmpdir = tempfile.mkdtemp( dir = os.getcwd() )

    track = infile[:-len(".bed")]

    if kwargs:
        segments = buildAnnotatorSegmentsFromDatabase( tmpdir,
                                                                 track, outfile,
                                                                 **kwargs )
    else:
        segments = buildAnnotatorSegments( tmpdir, infile, outfile )


    workspaces, synonyms = buildAnnotatorWorkSpace( tmpdir, outfile,
                                                    workspaces = ("mappable","genomic"),
                                                    gc_control=True )
    
    
    annotations = buildAnnotatorAnnotations( tmpdir, outfile, annotations="architecture" )
    
    runAnnotator( tmpdir, outfile, annotations, segments, workspaces, synonyms )

    shutil.rmtree( tmpdir )

############################################################
############################################################
############################################################
##
############################################################
def makeAnnotatorTracks( infiles, outfile, **kwargs ):
    '''check statistical overlap between intervals and selected ucsc tracks

    Annotator is run with the following parameters:
    1. Segments: the interval track
    2. Annotations:
       1. ucsc encode features
       2. disease intervals (regions of interest)
    3. Workspace: the full genome
    '''
    
    infile, infile_annotations = infiles

    track = infile[:-len(".bed")]

    tmpdir = tempfile.mkdtemp( dir = os.getcwd() )

    if kwargs:
        segments = buildAnnotatorSegmentsFromDatabase( tmpdir,
                                                                 track, outfile,
                                                                 **kwargs )
    else:
        segments = buildAnnotatorSegments( tmpdir, infile, outfile )
    
    workspaces, synonyms = buildAnnotatorWorkSpace( tmpdir, outfile,
                                                    workspaces = ("mappable","genomic"),
                                                    gc_control = True )

    annotations = buildAnnotatorAnnotations( tmpdir, outfile, bedfiles=(infile_annotations,) )


    runAnnotator( tmpdir, outfile, annotations, segments, workspaces, synonyms )

    shutil.rmtree( tmpdir )

############################################################
############################################################
############################################################
##
############################################################
def makeAnnotatorRegionsOfInterest( infiles, outfile, **kwargs ):
    '''check statistical overlap between intervals regions of interest.

    Annotator is run with the following parameters:
    1. Segments: the interval track
    2. Annotations:
       1. disease intervals (regions of interest)
    3. Workspace: mappable part of gene territories
    '''

    infile, infile_regions = infiles
    track = infile[:-len(".bed")]

    tmpdir = tempfile.mkdtemp( dir = os.getcwd() )

    annotations = buildAnnotatorAnnotations( tmpdir, outfile, bedfiles=(infile_regions,) )

    workspaces, synonyms = buildAnnotatorWorkSpace( tmpdir,
                                                    outfile,
                                                    workspaces = ("mappable", "gene-territories"),
                                                    gc_control = True )

    
    if kwargs:
        segments = buildAnnotatorSegmentsFromDatabase( tmpdir,
                                                       track, outfile,
                                                       **kwargs )
    else:
        segments = buildAnnotatorSegments( tmpdir, infile, outfile )



    runAnnotator( tmpdir, outfile, annotations, segments, workspaces, synonyms )

    shutil.rmtree( tmpdir )


############################################################
############################################################
############################################################
##
############################################################
def buildAnnotatorDistanceAnnotations( annotations = "expression" ):
    '''build an annotations file for annotator_distance.'''

    tmpfile = P.getTempFile( "." )
    tmpfilename = tmpfile.name 

    if annotations == "expression":
        dbhandle = sqlite3.connect( PARAMS["database"] )
        cc = dbhandle.cursor()

        statement = """
        SELECT gene_id,
        CASE WHEN %(annodist_master_expression_select)s THEN 'responsive' ELSE 'nonresponsive' END
        FROM probeset2transcript AS e,
        %(annodist_master_expression)s AS d 
        WHERE d.cluster_id = e.cluster_id
        """ % dict( locals().items() + PARAMS.items() )

        data = cc.execute( statement ).fetchall()
        tmpfile.write("gene_id\tlabel\n" )
        for gene_id, label in data: tmpfile.write("%s\t%s\n" % (gene_id, label ) )
        tmpfile.close()

    return tmpfilename

############################################################
############################################################
############################################################
##
############################################################
def makeAnnotatorDistance( infile, 
                           outfile, 
                           builder,
                           workspace, 
                           workspace_label="direction",
                           annotations = None ):
    '''check statistical association between intervals and 
    transcription start sites.

    '''

    to_cluster = True

    target_path = os.path.join( os.path.abspath( PARAMS["exportdir"] ), 
                                "annotator_distance", 
                                outfile )

    if os.path.exists( target_path): shutil.rmtree( target_path)
    try:
        os.makedirs( target_path )
    except OSError: 
        pass

    options = [] 

    if annotations: 
        options.append( "--filename-annotations=%s" % annotations )
        
    options = " ".join( options )

    statement = '''
	python %(scriptsdir)s/annotator_distance.py \
		--workspace=%(workspace)s \
		--segments=%(infile)s \
                --segments-format=bed \
		--counter=%(annodist_counter)s \
		--workspace-label=%(workspace_label)s \
		--sampler=permutation \
		--transform-counts=cumulative \
		--logscale=x \
                --remove-overhangs \
		--analysis=proximity \
		--num-samples=%(annodist_iterations)i \
		--num-bins=%(annodist_bins)i \
		--hardcopy=%(target_path)s/%%s.png \
		--output-filename-pattern=%(target_path)s/%%s.table \
                --workspace-builder=%(builder)s \
                --resolution=%(annodist_resolution_intergenic)s \
		--plot 
                %(options)s < /dev/null > %(outfile)s'''

    P.run()




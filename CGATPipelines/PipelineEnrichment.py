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
PipelineEnrichment.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

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

try:
    PARAMS = P.getParameters()
except IOError:
    pass

def outputSegments( outfile,
                    intervals,
                    section,
                    outfile_synonyms = None,
                    max_length = 100000,
                    remove_regex = None ):

    if section == "segments":
        prefix = "##Segs"
    elif section == "workspace":
        prefix = "##Work"
    else:
        raise ValueError("invalid section `%s`" % section )
    
    ninput, ncontigs, nsegments, ndiscarded = 0, 0, 0, 0
    for contig, gffs in intervals.items():
        ninput += 1
        if remove_regex and remove_regex.search(contig): continue

        if max_length:
            segments = [ x for x in gffs if x[1]-x[0] <= max_length ]
        else:
            segments = [ x for x in gffs ]

        nsegments += len(segments)
        ndiscarded += len(gffs) - len(segments)

        if outfile_synonyms:
            synonyms = collections.defaultdict( list )
            for x in segments:
                synonyms[x[2].source].append( (x[0], x[1]) )

            for key, segs in synonyms.iteritems():
                outfile.write( "%s\t%s\t%s\n" % (prefix, key, "\t".join(
                    [ "(%i,%i)" % x for x in segs ] )))
                outfile_synonyms.write("##Synonym\t%s\t%s\n" % (key, contig) )

        else:
            outfile.write( "%s\t%s\t%s\n" % (prefix, contig, "\t".join(
                [ "(%i,%i)" % x for x in segments ] )))

        ncontigs += 1
    return ninput, nsegments, ndiscarded, ncontigs

############################################################
############################################################
############################################################
## 
############################################################
def buildGenomeGCSegmentation( infile, outfile ):
    '''segment the genome into windows according to G+C content.'''

    to_cluster = True
    
    statement = '''
    python %(scriptsdir)s/fasta2bed.py 
        --method=fixed-width-windows 
        --window-size=%(enrichment_gc_window_size)i 
        --log=%(outfile)s.log 
    < %(genome)s.fasta > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
## 
############################################################
def buildAnnotatorGC( infile, outfile ):
    '''compute G+C regions.'''

    to_cluster = True
    statement = '''
    python %(scriptsdir)s/bed2bed.py 
        --method=bins
        --num-bins=%(enrichment_gc_bins)s
        --binning-method=%(enrichment_gc_method)s 
        --log=%(outfile)s.log 
    < %(infile)s > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
## 
############################################################
def buildIsochoresGC( infile, outfile ):
    '''compute isochores based on G+C content
    '''
    to_cluster = True
    
    statement = '''
    python %(scriptsdir)s/fasta2bed.py 
        --method=fixed-width-windows 
        --window-size=%(enrichment_gc_window_size)i 
        --log=%(outfile)s.log 
    < %(infile)s 
    | python %(scriptsdir)s/bed2bed.py 
        --method=bins 
        --num-bins=%(enrichment_gc_bins)s 
        --binning-method=%(enrichment_gc_method)s 
        --log=%(outfile)s.log 
    | gzip 
    > %(outfile)s'''

    P.run()

############################################################
############################################################
############################################################
## Annotator utility functions
############################################################
def buildWorkSpace( outfile, workspace):
    '''write genomic workspace.

    Available workspaces are:

    genomic
       the full genome
    intronic
       introns (requires annotator_regions to be set)
    exonic
       exonic (requires annotator_regions to be set)
    intergenic
       introns (requires annotator_regions to be set)
    geneterritories
       introns (requires annotator_geneterritories to be set)
    mappable
       mappable part of genome (requires annotator_mappability to be set )
    alignable
       only the alignable part of a genome (requires annotator_alignment to be set)


    If ``gc_control`` is True, the chromosomes will be divided into isochores
    (requiers the paramater ``annotator_gc_workspace`` to be set).
    '''

    to_cluster = True
    job_options = "-l mem_free=4000M"

    workspace = workspace.lower()

    if workspace == "genomic":
        P.checkParameter( "genome" )

        statement = '''
        python %(scriptsdir)s/index2bed.py 
                --genome=%(genome)s 
                --log=%(outfile)s.log 
                --remove-regex='%(enrichment_remove_pattern)s'
        > %(outfile)s
        '''

    elif workspace in ("intergenic", "intronic", "cds" ):

        P.checkParameter( "enrichment_regions" )

        workspace_upper = workspace.upper()

        statement = '''
        gunzip < %(enrichment_regions)s 
        | awk 'BEGIN {printf("track name=%(workspace)s\\n"); } 
               ($3 == "%(workspace)s" 
               || $3 == "%(workspace_upper)s") 
               && !( $1 ~ /%(enrichment_remove_pattern)s/)
               { printf("%%s\\t%%i\\t%%i\\n", $1, $4-1, $5); }'
        > %(outfile)s
        '''
    elif workspace == "unknown":

        P.checkParameter( "enrichment_regions" )
        statement = '''
        awk '($3 == "intronic" || $3 == "intergenic" )' 
        < %(enrichment_regions)s
        | python %(scriptsdir)s/gff2enrichment.py 
                --section=workspace 
                --max-length=0 
                --log=%(outfile)s.log 
                --remove-regex='%(enrichment_remove_pattern)s'
        > %(outfile)s
        '''

    elif workspace == "known":
        P.checkParameter( "enrichment_regions" )
        statement = '''
        awk '($3 == "CDS" || $3 ~ /UTR/ || $3 ~ /flank/)' 
        < %(enrichment_regions)s
        | python %(scriptsdir)s/gff2enrichment.py 
                --section=workspace 
                --max-length=0 
                --log=%(outfile)s.log 
                --remove-regex='%(enrichment_remove_pattern)s'
        > %(outfile)s
        '''

    elif workspace == "alignable":

        P.checkParameter( "enrichment_alignment" )
        statement = '''gunzip
        < %(enrichment_alignment)s 
        | sort -k10,10 
        | awk '$10 !~ /%(enrichment_remove_pattern)s/ \
            {if ($10!=l) {printf("\\n##Work\\t%%s", $10); l=$10;} \
            printf("\\t(%%i,%%i)", $12,$13); }; \
        END {printf ("\\n");}'\
        > %(outfile)s
        '''

    elif workspace == "gene-territories":

        P.checkParameter( "enrichment_geneterritories" )
        statement = '''
        python %(scriptsdir)s/gff2enrichment.py \
                --section=workspace \
                --max-length=0 \
                --log=%(outfile)s.log \
                --remove-regex='%(enrichment_remove_pattern)s'
        < %(enrichment_geneterritories)s > %(outfile)s
        '''

    elif workspace == "mappable":

        P.checkParameter( "enrichment_mappability" )
        statement = '''
        python %(scriptsdir)s/bed2gff.py < %(enrichment_mappability)s 
        | python %(scriptsdir)s/gff2enrichment.py 
                --section=workspace 
                --max-length=0 
                --log=%(outfile)s.log 
                --remove-regex='%(enrichment_remove_pattern)s'
        > %(outfile)s
        '''
    else:
        raise P.PipelineError("unknown workspace '%s'" % workspace )

    P.run()

############################################################
############################################################
############################################################
##
############################################################
def buildAnnotatorAnnotations( tmpdir, outfile,
                               annotations=None,
                               bedfiles = None,
                               gfffiles = None,
                               gofile = None ):
    '''write annotations in annotator format.
    '''

    tmpannotations = os.path.join( tmpdir, "annotations" )
    to_cluster = True
    job_options = "-l mem_free=4000M"
    
    if annotations == "architecture":
        statement = '''
         cat %(promotors)s %(annotation)s 
         | python %(scriptsdir)s/gff2annotator2tsv.py 
         	--section=annotations-gff 
         	--log=%(outfile)s.log 
                --remove-regex='%(annotator_remove_pattern)s'
         > %(tmpannotations)s
        '''
    elif annotations=="go":
        statement = '''
        python %(scriptsdir)s/gff2annotator2tsv.py 
        --section=annotations-go 
        --input-filename-map=<(cut -f 2,4 < %(gofile)s) 
        --log=%(outfile)s.log
        --remove-regex='%(annotator_remove_pattern)s'
        < %(annotator_geneterritories)s  
        > %(tmpannotations)s
        '''
    elif bedfiles:
        bedfiles = " ".join(bedfiles)
        statement = '''
        cat %(bedfiles)s 
        | python %(scriptsdir)s/bed2annotator2tsv.py 
        --max-length=0 
        --merge 
        --section=annotations 
        --log=%(outfile)s.log 
        > %(tmpannotations)s
        '''
    else:
        raise P.PipelineError("unknown annotations '%s'" % annotations )
    
    P.run()

    return tmpannotations

def buildGeneSetAnnotations( infiles, outfile, slice ):
    '''build annotations of all sets from database.

    ``slice`` can be any of the slices in the ``annotation`` 
    tables.'''

    statement = '''SELECT gene_id FROM %(track)s_annotation as a WHERE %(where)s'''

    if slice == "all": where = "'1'"
    else: where = "is_%(slice)s" % locals()

    dbhandle = sqlite3.connect( PARAMS["database"] )
    
    subsets = []
    
    for f in infiles:

        assert f.endswith( ".gtf.gz" )
        track = f[:-len(".gtf.gz")]
        key = "%s.%s" % (track,slice)

        cc = dbhandle.cursor()
        data = [x[0] for x in cc.execute( statement % locals() ).fetchall() ]
        E.info( "%s: adding %i genes" % (key, len(data)))

        filename = outfile + ".tmp.%s" % key
        outf = open( filename, "w" )
        outf.write( "gene_id\n%s\n" % "\n".join( map(str, data) ))
        outf.close()

        subsets.append( "--subset=%s" % ",".join( (track, key, filename) ) )
        
    infiles = " ".join(infiles)
    subsets = " ".join(subsets)

    to_cluster = True
    statement = '''
	python %(scriptsdir)s/gff2annotator2tsv.py 
		--section=annotations-genes 
		--log=%(outfile)s.log 
                --remove-regex='%(annotator_remove_pattern)s'
		%(subsets)s
		%(infiles)s
	> %(outfile)s
    '''

    P.run()

    statement = '''
    rm -f %(outfile)s.tmp* 
    '''

    P.run()

############################################################
############################################################
############################################################
##
############################################################
def buildAnnotatorSlicedSegments( tmpdir, outfile, track, slice ):
    '''slice segments.'''
    
    tmpsegments = os.path.join( tmpdir, "segments" )
    to_cluster = True

    if slice == "all": where = "'1'"
    else: where = "is_%(slice)s" % locals()
    
    statement = '''
        %(cmd-sql)s %(database)s 
        "SELECT g.* FROM %(track)s_gtf as g, %(track)s_annotation AS a WHERE a.gene_id = g.gene_id AND %(where)s"
        | python %(scriptsdir)s/gtf2tsv.py --invert 
	| python %(scriptsdir)s/gff2annotator2tsv.py 
               --remove-regex='%(annotator_remove_pattern)s'
               --log=%(outfile)s.log 
               --section=segments 
        > %(tmpsegments)s
    '''

    P.run()

    if os.path.getsize(tmpsegments) == 0:
        return None
    else:
        return tmpsegments

############################################################
############################################################
############################################################
##
############################################################
def buildAnnotatorSegments(tmpdir, infile, outfile ):
    '''convert segments in bed format to annotator format
    from infile to outfile.
    '''

    tmpsegments = os.path.join( tmpdir, "segments" )
    to_cluster = True
    
    statement = '''
        python %(scriptsdir)s/bed2gff.py < %(infile)s 
	| python %(scriptsdir)s/gff2annotator2tsv.py 
                --remove-regex='%(annotator_remove_pattern)s'
                --log=%(outfile)s.log --section=segments 
    > %(tmpsegments)s
    '''

    P.run()

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
    java -Xmx8000M -cp %(annotator_dir)s/commons-cli-1.0.jar:%(annotator_dir)s/Annotator.jar app.Annotator 
    -verbose 4 -iterations %(annotator_iterations)s 
    -annotation %(tmpannotations)s 
    -segments %(tmpsegments)s 
    -bucketsize %(annotator_bucketsize)i 
    %(workspace_options)s 
    %(options)s 
    > %(outfile)s '''
    
    P.run( **dict( locals().items() + PARAMS.items() ) )

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
	python %(scriptsdir)s/annotator2tsv.py \
		--method=fdr-table \
		--fdr-method=%(fdr_method)s \
		--log=%(outfile)s.log \
                --regex-id="(.*)%(suffix)s" \
                %(infile)s > %(tmpfilename)s
        '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

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

    P.run( **dict( locals().items() + PARAMS.items() ) )
    os.unlink( tmpfilename )
    os.unlink( tmpfilename2 )

############################################################
############################################################
############################################################
## import annotator GO results
############################################################
def importAnnotator( infiles, outfile, regex_id, table, 
                     fdr_method,
                     with_slice = False ):
    '''generic import of annotator results.
    
    If with-slice is true, the first id field is assumed
    to be a concatenation of track.slice.
    '''

    infile = " ".join(infiles)

    if with_slice:
        transform = '''perl -p -e "s/^id/track\\tslice/; s/\\./\\t/" '''
    else:
        transform = '''sed "s/^id/track/"'''

    statement = '''
	python %(scriptsdir)s/annotator2tsv.py 
		--method=fdr-table 
		--fdr-method=%(fdr_method)s 
		--log=%(outfile)s.log 
                --regex-id="%(regex_id)s" 
                %(infile)s 
        | %(transform)s
        |python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                --table=%(table)s 
        > %(outfile)s
    '''
    
    P.run()


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
def makeAnnotatorArchitecture( infile, outfile, 
                               workspaces = ("mappable","genomic"),
                               **kwargs ):
    '''check statistical overlap between intervals and and other genomic features
    defined in the file PARAMS["annotations"].

    Annotator is run with the following parameters:

      1 Segments: the interval track

      2 Annotations:

         1 genomic architecture (PARAMS["annotation"])
         2 promotors (PARAMS["promotors"])

      3 Workspace: the full genome
    '''
    
    tmpdir = tempfile.mkdtemp( dir = os.getcwd() )

    track = infile[:-len(".bed")]

    segments = buildAnnotatorSegments( tmpdir, infile, outfile )


    workspaces, synonyms = buildAnnotatorWorkSpace( tmpdir, outfile,
                                                    workspaces = workspaces,
                                                    gc_control=True )
    
    
    annotations = buildAnnotatorAnnotations( tmpdir, outfile, annotations="architecture" )
    
    runAnnotator( tmpdir, outfile, annotations, segments, workspaces, synonyms )

    shutil.rmtree( tmpdir )

############################################################
############################################################
############################################################
##
############################################################
def makeAnnotator( infile, outfile, 
                   segments,
                   annotations,
                   workspaces = ("mappable","genomic"),
                   gc_control = True ):
    '''check statistical overlap between intervals and and other genomic features
    '''
    
    tmpdir = tempfile.mkdtemp( dir = os.getcwd() )

    # build work spaces
    workspaces, synonyms = buildAnnotatorWorkSpace( tmpdir, outfile,
                                                    workspaces = workspaces,
                                                    gc_control=gc_control )
        
    annotations = buildAnnotatorAnnotations( tmpdir, outfile,  )
    
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
      1 Segments: the interval track

      2 Annotations:
         1 ucsc encode features
         2 disease intervals (regions of interest)

      3 Workspace: the full genome
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

      1 Segments: the interval track

      2 Annotations:
         1 disease intervals (regions of interest)

      3 Workspace: mappable part of gene territories

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



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
annotator2tsv.py - analyse annotator results
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes annotator results, parses them, 
filters by FDR and outputs the results 

Usage
-----

Example::

   python annotator2tsv.py --help

Type::

   python annotator2tsv.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import optparse
import collections
import bisect
import os

USAGE="""python %s [OPTIONS] input1 input2

parse annotator files and compute overlap.

Note: this currently does a brute force overlap. In the future
use bme or rtree indices.

Version: $Id: annotator2tsv.py 2861 2010-02-23 17:36:32Z andreas $
""" % sys.argv[0]

import CGAT.Experiment as E
import CGAT.GTF as GTF
import numpy
import CGAT.Stats as Stats

class Identifier:
    def __init__( self, line = None):
        if line: self.read(line)
            
    def read( self, line ):
        data = re.split("\s", line[:-1])
        assert data[0] == "##Id", "formatting error in line %s" % line
        self.mId = int(data[1])
        self.contig = data[2]
        self.mSegments = map( lambda x: tuple(map(int, x[1:-1].split(","))), data[3:] )

    def __str__(self):
        return "##Id\t%s\t%s\t%s" % (self.mId, 
                                     self.contig, 
                                     "\t".join( ["(%i,%i)" % x for x in self.mSegments] ) )

class Workspace:
    def __init__( self, line = None):
        if line: self.read(line)
            
    def read( self, line ):
        data = re.split("\s", line[:-1])
        assert data[0] == "##Work", "formatting error in line %s" % line
        self.contig = data[1]
        self.mSegments = map( lambda x: tuple(map(int, x[1:-1].split(","))), data[2:] )

    def __str__(self):
        return "##Work\t%s\t%s" % (self.contig, "\t".join( ["(%i,%i)" % x for x in self.mSegments] ) )

class Segments:
    def __init__( self, line = None):
        if line: self.read(line)
            
    def read( self, line ):
        data = re.split("\s", line[:-1])
        assert data[0] == "##Segs", "formatting error in line %s" % line
        self.contig = data[1]
        self.mSegments = map( lambda x: tuple(map(int, x[1:-1].split(","))), data[2:] )
    def __str__(self):
        return "##Segs\t%s\t%s" % (self.contig, "\t".join( ["(%i,%i)" % x for x in self.mSegments] ) )

class Annotation:
    def __init__( self, line = None):
        if line: self.read(line)
            
    def read( self, line ):
        data = re.split("\t", line[:-1])
        assert data[0] == "##Ann", "formatting error in line %s" % line
        self.mId = data[1]
        # deal with empty annotations
        if data[2] == "": 
            self.mIdentifiers = []
            return

        try:
            self.mIdentifiers = map(int, data[2:])
        except ValueError, msg:
            raise ValueError( "formatting error in line %s: %s" % (line, msg))

    def __str__(self):
        return "##Ann\t%s\t%s" % (self.mId, "\t".join(map(str, self.mIdentifiers)) )

def readSegments( infile ):
    segments = []
    for line in infile:
        if line.startswith( "##Segs" ):
            segments.append( Segments( line=line ) )


    return segments

def readWorkspace( infile ):
    workspace = []
    for line in infile:
        if line.startswith( "##Work" ):
            workspace.append( Workspace( line=line ) )

    return workspace

def readAnnotations( infile ):
    identifiers = []
    annotations = []
    for line in infile:
        if line.startswith( "##Id" ):
            identifiers.append( Identifier( line=line ) )
        elif line.startswith( "##Ann" ):
            annotations.append( Annotation( line = line ) )    

    map_query2annotation, map_annotation2identifier = {}, {}
    for identifier in identifiers:
        map_annotation2identifier[identifier.mId] = identifier
    for annotation in annotations:
        map_query2annotation[annotation.mId] = annotation

    return identifiers, annotations, map_query2annotation, map_annotation2identifier

def doQuery( options, args):
    """query annotator input for overlapping segments."""

    options.queries.extend( args )

    identifiers, annotations, map_query2annotation, map_annotation2identifier = readAnnotations( open( options.filename_annotations, "r") )
    segments = readSegments( open( options.filename_segments, "r") )
    workspace = readWorkspace( open( options.filename_workspace, "r") )

    if options.filename_segments_gtf:
        with open(options.filename_segments_gtf, "r") as infile:
            segments_gtfs = GTF.readAsIntervals( GTF.iterator( infile ), with_records = True )
    else:
        segments_gtfs = {}

    if options.filename_annotations_gtf:
        with open(options.filename_annotations_gtf, "r") as infile:
            annotations_gtfs = GTF.readAsIntervals( GTF.iterator( infile ), with_records = True )
    else:
        annotations_gtfs = {}

    E.info( "identifiers=%i, annotations=%i, segments=%i, workspace=%i" % (len(identifiers), len(annotations), len(segments), len(workspace)) )

    if options.loglevel >= 3:
        for x in segments: options.stdlog.write( str(x) + "\n" )
        for x in workspace: options.stdlog.write( str(x) + "\n" )
        for x in annotations: options.stdlog.write( str(x) + "\n" )
        for x in identifiers: options.stdlog.write( str(x) + "\n" )

    if options.output_format == "table":
        options.stdout.write("query\tcontig\tanno_start\tanno_end\tanno_id\tseg_start\tseg_end\tgene_id\n" )

    if "all" in options.queries:
        options.queries = map_query2annotation.keys()
        E.debug( "doing all queries: %i" % len(options.queries) )

    for query in options.queries:

        annotation = map_query2annotation[query] 

        # get segments belonging to annotation
        annotation_segments = []
        for id in annotation.mIdentifiers:
            identifiers = map_annotation2identifier[id]
            for start, end in identifiers.mSegments:
                annotation_segments.append( (identifiers.contig, start, end) )

        E.info( "annotation_segments=%s"% len(annotation_segments ) )
        E.debug( "annotation_segments=%s"% str(annotation_segments ) )
        # quick and dirty overlap with workspace
        workspace_overlaps = []
        for contig, start, end in annotation_segments:
            finished = False
            for wrk in workspace:
                if wrk.contig != contig: continue
                for wrk_start, wrk_end in wrk.mSegments:
                    if min(wrk_end, end) - max(wrk_start, start) > 0:
                        workspace_overlaps.append( (contig, start, end ) )
                        finished = True
                        break
                if finished: break

        E.info( "workspace_overlaps=%s"% len(workspace_overlaps ) )
        E.debug( "workspace_overlaps=%s" % str(workspace_overlaps) )

        # Quick and dirty overlap with segments
        overlaps = []
        for contig, start, end in workspace_overlaps:
            for segment in segments:
                if segment.contig != contig: continue
                for seg_start, seg_end in segment.mSegments:
                    if min(seg_end, end) - max(seg_start, start) > 0:
                        overlaps.append( (segment.contig, start, end, seg_start, seg_end) )

        for contig, anno_start, anno_end, seg_start, seg_end in overlaps:
            seg_gene_ids = []
            if contig in segments_gtfs:
                for start,end,gtf in segments_gtfs[contig]:
                    if seg_start == start and seg_end == end:
                        seg_gene_ids.append( gtf.gene_id )

            anno_gene_ids = []
            if contig in annotations_gtfs:
                for start,end, gtf in annotations_gtfs[contig]:
                    if anno_start == start and anno_end == end:
                        anno_gene_ids.append( gtf.gene_id )

            if options.output_format == "table":
                options.stdout.write( "%s\t%s\t%i\t%i\t%s\t%i\t%i\t%s\n" % \
                                          (query, contig, 
                                           anno_start, anno_end, 
                                           ";".join(anno_gene_ids),
                                           seg_start, seg_end, 
                                           ";".join(seg_gene_ids) ) )
            elif options.output_format == "gff":
                options.stdout.write( "\t".join( (contig, 
                                                  "anno", 
                                                  "anno",
                                                  str(seg_start + 1), 
                                                  str(seg_end), 
                                                  ".",
                                                  "+",
                                                  ".",
                                                  'annotation_id = "%s" % query') ) + "\n" )
                
class AnnotatorTest:
    def __init__(self):
        self.mFDR = 0
        self.mQValue = 1

class DataFDR:
    def __init__(self):
        pass

def readAnnotator( infile,
                   delims = "",
                   ignore = "" ):
    """parse annotator output.

    returns (data, fdr).
    """

    kw = "Iteration"
    data = []
    annotator_fdr = {}
    annotator_level = None

    synonyms = {}
    input_files = {}
    
    pvalues, qvalues = [], []

    for line in infile:

        # skip meaningless lines
        if line.startswith("Iteration"): continue       # skip iteration comment lines
        if line.startswith( "Z" ): continue            # skip header line of results
        if len(line) == 1: continue                  # skip trailing blank lines

        if line.startswith( "# Found duplicated annotation:" ):
            term = re.match( "# Found duplicated annotation: (.*)$", line[:-1] ).groups()[0]
            synonyms[term] = []
            continue
        
        if line.startswith( "#  Removing:" ):
            synonym = re.match( "#  Removing:\s+(.*)$", line[:-1] ).groups()[0]
            synonyms[term].append( synonym )
            continue

        if line.startswith( "# Reading file" ):
            filename, option = re.search( "# Reading file (\S+) for option (\S+)", line ).groups()
            filename = os.path.basename( filename )
            if option.startswith("workspace"):
                if "workspace" not in input_files:
                    input_files["workspace"] = filename
                else:
                    input_files["workspace"] += "-%s" % filename
            else:
                input_files[option] = filename

        # parse FDR
        # annotator only reports the FDR for fixed p-value intervals.
        if line.startswith( "--" ):
            if line.startswith("-- False"):
                annotator_level = float(re.search( "-- False Discovery summary for p-value (.+):", line ).groups()[0])
                annotator_fdr[annotator_level] = {}
            elif line.startswith("--  Category"):
                pass
            else:
                if re.search("insufficiently", line): continue
                dd = re.split("\s+", line[4:-1])

                d = DataFDR()
                d.mObserved, d.mAverage, d.mMedian, d.m95 = map(float,dd[1:])
                if d.mObserved > 0:
                    d.mQValue = d.mAverage / d.mObserved
                else:
                    d.mQValue = 1.0
                annotator_fdr[annotator_level][dd[0]] = d
                if dd[0] == "Significant":
                    pvalues.append( annotator_level )
                    qvalues.append( d.mQValue )

        if len(line[:-1].split('\t')) != 9: continue # HACK: accounts for a bug in Annotator output


        try:
            (z, percentchange, pvalue, observed, expected, low95, up95, stddev, description) = line[:-1].split('\t')[:9]
        except ValueError:
            raise ValueError("# parsing error in line: %s" % line[:-1])
        
        # remove leading .
        if description.startswith("."): description = description[1:]

        pvalue = float( pvalue )

        
        d = AnnotatorTest()
        d.mAnnotation = description
        d.mPValue = pvalue
        d.mFoldChange = 1.0 + float(percentchange) / 100.0
        d.mObserved = float(observed)
        d.mExpected = float(expected)
        d.mCI95 = float(low95), float(up95)
        d.mStdDev = float(stddev)
        d.mQValue = 1.0

        ## apply filters
        for c in delims:
            d.mAnnotation = d.mAnnotation.split(c)[0]
        for c in ignore:
            d.mAnnotation = d.mAnnotation.replace(c, '')
            
        data.append(d)

    # sort by increasing value
    pvalues.reverse()
    qvalues.reverse()
    # compute fdr
    # the q-value of each observation is set to the
    # qvalue for the next highest P-Value.
    if len(pvalues) == 0:
        return data, annotator_fdr, synonyms, input_files

    max_pvalue = pvalues[-1]
    
    for d in data:
        pvalue = d.mPValue
        if pvalue > max_pvalue: d.mQvalue = 1.0
        else: d.mQValue = qvalues[max(0, bisect.bisect(pvalues, pvalue) - 1)]

    return data, annotator_fdr, synonyms, input_files

def doOldFDR( options, args ):
    """apply fdr to output of annotator."""

    # read input
    annotators = []
    for filename in args:
        infile = open(filename,"r")
        annotators.append( readAnnotator( infile ))
        infile.close()
        
    # apply filters and create diagnostic plots    
    for filename, data in zip(args, annotators):
        ninput = len(data)
        pvalues = [x.mPValue for x in data]
        vlambda = numpy.arange( 0, max( pvalues ), 0.05 )
        try:
            qvalues = Stats.doFDR( pvalues, vlambda = vlambda, fdr_level = options.fdr )
        except ValueError, msg:
            E.warn( "%s: fdr could not be computed - no filtering: %s" % (filename, msg) )
            continue

        qvalues.plot( filename + "_diagnostics.png" )

        data = [ x[0] for x in zip( data, qvalues.mPassed ) if x[1] ]

def doFDR( options, args ):

    # read input
    annotators = []
    for filename in args:
        infile = open(filename,"r")
        annotators.append( readAnnotator( infile ))
        infile.close()

    do_filter = options.fdr_qvalue != None
    
    extra_headers = set()
    for data, fdr, synonyms, input_files in annotators:
        for key, value in input_files.iteritems():
            extra_headers.add( key )
    extra_headers=sorted(list(extra_headers))

    # note: id used to be file
    options.stdout.write( "id\tover\tcategory\tpvalue\tfold\tobserved\texpected\tci95low\tci95high\tstddev\tfdr\tqvalue\t%s\n" %\
                              "\t".join( extra_headers ))

    # apply filters and create diagnostic plots    
    for filename, vv in zip(args, annotators):
        data, fdr, synonyms, input_files = vv

        ninput = len(data)

        E.info( "processing %s with %i data points" % (filename, ninput) )
        no_fdr = False

        if options.fdr_method in ("annotator", "annotator-estimate"):
            pvalues = fdr.keys()
            pvalues.sort()
            pvalues.reverse()
            for pvalue in pvalues:
                try:
                    d = fdr[pvalue]["Significant"]
                except KeyError:
                    continue

                if d.mObserved == 0:
                    E.info("no data after fdr" )
                    break
                
                elif d.mAverage / d.mObserved < options.fdr_qvalue:
                    E.info("filtering with P-value of %f" % pvalue )
                    if do_filter: data = [ x for x in data if x.mPValue < pvalue ]
                    for d in data:
                        if d.mPValue < pvalue: 
                            d.mFDR = 1
                            d.mQValue = options.fdr_qvalue
                    break
                else:
                    E.warn( "fdr could not be computed - compute more samples (at P = %f, actual fdr=%f)" % (pvalue, d.mAverage / d.mObserved) )
                    no_fdr = True

        if options.fdr_method == "estimate" or (options.fdr_method == "annotator-estimate" and no_fdr):

            E.info( "estimating FDR from observed P-Values" )
            pvalues = [x.mPValue for x in data]
            vlambda = numpy.arange(0, max( pvalues), 0.05 )
            try:
                qvalues = Stats.doFDR( pvalues, vlambda = vlambda, fdr_level = options.fdr_qvalue )
            except ValueError, msg:
                E.warn( "fdr could not be computed - no output: %s" % msg)
                no_fdr = True
            else:
                for d, p, q in zip( data, qvalues.mPassed, qvalues.mQValues ):
                    if p: d.mFDR = 1
                    d.mQValue = q

                if do_filter:
                    data = [ x[0] for x in zip( data, qvalues.mPassed ) if x[1] ]

        if do_filter and no_fdr: data = []

        nremoved = ninput - len(data)

        E.info( "%s: %i data points left, %i removed" % (filename, len(data), nremoved) )

        extra_values = []
        for key in extra_headers:
            if key in input_files:
                extra_values.append( input_files[key] )
            else:
                extra_values.append( "" )

        extra_values = "\t".join( map(str, extra_values ))

        for d in data:
            if d.mFoldChange < 1:
                code = "-"
            else:
                code = "+"

            try:
                id = re.search( options.regex_id, filename ).groups()[0]
            except AttributeError:
                id = filename
            options.stdout.write( "%s\t%s\t%s\t%e\t%6.4f\t%f\t%f\t%f\t%f\t%f\t%i\t%e\t%s\n" % \
                                      (id,
                                       code,
                                       d.mAnnotation,
                                       d.mPValue,
                                       d.mFoldChange,
                                       d.mObserved,
                                       d.mExpected,
                                       d.mCI95[0],
                                       d.mCI95[1],
                                       d.mStdDev,
                                       d.mFDR,
                                       d.mQValue,
                                       extra_values) )

def main():
    parser = E.OptionParser( version = "%prog version: $Id: annotator2tsv.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-g", "--filename-gtf", dest="filename_gtf", type="string",
                      help="filename with gtf information. Used to map a segment to a gene name [default=%default]." )

    parser.add_option("-q", "--query", dest="queries", type="string", action="append",
                      help="term to query with [default=%default]." )

    parser.add_option("-m", "--method", dest="method", type="choice", 
                      choices=("query", "fdr", "fdr-table", "table" ),
                      help="methods to apply [default=%default]." )

    parser.add_option( "--output-format", dest="output_format", type="choice",
                       choices=("table", "gff"),
                       help="output format for query option [default=%default]" )

    parser.add_option( "--fdr-method", dest="fdr_method", type="choice",
                       choices=("annotator", "estimate", "annotator-estimate" ),
                       help="use fdr method. Use annotator, estimate or use annotator and estimate if not enough samples [default=%default].")

    parser.add_option( "--fdr-qvalue", dest="fdr_qvalue", type="float",
                       help="filter output by FDR qvalue (requires annotator output). [default=%default]")

    parser.add_option("-w", "--workspace", dest="filename_workspace", type="string",
                      help="filename with workspace information [default=%default]." )

    parser.add_option("-a", "--annotations", dest="filename_annotations", type="string",
                      help="filename with annotation information [default=%default]." )

    parser.add_option("--annotations-gtf", dest="filename_annotations_gtf", type="string",
                      help="filename with annotation information [default=%default]." )

    parser.add_option("-s", "--segments", dest="filename_segments", type="string",
                      help="filename with segment information [default=%default]." )

    parser.add_option( "--segments-gtf", dest="filename_segments_gtf", type="string",
                      help="filename with segment information [default=%default]." )

    parser.add_option( "--regex-id", dest="regex_id", type="string",
                      help="regular expression to extract id from filename [default=%default]." )

    parser.set_defaults(
        ignore_strand = False,
        filename_workspace = "intronic.workspace",
        filename_annotations = "go_territories.annotations",
        filename_segments = "test.segments",
        filename_segments_gtf = "FastDown.gtf",
        filename_annotations_gtf = "../data/tg1_territories.gff",
        regex_id = "(.*)",
        queries = [],
        method = "fdr",
        fdr_qvalue = None,
        fdr_method = "annotator",
        output_format = "table",
        )

    (options, args) = E.Start( parser )

    if options.method == "query":
        doQuery( options, args )
    elif options.method in ( "fdr", "fdr-table", "table"):
        doFDR( options, args )

    E.Stop()    

if __name__ == "__main__":
    sys.exit(main())


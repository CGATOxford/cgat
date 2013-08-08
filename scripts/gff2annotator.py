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
gff2annotator2tsv.py - convert from gff to annotator format
=======================================================

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

   python gff2annotator2tsv.py --help

Type::

   python gff2annotator2tsv.py --help

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
import time
import os
import shutil
import tempfile
import math
import itertools
import glob
import collections

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools
import CGATPipelines.PipelineEnrichment as PipelineEnrichment 

USAGE="""python %s [OPTIONS] < stdin > stdout

convert a gff file into annotator compatible regions. Depending on the option --section
this script will create:

   segments
      a segments file

   annotations
      a file with annotations. Multiple gff files can be provided. 

   annotations-genes
      if gtf files are provided (ending in .gtf), multiple subsets can be created using the 
      --subsets options. In this case a list of gene_ids for each subset is required.

   annotations-go 
      a file with annotations. Input is a gtf file and a map of gene identifiers
      to categories. Requires --input-filename-map.

   annotations-gff
      take all annotations in a gff file and create individual annotations for each feature 
      encountered. Multiple files will be aggregated if they contain the same feature.

   workspace
      a file with a workspace

""" % sys.argv[0]

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: gff2annotator2tsv.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"])

        
    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome."  )

    parser.add_option( "-f", "--features", dest="features", type="string", 
                       help="feature to collect [default=None]."  )

    parser.add_option( "-i", "--files", dest="files", action="append",
                       help="use multiple annotations [default=None]."  )

    parser.add_option(  "-a", "--annotations", dest="annotations", type="string", 
                       help="aggregate name for annotations if only single file is provided from STDIN [default=None]."  )

    parser.add_option(  "--input-filename-map", dest="input_filename_map", type="string", 
                       help="filename with a map of gene_ids to categories [default=None]."  )

    parser.add_option(  "--output-filename-synonyms", dest="output_filename_synonyms", type="string", 
                       help="output filename for synonyms. For workspace building, the gff source will be used as the id (instead of the contig) [default=None]."  )

    parser.add_option( "-m", "--max-length", dest="max_length", type="string", 
                       help="maximum segment length [default=None]."  )

    parser.add_option( "-s", "--section", dest="section", type="choice", 
                       choices=("segments", "annotations", "annotations-genes", "annotations-go", "workspace", "annotations-gff" ),
                       help="annotator section [default=None]."  )

    parser.add_option( "--subset", dest="subsets", type="string", action="append",
                       help="add filenames to delimit subsets within the gff files. The syntax is filename.gff,label,filename.ids [default=None]."  )

    parser.add_option( "--remove-regex", dest="remove_regex", type="string", 
                       help="regular expression of contigs to remove [default=None]."  )

    parser.set_defaults(
        genome_file = None,
        feature = None,
        section = "segments",
        annotations = "annotations",
        max_length = 100000,
        files = [],
        subsets = [],
        input_filename_map = None,
        output_filename_synonyms = None,
        input_format = "gff",
        remove_regex = None,
        )

    (options, args) = E.Start( parser )

    options.files += args
    if len(options.files) == 0: options.files.append("-")
    options.files = list( itertools.chain( *[ re.split( "[,; ]+", x) for x in options.files ] ) )

    if options.subsets:
        subsets = collections.defaultdict( list )
        for s in options.subsets: 
            filename_gff,label,filename_ids = s.split( "," )
            subsets[filename_gff].append( (label,filename_ids) )
        options.subsets = subsets

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    if options.section == "segments":
        prefix = "##Segs"
    elif options.section.startswith( "annotations" ):
        prefix = "##Id"
    elif options.section == "workspace":
        prefix = "##Work"
    else:
        raise ValueError("unknown section %s" % options.section)
        
    ninput, ncontigs, nsegments, ndiscarded = 0, 0, 0, 0

    if options.remove_regex:
        options.remove_regex = re.compile( options.remove_regex )

    if options.section in ("segments", "workspace"):

        iterator = GTF.iterator_filtered( GFF.iterator( options.stdin ),
                                          feature=options.feature )

        if options.output_filename_synonyms:
            outfile_synonyms = open(options.output_filename_synonyms, "w")
            with_records = True
        else:
            outfile_synonyms = None
            with_records = False

        intervals =GTF.readAsIntervals( iterator, with_records = with_records )
        ninput, nsegments, ndiscarded, ncontigs = \
            PipelineEnrichment.outputSegments( options.stdout,
                                               intervals,
                                               options.section,
                                               outfile_synonyms = outfile_synonyms,
                                               max_length = options.max_length,
                                               remove_regex = options.remove_regex )
            
        if outfile_synonyms:
            outfile_synonyms.close()

    elif options.section == "annotations-go":

        assert options.input_filename_map, "please supply option --input-filename-map" 

        iterator = GTF.iterator_filtered( GTF.iterator( options.stdin ),
                                          feature=options.feature )

        geneid2categories = IOTools.readMultiMap( open( options.input_filename_map, "r") )

        category2segments = collections.defaultdict( list )

        for contig, gffs in GTF.readAsIntervals( iterator, with_gene_id = True ).items():
            if options.remove_regex and options.remove_regex.search( contig ): continue
            
            for start, end, geneid in gffs:
                if geneid not in geneid2categories: continue
                for category in geneid2categories[geneid]:
                    category2segments[category].append(nsegments)

                options.stdout.write( "%s\t%i\t%s\t(%i,%i)\n" % (prefix, nsegments, contig, start, end ) )
                nsegments += 1                        
            
        for category, segments in category2segments.iteritems():
            options.stdout.write("##Ann\t%s\t%s\n" % (category, "\t".join( ["%i" % x for x in segments ] ) ) )
            E.info( "set %s annotated with %i segments" % (category, len(segments)) )

    elif options.section == "annotations":

        for filename in options.files:

            E.info( "adding filename %s" % filename )

            start = nsegments
            is_gtf = False

            if filename == "-":
                iterator = GTF.iterator_filtered( GFF.iterator( sys.stdin ),
                                                  feature=options.feature )
                filename = options.annotations
            elif filename.endswith(".gtf"):
                is_gtf = True
                with open( filename, "r") as infile:
                    iterator = GTF.iterator_filtered( GTF.iterator( infile ),
                                                      feature=options.feature )
                
            else:
                with open( filename, "r") as infile:
                    iterator = GTF.iterator_filtered( GFF.iterator( infile ),
                                                      feature=options.feature )
           
            E.debug("processing %s" % (filename))

            if not options.subsets or filename not in options.subsets:
                for contig, gffs in GTF.readAsIntervals( iterator ).items():
                    if options.remove_regex and options.remove_regex.search( contig ): continue

                    for x in gffs:
                        options.stdout.write( "%s\t%i\t%s\t(%i,%i)\n" % (prefix, nsegments, contig, x[0], x[1] ) )
                        nsegments += 1

                options.stdout.write("##Ann\t%s\t%s\n" % (filename, "\t".join( ["%i" % x for x in range(start, nsegments) ] ) ) )
                E.info( "set %s annotated with %i segments" % (filename, nsegments - start) )

            else:
                raise ValueError("don't know how to filter %s" % filename )

    elif options.section == "annotations-gff":

        for filename in options.files:
            if filename == "-":
                iterator = GTF.iterator( sys.stdin )
            else:
                iterator = GTF.iterator_filtered( GFF.iterator( open( filename, "r") ) )

            segments = collections.defaultdict( list )
            for gff in iterator:
                segments[":".join((gff.source,gff.feature))].append( (gff.contig,gff.start, gff.end) )
  
            feature2segments = {}

            for feature, s in segments.iteritems():
                s.sort()

                s1 = nsegments

                for contig, start, end in s:
                    if options.remove_regex and options.remove_regex.search( contig ): continue

                    options.stdout.write( "%s\t%i\t%s\t(%i,%i)\n" % (prefix, nsegments, contig, start, end ) )
                    nsegments += 1

                feature2segments[feature] = (s1, nsegments)
            
        for feature, id_range in feature2segments.iteritems():
            start, end = id_range
            options.stdout.write("##Ann\t%s\t%s\n" % (feature, "\t".join( ["%i" % x for x in xrange( start,end) ] ) ) )
            E.info( "set %s annotated with %i segments" % (feature, end-start) )

    elif options.section == "annotations-genes":

        for filename in options.files:

            E.info( "adding filename %s" % filename )

            start = nsegments

            assert filename.endswith(".gtf") or filename.endswith(".gtf.gz"), \
                "requiring .gtf files for gene list filtering, received %s" % filename

            infile = IOTools.openFile( filename )
            iterator = GTF.iterator_filtered( GTF.iterator( infile ),
                                              feature=options.feature )
                
            E.debug("processing %s" % (filename))
            
            if not options.subsets or filename not in options.subsets:
                ## output all
                for contig, gffs in GTF.readAsIntervals( iterator ).items():
                    if options.remove_regex and options.remove_regex.search( contig ): continue

                    for x in gffs:
                        options.stdout.write( "%s\t%i\t%s\t(%i,%i)\n" % (prefix, nsegments, contig, x[0],x[1] ) )
                        nsegments += 1

                options.stdout.write("##Ann\t%s\t%s\n" % (filename, "\t".join( ["%i" % x for x in range(start, nsegments) ] ) ) )
                E.info( "set %s annotated with %i segments" % (filename, nsegments - start) )

            else:
                ## create subsets
                E.debug("applying subsets for %s" % filename )
                geneid2label, label2segments = collections.defaultdict(list) , {}
                for label, filename_ids in options.subsets[filename]:
                    gene_ids = IOTools.readList( open(filename_ids, "r") )
                    for gene_id in gene_ids: geneid2label[gene_id].append( label )
                    label2segments[label] = []

                for contig, gffs in GTF.readAsIntervals( iterator, with_gene_id = True ).items():

                    if options.remove_regex and options.remove_regex.search( contig ): continue

                    for start, end, gene_id in gffs:
                        if gene_id not in geneid2label: continue
                        for label in geneid2label[gene_id]:
                            label2segments[label].append(nsegments)
                            
                        options.stdout.write( "%s\t%i\t%s\t(%i,%i)\n" % (prefix, nsegments, contig, start, end ) )
                        nsegments += 1                        
                        
                for label, segments in label2segments.iteritems():
                    options.stdout.write("##Ann\t%s\t%s\n" % (label, "\t".join( ["%i" % x for x in segments ] ) ) )
                    E.info( "set %s (%s) annotated with %i segments" % (label, filename, len(segments)) )

    E.info( "ninput=%i, ncontigs=%i, nsegments=%i, ndiscarded=%i" % (ninput, ncontigs, nsegments, ndiscarded))

    E.Stop()




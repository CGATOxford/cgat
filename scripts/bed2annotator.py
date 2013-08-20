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
bed2annotator.py - convert bed to annotator format
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script converts a bed file into annotator compatible regions. Depending on the option --section
this script will create:

   segments
      a segments file

   annotations
      a file with annotations. Each bed track is a separate annotation.

   workspace
      a file with a workspace

Usage
-----

Example::

   python bed2annotator2tsv.py --help

Type::

   python bed2annotator2tsv.py --help

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
import CGAT.Bed as Bed
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: bed2annotator2tsv.py 2885 2010-04-07 08:46:50Z andreas $", 
                                    usage = globals()["__doc__"])
        
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

    parser.add_option( "-l", "--max-length", dest="max_length", type="string", 
                       help="maximum segment length [default=None]."  )

    parser.add_option( "-m", "--merge", dest="merge", action="store_true",
                       help="merge overlapping bed segments [default=%default]."  )

    parser.add_option( "-s", "--section", dest="section", type="choice", 
                       choices=("segments", "annotations", "workspace" ),
                       help="annotator section [default=None]."  )

    parser.add_option( "--subset", dest="subsets", type="string", action="append",
                       help="add filenames to delimit subsets within the gff files. The syntax is filename.gff,label,filename.ids [default=None]."  )

    parser.set_defaults(
        genome_file = None,
        feature = None,
        remove_random = True,
        section = "segments",
        annotations = "annotations",
        max_length = 100000,
        files = [],
        subsets = [],
        input_filename_map = None,
        merge = False,
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
    elif options.section == "annotations":
        prefix = "##Id"
    elif options.section == "workspace":
        prefix = "##Work"
    else:
        raise ValueError("unknown section %s" % options.section)

    if options.max_length:
        max_length = options.max_length
    else:
        max_length = 0

    ninput, ntracks, ncontigs, nsegments, ndiscarded = 0, 0, 0, 0, 0

    if options.section in ("annotations"):
        contigs = set()
        it = itertools.groupby( Bed.iterator( options.stdin ), key=lambda x: x.track["name"])

        map_track2segments = {}
        for track, beds in it:
            ntracks += 1
            map_track2segments[track] = []
            first_segment = nsegments

            beds = list(beds)

            if options.merge: beds = Bed.merge( beds )

            for bed in beds:
                contig, start, end = bed.contig, bed.start, bed.end

                if options.remove_random and "random" in contig: continue
                
                if max_length > 0 and end - start > max_length:
                    ndiscarded += 1
                    continue

                contigs.add(contig)
                map_track2segments[track].append( nsegments )
                options.stdout.write( "%s\t%i\t%s\t(%i,%i)\n" % (prefix, nsegments, contig, start, end ) )
                nsegments += 1
                
            options.stdout.write("##Ann\t%s\t%s\n" % (track, "\t".join( ["%i" % x for x in range(first_segment, nsegments) ] ) ) )
            E.info( "track %s: annotated with %i segments" % (track, nsegments - first_segment) )

        ncontigs = len(contigs)
        E.info( "ninput=%i, ntracks=%i, ncontigs=%i, nsegments=%i, ndiscarded=%i" % (ninput, ntracks, ncontigs, nsegments, ndiscarded))

    E.Stop()




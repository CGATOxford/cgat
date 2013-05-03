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
windows2gff.py - create genomic windows/tiles
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Given indexed genome sequences, output a collection of windows.

Usage
-----

Example::

   python windows2gff.py --fixed-windows=500,250

creates windows of size 500 bp in steps of 250 bp.

Type::

   python windows2gff.py --help

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

import CGAT.Experiment as E
import CGAT.GFF as GFF
import CGAT.IndexedFasta as IndexedFasta

def getFixedWidthWindows( map_contig2size, options ):
    """return a list of fixed contig sizes."""
    
    assert( options.fixed_width_windows )
    v = map(int, options.fixed_width_windows.split(","))
    if len(v) == 2:
        window_size, window_increment = v
    elif len(v) == 1:
        window_size, window_increment = v[0], v[0]
    else:
        raise ValueError( "could not parse window size '%s': should be size[,increment]" % options.fixed_width_windows )

    for contig, size in map_contig2size.items():
        E.info("processing %s" % contig)
        for x in range(0, size, window_increment):
            if x + window_size > size: continue
            gff = GFF.Entry()
            gff.feature = "window"
            gff.source = "window"
            gff.contig = contig
            gff.start = x
            gff.end = min(size, x + window_size)
            yield gff

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: windows2gff.py 2861 2010-02-23 17:36:32Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome (indexed) for getting contig sizes."  )

    parser.add_option( "-e", "--gff-file", dest="gff_file", type="string",
                       help="gff file to use for getting contig sizes."  )

    parser.add_option( "-f", "--fixed-width-windows=", dest="fixed_width_windows", type="string",
                       help="fixed width windows. Supply the window size as a parameter. Optionally supply an offset." )

    parser.add_option( "-o", "--output-format=", dest="output_format", type="choice",
                       choices=("gff", "bed"),
                       help="output format [%default]" )

    parser.set_defaults(
        genome_file = None,
        fixed_windows = None,
        features = [],
        output_format = "gff",
        )

    (options, args) = E.Start( parser )

    map_contig2size = {}

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
        map_contig2size = fasta.getContigSizes( with_synonyms = False )
    else:
        fasta = None

    if options.gff_file:
        infile = open( options.gff_file, "r" )
        gff = GFF.readFromFile( infile )   
        infile.close()
        for g in gff:
            try:
                map_contig2size[g.mName] = max( map_contig2size[g.mName], g.end )
            except ValueError:
                map_contig2size[g.mName] = g.end 
            
    else:
        gff = None

    if map_contig2size == None:
        raise ValueError( "no source of contig sizes supplied" )

    if options.fixed_width_windows:
        windows = getFixedWidthWindows( map_contig2size, options )
    
    noutput = 0
    if options.output_format == "gff":
        for g in windows:
            options.stdout.write( str(g) + "\n" )
    
    elif options.output_format == "bed":
        for g in windows:
            noutput += 1
            options.stdout.write("\t".join(map(str, 
                                               (g.contig, 
                                               g.start,
                                               g.end )) ) + "\n" )

    E.info( "noutput=%i\n" % (noutput ))

    E.Stop()

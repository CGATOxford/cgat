################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
beds2counts.py - compute overlap stats between multiple bed files
==================================================================

:Author: Nick Ilott 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes multiple bed files e.g. from multiple samples from the same experiment. It 
asseses the overlap between samples and outputs a count for each merged interval corresponding
to the number of samples that a particular interval was found in.

Writes the output to stdout.

Usage
-----

Example::

   python beds2counts.py file1.bed file2.bed > output.bed

Type::

   python beds2counts.py --help

for command line help.

Documentation
-------------

Code
----

'''
import tempfile
import os
import sys
import re
import optparse
import pybedtools
import CGAT.Experiment as E
import CGAT.Bed as Bed
import collections
import CGAT.IOTools as IOTools
import CGAT.IndexedGenome as IndexedGenome

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-i", "--infiles", dest="infiles", type="string",
                      metavar = "bed",
                      action="append", help="supply list of bed files"  )

    parser.set_defaults( infiles = [] )
    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    options.infiles.extend( args )
    if len(options.infiles) == 0:
        raise ValueError( 'please provide at least 1 bed file' )

    E.info("concatenating bed files")
    # concatenate the list of files
    tmp = tempfile.NamedTemporaryFile(delete = False)
    tmp_merge = tempfile.NamedTemporaryFile(delete = False)
    infs = options.infiles
    for inf in infs:
        for bed in Bed.iterator(IOTools.openFile(inf)):
            tmp.write("%s\n" % bed)
    tmp.close()
    
    E.info("merging bed entries")
    # merge the bed entries in the file
    name = tmp.name
    tmp_bed = pybedtools.BedTool(name)
    tmp_bed.merge().saveas(tmp_merge.name)
    tmp_merge.close()

    E.info("indexing bed entries")
    # index the bed entries                                                                                                                                                                                                                  
    merged = IndexedGenome.Simple()
    for bed in Bed.iterator(IOTools.openFile(tmp_merge.name)):
        merged.add(bed.contig, bed.start, bed.end)

    counts = collections.defaultdict(int)
    # list of samples                                                                                                                                                                                                                        
    samples = options.infiles
    
    E.info("counting no. samples overlapping each interval")
    for sample in samples:
        found = set()
        for bed in Bed.iterator(IOTools.openFile(sample)):
            if merged.contains(bed.contig, bed.start, bed.end):
                key =  [bed.contig] + [x for x in merged.get(bed.contig, bed.start, bed.end)]
                key = (key[0], key[1][0], key[1][1])
                if key in found: continue
                found.add(key)

                # tuple of interval description as key - (contig, start, end)                                                                                                                                                                
                counts[key] += 1

    # open outfile                                                                                                                                                                                                                           
    options.stdout.write("contig\tstart\tend\tcount\n")
    
    E.info("outputting result")
    for interval, count in counts.iteritems():
        options.stdout.write("\t".join(map(str,interval)) + "\t" + str(count) + "\n")

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


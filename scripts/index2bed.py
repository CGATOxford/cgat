################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: index2gff.py 2880 2010-04-07 08:44:13Z andreas $
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
index2bed.py - convert indexed fasta file to bed file
=====================================================

:Author: Andreas Heger
:Release: $Id: index2gff.py 2880 2010-04-07 08:44:13Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os
import sys
import re
import optparse

import CGAT.IndexedFasta as IndexedFasta
import CGAT.Bed as Bed
import CGAT.Experiment as E

def getFixedWidthWindows( map_contig2size, window_size ):
    """return a list of fixed contig sizes."""
    

    for contig, size in map_contig2size.items():
        E.info("processing %s" % contig)
        for x in range(0, size, window_increment):
            if x + window_size > size: continue
            gff = GTF.Entry()
            gff.feature = "window"
            gff.source = "window"
            gff.contig = contig
            gff.start = x
            gff.end = min(size, x + window_size)
            yield gff


def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: index2gff.py 2880 2010-04-07 08:44:13Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option( "--remove-regex", dest="remove_regex", type="string", 
                       help="regular expression of contigs to remove [default=None]."  )

    parser.add_option( "-e", "--gff-file", dest="gff_file", type="string",
                       help="gff file to use for getting contig sizes."  )

    parser.add_option( "-f", "--fixed-width-windows=", dest="fixed_width_windows", type="string",
                       help="fixed width windows. Supply the window size as a parameter. Optionally supply an offset." )

    parser.set_defaults(
        genome_file = None,
        remove_regex = None,
        fixed_windows = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if options.remove_regex:
        remove_regex = re.compile( options.remove_regex )
    else:
        remove_regex = None

    if options.fixed_width_windows:
        v = map(int, options.fixed_width_windows.split(","))
        if len(v) == 2:
            window_size, window_increment = v
        elif len(v) == 1:
            window_size, window_increment = v[0], v[0]
        else:
            raise ValueError( "could not parse window size '%s': should be size[,increment]" % options.fixed_width_windows )

    if options.gff_file:
        infile = open( options.gff_file, "r" )
        gff = GTF.readFromFile( infile )   
        infile.close()
        for g in gff:
            try:
                map_contig2size[g.mName] = max( map_contig2size[g.mName], g.end )
            except ValueError:
                map_contig2size[g.mName] = g.end 
            
    else:
        gff = None

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
        map_contig2size = fasta.getContigSizes( with_synonyms = False )
    else:
        fasta = None

    if map_contig2size == None:
        raise ValueError( "no source of contig sizes supplied" )

    ## do sth
    counter = E.Counter()

    for contig, size in map_contig2size.items():
        size = int(size)
        counter.input += 1

        if remove_regex and remove_regex.search(contig): 
            counter.skipped += 1
            continue

        if options.fixed_width_windows:
            for x in range(0, size, window_increment):
                if x + window_size > size: continue
                options.stdout.write( "%s\t%i\t%i" % (contig, x, min(size, x + window_size)))
                counter.windows += 1
        else:
            options.stdout.write( "%s\t%i\t%i" % (contig, 0, size))
            counter.windows += 1

        counter.output += 1

    E.info( str(counter) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


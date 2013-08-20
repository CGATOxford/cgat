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
bed2stats.py - summary of bed file contents
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes a bed-formatted file as input and outputs the number
of intervals and bases in the bed file. Counts can be computed per
contig or per track in the bed file.

Usage
-----

Example::

   python bed2table.py --help

Type::

   python bed2table.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import collections
import CGAT.Bed as Bed
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta

class Counter:

    headers = ["ncontigs", "nintervals", "nbases" ]

    def __init__(self):
        self.intervals_per_contig = collections.defaultdict(int)
        self.bases_per_contig = collections.defaultdict(int)
        
    def add( self, bed ):
        self.intervals_per_contig[bed.contig] += 1
        self.bases_per_contig[bed.contig] += bed.end - bed.start

    def __str__(self):
        return "%i\t%i\t%i" % (len( self.intervals_per_contig),
                               sum( self.intervals_per_contig.values()),
                               sum( self.bases_per_contig.values())) 


##------------------------------------------------------------
def main( argv = None ):

    if argv == None: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: gtf2table.py 2888 2010-04-07 08:48:36Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.add_option("-n", "--per-name", dest="per_name", action="store_true",
                      help="compute counts per name [default=%default]."  )

    parser.add_option("-c", "--per-contig", dest="per_contig", action="store_true",
                      help="compute counts per contig [default=%default]."  )

    parser.add_option("-t", "--per-track", dest="per_track", action="store_true",
                      help="compute counts per track [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        per_name = False,
        per_track = False,
        )

    (options, args) = E.Start( parser, argv )

    # get files
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    counts = collections.defaultdict( Counter )

    if options.per_track:
        keyf = lambda x: x.track
    elif options.per_name:
        keyf = lambda x: x.name
    elif options.per_contig:
        keyf = lambda x: x.contig
    else:
        keyf = lambda x: "all"

    for bed in Bed.iterator(options.stdin):
        counts[keyf(bed)].add( bed )

    outf = options.stdout

    key = "track"
    outf.write( "%s\t%s\n" % (key, "\t".join( Counter.headers) ))

    for key, count in counts.iteritems():
        outf.write( "%s\t%s\n" % ( key, str(count)) )
        
    E.Stop()

if __name__ == '__main__':
    sys.exit(main(sys.argv))

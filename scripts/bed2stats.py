'''
bed2stats.py - summary of bed file contents
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals Summary BED

Purpose
-------

This script takes a bed-formatted file as input and outputs the number
of intervals and bases in the bed file. Counts can be computed per
contig or per track in the bed file.

Note that a count of bases usually makes only sense if the intervals
submitted are non-overlapping.

Usage
-----

Example::

   python bed2table.py --help

Type::

   python bed2table.py --help

for command line help.

Command line options
--------------------

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
    headers_percent = ["ncontigs", "nintervals", "nbases", "pbases" ]

    def __init__(self):
        self.intervals_per_contig = collections.defaultdict(int)
        self.bases_per_contig = collections.defaultdict(int)
        self.size = None

    def setSize( self, size ):
        self.size = size
        
    def add( self, bed ):
        self.intervals_per_contig[bed.contig] += 1
        self.bases_per_contig[bed.contig] += bed.end - bed.start

    def __str__(self):
        bases = sum( self.bases_per_contig.values())
        if self.size == None:
            return "%i\t%i\t%i" % (len( self.intervals_per_contig),
                                   sum( self.intervals_per_contig.values()),
                                   bases,
                                   ) 
        else:
            return "%i\t%i\t%i\t%5.2f" % (len( self.intervals_per_contig),
                                          sum( self.intervals_per_contig.values()),
                                          sum( self.bases_per_contig.values()),
                                          100.0 * bases / self.size
                                          ) 


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

    parser.add_option("-p", "--add-percent", dest="add_percent", action="store_true",
                      help="add percentages [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        per_name = False,
        per_track = False,
        add_percent = False,
        )

    (options, args) = E.Start( parser, argv )

    # get files
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        if options.add_percent:
            raise ValueError("--add-percent option requires --genome-file")
        fasta = None

    if options.add_percent and not options.per_contig:
        raise NotImplementedError("--add-percent option requires --per-contig")

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
    if options.add_percent:
        outf.write( "%s\t%s\n" % (key, "\t".join( Counter.headers_percent) ))
    else:
        outf.write( "%s\t%s\n" % (key, "\t".join( Counter.headers) ))

    total_bases = 0
    for key, count in counts.iteritems():
        if options.add_percent:
            total_bases += fasta.getLength( key )
            count.setSize( fasta.getLength( key ) )

        outf.write( "%s\t%s\n" % ( key, str(count)) )
        
        
    

    E.Stop()

if __name__ == '__main__':
    sys.exit(main(sys.argv))

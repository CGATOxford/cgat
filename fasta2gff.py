################################################################################
#   Gene prediction pipeline 
#
#   $Id: fasta2gff.py 2861 2010-02-23 17:36:32Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
import os, sys, string, re, getopt, tempfile, time, optparse, math, glob

import Experiment as E
import IndexedFasta
import GFF, GTF

USAGE="""python %s [OPTIONS] 

extract fragments from a genome.

Version: $Id: fasta2gff.py 2861 2010-02-23 17:36:32Z andreas $
""" % sys.argv[0]

def writeHeader( outfile ):
    outfile.write( "\t".join( ("contig",
                               "nresidues",
                               "ngaps",
                               "nseqregions",                                      
                               "ngapregions", 
                               "nA", "nC", "nG", "nT",
                               "nN", "nX", "nO" ) ) + "\n" )
    
    

def main():

    parser = optparse.OptionParser( version = "%prog version: $Id: fasta2gff.py 2861 2010-02-23 17:36:32Z andreas $")

    parser.add_option( "-g", "--genome-file", dest="genome_file", type="string",
                       help="filename with genome."  )

    parser.add_option("-a", "--as-gtf", dest="as_gtf", action="store_true",
                      help="output as gtf."  )

    parser.add_option( "-f", "--fragment-size", dest="fragment_size", type="int",
                       help="fixed size of fragments [default=%default]."  )

    parser.add_option( "-s", "--sample-size", dest="sample_size", type="int",
                       help="fixed size of fragments."  )

    parser.set_defaults(
        as_gtf = False,
        genome_file = None,
        fragment_size = 1000,
        sample_size = 10000,
        pattern_id = "%08i",
        )

    (options, args) = E.Start( parser )

    fasta = IndexedFasta.IndexedFasta( options.genome_file )
    contigs = fasta.getContigSizes()

    if options.as_gtf:
        entry = GTF.Entry()
    else:
        entry = GFF.Entry()

    n = 0
    entry.feature = "exon"
    entry.source = "random"

    for x in range( options.sample_size):

        entry.contig, entry.strand, entry.start, entry.end = fasta.getRandomCoordinates( options.fragment_size )
        if entry.strand == "-":
            l = contigs[entry.contig]
            entry.start, entry.end = l-entry.end, l-entry.start

        if options.as_gtf:
            entry.gene_id = options.pattern_id % n
            entry.transcript_id = entry.gene_id

        options.stdout.write( str(entry) + "\n" )
        n += 1

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

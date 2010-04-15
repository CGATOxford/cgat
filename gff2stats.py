################################################################################
#   Gene prediction pipeline 
#
#   $Id: gff2stats.py 2781 2009-09-10 11:33:14Z andreas $
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
import sys, string, re, optparse, collections

USAGE="""python %s [OPTIONS] input1 input2

compute simple statistics on the size of features and the 
distance between features in a gff file.

Methods:
  hist: histogram of distance/sizes
  stats: descriptive statistics of distance/sizes
  overlaps: output overlapping features

Version: $Id: gff2stats.py 2781 2009-09-10 11:33:14Z andreas $
""" % sys.argv[0]

import Experiment as E
import GFF, GTF

##------------------------------------------------------------------------
def main( argv = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id: gff2stats.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("--is-gtf", dest="is_gtf", action="store_true",
                      help="input is gtf.")

    parser.set_defaults(
        is_gtf = False,
        )

    (options, args) = E.Start( parser, add_output_options = True )

    is_gtf = options.is_gtf

    if is_gtf:
        iterator = GTF.iterator
        counts_gene_ids = collections.defaultdict( int )
        counts_transcript_ids = collections.defaultdict( int )
    else:
        iterator = GFF.iterator

    counts_contigs = collections.defaultdict( int )
    counts_strands = collections.defaultdict( int )
    counts_features = collections.defaultdict( int )
    counts_sources = collections.defaultdict( int )
        
    for entry in iterator( sys.stdin ):
        counts_contigs[entry.contig] += 1
        counts_features[entry.feature] += 1
        counts_sources[entry.source] += 1
        counts_strands[entry.strand] += 1
        if is_gtf:
            counts_gene_ids[entry.gene_id] += 1
            counts_transcript_ids[entry.transcript_id] += 1
    
    def printStats( title, stats ):
        outfile = options.stdout
        outfile.write( "%s: %i\n" % (title, len(stats) ) )

    printStats( "contigs", counts_contigs )
    printStats( "features", counts_features )
    printStats( "sources", counts_sources )
    printStats( "strands", counts_strands )
    if is_gtf:
        printStats( "gene_ids", counts_gene_ids )
        printStats( "transcript_ids", counts_transcript_ids )

    E.Stop()

if __name__=="__main__":
    sys.exit( main( sys.argv ) )
    

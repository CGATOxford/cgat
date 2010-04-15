################################################################################
#   Gene prediction pipeline 
#
#   $Id: fasta2spliced.py 2861 2010-02-23 17:36:32Z andreas $
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
import sys, string, re, optparse

USAGE="""python %s [OPTIONS] < input.fasta > output.fasta

Version: $Id: fasta2spliced.py 2861 2010-02-23 17:36:32Z andreas $
""" % sys.argv[0]

import Experiment
import IndexedFasta
import Genomics
import GFF

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: fasta2spliced.py 2861 2010-02-23 17:36:32Z andreas $")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("-r", "--filename-regions", dest="filename_regions", type="string",
                      help="filename with region information in GFF format."  )

    parser.add_option( "-p", "--output-filename-pattern", dest="output_filename_pattern", type="string" ,
                       help="OUTPUT filename pattern for additional data [%default].")

    parser.add_option( "--joined", dest="joined", action="store_true",
                       help="output mode. If joined, all is output in one big chromosome. Otherwise, each are single fragments [%default].")

    parser.add_option( "--only-first", dest="only_first", action="store_true",
                       help="only output the first possible splice site [%default].")

    parser.set_defaults(
        genome_file = "genome",
        filename_regions = None,
        output_format = "%08i",
        output_filename_pattern = "%s",
        methods = [],
        splice_pairs = ( ("GT", "AG"), ),
        min_intron_size = 30,
        max_intron_size = 25000,
        search_area = 5, # 10
        read_length = 32,
        only_first = False,
        joined = False,
        max_join_length = 1000000, # 100000000
        format_id = "seg%05i",
        )

    (options, args) = Experiment.Start( parser )

    genome = IndexedFasta.IndexedFasta( options.genome_file )

    assert options.filename_regions != None, "please supply a gff formatted filename with regions"

    regions = GFF.readAsIntervals( GFF.iterator( open(options.filename_regions, "r" ) ) )

    # build pairs for complement
    reverse_splice_pairs = []
    forward_splice_pairs = options.splice_pairs
    left_tokens, right_tokens = {}, {}
    x = 0
    for a,b in forward_splice_pairs:
        assert len(a) == 2, "only two-residue patterns allowed"
        assert len(b) == 2, "only two-residue patterns allowed"

        ca, cb = Genomics.complement( a ), Genomics.complement( b ) 
        reverse_splice_pairs.append( (b,a) )
        left_tokens[a] = x
        left_tokens[cb] = x+1
        right_tokens[b] = x
        right_tokens[ca] = x+1
        x += 2

    search_area = options.search_area
    read_length = options.read_length
    joined = options.joined

    ninput, noutput = 0, 0

    if joined:
        outfile_coordinates = open( options.output_filename_pattern % "coords", "w" )
        outfile_coordinates.write( "segment\tpos\tcontig\t5start\t3start\n" )
        out_contig = 1
        options.stdout.write( ">%s\n" % (options.format_id % out_contig ))
        nbases = 0
        separator = "N" * read_length
        lseparator = len(separator)

    contig_sizes = genome.getContigSizes()
    # collect possible start/end points of introns
    for contig, lcontig in contig_sizes.items():

        ninput += 1

        nintrons = 0
        if contig not in regions:
            if options.loglevel >= 2:
                options.stdlog.write("# skipped %s - no intervals defined\n" % (contig))
            continue

        sequence = genome.getSequence( contig, as_array = True )

        if options.loglevel >= 2:
            options.stdlog.write("# processing %s of length %i\n" % (contig, len(sequence)))

        regions[contig].sort()
        
        left_positions, right_positions = [], []

        def addPositions( start, end, tokens, positions, forward = True, first = False ):

            area = sequence[start:end].upper()
            if forward:
                for x in range(len(area)-1):
                    t = area[x:x+2]
                    if t in tokens: 
                        positions.append( (start+x,tokens[t] ) )
                        if first: return True
                    
            else:
                for x in range(len(area)-2,-1,-1):
                    t = area[x:x+2]
                    if t in tokens: 
                        positions.append( (start+x,tokens[t] ) )
                        if first: return True
            return False

        intron_start = regions[contig][0][1]
        for exon_start,exon_end in regions[contig][1:]:
            
            intron_end = exon_start
            if options.only_first:
                if not addPositions( intron_start, intron_start+search_area, left_tokens, left_positions, forward=True, first = True ):
                    addPositions( intron_start-search_area, intron_start, left_tokens, left_positions, forward=False, first = True )

                if not addPositions( intron_end-search_area, intron_end, right_tokens, right_positions, forward=False, first = True ):
                    addPositions( intron_end, intron_end+search_area, right_tokens, right_positions, forward=True, first = True )

            else:
                addPositions( intron_start-search_area, intron_start+search_area, left_tokens, left_positions, forward=True, first = False )
                addPositions( intron_end-search_area, intron_end+search_area, right_tokens, right_positions, forward=True, first = False )
            intron_start = exon_end

        if options.loglevel >= 2:
            options.stdlog.write("# %s: left=%i, right=%i\n" % (contig, len(left_positions), len(right_positions) ))
        
        # build possible introns
        #
        # iterate over left positions and collect right positions within a radius
        # given by min_intron_size and max_intron_size.
        # left_positions and right_positions are sorted
        ri, mr = 0, len(right_positions)

        for l,t in left_positions:
            lower_bound, upper_bound = l + options.min_intron_size, l + options.max_intron_size
            while ri < mr and right_positions[ri][0] < lower_bound: ri += 1
            rri = ri

            while rri < mr and right_positions[rri][0] < upper_bound:
                if right_positions[rri][1] == t:
                    # positions are start/end of splice motif
                    # so add two on the right side
                    r = right_positions[rri][0]+2
                    lmin = max(0, l-read_length )
                    rmax = min( lcontig, r + read_length )

                    if options.loglevel >= 3:
                        options.stdlog.write("# adding intron on %s: l=%i, r=%i, t=%i, %s %s %s %s\n" %\
                                                 (contig, l, r, t, 
                                                  sequence[lmin:l], 
                                                  sequence[l:l+2], 
                                                  sequence[r-2:r], 
                                                  sequence[r:rmax] ) )
                        
                    if joined:
                        outfile_coordinates.write("%s\t%i\t%s\t%i\t%i\n" % (options.format_id % out_contig, nbases, contig, lmin, r ) )

                        s = sequence[lmin:l] + sequence[r:rmax]
                        options.stdout.write( "%s\n%s\n" % (s, separator ) )
                                              
                        nbases += len(s) + lseparator

                        if nbases > options.max_join_length:
                            nbases = 0
                            out_contig += 1
                            options.stdout.write( ">%s\n" % (options.format_id % out_contig ) )

                    else:
                        options.stdout.write( ">%s_%i_%i\n%s%s\n" % (contig, lmin, r,     
                                                                     sequence[lmin:l], 
                                                                     sequence[r:rmax] ) )

                    
                    nintrons += 1
                    noutput += 1
                rri += 1

        if options.loglevel >= 1:
            options.stdlog.write( "# contig %s: %i introns\n" % (contig, nintrons))
        
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i\n" % (ninput, noutput) )

    Experiment.Stop()

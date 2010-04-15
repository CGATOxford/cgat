####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <andreas.heger@helsinki.fi>
##
## $Id: blat2gff.py 2781 2009-09-10 11:33:14Z andreas $
##
##
####
####

USAGE="""python blat2gff.py [OPTIONS] < input > output

convert blat output to gff output.
"""

import sys, re, string, optparse, time, os, glob

import Experiment as E
import IOTools
import Blat
import GFF, GTF

def main():

    parser = optparse.OptionParser( version = "%prog version: $Id: blat2gff.py 2781 2009-09-10 11:33:14Z andreas $", usage=USAGE )

    parser.add_option("-a", "--as-gtf", dest="as_gtf", action="store_true",
                      help="output as gtf."  )

    parser.add_option("-s", "--filename-strand", dest="filename_strand", type="string",
                      help="set strand information according to file [default=%DEFAULT]."  )

    parser.set_defaults( as_gtf = False,
                         filename_strand = None,
                         test = None )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    ####################################
    if options.filename_strand:
        map_id2strand = IOTools.readMap( open(options.filename_strand,"r")) 
    else:
        map_id2strand = {}

    iterator = Blat.BlatIterator( sys.stdin )

    ninput, noutput, nskipped = 0, 0, 0

    if options.as_gtf:
        gff = GTF.Entry()
    else:
        gff = GFF.Entry()

    gff.source = "blat"
    gff.feature = "exon"

    ids = {}

    while 1:
        
        if options.test and ninput >= options.test:
            break

        match = iterator.next()
        
        if match == None: break

        ninput += 1

        if match.mQueryId not in ids:
            ids[match.mQueryId] = 1
            id = match.mQueryId
        else:
            id = match.mQueryId + ":%i" % ids[match.mQueryId]
            ids[match.mQueryId] += 1

        if options.as_gtf:
            gff.contig = match.mSbjctId
            gff.gene_id = id
            gff.transcript_id = id
        else:
            gff.contig = match.mSbjctId
            gff.clearAttributes()
            gff.addAttribute( "gene_id", id )

        if id in map_id2strand:
            gff.strand = map_id2strand[id]
        else:
            gff.strand = match.strand

        for qstart, sstart, size in match.getBlocks():
            
            gff.start = sstart
            gff.end = sstart + size
            options.stdout.write( str( gff ) + "\n" )
        
        noutput += 1

    E.info("ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped) )

    E.Stop()

if __name__ == '__main__':
    sys.exit(main())

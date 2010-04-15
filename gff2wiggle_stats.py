####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: gff2wiggle_stats.py 2781 2009-09-10 11:33:14Z andreas $
##
##
####
####

USAGE="""python blat2assembly.py [OPTIONS] > output

Assemble blat matches overlapping on sbjct.
"""

import sys, re, string, optparse, time, os, glob

import Experiment
import Stats
import GFF, GTF

import Wiggle
import alignlib

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: gff2wiggle_stats.py 2781 2009-09-10 11:33:14Z andreas $", usage=USAGE )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option("--wiggle-files", dest="wiggle_files", type="string",
                      help="glob expression for wiggle files [%default]."  )

    parser.add_option("--prefix", dest="prefix", type="string",
                      help="prefix to add to contig names before lookup [%default]."  )

    parser.add_option( "-a", "--as-gtf", dest="as_gtf", action="store_true",
                       help="input is in gtf format.")

    parser.add_option( "--test", dest="test", type="int",
                       help="test - stop after # rows of parsing [%default]." )

    parser.add_option( "--with-values", dest="with_values", action="store_true",
                       help="output values in last column [%default]." )

    parser.set_defaults( wiggle_files = "*.data.bz2",
                         prefix = "",
                         with_values = False,
                         test = None )

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    # open indexed access to wiggles
    wiggle_files = glob.glob( options.wiggle_files )
    if not wiggle_files:
        raise IOError( "could not find wiggle files with '%s'" % options.wiggle_files )

    index = Wiggle.WiggleMultiIndexedAccess( wiggle_files,
                                             keep_open = True,
                                             use_cache = False )

    if options.as_gtf:
        iterator = GTF.flat_gene_iterator( GTF.iterator( sys.stdin) )
        id = "gene_id"
    else:
        iterator = GFF.chunk_iterator( GTF.iterator( sys.stdin ) )
        id = "query"

    ninput, noutput, nskipped = 0, 0, 0

    options.stdout.write( "%s\tnali\t%s" % (id, "\t".join( Stats.DistributionalParameters().getHeaders() )))
    if options.with_values:
        options.stdout.write( "\tvalues" )
    options.stdout.write("\n" )

    for gffs in iterator:

        if options.test and ninput >= options.test:
            break

        ninput += 1
        
        if options.loglevel >= 2:
            for gff in gffs:
                options.stdlog.write( str(gff) + "\n" )

        map_genome2query = alignlib.makeAlignmentBlocks()
        for gff in gffs:
            map_genome2query.addDiagonal( gff.start, gff.end, 0 )

        data = index.get( options.prefix + gff.contig,
                          map_genome2query.getRowFrom(),
                          map_genome2query.getRowTo() )

        values = []
        for x, vv in data:
            for v in vv:
                if map_genome2query.mapRowToCol( x ) >= 0:
                    values.append( v )
                x += 1
        if len(values) == 0:
            nskipped += 1
            continue
            
        noutput += 1

        if options.loglevel >= 2:
            options.stdlog.write( "# %s\n" % ",".join( [ "%5.3f" % v for v in values] ))

        s = Stats.DistributionalParameters( values )
        if options.as_gtf:
            id = gff.gene_id
        else:
            id = gff.mAttributes
        options.stdout.write( "%s\t%i\t%s" % ( id,
                                               map_genome2query.getNumAligned(),
                                               str( s ) ) )
        
        if options.with_values:
            options.stdout.write( "\t%s" % ",".join( [ "%5.3f" % v for v in values] ))
            
        options.stdout.write("\n" )    

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped) )

    Experiment.Stop()

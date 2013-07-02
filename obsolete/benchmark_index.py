################################################################################
#
#   Gene prediction pipeline 
#
#   $Id: benchmark_index.py 2782 2009-09-10 11:40:29Z andreas $
#
#   Copyright (C) 2007 Andreas Heger
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
#
#################################################################################
"""Subroutines for working on I/O of large genomic files.
"""

USAGE="""python bnechmark_index.py [options] file [ files ]

benchmarking of indexing methods for fasta files.
"""

import tempfile
import timeit
import sys
import os
from CGAT.IndexedFasta import *
import CGAT.Stats as Stats

if __name__ == "__main__":

    import Experiment

    parser = optparse.OptionParser( version = "%prog version: $Id: benchmark_index.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.set_defaults(
        methods = [ "lzo", "dictzip", "zlib", "gzip" ],
        num_iterations = 10,
        benchmark_num_iterations = 100000,
#        fragment_sizes = (100, 1000, 10000, 100000, 1000000),
        fragment_sizes = (100000, 1000000),        
        random_access_points = 1024,
        verify_fragment_size = 100,
        verify_num_iterations = 0,
        )
    
    (options, args) = Experiment.Start( parser )

    tempdir = tempfile.mkdtemp()

    ##############################################################################
    ##############################################################################
    ##############################################################################
    ## Compression of data
    ##############################################################################    
    options.stdout.write("method\ttime\tnerrors1\tnerrors2\n" )

    if options.stdlog >= 1:
        options.stdlog.write("# building uncompressed database.\n" )
        options.stdlog.flush()
                             
    dbfile_uncompressed = tempdir + "/uncompressed"
    timer = timeit.Timer( stmt="createDatabase( '%s', %s, random_access_points=%i )" % (dbfile_uncompressed, str(args),
                                                                                        options.random_access_points ),
                          setup="from __main__ import createDatabase" )
                                                                                                        
    t = timer.timeit( number = 1 )    
    options.stdout.write("uncompressed\t%i\t0\t0\n" % (t ) )
    dbfiles = []
    
    for compression in options.methods:

        if options.stdlog >= 1:
            options.stdlog.write("# building compressed database with method %s.\n" % (compression) )
            options.stdlog.flush()
        
        dbfile = tempdir + "/" + compression
        
        timer = timeit.Timer( stmt="createDatabase( '%s', %s, random_access_points=%i, compression='%s')" % (dbfile, str(args),
                                                                                                              options.random_access_points,
                                                                                                              compression ),
                              setup="from __main__ import createDatabase" )

        t = timer.timeit( number = 1 )

        if options.stdlog >= 1:
            options.stdlog.write("# verifying compressed database with method %s.\n" % (compression) )
            options.stdlog.flush()

        fasta = IndexedFasta( dbfile )
        nerrors1 = verify( IndexedFasta(dbfile), IndexedFasta(dbfile_uncompressed),
                           options.verify_num_iterations, options.verify_fragment_size,
                           quiet = True )
        nerrors2 = verify( IndexedFasta(dbfile_uncompressed), IndexedFasta(dbfile),
                           options.verify_num_iterations, options.verify_fragment_size,                           
                           quiet = True )        

        options.stdout.write("%s\t%i\t%i\t%i\n" % (compression, t, nerrors1, nerrors2 ))
        options.stdout.flush()
        
        dbfiles.append( dbfile )

    ##############################################################################
    ##############################################################################
    ##############################################################################
    ## random sampling of data points
    ##############################################################################    
    options.stdout.write("//\n")
    
    options.stdout.write( "method\tsize\t%s\tvalues\n" % ("\t".join(Stats.DistributionalParameters().getHeaders())))
    options.stdout.flush()
        
    for fragment_size in options.fragment_sizes:

        times = [ [] for x in range(len(options.methods)+1)] 

        for iteration in range(options.num_iterations):

            for x in range(len(options.methods)):

                if options.stdlog >= 1:
                    options.stdlog.write("# fragment_size=%i, iteration=%i/%i, method=%s.\n" % (fragment_size, iteration, options.num_iterations,options.methods[x]) )
                    options.stdlog.flush()

                timer = timeit.Timer( stmt="benchmarkRandomFragment( fasta = fasta, size = %i)" % (fragment_size),
                                      setup="""from __main__ import benchmarkRandomFragment,IndexedFasta\nfasta=IndexedFasta( "%s" )""" % (dbfiles[x]) )

                t = timer.timeit( number = options.benchmark_num_iterations )

                times[x].append( t )

            if options.stdlog >= 1:
                options.stdlog.write("# fragment_size=%i, iteration=%i/%i, method=%s.\n" % (fragment_size, iteration, options.num_iterations, "uncompressed") )
                options.stdlog.flush()

            timer = timeit.Timer( stmt="benchmarkRandomFragment( fasta = fasta, size = %i)" % (fragment_size),
                                  setup="""from __main__ import benchmarkRandomFragment,IndexedFasta\nfasta=IndexedFasta( "%s" )""" % (dbfile_uncompressed) )

            t = timer.timeit( number = options.benchmark_num_iterations )
            
            times[-1].append( t )


        for x in range(len(options.methods)):
            values = times[x]
            options.stdout.write( "%s\t%i\t%s\t%s\n" % (options.methods[x], fragment_size, str(Stats.DistributionalParameters(values)), ",".join(map(str,values))))

        values = times[-1]
        options.stdout.write( "%s\t%i\t%s\t%s\n" % ("uncompressed", fragment_size, str(Stats.DistributionalParameters(values)), ",".join(map(str,values))))
        options.stdout.flush()


    Experiment.Stop()

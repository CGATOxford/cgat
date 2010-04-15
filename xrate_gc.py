################################################################################
#   Gene prediction pipeline 
#
#   $Id: xrate_gc.py 2781 2009-09-10 11:33:14Z andreas $
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

import Experiment

USAGE="""python %s [OPTIONS] 

Version: $Id: xrate_gc.py 2781 2009-09-10 11:33:14Z andreas $


Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
""" % sys.argv[0]

import MatlabTools
import XGram
import XGram.Exceptions

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: xrate_gc.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-p", "--filename-input-pattern", dest="filename_input_pattern", type="string",
                      help="shell like glob pattern for collecting all files to be parsed." )

    parser.add_option("-e", "--pattern-id", dest="pattern_id", type="string",
                      help="pattern for extracting identifier from filename." )

    parser.set_defaults(
        filename_input_pattern = None,
        pattern_id = None,
        report_step = 100,
        counts_format = "%6.4f",
        set_missing_to_zero = True,
        )
    
    (options, args) = Experiment.Start( parser )

    ninput, noutput, nerrors = 0, 0, 0


    keys = None
    table = []
    
    if not options.filename_input_pattern:
        lines = sys.stdin.readlines()
        model = XGram.Parser.parseGrammar( lines )
        counts = model.mGrammar.mObservedCounts
        keys = set(counts.keys())
        table.append( ("stdin", counts) )
        
    else:
        files = glob.glob( options.filename_input_pattern )

        if options.loglevel >= 1:
            options.stdlog.write("# collecting data from %i files\n" % len(files) )
        
        for file in files:
            ninput += 1

            if options.loglevel >= 1 and (ninput % options.report_step == 0):
                options.stdlog.write( "# iteration: %i/%i %5.2f%%\n" % (ninput, len(files), 100.0 * ninput/len(files)) ) 
            
            infile = open(file, "r")
            lines = infile.readlines()
            infile.close()
            try:
                model = XGram.Parser.parseGrammar( lines )
            except XGram.Exceptions.ParsingError:
                nerrors += 1
                if options.loglevel >= 1:
                    options.stdlog.write( "# parsing error in file %s\n" % file )
                    options.stdlog.flush()
                continue

            counts = model.mGrammar.mObservedCounts
            skeys = set(counts.keys())

            if keys == None:
                keys = skeys
            else:
                if len(skeys.symmetric_difference( keys )) > 0:
                    if options.set_missing_to_zero:
                        for k in keys.difference(skeys):
                            counts[k] = { 'counts' : 0 }
                    else:
                        nerrors += 1
                        if options.loglevel >= 1:
                            options.stdlog.write( "# missing columns in file %s: %s\n" % (file, str(skeys.symmetric_difference(keys)) ) )
                            options.stdlog.flush()
                        continue

            table.append( (file, counts) )
            
    keys = list(keys)
    keys.sort()

    if options.pattern_id:
        rx = re.compile(options.pattern_id)
        extract_id = lambda x: rx.search(x).groups()[0]
    else:
        extract_id = lambda x: x

    ## get columns for percent_GC computation
    nodes = map( lambda x: x[len("RATE_"):], filter( lambda x: "RATE_" in x, keys))        

    ## print header
    options.stdout.write("id\t" + "\t".join( keys ) )
    for node in nodes:
        options.stdout.write("\tfGC_%s" % node )
        options.stdout.write("\tvGC_%s" % node )            
    options.stdout.write("\n")

    ## get nodes
    for title, row in table:
        noutput += 1
        id = extract_id( title )
        options.stdout.write( "%s\t%s" % (id, "\t".join([ options.counts_format % (row[x]['counts']) for x in keys ] )))

        for node in nodes:
            a = row["pGC_%s" % node]['counts']
            b = row["pAT_%s" % node]['counts']
            t = a + b
            if t > 0:
                r = 100.0 * a / t
                ## variance of beta-distributed random variable X
                v = a * b / ( t * t * (t + 1) )
            else:
                r = 0
                v = 0
            options.stdout.write( "\t%6.4f\t%6.4f" % (r, v) )

        options.stdout.write("\n")

    if options.loglevel >= 1:
        options.stdout.write("# ninput=%i, noutput=%i, nerrors=%i\n" % (ninput, noutput, nerrors ) )
        
    Experiment.Stop()





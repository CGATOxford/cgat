import sys, re, string, os, optparse

USAGE = """python %s < stdin > stdout

read in data and append columns to a density histogram

-> relative frequencies
-> cumulative counts and frequencies in both directions

'#' at start of line is a comment

""" % sys.argv[0]

import Experiment
import IOTools
import numpy

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: histogram2histogram.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-i", "--is-int", dest="is_ints", action="store_true",
                      help="categories are integers."  )
    parser.add_option("-m", "--method", dest="method", type="string",
                       help="method(s) to apply." )
    parser.add_option("--no-headers", dest="headers", action="store_false",
                      help="histogram has no headers.")
    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to use for plotting.")
    parser.add_option("", "--truncate", dest="truncate", type="string",
                      help="truncate at range."  )
    parser.add_option("", "--no-out-of-range", dest="cumulate_out_of_range", action = "store_false",
                      help="add up bins out of range."  )
    parser.add_option( "--format-bin", dest="format_bin", type="string",
                      help="format for bins."  )
    parser.add_option( "--format-val", dest="format_val", type="string",
                      help="format for vals."  )

    parser.set_defaults(
        is_ints = False,
        method = "cumul",
        columns = "all",
        headers = True,
        truncate = None,
        cumulate_out_of_range = True,
        format_bin="%6.4f",
        format_val="%6.4f",
        )

    (options, args) = Experiment.Start( parser )

    if options.truncate: options.truncate = map(float, options.truncate.split(","))
    options.method = options.method.split(",")
    data, legend = IOTools.readTable( sys.stdin,
                                      numeric_type=numpy.float32,
                                      take=options.columns,
                                      headers = options.headers,
                                      truncate= options.truncate,
                                      cumulate_out_of_range = options.cumulate_out_of_range )

    nfields = len(legend)

    ## note: because of MA, iteration makes copy of slices
    ## Solution: inplace edits.
    nrows, ncols = data.shape

    for method in options.method:
        if method == "cumul":
            l = [0] * ncols
            for x in range(nrows):
                for y in range(1, ncols):
                    data[x,y] += l[y]
                    l[y] = data[x,y]

        elif method == "rcumul":
            l = [0] * ncols
            for x in range(nrows-1,0,-1):
                for y in range(1, ncols):
                    data[x,y] += l[y]
                    l[y] = data[x,y]

        elif method == "normalize":
            m = [0] * ncols
            for x in range(nrows):
                for y in range(1, ncols):
                    ## the conversion to float is necessary
                    m[y] = max( m[y], float(data[x,y]) )

            for y in range(1, ncols):
                if m[y] == 0: m[y] = 1.0

            for x in range(nrows):
                for y in range(1, ncols):
                    data[x,y] = data[x,y] / m[y]
        else:
            raise "unknown method %s" % method

    print "\t".join(legend)

    format = options.format_bin + "\t" + "\t".join( [options.format_val] * (nfields-1)) 
    
    for d in data:
        print format % tuple(d)
    
    Experiment.Stop()










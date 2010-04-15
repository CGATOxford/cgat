import sys, re, string, os, optparse

USAGE = """python %s < stdin > stdout

read in data and append columns to a density histogram

-> relative frequencies
-> cumulative counts and frequencies in both directions

'#' at start of line is a comment

""" % sys.argv[0]

import Experiment
import Histogram

if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id: append_histogram.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-i", "--is-int", dest="is_ints", action="store_true",
                      help="categories are integers."  )

    parser.set_defaults(
        is_ints = False
        )

    (options, args) = Experiment.Start( parser )

    vals = []
    
    # retrieve histogram
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    # check if first line contains a header
    d = string.split(lines[0][:-1], "\t")[0]
    try:
        if options.is_ints:
            value = int(d)
        else:
            value = float( d )
    except ValueError:
        print string.join( (d, "counts", "frequency",
                            "cumulative counts", "increasing cumulative frequency",
                            "cumulative counts", "decreasing cumulative frequency"), "\t")
        del lines[0]
        
    data = map( lambda x: map(float, string.split(x[:-1], "\t")), lines)

    if len(data) == 0:
        raise "No data found."
        
    total = float(reduce( lambda x,y: x+y, map( lambda x: x[1], data)))

    cumul_down = int(total)
    cumul_up = 0

    if options.is_ints:
        form = "%i\t%i\t%6.4f\t%i\t%6.4f\t%i\t%6.4f"         
    else:
        form = "%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%6.4f" 
        
    for bin,val in data:
        percent     = float(val) / total
        cumul_up   += val
        percent_cumul_up = float(cumul_up) / total
        percent_cumul_down = float(cumul_down) / total        
        
        print form % \
              (bin, val, percent, cumul_up, percent_cumul_up, cumul_down, percent_cumul_down)

        cumul_down -= val

    Experiment.Stop()










####
####
##
## Project PairsDBTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: calculate_histogram.py 2782 2009-09-10 11:40:29Z andreas $
##
##
####
####


import sys, re, string, os, getopt, time

import Experiment
import Histogram

USAGE = """python calculate_histogram.py < stdin > stdout

read in data and build histogram of column

-c, --column            column to take [default = 0]
-a, --append=           append columns [normalize]
-n, --normalize         normalize column
--cumulative            cumulative histogram
--reverse-cumulative    reverse cumulative histogram
-i, --titles            use supplied titles
# at start of line is a comment
"""

param_loglevel = 1
param_separator = "//"
param_take = None
param_fill = 0
param_nonull = None
param_columns = [0,]
param_empty_bins = 1
param_titles = False

param_lower_limit = None
param_upper_limit = None
param_bin_size = None
param_scale = None
param_normalize = None
param_append = []

param_cumulative = False
param_reverse_cumulative = False

param_long_options = ["Verbose=", "nonull", "fill","take=", "column=", "show_empty",
                      "upper=", "lower=", "bin-size=", "scale=", "normalize", "append=", "titles",
                      "cumulative", "reverse-cumulative"]

param_short_options = "v:nft:c:eu:l:b:a:io"

##---------------------------------------------------------------------------------------------------------        
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      param_short_options,
                                      param_long_options)

    except getopt.error, msg:
        print USAGE
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)
        elif o in ("-t", "--take"):
            param_take = map(string.atoi, string.split(a, ","))
        elif o in ("-f", "--fill"):
            param_fill = 1
        elif o in ("-n", "--nonull"):
            param_nonull = ""
        elif o in ("-c", "--column"):
            if a == "all":
                param_columns = "all"
            else:
                param_columns = map(lambda x: int(x) - 1, string.split( a, ","))
        elif o in ("-e", "--show_empty"):
            param_empty_bins = 0
        elif o in ("-u", "--upper"):
            param_upper_limit = string.atof(a)
        elif o in ("-l", "--lower"):
            param_lower_limit = string.atof(a)
        elif o in ("-b", "--bin-size"):
            param_bin_size = string.atof(a)
        elif o in ("-s", "--scale"):
            param_scale = float(a)
        elif o in ("-o", "--normalize"):
            param_normalize = True
        elif o in ("-a", "--append"):
            param_append = string.split( a, ",")
        elif o in ("-i", "--titles"):
            param_titles = True
        elif o == "--cumulative":
            param_cumulative = True
        elif o == "--reverse-cumulative":
            param_reverse_cumulative = True
        
    histograms = []

    if param_loglevel >= 1:
        print Experiment.GetHeader()
        print Experiment.GetParams()    
    
    vals = []
    
    for x in param_columns: vals.append( [] )
    
    # retrieve histogram
    lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())

    ncols = len(string.split(lines[0][:-1], "\t"))
    if param_columns == "all":
        param_columns = range(ncols)
        for x in param_columns: vals.append( [] )

    if param_titles:
        data = lines[0][:-1].split("\t")
        del lines[0]
        param_titles = map( lambda x: data[x], param_columns)
        
    for l in lines:
        data = string.split(l[:-1], "\t")
            
        for x in range(len(param_columns)):
            try:
                v = string.atof(data[param_columns[x]])
            except IndexError:
                print "# IndexError in line:", l[:-1]
                continue
            except ValueError:
                continue

            if param_scale:
                v *= param_scale

            if param_upper_limit != None and v > param_upper_limit:
                v = param_upper_limit

            if param_lower_limit != None and v < param_lower_limit:
                v = param_lower_limit

            vals[x].append( v )

    lines = None

    hists = []
    titles = []
    
    for x in range(len(param_columns)):
        if param_loglevel >= 1:
            print "# column=%i, num_values=%i" % (param_columns[x], len(vals[x]))

        if len(vals[x]) == 0: continue
        
        h = Histogram.Calculate( vals[x], no_empty_bins = param_empty_bins, increment = param_bin_size)
        if param_scale: h = Histogram.Scale( h, 1.0 / param_scale )

        if param_normalize: h = Histogram.Normalize( h )
        if param_cumulative: h = Histogram.Cumulate( h )
        if param_reverse_cumulative: h = Histogram.Cumulate( h, direction = 0 )
        
        hists.append(h)

        for m in param_append:
            if m == "normalize":
                hists.append( Histogram.Normalize( h ) )

        if param_titles:
            titles.append( param_titles[x] )

    if titles:
        print "bin\t" + "\t".join(titles)

    if len(hists) == 1:
        Histogram.Print( hists[0], nonull = param_nonull )
    else:
        combined_histogram = Histogram.Combine( hists )
        Histogram.Print( combined_histogram, nonull = param_nonull )        

    if param_loglevel >= 1:        
        print Experiment.GetFooter()
    











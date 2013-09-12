################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
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
'''
r_table2scatter.py - R based plots and stats
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a table from a file or stdin. 
It can compute various stats (correlations, ...)
and/or plot the data using R.

Usage
-----

Example::

   python r_table2scatter.py --help

Type::

   python r_table2scatter.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import tempfile
import subprocess
import optparse
import time
import math
import code
import CGAT.Experiment as E
import CGAT.MatlabTools as MatlabTools
import numpy
import CGAT.Stats as Stats

from rpy2.robjects import r as R
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri

def readTable( lines,
               assign,
               separator = "\t",
               numeric_type = numpy.float,
               take_columns = "all",
               headers = True,
               row_names = None,
               colours = None,
               ):
    """read a matrix. There probably is a routine for this in Numpy, which
    I haven't found yet. 

    The advantage of using rpy's mechanism is to allow for 
    missing values.

    row_names: column with row names
    """

    # now import R as stdin has been read.
    take = take_columns
    
    handle, name = tempfile.mkstemp()
    os.close(handle)

    if take_columns == "all":
        num_cols = len(string.split(lines[0][:-1], "\t"))
        take = range( 0, num_cols)
    elif take_columns == "all-but-first":
        num_cols = len(string.split(lines[0][:-1], "\t"))        
        take = range( 1, num_cols)
        
    outfile = open(name, "w")
    c = []

    ## delete columns with colour information/legend
    to_delete = []
    if row_names != None:
        legend = []
        if take_columns == "all":
            to_delete.append( row_names )
    else:
        legend = None

    if colours != None:
        if take_columns == "all":
            to_delete.append( colours )

    to_delete.sort()
    to_delete.reverse()
    for x in to_delete: del take[x]

    ## get column headers
    if headers:
        headers = lines[0][:-1].split("\t")
        headers = map( lambda x: headers[x], take )        
        del lines[0]
    
    for l in lines:
        data = [x.strip() for x in l[:-1].split("\t") ]
        if not data or not [ x for x in data if x != ""] : continue
        outfile.write(string.join( map( lambda x: data[x], take ), "\t") + "\n")
        if row_names != None:
            legend.append( data[row_names] )
            
        if colours != None:
            c.append( data[colours] )
        
    outfile.close()

    # rpy.set_default_mode(rpy.NO_CONVERSION)
    # note that the conversion is not perfect. Some missing values are assigned to "nan", while some
    # are -2147483648. They seem to treated correctly, though, within R, but note that when computing
    # something like sum(), the result in python after conversion might be -2147483648.
    if headers:
        matrix = R("""%s <- read.table( '%s', na.string = c("NA", "na", 'nan', 'NaN'), col.names=c('%s'), sep="\t" )""" % (assign, name, "','".join(headers)) )
    else:
        matrix = R("""%s <- read.table( '%s', na.string = c("NA", "na", 'nan', 'NaN'), col.names=headers, sep="\t" )""" % (assign, name) )

    # rpy.set_default_mode(rpy.BASIC_CONVERSION)
    os.remove( name )

    return matrix, headers, c, legend

def writeMatrix( file,
                 matrix,
                 separator = "\t",
                 headers = [],
                 format = "%f" ):
    """write a matrix to file.

    if headers are given, add them to columns and rows.
    """
    
    if headers:
        file.write("\t" + string.join(headers, "\t") + "\n")

    nrows, ncols = matrix.shape

    for x in range(nrows):
        if headers: file.write(headers[x] + "\t")
        file.write(string.join (map(lambda x: format % x, matrix[x]), "\t") + "\n")

def FuncScatterDiagonal( data ):

    R.points( data )
    R.abline( 0, 1 )

def main( argv = None ):
    
    if argv == None: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: r_table2scatter.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="columns to take from table. Choices are 'all', 'all-but-first' or a ','-separated list of columns." )
    
    parser.add_option( "--logscale", dest="logscale", type="string",
                      help="log-transform one or both axes [default=%Default]."  )

    parser.add_option("-a", "--hardcopy", dest="hardcopy", type="string",
                      help="write hardcopy to file [default=%default].", 
                      metavar = "FILE" )

    parser.add_option("-f", "--file", dest="input_filename", type="string",
                      help="filename with table data [default=%default].",
                      metavar = "FILE")

    parser.add_option("-2", "--file2", dest="input_filename2", type="string",
                      help="additional data file [default=%default].",
                      metavar = "FILE")

    parser.add_option("-s", "--stats", dest="statistics", type="choice",
                      choices=("correlation", "spearman", "pearson", "count"),
                      help="statistical quantities to compute [default=%default]",
                      action = "append")
    
    parser.add_option("-p", "--plot", dest="plot", type="choice",
                      choices=("scatter", "pairs", "panel", "bar", "bar-stacked", 
                               "bar-besides", "1_vs_x", "matched", "boxplot", "scatter+marginal",
                               "scatter-regression" ),
                      help="plots to plot [default=%default]",
                      action = "append")

    parser.add_option("-t", "--threshold", dest="threshold", type="float",
                      help="min threshold to use for counting method [default=%default].")

    parser.add_option("-o", "--colours", dest="colours", type="int",
                      help="column with colour information [default=%default].")

    parser.add_option("-l", "--labels", dest="labels", type="string",
                      help="column labels for x and y in matched plots [default=%default].")

    parser.add_option("-d", "--add-diagonal", dest="add_diagonal", action="store_true",
                      help="add diagonal to plot [default=%default].")
    
    parser.add_option("-e", "--legend", dest="legend", type="int",
                      help="column with legend [default=%default].")

    parser.add_option("-r", "--options", dest="r_options", type="string",
                      help="R plotting options [default=%default].")
    
    parser.add_option("--format", dest="format", type="choice",
                      choices=("full", "sparse"),
                      help="output format [default=%default]." )

    parser.add_option( "--title", dest="title", type="string",
                       help="""plot title [default=%default].""")

    parser.add_option("", "--xrange", dest="xrange", type="string",
                      help="x viewing range of plot [default=%default]."  )

    parser.add_option("", "--yrange", dest="yrange", type="string",
                      help="y viewing range of plot[default=%default]."  )

    parser.add_option( "--allow-empty", dest="fail_on_empty", action="store_false",
                      help="do not fail on empty input [default=%default].")

    parser.add_option( "--fail-on-empty", dest="fail_on_empty", action="store_true",
                      help="fail on empty input [default=%default].")

    parser.set_defaults( \
        hardcopy = None,
        input_filename = "",
        input_filename2 = None,
        columns = "all",
        logscale = None,
        statistics = [],
        plot=[],
        threshold=0.0,
        labels = "x,y",
        colours= None,
        diagonal = False,
        legend = None,
        title = None,
        xrange = None,
        yrange = None,        
        r_options = "",
        fail_on_empty = True,
        format = "full")

    (options, args) = E.Start( parser )

    if len(args) == 1 and not options.input_filename:
        options.input_filename = args[0]

    if options.columns not in ("all", "all-but-first"):
        options.columns = map( lambda x: int(x)-1, options.columns.split(","))
        
    if options.colours: options.colours -= 1
    if options.legend: options.legend -= 1
    
    table ={}
    headers = []

    ## read data matrix
    if options.input_filename:
        lines = open(options.input_filename, "r").readlines()
    else:
        ## note: this will not work for interactive viewing, but
        ## creating hardcopy plots works.
        lines = sys.stdin.readlines()

    lines = filter( lambda x: x[0] != "#", lines)
    
    if len(lines) == 0:
        if options.fail_on_empty:
            raise IOError ( "no input" )
        E.warn( "empty input" )
        E.Stop()
        return

    matrix, headers, colours, legend = readTable( lines,
                                                  "matrix",
                                                  take_columns = options.columns,
                                                  headers=True,
                                                  colours=options.colours,
                                                  row_names = options.legend )

    if options.input_filename2:
        ## read another matrix (should be of the same format.
        matrix2, headers2, colours2, legend2 = readTable( lines,
                                                 "matrix2",
                                                 take_columns = options.columns,
                                                 headers=True,
                                                 colours=options.colours,
                                                 row_names = options.legend )
    R.assign("headers", headers)

    ndata = R( """length( matrix[,1] )""" )[0]

    if options.loglevel >=1:
        options.stdlog.write("# read matrix: %ix%i\n" % (len(headers),ndata) )

    if colours:
        R.assign("colours", colours)

    for method in options.statistics:

        if method == "correlation":
            cor = R.cor(matrix, use="pairwise.complete.obs" )
            writeMatrix( sys.stdout, cor, headers=headers, format = "%5.2f" )

        elif method == "pearson":
            options.stdout.write( "\t".join( ("var1", 
                                              "var2",
                                              "coeff",
                                              "passed",
                                              "pvalue",
                                              "n",
                                              "method",
                                              "alternative" )) + "\n" )
            for x in range(len(headers)-1):
                for y in range( x+1, len(headers)):
                    try:
                        result = R("""cor.test( matrix[,%i], matrix[,%i] )""" % (x + 1, y + 1))
                    except rpy.RPyException, msg:
                        E.warn( "correlation not computed for columns %i(%s) and %i(%s): %s" % (x, headers[x], y, headers[y], msg) )
                        options.stdout.write( "%s\t%s\t%s\t%s\t%s\t%i\t%s\t%s\n" % \
                                                  (headers[x], headers[y],
                                                   "na",
                                                   "na",
                                                   "na",
                                                   0,
                                                   "na",
                                                   "na" ))

                    else:
                        options.stdout.write( "%s\t%s\t%6.4f\t%s\t%e\t%i\t%s\t%s\n" % \
                                                  (headers[x], headers[y],
                                                   result.rx2('estimate').rx2('cor')[0], 
                                                   Stats.getSignificance( float(result.rx2('p.value')[0]) ),
                                                   result.rx2('p.value')[0],
                                                   result.rx2('parameter').rx2('df')[0],
                                                   result.rx2('method')[0], 
                                                   result.rx2('alternative')[0]) )

        elif method == "spearman":
            options.stdout.write( "\t".join( ("var1", "var2",
                                              "coeff",
                                              "passed",
                                              "pvalue",
                                              "method",
                                              "alternative" )) + "\n" )
            for x in range(len(headers)-1):
                for y in range( x+1, len(headers)):
                    result = R("""cor.test( matrix[,%i], matrix[,%i], method='spearman' )""" % (x + 1, y + 1))
                    options.stdout.write( "%s\t%s\t%6.4f\t%s\t%e\t%i\t%s\t%s\n" % \
                                              (headers[x], headers[y],
                                               result['estimate']['rho'], 
                                               Stats.getSignificance( float(result['p.value']) ),
                                               result['p.value'],
                                               result['parameter']['df'],
                                               result['method'], 
                                               result['alternative']))
                    
        elif method == "count":
            ## number of shared elements > threshold
            m, r, c = MatlabTools.ReadMatrix( open(options.input_filename, "r"),
                                              take = options.columns, 
                                              headers=True)
            mask = numpy.greater(m, options.threshold)
            counts = numpy.dot( numpy.transpose(mask), mask)
            writeMatrix( options.stdout, counts, headers=c, format = "%i" )

    if options.plot:

        # remove columns that are completely empty
        if "pairs" in options.plot:
            colsums = R('''colSums( is.na(matrix ))''')
            take = [ x for x in range(len(colsums)) if colsums[x] != ndata ]
            if take:
                E.warn( "removing empty columns %s before plotting" % str(take))
                matrix = R.subset( matrix, select=[ x+1 for x in take] )
                R.assign("""matrix""", matrix )
                headers = [ headers[x] for x in take ]
                if legend:
                    legend = [ headers[x] for x in take ]
        
        if options.r_options:
            extra_options = ", %s" % options.r_options
        else:
            extra_options = ""

        if options.legend != None and len(legend):
            extra_options += ", legend=c('%s')" % "','".join(legend)
            
        if options.labels:
            xlabel, ylabel = options.labels.split(",")
            extra_options += ", xlab='%s', ylab='%s'" % (xlabel, ylabel)
        else:
            xlabel, ylabel = "", ""

        if options.colours:
            extra_options += ", col=colours"
            
        if options.logscale:
            extra_options += ", log='%s'" % options.logscale

        if options.xrange:
            extra_options += ", xlim=c(%f,%f)" % tuple(map(float, options.xrange.split(",") ) )

        if options.yrange:
            extra_options += ", ylim=c(%f,%f)" % tuple(map(float, options.yrange.split(",") ) )

        if options.hardcopy:
            if options.hardcopy.endswith(".eps"): 
                R.postscript(options.hardcopy)
            elif options.hardcopy.endswith(".png"): 
                R.png(options.hardcopy, width=1024, height=768, type="cairo")
            elif options.hardcopy.endswith(".jpg"): 
                R.jpg(options.hardcopy, width=1024, height=768, type="cairo")

        for method in options.plot:

            if ndata < 100:
                point_size = "1"
                pch = "o"
            elif ndata < 1000:
                point_size = "1"
                pch = "o"
            else:
                point_size = "0.5"
                pch = "."

            if method == "scatter":
                R("""plot( matrix[,1], matrix[,2], cex=%s, pch="o" %s)""" % (point_size, extra_options) )

            if method == "scatter-regression":
                R("""plot( matrix[,1], matrix[,2], cex=%s, pch="o" %s)""" % (point_size, extra_options) )
                dat = R("""dat <- data.frame(x = matrix[,1], y = matrix[,2])""")
                R("""new <- data.frame(x = seq( min(matrix[,1]), max(matrix[,1]), (max(matrix[,1]) - min(matrix[,1])) / 100))""")
                mod = R("""mod <- lm( y ~ x, dat)""")
                R("""predict(mod, new, se.fit = TRUE)""")
                R("""pred.w.plim <- predict(mod, new, interval="prediction")""")
                R("""pred.w.clim <- predict(mod, new, interval="confidence")""")
                R("""matpoints(new$x,cbind(pred.w.clim, pred.w.plim[,-1]), lty=c(1,2,2,3,3), type="l")""")
                R.mtext(
                    "y = %f * x + %f, r=%6.4f, n=%i" % (mod["coefficients"]["x"], 
                                                        mod["coefficients"]["(Intercept)"], 
                                                        R("""cor( dat )[2]"""), 
                                                        ndata ),
                    3,
                    cex = 1.0)

            elif method == "pairs":
                if options.add_diagonal:
                    R( """panel.hist <- function( x,y,...  ) { points(x,y,...); abline(0,1); }""" )
                else:
                    R( """panel.hist <- function( x,y,...  ) { points(x,y,...); }""" )
                    
                # There used to be a argument na_action="na.omit", but removed this
                # as there appeared error messages saying "na.action is not a graphical parameter"
                # and the plots showed occasionally the wrong scale.
                # cex=point_size also caused trouble (error message: "X11 used font size 8 when 2 was requested" or similar)
                if options.colours:
                    R.pairs( matrix, 
                             pch=pch, 
                             col=colours,
                             main=options.title,
                             panel="panel.hist",
                             labels=headers, 
                             cex_labels=2.0 )
                else:
                    R.pairs( matrix, 
                             pch=pch,
                             panel="panel.hist", 
                             main=options.title,
                             labels=headers, 
                             cex_labels=2.0 )

            elif method == "boxplot":
                extra_options += ",main='%s'" % options.title

                # set vertical orientation
                if max( [len(x) for x in headers] ) > 40 / len(headers):
                    # remove xlabel:
                    extra_options = re.sub( ", xlab='[^']+'", "", extra_options )
                    extra_options += ", names.arg=headers, las=2"
                    R("""op <- par(mar=c(11,4,4,2))""") # the 10 allows the names.arg below the barplot

                R("""boxplot( matrix %s)""" % extra_options)
                
            elif method == "bar" or method == "bar-stacked":
                if not options.colours:
                    extra_options += ", col=rainbow(5)"

                # set vertical orientation
                if max( [len(x) for x in headers] ) > 40 / len(headers):
                    # remove xlabel:
                    extra_options = re.sub( ", xlab='[^']+'", "", extra_options )
                    extra_options += ", names.arg=headers, las=2"
                    R("""op <- par(mar=c(11,4,4,2))""") # the 10 allows the names.arg below the barplot
                
                R("""barplot(as.matrix(matrix), %s)""" % extra_options)


            elif method == "bar-besides":
                if not options.colours:
                    extra_options += ", col=rainbow(%i)" % ndata

                # set vertical orientation
                if max( [len(x) for x in headers] ) > 40 / len(headers):
                    # remove xlabel:
                    extra_options = re.sub( ", xlab='[^']+'", "", extra_options )
                    extra_options += ", names.arg=headers, las=2"
                    R("""op <- par(mar=c(11,4,4,2))""") # the 10 allows the names.arg below the barplot

                R("""barplot(as.matrix(matrix), beside=TRUE %s)""" % extra_options)

            elif method == "scatter+marginal":

                if options.title:
                    # set the size of the outer margins - the title needs to be added at the end
                    # after plots have been created
                    R.par(oma=R.c(0,0,4,0) )                     

                R( """matrix""" )
                R( """
x <- matrix[,1];
y <- matrix[,2];
xhist <- hist(x, breaks=20, plot=FALSE);
yhist <- hist(y, breaks=20, plot=FALSE);
top <- max(c(xhist$counts, yhist$counts));
nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), respect=TRUE );
par(mar=c(3,3,1,1)) ;
plot(x, y, cex=%s, pch="o" %s) ;
par(mar=c(0,3,1,1)) ;
barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0 ) ;
par(mar=c(3,0,1,1)) ;
title(main='%s');
barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE ) ;
title(main='%s');
""" % (point_size, extra_options, xlabel, ylabel))

                if options.title:
                    R.mtext( options.title,3,outer=True,line=1,cex=1.5)  

            elif method in ("panel", "1_vs_x", "matched") :

                if method == "panel":
                    pairs = []
                    for x in range(len(headers)-1):
                        for y in range( x+1, len(headers)):
                            pairs.append( (x, y) )
                        
                elif method == "1_vs_x":
                    pairs = []
                    for x in range(1, len(headers)):
                        pairs.append( (0, x) )
                    
                # print matching columns
                elif method == "matched":
                    pairs = []
                    for x in range(len(headers)-1):
                        for y in range( x+1, len(headers)):
                            if headers[x] == headers[y]:
                                pairs.append( (x, y) )
                                break

                w = int(math.ceil(math.sqrt(len(pairs))))
                h = int(math.ceil(float(len(pairs)) / w))

                PosInf = 1e300000
                NegInf = -1e300000             

                xlabel, ylabel = options.labels.split(",")

                R( """layout(matrix(seq(1,%i), %i, %i, byrow = TRUE))""" % (w*h, w, h))
                for a,b in pairs:
                    new_matrix = filter( lambda x: \
                                         x[0] not in (float("nan"), PosInf, NegInf) and \
                                         x[1] not in (float("nan"), PosInf, NegInf), \
                                         zip(matrix[a].values()[0], matrix[b].values()[0]) )
                    try:
                        R( """plot(matrix[,%i], matrix[,%i], main='%s versus %s', cex=0.5, pch=".", xlab='%s', ylab='%s' )""" % (\
                            a+1,b+1,headers[b], headers[a], xlabel, ylabel))
                    except rpy.RException, msg:
                        print "could not plot %s versus %s: %s" % (headers[b], headers[a], msg)

        if options.hardcopy:
            R['dev.off']()

    E.info( "matrix added as >matrix< in R.")
        
    if not options.hardcopy:
        if options.input_filename:
            interpreter = code.InteractiveConsole(globals())
            interpreter.interact()
        else:
            E.info( "can not start new interactive session as input has come from stdin.")

    E.Stop()

if __name__ == "__main__":
    main() 

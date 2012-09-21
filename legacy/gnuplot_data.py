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
gnuplot_data.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python gnuplot_data.py --help

Type::

   python gnuplot_data.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys, re, string, os, getopt, time, tempfile

USAGE = """python %s < stdin > stdout

plot a histogram (or a set of histograms).

-v, --verbose   verbosity
-l, --legend    set legend
-t, --title     set title
-c, --hardcopy  create hardcopy of picture
-w, --with      with command to gnuplot (for example, "points")
-f, --fit       do linear fit to data
-i, --files     plot a set of files
-o, --logscale=#    set logscale
-u, --function= plot a function
'#' at start of line is a comment

""" % sys.argv[0]

param_long_options = ["help", "verbose=", "legend=", "title=", "hardcopy=", "blackwhite", "with=", "fit", "files=",
                      "logscale=", "function="]
param_short_options = "hv:l:t:c:bw:fi:u:"
                      
param_legend = None
param_title   = None
param_hardcopy  = None
param_blackwhite = None
param_terminal = "postscript"
param_with="points"
param_fit = None

param_logscale = None

param_filenames = None

param_function = None

import Experiment
import Histogram
import Gnuplot

def PlotFit( g, data, cols=(0,1) ):

    fh1, fn1 = tempfile.mkstemp()
    a,b = cols
    os.close(fh1)
    outfile = open(fn1, "w")
    for d in data: outfile.write("%f\t%f\n" % (d[a], d[b]))
    outfile.close()
    
    parameters = {}
    fh2, fn2 = tempfile.mkstemp()
    fh3, fn3 = tempfile.mkstemp()    
    os.close(fh2)
    os.close(fh3)
    open(fn2, 'w').write('m=0\nb=0\n')
    g("f%i(x) = m * x + y0" % b) 
    g("fit f%i(x) '%s' using 1:2 via y0, m" % (b, fn1))
    g("replot f%i(x)" % (b))
    
##     g('fit m*x+b "%s" via "%s"' % (fn1, fn2) )    
##     g('update "%s" "%s"' % (fn2, fn3))
##     execfile( fn3, globals(), parameters )
##     g.replot( Gnuplot.Func( "%f*x + %f" % (parameters['m'], parameters['b']) ) )
        
    return [fn1, fn2, fn3]


def Fit( data, cols=(0,1) ):

    import scipy.linalg

    a,b = cols
    
    matrix = []
    imatrix = []
    
    for d in data:
        matrix.append([1.0, d[a]]) # for y = a + bx
        imatrix.append([1.0, d[b]]) # for x = a + by        
        
    coeffs = scipy.linalg.lstsq(matrix, map(lambda x: x[b], data))[0]
    icoeffs = scipy.linalg.lstsq(imatrix, map(lambda x: x[a], data))[0]    

    f = "%f + %f*x" % (coeffs[0], coeffs[1])

    r2 = coeffs[1] * icoeffs[1]

    return Gnuplot.Func( f, title="%s (r^2=%f)" % (f, r2))
##---------------------------------------------------------------------------------------------------------        
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      param_short_options,
                                      param_long_options)
                                      

    except getopt.error, msg:
        print USAGE, msg
        sys.exit(1)

    for o,a in optlist:
        if o in ("-l", "--legend"):
            param_legend = string.split( a, "," )
        elif o in ("-t", "--title"):
            param_title = a
        elif o in ("-c", "--hardcopy"):
            param_hardcopy = a
            if re.search( "\.png$", a):
                param_terminal = "png"
        elif o in ("-b", "--blackwhite"):
            param_blackwhite = 1
        elif o in ("-w", "--with"):
            param_with = a
        elif o in ("-f", "--fit"):
            param_fit = 1
        elif o in ("-u", "--function"):
            param_function = string.split(a,",")
        elif o in ("-i", "--files"):
            param_filenames = string.split(a, ",")
        elif o in ("-o", "--logscale"):
            param_logscale = a
            
    if len(args) > 0:
        print USAGE, "no arguments needed."
        sys.exit(1)
        
    print Experiment.GetHeader()
    print Experiment.GetParams()    

    if not param_hardcopy:
        g = Gnuplot.Gnuplot(debug=1, persist=1)
    else:
        g = Gnuplot.Gnuplot(debug=1)

    if param_filenames:
        filenames = param_filenames
    else:
        filenames = ["-"]

    if param_logscale:
        g("set logscale %s" % param_logscale)

    for filename in filenames:

        if filename == "-":
            lines = filter( lambda x: x[0] <> "#", sys.stdin.readlines())
        else:
            lines = filter( lambda x: x[0] <> "#", open(filename).readlines())
            
        if param_legend:
            data = map( lambda x: map(string.atof, string.split(x[:-1], "\t")), lines)
            legend = [param_legend[0]] +  param_legend[1:len(data[0])]
            del param_legend[1:len(data[0])]
        else:
            legend = string.split(lines[0][:-1], "\t")
            data = map( lambda x: map(string.atof, string.split(x[:-1], "\t")), lines[1:])        

        g.clear()
        
        if param_title:
            g.title( param_title )

        g.xlabel( legend[0] )

        for x in range(1, len(legend)):

            g.replot( Gnuplot.Data( data, cols=(0, x),
                                    xwith = param_with,
                                    title =legend[x]) )
            
            if param_fit:
                g.replot(Fit( data, cols=(0,x)))

        if param_function:
            for f in param_function:
                g.replot( Gnuplot.Func( f ) )

    g.refresh()

    if param_hardcopy:
        g.hardcopy( param_hardcopy,
                    terminal = param_terminal )




        
##         g.replot( Gnuplot.File( fn1,
##                                 using="1:2",
##                                 with = param_with,
##                                 title=legend[x]) )
#        print d.filename

##         temps += PlotFit( g, data, cols=(0, x) )

#    print d.filename

                  

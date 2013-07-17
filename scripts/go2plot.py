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
go2plot.py - 
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

   python go2plot.py --help

Type::

   python go2plot.py --help

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
import getopt
import time
import optparse
import math
import tempfile

""" program $Id: go2plot.py 2782 2009-09-10 11:40:29Z andreas $

"""
import CGAT.Experiment as E

doEps = True                 # generate eps instead of ps (but make sure it's just one page!)
largepvalues = True          # if true smallest p value is 10^-4 otherwise 10^-5.
annsize = 90                 # width of annotation text column
colsize = 9                  # width of column
colsep = 2                   # separation between (small) columns
rowheight = 7                # height of line
rowwidth = rowheight         # used for heading text only (?)
textangle = 45               # angle of heading text
labelsize = 70               # upper margin used for printing column labels
maincolsep = 15              # whitespace between (big) columns
greylines = 5                # number of lines per grey divider line
maxp = 1.0                   # only include annotations with p value < maxp in some input file
vertical = 0                 # if nonzero, prints a vertical grey line at (small) column 'vertical'

maxsize = rowheight-1
maxheight = 815 - labelsize - 50
maxwidth = 585 - labelsize/3 - 10

if largepvalues:
    # minimum 1/10,000
    pval0 = 0.0001
    pval1 = 0.0003
    pval2 = 0.001
    pval3 = 0.003
    pval4 = 0.01
    pval5 = 0.03
else:
    # minimum 1/100,000
    pval0 = 0.00001
    pval1 = 0.0001
    pval2 = 0.001
    pval3 = 0.01
    pval4 = 0.1
    pval5 = 0.3


def posx(x):
    return x*(colsize+colsep)+annsize

def posy(y):
    return (collen-y)*rowheight

def makeheader( input, annlist, numcols, numpages ):
    if doEps:
        print "%!PS-Adobe-2.0 EPSF-2.0"
    else:
        print "%!PS-Adobe-2.0"
        print "%%Pages: ",numpages
    print "%%BoundingBox: ", \
          10, \
          posy(collen+4)+60, \
          int(10+(numcols-1)*maincolsep + labelsize/3 + (numcols)*posx(len(input))),\
          int(posy(0)+labelsize)+60
    print """
/makebox {
    % x y size color
    gsave
    setrgbcolor
    /size exch def
    /y exch def
    /x exch def
     x y translate
     size size scale
     newpath
     -0.5 -0.5 moveto
     1 0 rlineto
     0 1 rlineto
     -1 0 rlineto
     0 -1 rlineto
     closepath
     fill
     grestore
} def
/min {/minarg1 exch def /minarg2 exch def
     minarg1 minarg2 lt {minarg1} {minarg2} ifelse}def
/max {/maxarg1 exch def /maxarg2 exch def
     maxarg1 maxarg2 gt {maxarg1} {maxarg2} ifelse}def
"""
    

def sameannotation(ann1,ann2):
    return ann1.replace('+','').replace('-','') == ann2.replace('+','').replace('-','')


def makecolheader(input, annlist, title = ""):
    print "/Times-Roman findfont ",colsize," scalefont setfont"
    for i in range(len(input)):
        print posx(i),posy(-1)," moveto gsave ",textangle," rotate"
        print "(",input[i],") show grestore"
    print "/Times-Roman findfont [",rowwidth,"0 0",rowheight,"0 0] makefont setfont"

    if title:
        print "gsave ", 10, posy(-8), "moveto ( %s ) show grestore" % (title)
        
    for i in range(len(annlist)):
        print "gsave 0 ",posy(i)," translate"
        print "newpath -5 -5 moveto ",annsize-5,0,"lineto",
        print annsize-5,5+rowheight,"lineto",0,5+rowheight,"lineto"
        print "closepath clip (",annlist[i],") 0 0 moveto show grestore"
        # helper lines
        if (greylines > 0 and (i+1)%greylines == 0) or \
           (greylines==0 and i<len(annlist)-1 and not sameannotation(annlist[i],annlist[i+1])):
            print "gsave 0.5 setgray"
            print "0.25 setlinewidth"
            print "[0.25 0.25] 0 setdash"
            print 0,posy(i),"moveto"
            print "(",annlist[i],") stringwidth pop ",annsize," min 0 rmoveto"
            print posx(len(input)),posy(i),"lineto stroke grestore"


def makecolor(color):
    if color == 0.5:
        print "0 0 0"
        return
    color = 1-color
    if color > 0.5:
        #red
        color = min(1,2*color-1)
        col0 = [0.1,0,0.5]
        col1 = [1,0,0]
    else:
        #green
        color = min(1,1-2*color)
        col0 = [0.1,0.2,0.05]
        col1 = [0.9,1,0.2]
    print col0[0]+color*(col1[0]-col0[0]),\
          col0[1]+color*(col1[1]-col0[1]),\
          col0[2]+color*(col1[2]-col0[2])


def makeps(i,j,element):
    (relrep,experp,pexp,ann) = element
    xpos = posx(i)
    if i==-1:
        xpos = 0
    logrelrep = math.log(1.0001 + (relrep / 100.0))
    color = 0.5 + (logrelrep / math.log(4.0))/2  # red=low, green=high
    if experp <= pval0:
        size = maxsize
    elif experp <= pval1:
        size = 0.8*maxsize
    elif experp <= pval2:
        size = 0.6*maxsize
    elif experp <= pval3:
        size = 0.45*maxsize
    elif experp <= pval4:
        size = 0.3*maxsize
    elif experp <= pval5:
        size = 0.2*maxsize
    else:
        size = 0
    if size > 0:
        print xpos,rowheight/2.0 + posy(j),size
        makecolor(color)
        print "makebox"
    else:
        print xpos,rowheight/2.0 + posy(j),1
        makecolor(0.5)
        print "makebox"
    if i==vertical:
        print "gsave 0.5 setgray"
        print "0.25 setlinewidth"
        print "[] 0 setdash"
        print posx(i)-colsize/2,posy(j),"moveto"
        print posx(i)-colsize/2,posy(j+1),"lineto stroke grestore"
    
def makefooter(lineno, footer = ""):
    print "/prtstr {dup stringwidth pop exch gsave 7",posy(lineno),"2 add rmoveto show grestore ",3+colsize," add 0 rmoveto currentpoint translate} def"
    print "/Times-Roman findfont ",colsize," scalefont setfont gsave"
    print "5 0 moveto"
    print "(P values: ) prtstr "
    makeps(-1,lineno,[-99,pval0,0,0])
    print "(<",pval0,") prtstr "
    makeps(-1,lineno,[-99,pval1,0,0])
    print "(<",pval1,") prtstr "
    makeps(-1,lineno,[-99,pval2,0,0])
    print "(<",pval2,") prtstr "
    makeps(-1,lineno,[-99,pval3,0,0])
    print "(<",pval3,") prtstr "
    makeps(-1,lineno,[-99,pval4,0,0])
    print "(<",pval4,") prtstr "
    makeps(-1,lineno,[-99,pval5,0,0])
    print "(<",pval5,") prtstr "
    print "grestore"
    # fold change
    print "gsave 0 ",-colsize,"translate"
    print "5 0 moveto"
    print "(Fold change:) prtstr "
    makeps(-1,lineno,[300,pval0,0,0])
    print "(>+4) prtstr"
    makeps(-1,lineno,[200,pval0,0,0])
    print "(+3) prtstr"
    makeps(-1,lineno,[100,pval0,0,0])
    print "(+2) prtstr"
    makeps(-1,lineno,[50,pval0,0,0])
    print ("(+1.5) prtstr")
    makeps(-1,lineno,[-33,pval0,0,0])
    print ("(-1.5) prtstr")
    makeps(-1,lineno,[-50,pval0,0,0])
    print "(-2) prtstr"
    makeps(-1,lineno,[-66,pval0,0,0])
    print "(-3) prtstr"
    makeps(-1,lineno,[-75,pval0,0,0])
    print "(<-4) prtstr"
    print "grestore"
    if footer:    
        print "gsave ", 0, posy(lineno + 3), "moveto ( %s ) show grestore" % (footer)
        
    
def makepageheader(page):
    print "%%Page: ",page,page
    print "gsave 11 60 translate"


def makepagetrailer(page):
    print "grestore"

def collect( infile, with_headers = False):
    """read input table."""

    data = []

    lines = filter(lambda x: x[0] != "#", infile.readlines())

    if len(lines) == 0: return data
    
    if with_headers:
        del lines[0]
        
    for line in lines:
        if line[0] == "#": continue

        try:
            (code, goid, scount, stotal, spercent, bcount, btotal, bpercent, pover, punder, goid, category, description) = line[:-1].split("\t")[:13]
        except ValueError:
            raise "# parsing error in line: %s" % line[:-1]
            sys.exit(1)
        
        if code == "+":
            p = pover
        else:
            p = punder
            
        data.append( (100.0 * ( float(spercent) / float(bpercent) - 1.0) , abs(float(p)), float(bpercent), description) )
        
    return data

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: go2plot.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-e", "--headers", dest="headers", action="store_true",
                      help="first row is a header [ignored]."  )
    parser.add_option("-t", "--title", dest="title", type="string",
                      help="page title.")
    parser.add_option("-f", "--footer", dest="footer", type="string",
                      help="page footer.")
    parser.add_option("-c", "--column-titles", dest="titles", type="string",
                      help="comma separated list of column titles [default: use filenames]."  )
    parser.add_option("-p", "--pattern-filename", dest="pattern_filename", type="string",
                      help="pattern to map columns to filename."  )

    parser.set_defaults(
        sortAlphabetically = True,
        headers = False,
        titles = "",
        pattern_filename = None,
        title = "",
        footer = "",        
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if len(args) == 0:
        raise "Please supply at least one input file."

    if options.pattern_filename:
        input = []
        titles = args
        for x in titles: input.append( options.pattern_filename % x )
    else:
        input = args
        
        if options.titles:
            titles=options.titles.split(",")
            if len(titles) != len(input):
                raise "Number of titles and files different: %i != %i" % (len(titles), len(input))            
        else:
            titles=input

    if options.loglevel >= 2:
        options.stderr.write( "# reading data from:\n" )
        for x in range(len(input)):
            options.stderr.write( "# %s: %s\n" % ( titles[x], input[x]) )
    
    data = []
    for filename in input:
        # collect all data
        try:
            values = collect( open(filename,"r"), with_headers = options.headers )
        except IOError:
            values = []
        data.append ( values )
        
    # collect all annotations
    # Also filter for max pvalue
    annotations = {}
    for i in range(len(input)):
        for line in data[i]:
            if line[1] < maxp:
                annotations[ line[3] ] = 0

    # Filter for pvalues and relreps in first two columns
    for ann in annotations:
        annotations[ann] = [1.0,1.0,1.0,1.0]

    # sort and filter annotations
    # (Code removed which did some filtering; the annotations data is not used)
    # By removing labels from annlist you can select the annotations you want to display
    annlist = [ ann for ann in annotations ]
    if options.sortAlphabetically:
        annlist.sort()

    # now make table

    # Number of columns
    columns = 1+int(len(annlist)*rowheight / maxheight)
    # Length of one column
    collen = int(maxheight / rowheight)
    # Max number of columns on a single page
    pagecols = int((maxwidth + maincolsep) / (posx(len(titles)) + maincolsep))
    if doEps:
        pagecols = columns
    # Number of pages needed
    numpages = 1 + int((columns-1)/pagecols)

    makeheader( titles, annlist, min(pagecols,columns), numpages )

    for page in range(numpages):
        if page == numpages - 1:
            colsonpage = columns - page*pagecols
        else:
            colsonpage = pagecols

        makepageheader(page+1)

        for column in range(colsonpage):
            globcolumn = page*pagecols + column
            colannlist = annlist[ globcolumn*collen:
                                  min((globcolumn+1)*collen, len(annlist)) ]
            print "gsave"
            print column * ( posx(len(titles)) + maincolsep )," 0 translate"

            makecolheader( titles, colannlist,
                           title = options.title )

            for i in range(len(titles)):
                for j in range(len(colannlist)):
                    # find element
                    element = 0
                    for line in data[i]:
                        if line[3] == colannlist[j]:
                            element = line
                    if element != 0:
                        # show
                        makeps(i,j,element)

            print "grestore"

        if colsonpage > 1:
            makefooter( collen+2, footer = options.footer)
        else:
            makefooter( len(colannlist)+2, footer = options.footer )

        makepagetrailer(page+1)

        if not doEps:
            print "showpage"





    E.Stop()

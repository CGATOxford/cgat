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
go2svg.py - plot results of a GO analysis
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script creates a colourful plot in :term:`svg` format from
the results of a GO analysis.

Usage
-----

Example::

   python go2svg.py --help

Type::

   python go2svg.py --help

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


import CGAT.Experiment as E
import CGAT.SVGdraw as SVGdraw
import bisect
import numpy
import CGAT.Stats as Stats

class DataPoint:
    def __init__(self):
        pass

class DataFDR:
    def __init__(self):
        pass

def Collect( infile, 
             with_headers = False, 
             annotator_format = False, 
             use_annotator_fdr = False, 
             delims = "",
             ignore = "",
             max_pvalue = 1.0,
             max_qvalue = None):
    """read input table."""

    data = []

    lines = filter(lambda x: x[0] != "#", infile.readlines())

    if len(lines) == 0: return data
    
    if with_headers:
        del lines[0]

    if annotator_format:        

        lines = [ line for line in lines if not line.startswith( "Iteration" ) ]
        annotator_fdr = {}
        annotator_level = None
        for line in lines:
            if len(line) == 1: continue                  # skip trailing blank lines
            
            if line.startswith( "--" ):
                if line.startswith("-- False"):
                    annotator_level = float(re.search( "-- False Discovery summary for p-value (.+):", line ).groups()[0])
                    annotator_fdr[annotator_level] = {}
                elif line.startswith("--  Category"):
                    pass
                else:
                    if re.search("insufficiently", line): continue
                    dd = re.split("\s+", line[4:-1])
                    d = DataFDR()
                    d.mObserved, d.mAverage, d.mMedian, d.m95 = map(float,dd[1:])
                    annotator_fdr[annotator_level][dd[0]] = d
                continue
            else:
                if line[0] == "Z": continue # skip header
                if len(line[:-1].split('\t')) != 9: continue # HACK: accounts for a bug in Annotator output

                try:
                    (z, percentchange, pvalue, observed, expected, low95, up95, stddev, description) = line[:-1].split('\t')[:9]
                except ValueError:
                    raise ValueError("# parsing error in line: %s" % line[:-1])

            d = DataPoint()
            d.mAnnotation = description
            d.mPValue     = float(pvalue)
            d.mFoldChange = 1.0 + float(percentchange) / 100.0
            data.append(d)
    else:
        
        for line in lines:
            try:
                (code, goid, scount, stotal, spercent, bcount, btotal, bpercent, ratio, pover, punder, goid, category, description) = line[:-1].split("\t")[:14]
            except ValueError:
                raise ValueError("# parsing error in line: %s" % line[:-1])
            
            if code == "+":
                p = pover
            else:
                p = punder

            d = DataPoint()
            d.mAnnotation = description        
            d.mPValue     = float(p)
            d.mFoldChange = float(spercent) / float(bpercent)
            data.append(d)

    # apply filters
    for c in delims:
        for d in data:
            d.mAnnotation = d.mAnnotation.split(c)[0]
    for c in ignore:
        for d in data:
            d.mAnnotation = d.mAnnotation.replace(c, '')

    ninput = len(data)
    no_fdr = False
    ## apply filters

    if ninput > 0:
        if max_qvalue != None:
            if use_annotator_fdr:
                pvalues = annotator_fdr.keys()
                pvalues.sort()
                pvalues.reverse()
                for pvalue in pvalues:
                    try:
                        d = annotator_fdr[pvalue]["Significant"]
                    except KeyError:
                        continue
                    if d.mObserved == 0:
                        E.info("no data remaining after fdr filtering" )
                        data = []
                        break
                    elif d.mAverage / d.mObserved < max_qvalue:
                        E.info("filtering with P-value of %f" % pvalue )
                        data = [ x for x in data if x.mPValue < pvalue ]
                        break
                else:
                    E.warn( "fdr could not be computed - compute more samples (at P = %f, actual fdr=%f)" % (pvalue, d.mAverage / d.mObserved) )
                    no_fdr = True

            if no_fdr:
                if use_annotator_fdr: E.info( "estimating FDR from observed P-Values" )

                pvalues = [x.mPValue for x in data]
                vlambda = numpy.arange(0, max( pvalues), 0.05 )
                try:
                    qvalues = Stats.doFDR( pvalues, vlambda = vlambda, fdr_level = max_qvalue )
                except ValueError, msg:
                    E.warn( "fdr could not be computed - no filtering: %s" % msg)
                    no_fdr = True
                else:
                    data = [ x[0] for x in zip( data, qvalues.mPassed ) if x[1] ]
        elif max_pvalue != None:
            data = [ x for x in data if x.mPValue < max_pvalue ]

    if no_fdr: data = []

    nremoved = ninput - len(data)

    return data, nremoved, no_fdr

# some definitions for the layout of the picture
BLACK = (0,0,0)
WHITE = (255,255,255)
RED   = (255,0,0)
GREEN = (0,255,0)
BLUE  = (0,0,255)
YELLOW = (255,255,0)
CYAN   = (0,255,255)
PURPLE = (255,0,255)
GREY = (128,128,128)

MAX_GREY = 240
GREY_COLORS = map( lambda x: (x,x,x), range(0,MAX_GREY))

DEFAULT_XWIDTH = 500
DEFAULT_YWIDTH = 500
DEFAULT_LINE_DISTANCE = 10

DEFAULT_OFFSET_Y = 20
DEFAULT_BOX_WIDTH = 15
DEFAULT_LINE_SPACING = 40
DEFAULT_TEXT_DISTANCE = 40
DEFAULT_ANNOTATION_DISTANCE = 20
DEFAULT_LINE_WIDTH = 2

DEFAULT_SCALE_FACTOR = 1.0

DEFAULT_COLOR_SEQUENCE   = BLACK
DEFAULT_COLOR_ANNOTATION = BLACK

class GoPlot:

    def __init__(self,
                 row_names,
                 col_names,
                 thresholds_size,
                 thresholds_colour,
                 alt_colours,
                 template = "screen",
                 max_pvalue = None,
                 max_qvalue = None,
                 mark_columns = None,
                 ):
        """If max_pvalue and max_qvalue are given, they are
        added to the footer.

        If mark_columns is given, the columns for which
        mark_columns is True are marked.

        """

        self.mElements = []
        self.mTemplate = template
        
        ## save row names
        self.mRowNames = row_names
        self.mColNames = col_names

        ## colours
        self.startColour  = RED
        self.mMiddleColour = YELLOW
        self.mStopColour   = GREEN
        self.mAltColours = alt_colours

        ## a space
        self.mSeparator = 10
        
        ## info
        self.mMaxPValue = max_pvalue
        self.mMaxQValue = max_qvalue

        ## width and height of a row/column
        self.mRowHeight = 150                               # GAL changed; was 200
        self.mColWidth  = 200

        ## maximum size of a box
        self.mMaxBoxSize = 140
        self.mRevertSize = True

        ## Height and width get set
        self.mHeaderFontSize  = self.mRowHeight / 1.5       # GAL changed; was 2.0
        self.mHeaderFont      = "Verdana"
        
        self.mThresholdsSize = thresholds_size
        self.mThresholdsSizeTitle = "P-Values"
        self.mThresholdsColour = thresholds_colour
        self.mThresholdsColourTitle = "Fold change"

        ## write a grid line every five rows
        self.mRowTicks = 5
        self.mColTicks = 5

        ## footer
        self.mFooterFrom = 10        
        self.mFooterFontSize  = self.mRowHeight / 1.5
        self.mFooterFont      = "Verdana"
        self.mFooter = None

        ## page margin
        self.mBottomMargin = 300
        
        ## Title
        self.mTitleFontSize  = self.mRowHeight
        self.mTitleFont      = "Verdana"
        self.mTitle = None

        if self.mTemplate == "screen":
            ## screen is default
            pass
        elif self.mTemplate == "publication":
            self.startColour  = RED
            self.mMiddleColour = WHITE
            self.mStopColour   = BLUE

        self.mMarkColumns = mark_columns
        if self.mMarkColumns:
            assert len(self.mMarkColumns) == len(self.mColNames), "length of mark_columns must equal length of columns"
            
        self.initializePlot()
        
    #####################################################################                
    def initializePlot( self ):
        """set various coordinates in the plot."""

        ## height of col header
        self.mHeaderHeight = max( map( len, self.mColNames) ) * self.mHeaderFontSize / 2

        if self.mTitle:
            self.mHeaderHeight += self.mSeparator + self.mTitleFontSize
            
        ## width of row header
        self.mHeaderWidth = max( map( len, self.mRowNames) ) * self.mHeaderFontSize / 2
        
        ## height of footer:
        self.mFooterHeight =  2 * max( self.mFooterFontSize, self.mMaxBoxSize) + self.mSeparator

        if self.mFooter:
            self.mFooterHeight += self.mSeparator + self.mFooterFontSize
        
        ## height and width of data section
        self.mDataWidth  = len(col_names) * self.mColWidth 
        self.mDataHeight = len(row_names) * self.mRowHeight

        ## build coordinates for rows and columns
        self.buildMapRow2Position()
        self.buildMapCol2Position()

        ## build colour map
        self.buildColourMap()

        self.mPageWidth = self.mHeaderWidth + self.mDataWidth + 1 * self.mSeparator
        self.mPageHeight = self.mHeaderHeight + self.mDataHeight + self.mFooterHeight + 2 * self.mSeparator + self.mBottomMargin

    #####################################################################
    def setTitle( self, title ):
        """set title."""
        self.mTitle = title

    #####################################################################
    def setFooter( self, footer ):
        """set footer."""
        self.mFooter = footer
        
    #####################################################################
    def writeTitle( self ):
        """write title into plot."""
        
        if self.mTitle:
            e = SVGdraw.text( self.mPageWidth / 2,
                              self.mTitleFontSize ,
                              self.mTitle,
                              self.mTitleFontSize,
                              self.mTitleFont,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              text_anchor = "middle" )

            self.mElements.append(e)

    #####################################################################        
    def buildColourMap( self ):
        """build map of thresholds to colours.

        This is two gradients:
        1: first half: start to middle
        2: second half: middle to end
        """

        self.mColours = []

        if self.mAltColours:

            for colidx in range(len(self.mThresholdsColour)+1):

                col = colidx / (len(self.mThresholdsColour) - 0.0)

                if col == 0.5:
                    x = 0
                    col0 = [0,0,0]
                    col1 = [0,0,0]
                elif col < 0.5:
                    x = min(1,1-2*col)
                    col0 = [26,0,128]
                    col1 = [255,0,0]
                else:
                    x = min(1,2*col-1)
                    col0 = [26,52,13]
                    col1 = [230,255,52]
                self.mColours.append( (col0[0] + x*(col1[0]-col0[0]),
                                       col0[1] + x*(col1[1]-col0[1]),
                                       col0[2] + x*(col1[2]-col0[2]) ) )

        else:

            num_steps = int(math.floor( (len(self.mThresholdsColour) + 1) / 2.0)) 
            
            d = map( lambda x, y: (x - y) / num_steps, self.mMiddleColour, self.startColour)
            for x in range(num_steps):
                self.mColours.append ( (self.startColour[0] + x * d[0],
                                        self.startColour[1] + x * d[1],
                                        self.startColour[2] + x * d[2] ) )

            # self.mColours.append( self.mMiddleColour )
            
            d = map( lambda x, y: (x - y) / num_steps, self.mStopColour, self.mMiddleColour)
            for x in range(1, num_steps):
                self.mColours.append ( (self.mMiddleColour[0] + x * d[0],
                                        self.mMiddleColour[1] + x * d[1],
                                        self.mMiddleColour[2] + x * d[2] ) )

            self.mColours.append( self.mStopColour )

    #####################################################################        
    def buildMapRow2Position( self ):
        
        ## build map of row_name to row
        self.mMapRow2Position = {}

        offset = self.mHeaderHeight + self.mSeparator
        for x in range(len(self.mRowNames)):
            self.mMapRow2Position[self.mRowNames[x]] = offset + x * self.mRowHeight

    #####################################################################        
    def buildMapCol2Position( self ):

        ## build map of row_name to row
        self.mMapCol2Position = {}

        for x in range(len(self.mColNames)):
            self.mMapCol2Position[self.mColNames[x]] = x * self.mColWidth

    #####################################################################
    def addValue( self, row, col, size, colour_value ):
        """add a dot in row/col.
        """

        ## decide the size of the box
        pos = bisect.bisect( self.mThresholdsSize, size )

        if self.mRevertSize: 
            size = self.mMaxBoxSize * (1.0 - float(pos) / len(self.mThresholdsSize))
        else:
            size = self.mMaxBoxSize * float(pos) / len(self.mThresholdsSize)

        d = (self.mMaxBoxSize - size) / 2
        
        x = self.mMapCol2Position[col] + d
        try: 
            y = self.mMapRow2Position[row] + d
        except KeyError:
            return

        # determine the colour of the box
        pos = bisect.bisect( self.mThresholdsColour, colour_value )
        colour = self.mColours[pos]

        e = SVGdraw.rect( x, y,
                          size, size, 
                          stroke = "black",
                          fill = "rgb(%i,%i,%i)" % colour )

        self.mElements.append(e)
        
    #####################################################################        
    def writeColHeaders( self ):
        """write row headers."""

        current_x = self.mColWidth / 2
        current_y = self.mHeaderHeight
        for i in range( len(self.mColNames) ):
            if self.mMarkColumns and self.mMarkColumns[i]:
                color = BLUE
                name = self.mColNames[i] + "*"
            else:
                color = BLACK
                name = self.mColNames[i]
            
            e = SVGdraw.text( current_x,
                              current_y,
                              name,
                              self.mHeaderFontSize,
                              self.mHeaderFont,
                              stroke = "rgb(%i,%i,%i)" % color,
                              text_anchor = "start",
                              transform="rotate(-45,%i,%i)" %( current_x, current_y ))

            self.mElements.append(e)
            
            current_x += self.mColWidth
            # current_y -= self.mColWidth / 2        # GAL added # AH removed?

    #####################################################################
    def writeGrid( self ):
        """add grid lines."""

        if self.mRowTicks:

            start_x = 0,
            end_x = self.mDataWidth + self.mSeparator + self.mHeaderWidth
            current_y = self.mHeaderHeight + self.mSeparator / 2 + self.mRowTicks * self.mRowHeight
            
            for x in range(self.mRowTicks, len(self.mRowNames), self.mRowTicks):

                e = SVGdraw.line( start_x,
                                  current_y,
                                  end_x,
                                  current_y,
                                  stroke = "rgb(%i,%i,%i)" % BLACK )
                
                self.mElements.append(e)
                
                current_y += self.mRowTicks * self.mRowHeight

        if self.mColTicks:            

            start_y = self.mHeaderHeight + self.mSeparator / 2
            end_y = start_y + self.mDataHeight
            
            current_x = self.mColTicks * self.mColWidth - self.mColWidth / 2
            
            for x in range(self.mColTicks, len(self.mColNames), self.mColTicks):

                e = SVGdraw.line( current_x,
                                  start_y,
                                  current_x,
                                  end_y,
                                  stroke = "rgb(%i,%i,%i)" % BLACK )
                
                self.mElements.append(e)
                
                current_x += self.mColTicks * self.mColWidth
            
    #####################################################################
    def writeRowHeaders( self ):
        """write row headers."""

        current_x = self.mDataWidth + self.mSeparator
        
        current_y = self.mHeaderHeight + self.mSeparator + self.mHeaderFontSize

        for name in self.mRowNames:
            e = SVGdraw.text( current_x,
                              current_y,
                              name,
                              self.mHeaderFontSize,
                              self.mHeaderFont,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              text_anchor = "start")
            
            self.mElements.append(e)
            
            current_y += self.mRowHeight

        self.mHeaderWidth = max( map( len, self.mRowNames) ) * self.mHeaderFontSize / 2
            
    #####################################################################
    def writeFooter( self ):
        """write footer.

        The footer contains the legend.
        """
        current_x = self.mFooterFrom
        current_y = self.mHeaderHeight + self.mDataHeight + 2 * self.mSeparator

        ###########################################################
        ## Draw legend 1: size of boxes
        e = SVGdraw.text( current_x,
                          current_y + self.mFooterFontSize,
                          self.mThresholdsSizeTitle,
                          self.mFooterFontSize,
                          self.mFooterFont,
                          stroke = "rgb(%i,%i,%i)" % BLACK,
                          text_anchor = "start")

        current_x += len(self.mThresholdsSizeTitle) * self.mFooterFontSize / 1.5 + self.mSeparator
        
        self.mElements.append(e)
        
        l = len(self.mThresholdsSize)
        for x in range(l):
            if self.mRevertSize:
                p = int(self.mMaxBoxSize * (1.0 - float(x) / l))
            else:
                p = int(self.mMaxBoxSize * (float(x) / l))
                
            e = SVGdraw.rect( current_x,
                              current_y + (self.mMaxBoxSize - p) / 2,
                              p, p,
                              stroke = "black",
                              fill = "rgb(%i,%i,%i)" % self.startColour )

            self.mElements.append(e)            
            current_x += self.mMaxBoxSize + self.mSeparator 

            t = "< %g" % (self.mThresholdsSize[x])
            e = SVGdraw.text( current_x,
                              current_y + self.mFooterFontSize,
                              t,
                              self.mFooterFontSize,
                              self.mFooterFont,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              text_anchor = "start")                              

            current_x += len(t) * self.mFooterFontSize / 1.5 + self.mSeparator 
            self.mElements.append(e)

        ###########################################################
        ## Draw legend 2: colour of boxes
        current_x = self.mFooterFrom            
        current_y += max( self.mFooterFontSize, self.mMaxBoxSize) + self.mSeparator

        e = SVGdraw.text( current_x,
                          current_y + self.mFooterFontSize,
                          self.mThresholdsColourTitle,
                          self.mFooterFontSize,
                          self.mFooterFont,
                          stroke = "rgb(%i,%i,%i)" % BLACK,
                          text_anchor = "start")

        current_x += len(self.mThresholdsColourTitle) * self.mFooterFontSize / 1.5 + self.mSeparator
        
        self.mElements.append(e)
        
        l = len(self.mThresholdsColour)

        for x in range(l+1):
            
            p = self.mMaxBoxSize

            if x < l:
                t = "< %g" % (self.mThresholdsColour[x])
            else:
                t = "> %g" % (self.mThresholdsColour[x-1])
                
            e = SVGdraw.rect( current_x,
                              current_y,
                              p, p,
                              stroke = "black",
                              fill = "rgb(%i,%i,%i)" % self.mColours[x] )

            self.mElements.append(e)            
            current_x += self.mMaxBoxSize + self.mSeparator 

            e = SVGdraw.text( current_x,
                              current_y + self.mFooterFontSize,
                              t,
                              self.mFooterFontSize,
                              self.mFooterFont,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              text_anchor = "start")                              


            current_x += len(t) * self.mFooterFontSize / 1.5 + self.mSeparator
            self.mElements.append(e)

        ###########################################################
        if self.mMaxPValue != None or self.mMaxQValue != None:
            current_y += max( self.mFooterFontSize / 1.5 , self.mMaxBoxSize) + self.mSeparator
            
            a = []
            if self.mMaxPValue: a.append( "P < %6.4f" % self.mMaxPValue )
            if self.mMaxQValue: a.append( "FDR = %6.4f" % self.mMaxQValue )

            e = SVGdraw.text( self.mPageWidth / 2,
                              current_y + self.mFooterFontSize,
                              " ".join( a ),
                              self.mFooterFontSize,
                              self.mFooterFont,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              text_anchor = "middle")

        ###########################################################
        if self.mFooter:

            current_y += max( self.mFooterFontSize / 1.5 , self.mMaxBoxSize) + self.mSeparator
            
            e = SVGdraw.text( self.mPageWidth / 2,
                              current_y + self.mFooterFontSize,
                              self.mFooter,
                              self.mFooterFontSize,
                              self.mFooterFont,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              text_anchor = "middle")
            
            self.mElements.append(e)
            
    #####################################################################        
    def finalizePlot( self ):
        """write remaining parts of the plot."""
        
        self.writeGrid()
        self.writeTitle()
        self.writeFooter()
        self.writeRowHeaders()
        self.writeColHeaders()        

    #####################################################################        
    def writeToFile(self, outfile):
        """write svg image to file.
        """
        self.finalizePlot()

        self.mRoot = SVGdraw.drawing()
        self.mDraw = SVGdraw.svg( (0, 0, self.mPageWidth, self.mPageHeight ) , "100%", "100%" )
        
        for e in self.mElements:
            self.mDraw.addElement( e )
            
        self.mRoot.setSVG(self.mDraw)

        tfile = tempfile.mktemp()
        
        self.mRoot.toXml( tfile )

        lines = open(tfile,"r").readlines()
        
        outfile.write(string.join(lines,""))
        outfile.write("\n")
        
        os.remove(tfile)

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: go2svg.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-e", "--headers", dest="headers", action="store_true",
                      help="first row is a header [ignored]."  )
    parser.add_option("-t", "--title", dest="title", type="string",
                      help="page title.")
    parser.add_option("-f", "--footer", dest="footer", type="string",
                      help="page footer.")
    parser.add_option( "--maxP", dest="max_pvalue", type="float",
                      help="maximum P-value displayed [default=%default].")
    parser.add_option( "--maxQ", dest="max_qvalue", type="float",
                      help="maximum Q-value for controlling for FDR [default=%default].")
    parser.add_option("-c", "--column-titles", dest="col_names", type="string",
                      help="comma separated list of column titles [default: use filenames]."  )
    parser.add_option("-p", "--pattern-filename", dest="pattern_filename", type="string",
                      help="pattern to map columns to filename."  )
    parser.add_option("-A", "--Annotator", dest="annotator", action="store_true",
                      help="use Annotator-style input files.")
    parser.add_option( "--annotator-fdr", dest="annotator_fdr", action="store_true",
                      help="use fdr computed from annotator [default=%default].")
    parser.add_option("-T", "--thresholds", dest="thresholds", type="string",
                      help="7 comma-separated fold-change threshold values")
    parser.add_option("-P", "--pvalues", dest="pvalues", type="string",
                      help="6 comma-separated p value threshold values"),
    parser.add_option("-C", "--altcolours", dest="altcolours", action="store_true",
                      help="Use alternative colour palette")
    parser.add_option("-X", "--delimiters", dest="delims", type="string",
                      help="Delimiter characters for annotation label")
    parser.add_option("-Z", "--ignore", dest="ignore", type="string",
                      help="Ignored characters in annotation label")
    parser.add_option( "--fdr", dest="fdr", type="float",
                      help="filter output by FDR (requires annotator output). [default=%default]")
    parser.add_option("-a", "--template", dest="template", type="choice",
                      choices=("screen", "publication" ),
                      help="layout template to choose - affects colours.")
    parser.add_option( "--sort-columns", dest="sort_columns", type="choice",
                       choices=( "unsorted", "similarity", "alphabetical", ),
                       help="sort columns. The default, unsorted, list columns in the order that they are supplied on the command line [default=%default]")

    parser.set_defaults(
        sortAlphabetically = True,
        headers = False,
        col_names = "",
        pattern_filename = None,
        title = "",
        footer = "",
        max_pvalue = None,
        max_qvalue = None,
        annotator = False,
        thresholds = "0.25,0.33,0.5,1.0,2.0,3.0,4.0",
        pvalues = "0.00001,0.0001,0.001,0.01,0.1",
        altcolours = False,
        delims = "",
        ignore = "",
        template = "screen",
        annotator_fdr = False,
        fdr=None,
        sort_columns = "unsorted",
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if len(args) == 0:
        raise IOError("Please supply at least one input file.")

    if options.pattern_filename:
        input = []
        col_names = args
        for x in col_names: input.append( options.pattern_filename % x )
    else:
        input = args
        
        if options.col_names:
            col_names=options.col_names.split(",")
            if len(col_names) != len(input):
                raise ValueError("Number of col_names and files different: %i != %i" % (len(col_names), len(input)))
        else:
            col_names=input

    E.info( "reading data for %i columns" % len(input))

    columns = []
    errors = []
    
    for col_name, filename in zip( col_names, input ):

        E.debug( "reading data for column %s from %s " % (col_name, filename ))
        # collect all columns
        try:
            values, nremoved, no_fdr = Collect( open(filename,"r"),
                                                with_headers = options.headers,
                                                annotator_format = options.annotator,
                                                delims = options.delims,
                                                ignore = options.ignore,
                                                use_annotator_fdr = options.annotator_fdr,
                                                max_pvalue = options.max_pvalue,
                                                max_qvalue = options.max_qvalue )
        except IOError:
            E.warn( "no data from %s" % filename )
            values = []
            no_fdr = False
            nremoved = 0

        E.info("read %i values from %s: %i significant, %i removed" % (len(values) + nremoved, filename, 
                                                                       len(values),
                                                                       nremoved))
        columns.append( (col_name, values) )
        errors.append( no_fdr )

    if sum( [ len(x) for x in columns]) == 0:
        raise IOError("no data read - please check supplied files.")

    # collect all annotations
    # Also filter for max pvalue
    annotations = set()
    for col_name, column in columns:
        for d in column:
            annotations.add( d.mAnnotation )
                
    E.info( "There are %i rows" % len(annotations) )

    # sort and filter annotations
    # (Code removed which did some filtering; the annotations data is not used)
    # By removing labels from annlist you can select the annotations you want to display
    row_names = list(annotations)
    if options.sortAlphabetically:
        row_names.sort()


    if options.sort_columns == "unsorted":
        pass
    elif options.sort_columns == "alphabetical":
        col_names.sort()
    elif options.sort_columns == "similarity":
        if len(row_names) * len(col_names) > 10000: 
            E.info( "no sorting as matrix too large" )
        else:
            import CorrespondenceAnalysis
            matrix = numpy.ones( (len(row_names), len(col_names)), numpy.float )
            map_rows = dict( zip( row_names, range(len(row_names) ) ) )
            x = 0
            for col_name, column in columns:
                for d in column:            
                    matrix[map_rows[d.mAnnotation],x] = d.mFoldChange
                x += 1
            row_indices, col_indices = CorrespondenceAnalysis.GetIndices( matrix )
            map_row_new2old = numpy.argsort(row_indices)
            map_col_new2old = numpy.argsort(col_indices)
            row_names = [ row_names[map_row_new2old[x]] for x in range(len(row_names))]
            col_names = [ col_names[map_col_new2old[x]] for x in range(len(col_names))]

    E.info( "columns have been sorted" )

    plot = GoPlot( row_names,
                   col_names,
                   thresholds_size = tuple( map( float, options.pvalues.split(',') ) ),
                   thresholds_colour = tuple( map(float, options.thresholds.split(',') ) ),
                   template = options.template, 
                   alt_colours = options.altcolours,
                   max_pvalue = options.max_pvalue, 
                   max_qvalue = options.max_qvalue,
                   mark_columns = errors )

    if options.title: plot.setTitle(options.title)
    if options.footer: plot.setFooter(options.footer)

    plot.initializePlot()
        
    for col_name, column in columns:
        for d in column:
            plot.addValue( d.mAnnotation, 
                           col_name,
                           d.mPValue, 
                           d.mFoldChange )
        
    plot.writeToFile(options.stdout)
    
    E.Stop()

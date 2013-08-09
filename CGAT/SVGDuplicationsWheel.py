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
SVGDuplicationsWheel.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, getopt, time, optparse, math, tempfile, bisect

USAGE="""python optic/plot_duplications.py [OPTIONS] < stdin > stdout

plots a wheel plot of duplications in SVG format.
The script requirse the following input files:

--contig-sizes: a filename with contig sizes to output. This is a tab-separated
    file with contig<tab>size on each line.

Example for contig sizes:
chr1    1000
chr2    1000
chr3    1000

The duplications are input via standard input as a tab-separated table containing
the following three fields:

   1. cluster_id: an identifier for the cluster
   2. locations: a ';' separated list of duplicated genes with their locations. Each
        entry is of the following format:
        gene_id:chr:strand:start:end
   3. tree: the tree of duplicated genes. This is used to color arcs by tree height.
   The trees are in newick format.

Example for standard input:

1       gene1:chr1:+:100:200;gene2:chr2:+:100:200       ((gene1:0.5,gene2:0.5):0.5);
2       gene3:chr1:+:300:400;gene4:chr3:+:300:400       ((gene3:0.2,gene4:0.2):0.2);

Output: the svg formatted file is writtent to standard output. Comment lines
start with a '#'. You can remove this by setting -v 0 or --verbose=0.
"""

import Experiment
import IOTools
import SVGdraw 
import TreeTools

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
ORANGE = (255,125,0)

class DuplicationPlot:

    def __init__(self,
                 contigs,
                 map_contig2size,
                 num_entries,
                 template = "screen",
                 ):

        self.mElements = {}
        
        ## save contig sizes
        self.contigs = contigs
        self.mMapContig2Size = map_contig2size

        ## get row size
        self.mNumEntries = num_entries

        self.mRadius = 3000
        self.mRadiusIncrement = 40
        self.mRadiusStart = self.mRadius

        self.mMinValue = 0.0
        self.mMaxValue = 1.0

        ## planes
        self.mPlaneJoins = 4
        self.mPlaneLinks = 5
        self.mPlaneGrid = 10

        ## set of previous points
        self.mPreviousPoints = set()

        ## whether or not coordinates are forward coordinates
        self.mForwardCoordinates = False

        self.mTemplate = template
        
    def setFormattingOptions( self):
        """set formatting options according to template.

        This function is called, after basic dimensions have been set
        like the size of the panel and radius increments.
        """

        ## a space
        self.mSeparator = 10
        
        ## colours for duplications
        ## coloured by distance
        self.mLinkStrokeWidthSymbol = 1

        ## radial connection (same contig)
        self.mLinkRadStrokeWidth = 10
        ## arc connection (different contigs)
        self.mLinkArcStrokeWidth = 10
        
        self.mLinkStrokeWidthMasked = 1
        self.mLinkColourMasked = GREY
        self.mLinkColourSymbolMasked = GREY
        
        self.startColour  = RED
        self.mMiddleColour = YELLOW
        self.mStopColour   = GREEN
        
        ## size of symbol
        self.mLinkSymbolSize = self.mRadiusIncrement / 4

        ## Parameters for the grid
        self.mGridStrokeWidth = 1            
        self.mGridColour = GREY

        ## maximum size of a box
        self.mMaxBoxSize = 100
        self.mTextSize   = 100

        ## Height and width get set
        self.mHeaderFontSize  = self.mTextSize
        self.mHeaderFont      = "Verdana"

        self.mGridFontSize    = self.mRadius / 10
        self.mGridFont        = "Verdana"

        ## Scale
        self.mScaleFontSize   = self.mRadius / 10
        self.mScaleFont       = "Verdana"
        self.mScaleNumTicks = 20
        self.mScaleBoxSizeY = 30
        
        ## Footer
        self.mFooterFontSize  = self.mTextSize
        self.mFooterFont      = "Verdana"
        self.mFooter = None

        ## Title
        self.mTitleFontSize  = self.mTextSize
        self.mTitleFont      = "Verdana"
        self.mTitle = None

        self.mAngleResolution = 1

        if self.mTemplate == "screen":
            ## screen is default
            pass
        elif self.mTemplate == "publication":
            ## for publication: increase line widths and
            ## symbol size by a factor
            factor = 4
                    
            self.mAngleResolution = 2

            self.mLinkSymbolSize  *= factor
            self.mLinkRadStrokeWidth *= factor
            self.mLinkArcStrokeWidth *= factor            
            self.mGridStrokeWidth = self.mLinkRadStrokeWidth / 2

            self.startColour  = RED
            self.mMiddleColour = ORANGE
            self.mStopColour   = YELLOW

        elif self.mTemplate == "poster":
            ## for publication: increase line widths and
            ## symbol size by a factor
            ## use red -> blue as links
            ## reduce thickness of arcs
            factor = 4
                    
            self.mAngleResolution = 2

            self.mLinkSymbolSize  *= factor
            self.mLinkRadStrokeWidth *= factor / 2
            self.mLinkArcStrokeWidth *= factor / 3
            self.mGridStrokeWidth = self.mLinkRadStrokeWidth / 2

            self.mGridColour   = BLACK
            self.startColour  = RED
            self.mMiddleColour = None
            self.mStopColour   = BLUE

        ## elements within wheel
        self.mWheelElements = []
        
    #####################################################################
    def addElement( self, element, plane = 0):
        """add element to list in plane.
        
        The list is later sorted by priority.
        """
        if plane not in self.mElements:
            self.mElements[plane] = []
        self.mElements[plane].append( element )

    #####################################################################
    def addWheelElement( self, element, plane = 0):
        self.mWheelElements.append( (element, plane ) )
        
    #####################################################################                        
    def initializePlot( self ):
        """set dimensions of the plot.

        Header is left and top, footer is right and bottom.
        """

        ## set formatting options according to user input
        self.setFormattingOptions( )

        ## build coordinates for rows and columns
        self.buildMap2Position()

        ## build colour map
        self.buildColourMap()

        self.mTotalSize = sum( [self.mMapContig2Size[x] for x in self.contigs] )

        ## smallest radius
        self.mRadiusStart = self.mRadius

        ## radius to fall back on
        self.mRadiusFallBack = self.mRadius
        
        self.mDataMiddleX, self.mDataMiddleY = 0, 0

        ## maximum coordinate of longest range for fallback
        self.mLastMax = 0

        ## maximum used radius in slice
        self.mRadiusMax = self.mRadius
        
        ## maximum coordinate of previous range
        self.mPreviousMax = 0
        
        self.mLastChr = None

        ## keeps track of written separators, so that they do not get written multiply
        self.mSeparators = {}

        ## minimum distance for resolving dots
        self.mMinDistance = self.mTotalSize / 360.0 

        ## whether or not last entry was cis
        self.mPreviousCis = False

        ## number format for legend
        self.mFormatNumberLegend = "%5.2f"

        
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

            self.addElement(e)

    #####################################################################        
    def buildColourMap( self ):
        """build map of heights to colours.

        If self.mMiddleColour is defined, the gradient will be built in two steps:
        1: first half: start to middle
        2: second half: middle to end

        Otherwise, only one gradient is built.
        """

        self.mColours = []
        self.mColourThresholds = []
        value = self.mMinValue
        num_intervalls = 200
        
        if self.mMiddleColour:
            num_steps = num_intervalls / 2
            
            self.mColourIncrement = float((self.mMaxValue - self.mMinValue)) / num_steps / 2

            d = map( lambda x, y: (x - y) / num_steps, self.mMiddleColour, self.startColour)
            for x in range(num_steps):
                self.mColours.append ( (self.startColour[0] + x * d[0],
                                        self.startColour[1] + x * d[1],
                                        self.startColour[2] + x * d[2] ) )
                value += self.mColourIncrement
                self.mColourThresholds.append( value )

            self.mColours.append( self.mMiddleColour )
            value += self.mColourIncrement
            self.mColourThresholds.append( value )

            start_colour = self.mMiddleColour
        else:
            num_steps = num_intervalls             
            start_colour = self.startColour
            self.mColourIncrement = float((self.mMaxValue - self.mMinValue)) / num_steps 

        d = map( lambda x, y: (x - y) / num_steps, self.mStopColour, start_colour)
        for x in range(num_steps):
            self.mColours.append ( (start_colour[0] + x * d[0],
                                    start_colour[1] + x * d[1],
                                    start_colour[2] + x * d[2] ) )

            value += self.mColourIncrement
            self.mColourThresholds.append( value )

        self.mColours.append( self.mStopColour )

    #####################################################################        
    def buildMap2Position( self ):
        """build map of contigs to start."""
        self.contig2Position = {}

        l = 0
        for x in self.contigs:
            self.contig2Position[x] = l
            l += self.mMapContig2Size[x]

    #####################################################################
    def addSeparator( self ):
        """add separator on circles."""
        if self.mRadius not in self.mSeparators:
            e = SVGdraw.circle( self.mDataMiddleX, self.mDataMiddleY, self.mRadius, fill="none",
                                stroke = "rgb(%i,%i,%i)" % self.mGridColour,
                                stroke_width = self.mGridStrokeWidth )
            self.addWheelElement( e, self.mPlaneGrid )
        self.mSeparators[self.mRadius] = 1
        
    #####################################################################
    def getPosition( self, chr, strand, residue ):
        """return position on arc where a residue lies."""
        if self.mForwardCoordinates:
            return self.contig2Position[chr] + residue
        else:
            if strand in ("-", 0, "0"):
                return self.contig2Position[chr] + self.mMapContig2Size[chr] - residue
            else:
                return self.contig2Position[chr] + residue            

    #####################################################################
    def getAngle( self, pos ):
        """return angle on arc for position pos."""
        
        return 360.0 * pos / self.mTotalSize

    #####################################################################
    def getPosOnArc( self, angle, radius):
        """return position on arc centered on x, y at angle and radius."""
        
        degPerRad = math.pi / 180.0
        
        dx = radius* math.cos( angle* degPerRad)
        dy = radius* math.sin( angle* degPerRad)
        
        return dx + self.mDataMiddleX, dy + self.mDataMiddleY

    #####################################################################
    def getPosOnArcForRadius( self, x, radius ):
        """get pos on arc for coordinate x."""
        return self.getPosOnArc( self.getAngle(x), radius)

    #####################################################################
    def addDuplication( self, entries, map_gene2pos, height,
                        url = None,
                        with_separator=True,
                        link_to_previous = False,
                        quality2symbol = {},
                        quality2mask = {}):
        """add a dot in row/col."""
        
        mi, ma = None, 0

        pos = bisect.bisect( self.mColourThresholds, height )
        master_colour = self.mColours[pos]
        
        chrs = {}
        points = []

        if not link_to_previous:
            self.mPreviousPoints = {}

        ########################################################
        ########################################################
        ########################################################            
        ## convert gene list to a set of points
        ########################################################            
        for species, transcript, gene, quality in entries:

            chr, strand, first_res, last_res = map_gene2pos[gene]
            chrs[chr] = 1
            pos1 = self.getPosition( chr, strand, first_res )
            pos2 = self.getPosition( chr, strand, last_res )
            
            a = min( pos1, pos2 )
            b = max( pos1, pos2 )

            if mi == None:
                mi = a
            else:
                mi = min(a, mi)
            ma = max(b, ma)
            
            points.append( (pos1, pos2, gene, quality, chr) )

        ########################################################
        ########################################################
        ########################################################            
        ## decide whether we need to increment the radius
        ########################################################                        
        cis = len(chrs) == 1

        old_radius = self.mRadius
        is_overlap = False
        if cis:
            if not self.mLastChr:
                self.mLastChr = chr

            if chr != self.mLastChr:
                self.mRadius = self.mRadiusFallBack
                self.mLastMax = ma
                self.mPreviousMax = ma
                self.mLastChr = chr
            else:
                if self.mPreviousMax + self.mMinDistance > mi:
                    ## overlap due to close proximitiy
                    self.mRadius += self.mRadiusIncrement
                    if with_separator: self.addSeparator()

                    ## true overlap
                    if self.mPreviousMax > mi:
                        is_overlap = True
                    
                elif self.mLastMax + self.mMinDistance > mi:
                    pass
                else:
                    self.mRadius = self.mRadiusFallBack
                    
                self.mLastMax = max(self.mLastMax, ma)
        else:
            if self.mLastMax > mi:
                self.mRadius += self.mRadiusIncrement
                if with_separator: self.addSeparator()

        self.mPreviousMin = mi
        self.mPreviousMax = ma
        self.mPreviousCis = cis
        
        self.mRadiusMax = max(self.mRadius, self.mRadiusMax)

        ########################################################
        ########################################################
        ########################################################            
        ## draw points 
        ########################################################                    
        link_colour = master_colour
        link_rad_width = self.mLinkRadStrokeWidth
        link_arc_width = self.mLinkArcStrokeWidth
        new_points = {}
        for pos1, pos2, gene, quality, chr in points:

            angle = self.getAngle( (pos1 + pos2) / 2 )
            
            x,y = self.getPosOnArc( angle, self.mRadius )

            try:
                symbol = quality2symbol[quality]
            except KeyError:
                symbol = "rect"

            if quality in quality2mask:
                colour = self.mLinkColourSymbolMasked
                link_colour = self.mLinkColourMasked
                link_rad_width = self.mLinkStrokeWidthMasked
                link_arc_width = self.mLinkStrokeWidthMasked
            else:
                colour = master_colour

            if gene in self.mPreviousPoints:
                continue
            
            new_points[gene] = (x, y, angle, quality, chr)
            
            if symbol == "circle":
                ee = SVGdraw.circle( x, y, self.mLinkSymbolSize,
                                     fill = "rgb(%i,%i,%i)" % colour,
                                     stroke="black",
                                     stroke_width = self.mLinkStrokeWidthSymbol )
            elif symbol == "rect":
                ee =  SVGdraw.rect( x-self.mLinkSymbolSize/2, y-self.mLinkSymbolSize/2,
                                    self.mLinkSymbolSize, self.mLinkSymbolSize,
                                    fill = "rgb(%i,%i,%i)" % colour,
                                    stroke="black",
                                    stroke_width = self.mLinkStrokeWidthSymbol )

            if url:
                e = SVGdraw.link( url % gene )
                e.addElement( ee )
            else:
                e = ee
                
            self.addWheelElement( e )                
            
        ########################################################
        ########################################################
        ########################################################
        ## write all arcs in between old points and new points
        ## cis:   circular arc
        ## trans: radial arc
        ########################################################   

        angles = []
        
        for x1,y1,angle1,quality1,chr1 in new_points.values():

            ## reduce clutter by not writing arc to the same angle
            for x2,y2,angle2,quality2,chr2 in self.mPreviousPoints.values():

                for a in angles:
                    if a - self.mAngleResolution < angle2 < a + self.mAngleResolution:
                        break
                else:
                    angles.append( angle2 )

                    d = SVGdraw.pathdata( x1, y1 )

                    if chr1 == chr2:
                        d.relellarc( self.mRadius, self.mRadius, 0, 0, 1, x2-x1, y2-y1 )
                        link_width = link_rad_width
                    else:
                        d.relellarc( self.mRadius * 2, self.mRadius * 2, 0, 0, 0, x2-x1, y2-y1 )
                        link_width = link_arc_width
                        
                    e = SVGdraw.path( d,
                                      fill = "none",
                                      stroke = "rgb(%i,%i,%i)" % link_colour,
                                      stroke_width = link_width )

                    self.addWheelElement(e, self.mPlaneLinks)


        ## plot lines between new points
        new_genes = new_points.keys()

        for g1 in range(len(new_genes)-1):
            
            x1,y1,angle1,quality1,chr1 = new_points[new_genes[g1]]
            
            for g2 in range(g1+1, len(new_genes)):
                
                x2,y2,angle2,quality2,chr2 = new_points[new_genes[g2]]

                for a in angles:
                    if a - self.mAngleResolution < angle2 < a + self.mAngleResolution:
                        break
                else:
                    angles.append( angle2 )

                    d = SVGdraw.pathdata( x1, y1 )

                    if chr1 == chr2:
                        d.relellarc( self.mRadius, self.mRadius, 0, 0, 1, x2-x1, y2-y1 )
                        link_width = link_rad_width                        
                    else:
                        d.relellarc( self.mRadius * 2, self.mRadius * 2, 0, 0, 0, x2-x1, y2-y1 )
                        link_width = link_arc_width                        

                    e = SVGdraw.path( d,
                                      fill = "none",
                                      stroke = "rgb(%i,%i,%i)" % link_colour,
                                      stroke_width = link_width )

                    self.addWheelElement(e, self.mPlaneLinks)
                
        ## add new points to old points
        for k, v in new_points.items():
            self.mPreviousPoints[k] = v

    #####################################################################
    def pushRadius( self ):
        """push radius as fallback."""
        
        self.mRadiusFallBack = self.mRadiusMax + self.mRadiusIncrement
        self.mRadius = self.mRadiusFallBack
        # self.addSeparator()
        self.mLastMax = 0
        self.mPreviousMax = 0
        
    #####################################################################
    def writeGrid( self ):
        """add grid lines."""

        middlex = self.mDataMiddleX
        middley = self.mDataMiddleY

        ## print separators
        for c in range(len(self.contigs)):

            contig = self.contigs[c]
            
            pos = self.getPosition( contig, "+", 0)
            angle = self.getAngle( pos )

            x,y = self.getPosOnArc( angle, self.mRadius )
            
            e = SVGdraw.line( middlex,
                              middley,
                              x, 
                              y,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              stroke_width = self.mGridStrokeWidth )

            self.addElement(e, self.mPlaneGrid )

            if c < len(self.contigs) - 1:
                next_angle = self.getAngle(self.getPosition( self.contigs[c+1], "+", 0 ) )
            else:
                next_angle = 360
            
            x,y = self.getPosOnArc( angle + float(next_angle - angle) / 2, self.mRadiusStart / 2)
            
            e = SVGdraw.text( x,
                              y,
                              contig,
                              self.mGridFontSize,
                              self.mGridFont,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              text_anchor = "start" )
            # do not rotate text:
            # transform="rotate(%i,%i,%i)" % ( angle, x, y ))

            self.addElement(e)

    #####################################################################
    def writeScale( self ):
        """write scales."""

        current_x = self.mScaleX
        current_y = self.mScaleY + self.mScaleHeight

        nboxes = len(self.mColourThresholds)
        # box size for legend in x-direction
        # subtract size of right-most axis label so that it takes the
        # same width as self.mDataWidth.
        box_size_x = math.ceil( (self.mDataWidth -
                                (self.mScaleFontSize * len(self.mFormatNumberLegend % self.mColourThresholds[-1])))
                                / nboxes )

        # change font size such that it labels will fit between tick-marks
        self.mScaleFontSize = min( self.mScaleFontSize,
                                   (box_size_x * self.mScaleNumTicks * 1.5) / len(self.mFormatNumberLegend % self.mColourThresholds[-1]) )
        
        for x in range(nboxes):

            e = SVGdraw.rect( current_x,
                              current_y,
                              box_size_x,
                              self.mScaleBoxSizeY,                                   
                              fill = "rgb(%i,%i,%i)" % self.mColours[x],
                              stroke = "rgb(%i,%i,%i)" % self.mColours[x])                                   
            
            self.addElement(e)
            
            if x % self.mScaleNumTicks == 0:

                e = SVGdraw.line( current_x,
                                  current_y,
                                  current_x,
                                  current_y + self.mScaleBoxSizeY,
                                  stroke = "rgb(%i,%i,%i)" % BLACK,
                                  stroke_width = 5)
                self.addElement(e)
                
                e = SVGdraw.text( current_x,
                                  current_y - self.mScaleBoxSizeY,
                                  self.mFormatNumberLegend % self.mColourThresholds[x],
                                  self.mScaleFontSize,
                                  self.mScaleFont,
                                  stroke = "rgb(%i,%i,%i)" % BLACK,
                                  text_anchor = "start")
                self.addElement(e)
            
            
            current_x += box_size_x

        
    #####################################################################
    def writeFooter( self ):
        """write footer."""
        
        current_x = self.mFooterX
        current_y = self.mFooterY

        if self.mFooter:

            current_y += max( self.mFooterFontSize, self.mMaxBoxSize) + self.mSeparator
            
            e = SVGdraw.text( self.mPageWidth / 2,
                              current_y + self.mFooterFontSize,
                              self.mFooter,
                              self.mFooterFontSize,
                              self.mFooterFont,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              text_anchor = "middle")
            
            self.addElement(e)

    #####################################################################
    def writeWheel( self ):
        """write duplications.

        shifts everything according to final canvas size.
        There must a way to do this in SVG using relative
        coordinates.
        """
        
        def shiftElement( element, shiftx, shifty):
            if 'x' in element.attributes:
                element.attributes['x'] += shiftx
            if 'y' in element.attributes:
                element.attributes['y'] += shifty
            if 'cx' in element.attributes:
                element.attributes['cx'] += shiftx
            if 'cy' in element.attributes:
                element.attributes['cy'] += shifty
            if 'x1' in element.attributes:
                element.attributes['x1'] += shiftx
            if 'y1' in element.attributes:
                element.attributes['y1'] += shiftx
            if 'd' in element.attributes:
                ## assumes that path starts with a move
                path = element.attributes['d'].split(" ")
                if path[0] != "M":
                    raise "don't know what to do with paths not starting with 'M': %s" % (" ".join(path))
                path[1] = str(float(path[1]) + shiftx)
                path[2] = str(float(path[2]) + shifty)
                element.attributes['d'] = " ".join(path)
                
            for e in element.elements:
                shiftElement( e, shiftx, shifty )

        shiftx = self.mDataX + self.mDataWidth / 2
        shifty = self.mDataY + self.mDataHeight / 2
        
        for element, plane in self.mWheelElements:

            shiftElement( element, shiftx, shifty )
            self.addElement( element, plane )
        
    #####################################################################        
    def finalizePlot( self ):
        """write remaining parts of the plot."""

        ## add last circle
        self.mRadius = self.mRadiusMax + self.mRadiusIncrement
        self.addSeparator()

        ## From user input, build canvas size
        ## height width of various sections
        self.mHeaderWidth, self.mHeaderHeight = 0, 0
        self.mFooterWidth, self.mFooterHeight = 0, 0
        self.mDataWidth, self.mDataHeight = 0, 0

        if self.mTitle:
            self.mHeaderHeight += self.mSeparator + self.mTitleFontSize
        
        ## corrected height and width of data section
        self.mDataX = self.mHeaderWidth
        self.mDataY = self.mHeaderHeight
        self.mDataWidth  = 2 * self.mRadiusMax + 2 * self.mSeparator
        self.mDataHeight = self.mDataWidth
        self.mDataMiddleX = self.mDataX + self.mDataWidth / 2
        self.mDataMiddleY = self.mDataY + self.mDataHeight / 2

        ## scale
        self.mScaleX = self.mHeaderWidth + self.mSeparator
        self.mScaleY = self.mHeaderHeight + self.mSeparator + self.mDataHeight + self.mSeparator
        self.mScaleHeight =  2 * max( self.mScaleFontSize, self.mMaxBoxSize) + self.mSeparator
        self.mScaleWidth = 0

        ## height of footer 
        if self.mFooter:
            self.mFooterHeight = 2 * self.mFooterFontSize
        else:
            self.mFooterHeight = 0

        self.mFooterX = self.mHeaderWidth + self.mDataWidth + self.mScaleWidth + self.mFooterWidth + 3 * self.mSeparator
        self.mFooterY = self.mHeaderHeight + self.mDataHeight + self.mScaleHeight + self.mFooterHeight + 3 * self.mSeparator 

        ## compute actual page width and height
        self.mPageWidth = self.mHeaderWidth + self.mDataWidth + self.mScaleWidth + self.mFooterWidth + 3 * self.mSeparator
        self.mPageHeight = self.mHeaderHeight + self.mDataHeight + self.mScaleHeight + self.mFooterHeight + 3 * self.mSeparator        

        self.writeWheel()
        self.writeGrid()
        self.writeTitle()
        self.writeScale()
        self.writeFooter()

    #####################################################################        
    def writeToFile(self, outfile):
        """write svg image to file.
        """
        self.finalizePlot()

        self.mRoot = SVGdraw.drawing()
        self.mDraw = SVGdraw.svg( (0, 0, self.mPageWidth, self.mPageHeight ) , "100%", "100%" )

        kk = self.mElements.keys()
        kk.sort()
        kk.reverse()
        for k in kk:
            for e in self.mElements[k]:
                self.mDraw.addElement( e )
            
        self.mRoot.setSVG(self.mDraw)

        tfile = tempfile.mktemp()
        
        self.mRoot.toXml( tfile )

        lines = open(tfile,"r").readlines()
        
        outfile.write(string.join(lines,""))
        outfile.write("\n")
        
        os.remove(tfile)

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: SVGDuplicationsWheel.py 2784 2009-09-10 11:41:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-e", "--headers", dest="headers", action="store_true",
                      help="first row is a header [ignored]."  )
    parser.add_option("-t", "--title", dest="title", type="string",
                      help="page title.")
    parser.add_option("-f", "--footer", dest="footer", type="string",
                      help="page footer.")
    parser.add_option("-c", "--contig-sizes", dest="filename_contig_sizes", type="string",
                      help="filname with contig sizes." )
    parser.add_option("-r", "--radius", dest="radius", type="int",
                      help="radius.")
    parser.add_option("-i", "--increment", dest="radius_increment", type="int",
                      help="radius increment.")
    parser.add_option("-u", "--url", dest="url", type="string",
                      help="string to build url for annotation.")
    parser.add_option( "--min-contig", dest="min_contig_size", type="string",
                      help="minimum contig size to delineate.")
    parser.add_option( "--min-value", dest="min_value", type="float",
                       help="minimum branch length.")
    
    parser.add_option( "--max-value", dest="max_value", type="float",
                       help="maximum branch length.")
    
    parser.set_defaults(
        filename_contig_sizes = None,
        headers = False,
        titles = "",
        pattern_filename = None,
        title = "",
        footer = "",
        radius = 3000,
        min_value = 0.0,
        max_value = 0.2,
        url=None,
        radius_increment = 40,
        min_contig_size = 10000,
        remove_empty_contigs = True,
        separator = "|",
        quality2symbol = { 'CG' : "circle", 'PG' : "circle", 'SG' : "circle" },
        quality2mask = ("RG", "CP", "PP", "SP", "RP", "CF", "PF", "SF", "UG", "UP", "UF", "BF", "UK" ),
        sort_by_size = True,
        input_format = "pairwise",
        )

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    if options.filename_contig_sizes:
        map_contig2size = IOTools.ReadMap( open(options.filename_contig_sizes, "r"),
                                           map_functions=(str, int) )

    ## read data and get contigs that are used (i.e.: remove empty contigs)
    chrs = {}
    lines = sys.stdin.readlines()

    if options.remove_empty_contigs:
        for line in lines:
            if line[0] == "#": continue

            d = line[:-1].split( "\t" )

            cluster_id, in_locations, in_tree = d[:3]

            for l in in_locations.split(";"):
                gene_id, chr, strand, sbjct_from, sbjct_to = l.split(":")
                if chr not in map_contig2size: continue
                chrs[chr] = 1
        for k in map_contig2size.keys():
            if k not in chrs:
                del map_contig2size[k]

    k = map_contig2size.keys()
    
    if len(k) == 0:
        Experiment.Stop()
        sys.exit(0)
    
    k.sort()
    
    if options.sort_by_size:
        k.sort( lambda x,y: cmp(map_contig2size[x], map_contig2size[y]))
        
    plot = DuplicationPlot( k,
                            map_contig2size,
                            num_entries = 0 )

    plot.mRadiusIncrement = options.radius_increment
    plot.mRadius = options.radius
    plot.mMaxValue = options.max_value
    plot.mMinValue = options.min_value
    
    if options.title: plot.setTitle(options.title)
    if options.footer: plot.setFooter(options.footer)    

    plot.initializePlot()
    
    data = []

    if options.input_format == "pairwise":
        
        ## read data from pairwise analysis
        ## format is: cluster_id, locations of duplications, tree of duplications
    
        for line in lines:
            if line[0] == "#": continue

            d = line[:-1].split( "\t" )

            cluster_id, in_locations, in_tree = d[:3]

            mi, ma = 0,0
            found = False
            n = 0
            chrs = {}
            for l in in_locations.split(";"):
                gene_id, chr, strand, sbjct_from, sbjct_to = l.split(":")
                if chr not in map_contig2size: continue
                chrs[chr] = 1
                sbjct_from, sbjct_to = int(sbjct_from), int(sbjct_to)

                xi = plot.getPosition( chr, strand, sbjct_from)
                xa = plot.getPosition( chr, strand, sbjct_to)

                if not mi:
                    mi = xi
                else:
                    mi = min(mi, xi)

                n += 1
                ma = max(ma, xa )
                found = True

            if not found: continue
            cis = len(chrs) == 1
            data.append( (cis, n, mi, ma, cluster_id, in_locations, in_tree) )
        
    data.sort()

    plot.mNumEntries = len(data)
    plot.initializePlot()
    
    last_ndups = 0

    for cis, ndups, mi, ma, cluster_id, in_locations, in_tree in data[:]:

        if ndups != last_ndups:
            plot.pushRadius()
            plot.addSeparator()
            
        last_ndups = ndups
        
        map_gene2location = {}
        for l in in_locations.split(";"):
            gene_id, chr, strand, sbjct_from, sbjct_to = l.split(":")
            if chr not in map_contig2size:
                continue
            sbjct_from, sbjct_to = int(sbjct_from), int(sbjct_to)
            map_gene2location[gene_id] = (chr, strand, sbjct_from, sbjct_to)

        if not map_gene2location: continue
        
        tree = TreeTools.Newick2Tree( in_tree )

        ## the last subset is all nodes again.
        s = TreeTools.GetSubsets( tree )
        is_first = True
        for children, height, branchlength in s[:-1]:
            if len(children)==1: continue
            data = []
            for c in map(lambda x: x.split(options.separator), children):
                if len(c) == 2:
                    data.append( (c[0], c[1], c[1], "CG" ) )
                elif len(c) == 1:
                    data.append( ("unk", c[0], c[0], "CG" ) )
                elif len(c) == 3:
                    data.append( (c[0], c[2], c[3], "CG" ) )

            for species, transcript, gene, quality in data:
                if not gene in  map_gene2location:
                    print >>sys.stderr, in_tree + "\n" +  in_locations + "\n", map_gene2location
                

            plot.addDuplication( data, map_gene2location, height,
                                 url=options.url,
                                 with_separator=is_first,
                                 link_to_previous=not is_first,
                                 quality2symbol = options.quality2symbol,
                                 quality2mask = options.quality2mask)
            is_first = False

    plot.writeToFile(sys.stdout)
    
    Experiment.Stop()

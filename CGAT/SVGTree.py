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
SVGTree.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, getopt, time, optparse, math, tempfile, bisect

""" program $Id: SVGTree.py 2784 2009-09-10 11:41:14Z andreas $
class to plot a phylogenetic tree as an SVG file.
"""
import Experiment
import IOTools
import SVGdraw 
import TreeTools
import Tree

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
ORANGE = (255,165,0)
PINK = (255,192,203)
PLUM= (221,160,221)
OLIVE = (128,128,0)
BROWN = (165,42,42)

COLOURS= (BLACK, RED, GREEN, BLUE, YELLOW, CYAN, PURPLE, GREY, ORANGE, PINK,
          PLUM, OLIVE, BROWN)

MAX_GREY = 240
GREY_COLORS = map( lambda x: (x,x,x), range(0,MAX_GREY))

###################################################################
class BorderDecorator:
    """class for decorating the border in a tree plot."""

    mFontSize = 36
    mFont = "Verdana"
    
    def __init__(self, tree):
        self.mTree = tree

class BorderDecoratorLocations(BorderDecorator):

    mTickWidth = 10
    mDefaultColour = BLACK
    
    def __init__(self, tree, map_id2location,
                 extract_species = None, map_species2colour = None,
                 max_separation = 0,
                 print_location = False,
                 *args, **kwargs):

        """
        provide extract_species and map_species2colour to colour
        by species.
        
        max_separation: maximum separation between syntenic locations.
        set 0, if any distance is ok and only synteny is enough.
        """
        BorderDecorator.__init__(self, tree, *args, **kwargs)

        self.mMaxSeparation = max_separation        
        self.mMapId2Location = map_id2location
        self.mExtractSpecies = extract_species
        self.mMapSpecies2Colour = {}

        if map_species2colour:
            for s, c in map_species2colour.items():
                self.mMapSpecies2Colour[s] = tuple(map(int, c.split(",")))
        
        self.mPrintLocation = print_location
        
    def getWidth( self ):

        t = self.mTree.get_terminals()

        ## check if we have terminals which need to be connected
        count = 0
        for i in range(len(t)-1):
            node_id1 = t[i]
            taxon1 = self.mTree.node(node_id1).data.taxon

            for j in range(i+1, len(t)):
                node_id2 = t[j]

                taxon2 = self.mTree.node(node_id2).data.taxon                
                
                l1 = self.mMapId2Location[taxon1]
                l2 = self.mMapId2Location[taxon2]                
                if l1.contig != l2.contig:
                    continue

                if self.mMaxSeparation:
                    s = min( abs(l1.mFrom - l2.mTo), abs(l1.mTo - l2.mFrom))
                    if s >= self.mMaxSeparation: continue

                count += 1

        if count == 0:
            return 0
        
        ## rule of thumb
        if self.mPrintLocation:
            return 400
        else:
            return 100
    
    def getElements(self, x, y, map_node2height ):

        t = self.mTree.get_terminals()
        
        elements = []

        ## print locations
        if self.mPrintLocation:
            for i in range(len(t)):
                node_id1 = t[i]
                taxon1 = self.mTree.node(node_id1).data.taxon
                y1 = map_node2height[node_id1] + y

                elements.append( SVGdraw.text( x, y1,
                                               str(self.mMapId2Location[taxon1]),
                                               self.mFontSize,
                                               self.mFont,
                                               stroke = "rgb(%i,%i,%i)" % BLACK,
                                               text_anchor = "left" ))

        
        ## print connectors
        for i in range(len(t)-1):
            node_id1 = t[i]
            taxon1 = self.mTree.node(node_id1).data.taxon
            y1 = map_node2height[node_id1] + y

            for j in range(i+1, len(t)):
                node_id2 = t[j]

                taxon2 = self.mTree.node(node_id2).data.taxon                

                if self.mExtractSpecies:
                    species1 = self.mExtractSpecies(taxon1)
                    species2 = self.mExtractSpecies(taxon2)

                    if species1 != species2: continue

                    if species1 not in self.mMapSpecies2Colour:
                        self.mMapSpecies2Colour[species1] = COLOURS[len(self.mMapSpecies2Colour) % len(COLOURS) ]

                    colour = self.mMapSpecies2Colour[species1]
                    
                else:
                    colour = self.mDefaultColour
                    
                l1 = self.mMapId2Location[taxon1]
                l2 = self.mMapId2Location[taxon2]                
                if l1.contig != l2.contig:
                    continue

                if self.mMaxSeparation:
                    s = min( abs(l1.mFrom - l2.mTo), abs(l1.mTo - l2.mFrom))
                    if s >= self.mMaxSeparation: continue
                    
                y2 = map_node2height[node_id2] + y

                distance = y2 - y1

                d = SVGdraw.pathdata( x, y1 )

                d.line( x + self.mTickWidth, y1 )
                d.ellarc( distance, distance, 0, 0, 1, x + self.mTickWidth, y2 )
                d.line( x, y2 )                

                
                e = SVGdraw.path( d,
                                  fill = "none",
                                  stroke = "rgb(%i,%i,%i)" % colour,
                                  stroke_width = 1 )

                elements.append( e )

        return elements

class BorderDecoratorClusters(BorderDecorator):

    mBoxWidth = 10
    mDefaultColour = BLACK
    
    def __init__(self, tree, map_id2cluster,
                 map_cluster2colour = None,
                 *args, **kwargs):

        """
        plot a box according to cluster.
        
        """
        BorderDecorator.__init__(self, tree, *args, **kwargs)

        self.mMapId2Cluster = map_id2cluster
        self.mMapCluster2Colour = {}

        if map_cluster2colour:
            for s, c in map_cluster2colour.items():
                self.mMapCluster2Colour[s] = tuple(map(int, c.split(",")))
        
    def getWidth( self ):
        return self.mBoxWidth
    
    def getElements(self, x, y, map_node2height ):

        elements = []
        for n in self.mTree.get_terminals():

            t = self.mTree.node(n).data.taxon

            if t not in self.mMapId2Cluster:
                continue
            cluster = self.mMapId2Cluster[t]
            if cluster not in self.mMapCluster2Colour:
                self.mMapCluster2Colour[cluster] = COLOURS[len(self.mMapCluster2Colour) % len(COLOURS) ]

            colour = self.mMapCluster2Colour[cluster]
            elements.append( SVGdraw.rect( x, y + map_node2height[n] - self.mBoxWidth / 2, self.mBoxWidth, self.mBoxWidth,
                                           stroke = "rgb(%i,%i,%i)" % colour,
                                           fill = "rgb(%i,%i,%i)" % colour) )
                
        return elements

###################################################################
class BranchDecorator:
    """class for decorating the border in a tree plot."""
    def __init__(self, tree):
        self.mTree = tree
        self.mFontSize = 36
        self.mFont = "Verdana"
        self.mFontColour = BLACK
        
    def getHeight( self, node_id ):
        return 1

    def getWidth( self, node_id ):
        return 1

    def setFontSize( self, size ):
        """set font size."""
        self.mFontSize = size

    def setFontColour( self, colour ):
        self.mFontColour = colour
        
class BranchDecoratorVertical( BranchDecorator ):

    
    def __init__(self, *args, **kwargs):
        BranchDecorator.__init__( self, *args, **kwargs)

    def getElements( self, node_id, x, y1, y2 ):

        e = []
        e.append(SVGdraw.line( x, y1, x, y2,
                                stroke = "rgb(%i,%i,%i)" % BLACK,
                                stroke_width = 1 ))
        return e

        
class BranchDecoratorHorizontal( BranchDecorator ):
    """default branch decorator - default is plane line."""
    
    def __init__(self, *args, **kwargs):
        BranchDecorator.__init__( self, *args, **kwargs)
        self.mNoBranch = False
        
    def setNoBranch( self ):
        self.mNoBranch = True
        
    def getElements( self, node_id, x1, x2, y ):

        e = []
        if not self.mNoBranch:
            e.append(SVGdraw.line( x1, y, x2, y,
                                   stroke = "rgb(%i,%i,%i)" % BLACK,
                                   stroke_width = 1 ) )
        return e

class BranchDecoratorHorizontalBranchLength( BranchDecoratorHorizontal ):
    """branch decorator - adds branch length to branch."""

    mBranchLengthFormat = "%5.2f"
    
    def __init__(self, *args, **kwargs):
        BranchDecoratorHorizontal.__init__( self, *args, **kwargs)
        self.mWritten = 0
        
    def getElements( self, node_id, x1, x2, y ):

        e = BranchDecoratorHorizontal.getElements( self, node_id, x1, x2, y )

        text = self.getText( node_id )

        d = x2 - x1
        df = len(text) * (self.mFontSize / 2)
        if df > d:
            x = x1 - df
        else:
            x = x1 + (x2 - x1 - df) / 2

        e.append( SVGdraw.text( x, y - 5,
                                text,
                                self.mFontSize,
                                self.mFont,
                                stroke = "rgb(%i,%i,%i)" % self.mFontColour,
                                text_anchor = "left" ))

        return e

    def getText( self, node_id ):
        """return text."""        
        return self.mBranchLengthFormat % self.mTree.node(node_id).data.branchlength

    def getHeight( self, node_id ):
        return 5 + self.mFontSize
        
class BranchDecoratorHorizontalBranchLengthError( BranchDecoratorHorizontalBranchLength ):
    """branch decorator - adds branch length to branch + error."""

    mBranchErrorFormat = "%5.2f"

    def __init__(self, tree, error_tree, *args, **kwargs):
        BranchDecoratorHorizontalBranchLength.__init__( self, tree, *args, **kwargs)
        self.mErrorTree = error_tree
        
    def getText( self, node_id ):
        """return text."""
        return BranchDecoratorHorizontalBranchLength.getText(self, node_id ) +\
               "+/-" + self.mBranchErrorFormat % self.mErrorTree.node(node_id).data.branchlength        

class BranchDecoratorHorizontalBranchLengthWithKaks( BranchDecoratorHorizontalBranchLength ):
    """branch decorator - adds branch length to branch + kaks."""

    mBranchKaksFormat = "%.2f"

    def __init__(self, tree, error_tree, *args, **kwargs):
        BranchDecoratorHorizontalBranchLength.__init__( self, tree, *args, **kwargs)
        self.mKaksTree = error_tree
        
    def getText( self, node_id ):
        """return text."""
        d = self.mKaksTree.node(node_id).data.branchlength
        if d:
            text = " (%s)" % (self.mBranchKaksFormat % self.mKaksTree.node(node_id).data.branchlength)
        else:
            text = ""

        return BranchDecoratorHorizontalBranchLength.getText(self, node_id ) + text

class BranchDecoratorHorizontalAboveBelow( BranchDecoratorHorizontalBranchLength ):
    """branch decorator - displays two branch values.

    One decorator is displayed above and the the other one below the branch.
    The branch itself is rendered as a horizontal line.
    """

    def __init__(self, tree, decorator1, decorator2, *args, **kwargs):

        BranchDecoratorHorizontalBranchLength.__init__( self, tree, *args, **kwargs)
        self.mDecorator1 = decorator1
        self.mDecorator2 = decorator2

        ## do not write branch
        self.mDecorator1.setNoBranch()
        self.mDecorator2.setNoBranch()        
        
    def getElements( self, node_id, x1, x2, y ):

        e = []
        e.append(SVGdraw.line( x1, y, x2, y,
                               stroke = "rgb(%i,%i,%i)" % BLACK,
                               stroke_width = 1 ) )
        
        e += self.mDecorator1.getElements( node_id, x1, x2, y )
        e += self.mDecorator2.getElements( node_id, x1, x2, y + self.mFontSize )        
        
        return e

###################################################################
class NodeDecorator:
    """class for decorating nodes.

    The default decorator does nothing. It implements the
    methods to plot symbols. Symbol choices are

    * square
    * circle
    
    """

    mFontSize  = 36
    mFont      = "Verdana"
    mFontStyle = None
    
    def __init__(self, tree,
                 plot_label = False,
                 plot_symbol = None):
        
        self.mTree = tree
        self.mPlotLabel = plot_label
        self.mPlotSymbol = plot_symbol

    def getColour( self, node_id, x, y ):
        return BLACK
    
    def getElements( self, node_id, x, y ):

        e = []
        colour = self.getColour( node_id, x, y )
        
        if self.mPlotSymbol == "circle":
            e.append( SVGdraw.circle( x + self.mFontSize / 2,
                                      y,
                                      self.mFontSize/2,
                                      stroke = "rgb(%i,%i,%i)" % BLACK,
                                      fill = "rgb(%i,%i,%i)" % colour) )

        elif self.mPlotSymbol == "square":
            e.append( SVGdraw.rect( x, y-self.mFontSize/2,
                                    self.mFontSize,
                                    self.mFontSize,
                                    stroke = "rgb(%i,%i,%i)" % BLACK,
                                    fill = "rgb(%i,%i,%i)" % colour) )
            
        return e

    def getHeight( self, node_id ):
        if self.mPlotSymbol != None:
            return self.mFontSize
        else:
            return 0
        
    def getWidth( self, node_id ):
        if self.mPlotSymbol != None:
            return self.mFontSize
        else:
            return 0

    def setPlotLabel( self, plot_label ):
        self.mPlotLabel = plot_label

    def getLegend( self, x, y ):
        """add legend.
        returns elements and maximum x and y output.
        """
        return [], x, y

    def getFooterHeight( self ):
        return 0

    def setFontStyle( self, style ):
        self.mFontStyle = style
    
###################################################################        
class NodeDecoratorExternal(NodeDecorator):
    """class for decorating external nodes. Simply write name.
    """
    
    def __init__(self, *args, **kwargs):
        NodeDecorator.__init__( self, *args, **kwargs)
    
    def getElements( self, node_id, x, y, x_label = None, y_label = None ):
        
        e = NodeDecorator.getElements( self, node_id, x, y )

        if x_label == None: x_label = x
        if y_label == None: y_label = y        
        e.append( SVGdraw.text( x_label, y_label,
                                self.mTree.node(node_id).data.taxon,
                                self.mFontSize,
                                self.mFont,
                                stroke = "rgb(%i,%i,%i)" % BLACK,
                                font_style = self.mFontStyle,
                                text_anchor = "left" ))
        return e

    def getHeight( self, node_id ):
        return max( NodeDecorator.getHeight(self, node_id), self.mFontSize)

    def getWidth( self, node_id ):
        m = NodeDecorator.getWidth( self, node_id )
        if self.mPlotLabel:
            l = len(self.mTree.node(node_id).data.taxon) * self.mFontSize
        else:
            l = 0
        return max(l, m)
    
###################################################################
class NodeDecoratorBySpecies(NodeDecorator):
    """class for decorating external nodes.

    Colour by species and rename.
    """

    def __init__(self, tree,
                 plot_label = None,
                 plot_symbol = None,
                 extract_species = str,
                 map_species2colour = None,
                 map_species2name = None,
                 map_taxon2url = None,
                 *args, **kwargs):
        
        """colours are a ','-separated tuples of RGB values.

        """
        ## plot labels unless explicitely set to false:
        if plot_label == None: plot_label = True
        
        NodeDecorator.__init__( self, tree,
                                plot_label = plot_label,
                                plot_symbol = plot_symbol,
                                *args, **kwargs )
        
        self.mExtractSpecies = extract_species
        self.mMapSpecies2Colour = {}

        if map_species2colour:
            for s, c in map_species2colour.items():
                self.mMapSpecies2Colour[s] = tuple(map(int, c.split(",")))

        self.mMapTaxon2URL = map_taxon2url
        if map_species2name:
            self.mMapSpecies2Name = map_species2name
        else:
            self.mMapSpecies2Name = {}

    def getColour( self, node_id, x, y ):

       species = self.mExtractSpecies(self.mTree.node(node_id).data.taxon)
        
       if species not in self.mMapSpecies2Colour:
            self.mMapSpecies2Colour[species] = COLOURS[len(self.mMapSpecies2Colour) % len(COLOURS) ]

       return self.mMapSpecies2Colour[species]
            
    def getElements( self, node_id, x, y, x_label = None, y_label = None):

        e = NodeDecorator.getElements( self, node_id, x, y )

        if x_label == None: x_label = x
        if y_label == None: y_label = y        
        
        t = self.mTree.node(node_id).data.taxon        
        species = self.mExtractSpecies(t)
        
        if species not in self.mMapSpecies2Colour:
            self.mMapSpecies2Colour[species] = COLOURS[len(self.mMapSpecies2Colour) % len(COLOURS) ]

        if species in self.mMapSpecies2Name:
            tx = re.sub( species, "%s" % self.mMapSpecies2Name[species], t)
        else:
            tx = t

        colour = self.getColour( node_id, x, y )
        
        if self.mPlotLabel:
            ee = SVGdraw.text( x_label, y_label,
                               tx,
                               self.mFontSize,
                               self.mFont,
                               stroke = "rgb(%i,%i,%i)" % colour,
                               text_anchor = "left" )

            if self.mMapTaxon2URL != None:
                url = self.mMapTaxon2URL(t)
                if url:
                    l = SVGdraw.link( url )
                    l.addElement( ee )
                    e.append( l )
                else:
                    e.append(ee)
            else:
                e.append( ee )
        return e

    def getHeight( self, node_id ):
        return max( NodeDecorator.getHeight(self, node_id), self.mFontSize)

    def getWidth( self, node_id ):
        m = NodeDecorator.getWidth( self, node_id)
        
        if self.mPlotLabel:
            t = self.mTree.node(node_id).data.taxon
            if t in self.mMapSpecies2Name:
                t = self.mMapSpecies2Name[t]
            l = len(t) * self.mFontSize * 2 / 3
        else:
            l = 0
            
        return max(m, l)
        

###################################################################
    
class SVGTree:

    def __init__(self, tree ):

        self.mElements = {}

        self.mTree = tree
        
        ## footer
        self.mFooterFrom = 10        
        self.mFooterFontSize  = 36
        self.mFooterFont      = "Verdana"
        self.mFooter = None
        
        ## Title
        self.mTitleFontSize  = 36
        self.mTitleFont      = "Verdana"
        self.mTitle = None

        ## Row/column labels
        self.mLabelFontSize  = 36
        self.mLabelFont      = "Verdana"

        ## default data area size (without terminal node labels)
        self.mDefaultHeight = 1000
        self.mDefaultWidth = 1000

        ## Row/column labels
        self.mTerminalLabelSeparator = 5
        self.mTerminalLabelOffset    = 5

        ## factor by which to scale branch lengths
        self.mBranchScaleFactor = 0
        ## factor by which to scale heights
        self.mHeightScaleFactor = 0
        
        self.mDecoratorExternalNodes = NodeDecoratorExternal(self.mTree)
        self.mDecoratorInternalNodes = NodeDecorator(self.mTree)

        self.mSeparatorWidth  = self.mLabelFontSize
        self.mSeparatorHeight = 2 * self.mLabelFontSize        

        ## increment of ruler in branch length
        self.mRulerIncrement = 0.1
        self.mRulerFormat = "%5.2f"
        self.mRulerTickSize = 5
        self.mRulerElements = ("ruler","left-ticks", "scale" )
        self.mRulerFontSize = 36
        self.mRulerFont      = "Verdana"

        ## whether or not terminal labels are left justified
        self.mLeftJustifiedExternalNodes = True

        ## decorators for the border
        self.mBorderDecorators= []

        self.mDecoratorVerticalBranches = BranchDecoratorVertical(self.mTree)
        self.mDecoratorHorizontalBranches = BranchDecoratorHorizontal(self.mTree)        


    #####################################################################
    def setRulerElements( self, ruler_elements = [] ):
        self.mRulerElements = ruler_elements
        
    #####################################################################
    def setBranchScale( self, scale ):
        self.mBranchScaleFactor = scale
        
    #####################################################################
    def setHeightScale( self, scale ):
        self.mHeightScaleFactor = scale

    #####################################################################
    def addBorderDecorator( self, decorator ):
        self.mBorderDecorators.append( decorator )
        
    #####################################################################
    def setDecoratorExternalNodes( self, decorator ):
        self.mDecoratorExternalNodes = decorator        

    #####################################################################
    def getDecoratorExternalNodes( self ):
        return self.mDecoratorExternalNodes

    #####################################################################
    def setDecoratorInternalNodes( self, decorator ):
        self.mDecoratorInternalNodes = decorator        

    #####################################################################
    def getDecoratorInternalNodes( self ):
        return self.mDecoratorInternalNodes
    
    #####################################################################
    def setDecoratorHorizontalBranches( self, decorator ):
        self.mDecoratorHorizontalBranches = decorator        

    #####################################################################
    def getDecoratorHorizontalBranches( self ):
        return self.mDecoratorHorizontalBranches

    #####################################################################
    def setDecoratorVerticalBranches( self, decorator ):
        self.mDecoratorVerticalBranches = decorator        

    #####################################################################
    def setDecoratorVerticalBranches( self ):
        return self.mDecoratorVerticalBranches

    ###################################################################
    def setFontSize( self, size ):
        """set font size in all elements."""
        self.mDecoratorInternalNodes.mFontSize = size
        self.mDecoratorExternalNodes.mFontSize = size
        self.mDecoratorVerticalBranches.mFontSize = size
        self.mDecoratorHorizontalBranches.mFontSize = size        
    
    #####################################################################
    def addElement( self, element, plane = 0):
        """add element to list in plane.
        
        The list is later sorted by priority.
        """
        if plane not in self.mElements:
            self.mElements[plane] = []
            
        self.mElements[plane].append( element )

    #####################################################################
    def addElements( self, elements, plane = 0):
        """add multiple elments to a plane."""
        for e in elements:
            self.addElement( e, plane )

    #####################################################################                        
    def initializePlot( self ):
        """set various coordinates in the plot.

        Note:

        Width  = X = coordinate 1
        Height = Y = coordinate 2
        
        """

        self.mNTaxa = len(self.mTree.get_taxa())
        self.mNNodes = max( self.mTree.chain.keys() ) + 1

        self.calculateCoordinates()
        
        self.calculateCanvasSize( )

    #####################################################################
    def getHeaderHeight( self ):
        """return height of header"""
        return 0

    #####################################################################
    def getHeaderWidth( self ):
        """return width of header"""
        return 0
    
    #####################################################################
    def calculateCanvasSize( self ):
        
        ## calculate size of data area
        self.mDataWidth  = max(self.mNodeWidthsEnd ) + max([ self.mDecoratorExternalNodes.getWidth( x ) for x in self.mTree.get_terminals() ])
        self.mDataHeight = self.mMaxNodeHeight

        ## calculate border width
        self.mBorderWidth = 0
        for x in self.mBorderDecorators:
            self.mBorderWidth += x.getWidth()
            
        ## set the page size
        self.mPageWidth  = self.getHeaderWidth() + self.mDataWidth + self.mBorderWidth
        self.mPageHeight = self.getHeaderHeight() + self.mDataHeight + self.getFooterHeight()

        ## patch: white background
##         self.addElement( SVGdraw.rect( 0, 0, self.mPageWidth, self.mPageHeight,
##                                        stroke = "rgb(%i,%i,%i)" % WHITE,
##                                        fill = "rgb(%i,%i,%i)" % WHITE) )

    #####################################################################
    def getFooterHeight( self ):

        footer_height = 0
        if "ruler" in self.mRulerElements:
            footer_height += 2 * self.mRulerTickSize + 1 + self.mSeparatorHeight
        if "scale" in self.mRulerElements:
            footer_height += 2 * self.mRulerTickSize + 1 + self.mSeparatorHeight + self.mRulerFontSize

        return footer_height + \
               self.mDecoratorExternalNodes.getFooterHeight() + \
               self.mDecoratorInternalNodes.getFooterHeight()

    #####################################################################
    def calculateCoordinates( self ):
        
        self.mNodeHeights = [0] * self.mNNodes
        self.mNodeWidthsStart = [0] * self.mNNodes
        self.mNodeWidthsEnd  = [0] * self.mNNodes        

        ## if no scales are given, try to do best fit
        if self.mHeightScaleFactor == 0:
            rescale_height = True
            self.mHeightScaleFactor = 1
        else:
            rescale_height = False

        if self.mBranchScaleFactor == 0:
            rescale_width = True
            self.mBranchScaleFactor = 100
        else:
            rescale_width = False

        ##########################################################
        ## Get Vertical coordinates
        ## Label nodes by their height. Terminal nodes have integer coordinates.
        ## Internal nodes have fractional coordinates (the average between the two
        ## children)
        counter = [0]        
        def updateHeights (node_id):
            l = len(self.mTree.node(node_id).succ)
            if l:
                ## set node height for internal node
                t = 0
                for x in self.mTree.node(node_id).succ:
                    t += self.mNodeHeights[x]
                    # used to use use the following to take into account
                    # the height of symbols. This is wrong and better done by
                    # pre-traversal of the tree
                    # self.mNodeHeights[node_id] = float(t) / float(l) + max( self.mDecoratorInternalNodes.getHeight( node_id ), self.mDecoratorHorizontalBranches.getHeight( node_id ))
                    # instead: use uncorrected heights.
                    self.mNodeHeights[node_id] = float(t) / float(l) 
            else:
                ## set node height for external node
                self.mNodeHeights[node_id] = counter[0]
                counter[0] += max( self.mDecoratorExternalNodes.getHeight( node_id ),
                                   self.mDecoratorHorizontalBranches.getHeight( node_id) ) \
                                   * self.mHeightScaleFactor + self.mTerminalLabelSeparator
                
        TreeTools.TreeDFS( self.mTree, self.mTree.root,
                           post_function = updateHeights )

        self.mMaxNodeHeight = counter[0]

        ##########################################################
        ## Get horizontal coordinates
        def updateWidths( node_id ):

            node = self.mTree.node( node_id )
            d = node.data.branchlength
            ## set default branchlength to 0.01 for empty branch lengths
            ## TODO: deal with trees without branch lengths later.
            if d <= 0.0: d = 0.01
            right = self.mNodeWidthsStart[node_id] + int(d * self.mBranchScaleFactor)
            self.mNodeWidthsEnd[node_id] = right
            for s in node.succ:
                self.mNodeWidthsStart[s] = right

        TreeTools.TreeDFS( self.mTree, self.mTree.root,
                           pre_function = updateWidths )

        if rescale_height:
            m = max( self.mNodeHeights )
            f = float(self.mDefaultHeight) / m
            if 100 * f < 1: f = 0.01
            self.mHeightScaleFactor = 100 * f
            self.mNodeHeights = map(lambda x: int(x * f), self.mNodeHeights)
            self.mMaxNodeHeight *= f
            
        if rescale_width:
            m = max( self.mNodeWidthsEnd )
            f = float(self.mDefaultWidth) / m
            self.mBranchScaleFactor = 100 * f            
            self.mNodeWidthsStart = map(lambda x: int(x * f), self.mNodeWidthsStart)
            self.mNodeWidthsEnd = map(lambda x: int(x * f), self.mNodeWidthsEnd)            

        ## add a safety margin for decorators writing above the line. This
        ## is a patch and should be changed such that decorators report
        ## their correct height

        for x in range(self.mNNodes):
            self.mNodeHeights[x] += 45

    #####################################################################
    def plotLines( self ):
        """plot lines of the tree"""
        
        ## plot tree in dfs manner
        def plotLines( node_id ):

            node = self.mTree.node( node_id )

            left = self.mNodeWidthsStart[node_id]
            right = self.mNodeWidthsEnd[node_id]
            height = self.mNodeHeights[node_id] 

            if right != left and node_id != self.mTree.root:
                self.addElements( self.mDecoratorHorizontalBranches.getElements(
                    node_id,
                    self.getHeaderWidth() + left,
                    self.getHeaderWidth() + right,
                    self.getHeaderHeight() + height ))
                                 

            for s in node.succ:

                new_height = self.mNodeHeights[s]
                self.addElements( self.mDecoratorVerticalBranches.getElements(
                    node_id,
                    self.getHeaderWidth() + right,
                    self.getHeaderHeight() + height,
                    self.getHeaderHeight() + new_height ))
                
        TreeTools.TreeDFS( self.mTree, self.mTree.root,
                           pre_function = plotLines )

        
    #####################################################################        
    def plotTree( self ):
        """plot the tree."""

        self.plotLines()
        
        self.plotExternalNodes()
        self.plotInternalNodes()        
        
    #####################################################################
    def plotExternalNodes( self ):
        """plot terminal node labels."""

        max_x = max(self.mNodeWidthsEnd)
        for node_id in self.mTree.get_terminals():

            node = self.mTree.node( node_id )

            x = self.mNodeWidthsEnd[node_id]
            y = self.mNodeHeights[node_id]
            
            if self.mLeftJustifiedExternalNodes:
                x_label = max_x
            else:
                x_label = x
                
            e = self.mDecoratorExternalNodes.getElements( node_id,
                                                          self.getHeaderWidth() + x,
                                                          self.getHeaderHeight() + y,
                                                          self.getHeaderWidth() + x_label,
                                                          self.getHeaderHeight() + y )
            
            self.addElements(e)

    #####################################################################
    def plotInternalNodes( self ):
        """plot terminal node labels."""

        
        for node_id in self.mTree.chain.keys():

            node = self.mTree.node( node_id )
            if node.succ == []: continue
            
            x = self.mNodeWidthsEnd[node_id]
            y = self.mNodeHeights[node_id]

            e = self.mDecoratorInternalNodes.getElements( node_id,
                                                          self.getHeaderWidth() + x,
                                                          self.getHeaderHeight() + y )
            
            self.addElements(e)
        
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
    def writeFooter( self ):
        """write footer.

        The footer contains the legend.
        """
        total_branch_length = (max(self.mNodeWidthsEnd) - min(self.mNodeWidthsStart)) / self.mBranchScaleFactor
        
        self.mFooterX = self.getHeaderWidth()
        self.mFooterY = self.getHeaderHeight() + self.mDataHeight + self.mSeparatorHeight

        ruler_start = self.mFooterX
        ruler_end = self.mFooterX + int(total_branch_length * self.mBranchScaleFactor)

        if "ruler" in self.mRulerElements:
            ## full length ruler with tick marks and labels
            e = SVGdraw.line( ruler_start,
                              self.mFooterY + self.mRulerTickSize + 1,
                              ruler_end,
                              self.mFooterY + self.mRulerTickSize + 1,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              stroke_width = 1 )
            self.addElement( e )

            ## get reasonable intervalls

        increment = self.mRulerIncrement * self.mBranchScaleFactor

        ## adjust increment for extremely long trees
        if (ruler_end - ruler_start) / increment > 1000:
            increment = (ruler_end - ruler_start) / 1000.0
            self.mRulerIncrement = increment / self.mBranchScaleFactor
        
        if "right-ticks" in self.mRulerElements:

            x = ruler_end
            while x >= ruler_start:
                e = SVGdraw.line( x,
                                  self.mFooterY,
                                  x,
                                  self.mFooterY + 2 * self.mRulerTickSize + 1,
                                  stroke = "rgb(%i,%i,%i)" % BLACK,
                                  stroke_width = 1 )
                self.addElement( e )
                x -= self.mRulerIncrement * self.mBranchScaleFactor
                
            self.mFooterY += 2 * self.mRulerTickSize + 1 + self.mSeparatorHeight

        if "left-ticks" in self.mRulerElements:
            
            x = ruler_start
            while x <= ruler_end:
                e = SVGdraw.line( x,
                                  self.mFooterY,
                                  x,
                                  self.mFooterY + 2 * self.mRulerTickSize + 1,
                                  stroke = "rgb(%i,%i,%i)" % BLACK,
                                  stroke_width = 1 )
                self.addElement( e )
                x += increment
                
            self.mFooterY += 2 * self.mRulerTickSize + 1 + self.mSeparatorHeight

        if "scale" in self.mRulerElements:

            w = int(self.mRulerIncrement * self.mBranchScaleFactor)
            
            e = SVGdraw.line( ruler_end,
                              self.mFooterY + self.mRulerTickSize + 1,
                              ruler_end - w,
                              self.mFooterY + self.mRulerTickSize + 1,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              stroke_width = 1 )
            
            self.addElement( e )
            
            e = SVGdraw.line( ruler_end,
                              self.mFooterY,
                              ruler_end,
                              self.mFooterY + 2 * self.mRulerTickSize + 1,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              stroke_width = 1 )
            
            self.addElement( e )

            e = SVGdraw.line( ruler_end - w,
                              self.mFooterY,
                              ruler_end - w,
                              self.mFooterY + 2 * self.mRulerTickSize + 1,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              stroke_width = 1 )
            
            self.addElement( e )

            e = SVGdraw.text( ruler_end - w / 2,
                              self.mFooterY + 2 * self.mRulerTickSize + 1 + self.mRulerFontSize ,
                              self.mRulerFormat % self.mRulerIncrement,
                              self.mRulerFontSize,
                              self.mRulerFont,
                              stroke = "rgb(%i,%i,%i)" % BLACK,
                              text_anchor = "middle" )

            self.addElement(e)
            self.mFooterY += 2 * self.mRulerTickSize + 1 + self.mRulerFontSize 
            
        e, self.mFooterX, self.mFooterY = self.mDecoratorExternalNodes.getLegend( self.mFooterX, self.mFooterY )
        self.addElements(e)
        e, self.mFooterX, self.mFooterY = self.mDecoratorInternalNodes.getLegend( self.mFooterX, self.mFooterY )
        self.addElements(e)

    #####################################################################        
    def finalizePlot( self ):
        """build plot."""

        self.plotTree()

        x = self.mDataWidth
        for b in self.mBorderDecorators:
            self.addElements( b.getElements( x, self.getHeaderHeight(), self.mNodeHeights ) )
            x += b.getWidth()
                              
        # self.writeTitle()
        
        self.writeFooter()

        
    #####################################################################        
    def writeToFile(self, outfile):
        """write svg image to file.
        """
        self.finalizePlot()

        kk = self.mElements.keys()
        kk.sort()
        kk.reverse()

        ## make sure the image size is ok
        min_x, min_y, max_x, max_y = 0, 0, 0, 0
        
        for k in kk:
            for e in self.mElements[k]:
                for x in ('x', 'x2', 'x1'):
                    if x in e.attributes:
                        v = e.attributes[x] 
                        min_x = min(min_x, v )
                        max_x = max(max_x, v )
                for y in ('y', 'y2', 'y1'):
                    if y in e.attributes:
                        v = e.attributes[y]
                        min_y = min(min_y, v )
                        max_y = max(max_y, v )

        min_x, min_y = int(math.floor(min_x)), int(math.floor(min_y))
        max_x, max_y = int(math.floor(max_x)), int(math.floor(max_y))        

        for k in kk:
            for e in self.mElements[k]:
                for x in ('x', 'x2', 'x1'):
                    if x in e.attributes:
                        e.attributes[x] -= min_x
                for x in ('y', 'y2', 'y1'):
                    if y in e.attributes:
                        e.attributes[y] -= min_y

        ## now add all the elements
        self.mRoot = SVGdraw.drawing()
        self.mDraw = SVGdraw.svg( (0, 0, self.mPageWidth - min_x, self.mPageHeight - min_y ) , "100%", "100%" )
                        
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

    parser = E.OptionParser( version = "%prog version: $Id: SVGTree.py 2784 2009-09-10 11:41:14Z andreas $")

    parser.add_option("-i", "--title", dest="title", type="string",
                      help="page title.")
    parser.add_option("-f", "--footer", dest="footer", type="string",
                      help="page footer.")
    parser.add_option("-s", "--filename-tree", dest="filename_tree", type="string",
                      help="filename with tree."  )
    parser.add_option("-t", "--tree", dest="tree", type="string",
                      help="tree."  )
    parser.add_option( "-r", "--species-regex", dest="species_regex", type="string" ,
                       help="regular expression to extract species from identifier.")
    parser.add_option( "--colour-by-species", dest="colour_by_species", action="store_true",
                      help="colour by species."  )

    parser.add_option( "--branch-scale", dest="branch_scale", type="float",
                      help="branch length scale factor."  )
    parser.add_option( "--height-scale", dest="height_scale", type="float",
                      help="height scale factor."  )
    
    parser.set_defaults(
        titles = "",
        title = "",
        footer = "",
        filename_tree = None,
        species_regex ="^([^|]+)\|",
        colour_by_species = None,
        tree = None,
        branch_scale = 0,
        height_scale = 0,
        )

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    if options.filename_tree:
        tree_lines = open(options.filename_tree, "r").readlines()
    elif options.tree:
        tree_lines = options.tree
    else:
        raise "please supply a species tree."

    nexus = TreeTools.Newick2Nexus( tree_lines )
    Tree.updateNexus( nexus )
    tree = nexus.trees[0]

    if options.loglevel >= 2:
        tree.display()
    
    plot = SVGTree( tree )

    plot.setBranchScale( options.branch_scale )
    plot.setHeightScale( options.height_scale )    
    
    if options.colour_by_species:
        rx = re.compile(options.species_regex)
        extract_species = lambda x: rx.search( x ).groups()[0]
        plot.setDecoratorExternalNodes( NodeDecoratorBySpecies( tree, extract_species= extract_species ) )
    
    plot.initializePlot()
    
    plot.writeToFile(sys.stdout)
    
    Experiment.Stop()

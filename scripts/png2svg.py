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
png2svg.py - arrange several bitmat images into a single image
==============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script takes a collection of .png files and outputs an .svg file
arranging the images side by side or in different layers.

Usage
-----

Example::

   python png2svg.py --help

Type::

   python png2svg.py --help

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
import subprocess
import random

from types import *
""" program $Id: png2svg.py 2782 2009-09-10 11:40:29Z andreas $

arrange png files in an svg image.
"""
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.SVGdraw as SVGdraw

###------------------------------------------------------------------------------------------------

###################################################################

class Image:

    def __init__(self, filename, x, y, width, height, opacity = 1.0 ):
        self.mFilename = filename
        self.mWidth = width
        self.mHeight = height
        self.mX = x
        self.mY = y

        self.mOpacity = opacity

    def getSVG( self, xoffset, yoffset ):
        return SVGdraw.image( self.mFilename,
                              self.mX + xoffset,
                              self.mY + yoffset,
                              self.mWidth,
                              self.mHeight,
                              opacity = self.mOpacity )
    
class PngPlot:

    def __init__(self, filenames ):

        self.mElements = {}

        self.mFilenames = filenames
        
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

        self.mSeparatorWidth  = 10
        self.mSeparatorHeight = 10

        self.mKeepAspectRatio = True

        self.mMapImage2Coords = []

    #####################################################################
    def setCanvasSize( self, width, height ):
        """set canvas width and size."""
        self.mDefaultWidth, self.mDefaultHeight = width, height

    #####################################################################
    def addElement( self, element, plane = 0):
        """add element to list in plane.
        
        The list is later sorted by priority.
        """
        if plane not in self.mElements:
            self.mElements[plane] = []
            
        self.mElements[plane].append( element )

    def getHeaderHeight(self):
        return 0
    def getHeaderWidth(self):
        return 0
    def getFooterHeight(self):
        return 0
    def getFooterWidth(self):
        return 0

    #####################################################################
    def addElements( self, elements, plane = 0):
        """add multiple elments to a plane."""
        for e in elements:
            self.addElement( e, plane )

    #####################################################################                        
    def calculateCanvasSize( self ):
        """calculate the size of the canvas."""

        ## set the page size
        self.mPageWidth  = self.getHeaderWidth() + self.mDataWidth + self.getFooterWidth()
        self.mPageHeight = self.getHeaderHeight() + self.mDataHeight + self.getFooterHeight()
    
    #####################################################################                        
    def initializePlot( self ):
        """set various coordinates in the plot.

        Note:

        Width  = X = coordinate 1
        Height = Y = coordinate 2
        
        """

        self.calculateCoordinates()
        
        self.calculateCanvasSize( )

    #####################################################################
    def setTitle( self, title ):
        """set title."""
        self.mTitle = title

    #####################################################################
    def setFooter( self, footer ):
        """set footer."""
        self.mFooter = footer

    #####################################################################
    def writeImages( self ):
        """write images."""

        xoffset = self.getHeaderWidth()
        yoffset = self.getHeaderHeight()
        
        for i in range(len(self.mFilenames) ):
            image = self.mMapImage2Coords[i]
            e = image.getSVG( xoffset, yoffset )
            self.addElement( e )
        
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
    def finalizePlot( self ):
        """build plot."""

        self.writeImages()
        self.writeTitle()
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

class Grid( PngPlot ):

    """grid layout."""
    def __init__(self, filenames, force_square = False, *args, **kwargs ):

        PngPlot.__init__(self, filenames, *args, **kwargs )

        self.mForceSquare = force_square
    
    #####################################################################                        
    def calculateCoordinates( self ):
        """calculate coordinates."""

        ## square image
        nimages = len(self.mFilenames)
        nimages_x = int(math.sqrt( float(nimages) * (float(self.mDefaultWidth) / self.mDefaultHeight) ) )
        if self.mForceSquare:
            nimages_y = nimages_x
            del self.mFilenames[nimages_y * nimages_x:]
        else:
            nimages_y = int(math.ceil(float(nimages) / nimages_x))

        image_width  = int(float((self.mDefaultWidth  - (nimages_x * self.mSeparatorWidth)) / nimages_x))
        ## add +1 for image_height, as the we use left/upper coordiantes
        image_height = int(float((self.mDefaultHeight - ((nimages_y) * self.mSeparatorHeight)) / (nimages_y)))

        if self.mKeepAspectRatio:
            m = min(image_width, image_height)
            image_width, image_height = m, m

        coord_y = -(self.mSeparatorHeight + image_height)
        coord_x = 0
        
        for x in range(len(self.mFilenames)):
            if x % nimages_x == 0:
                coord_x = 0
                coord_y += self.mSeparatorHeight + image_height
                
            self.mMapImage2Coords.append( Image(self.mFilenames[x], coord_x, coord_y, image_width, image_height) )
            
            coord_x += self.mSeparatorWidth + image_width 
            
        self.mDataWidth = self.mDefaultWidth
        self.mDataHeight = self.mDefaultHeight

class RandomLayers( PngPlot ):

    """elements are positioned randomly in several layers."""
    
    def __init__(self, filenames,
                 num_layers,
                 image_width = 100,
                 image_height = 100,
                 *args, **kwargs ):

        PngPlot.__init__(self, filenames, *args, **kwargs )

        self.mNumLayers = num_layers
        self.mImageHeight = image_height
        self.mImageWidth = image_width
        
    #####################################################################                        
    def calculateCoordinates( self ):
        """calculate coordinates."""

        nimages = len(self.mFilenames)
        images_per_layer = int(float( nimages ) / self.mNumLayers)
        opacity_intervall = 1.0 / self.mNumLayers

        opacity = 0 
        for x in range(nimages):
            
            if x % images_per_layer == 0:
                opacity += opacity_intervall

            image_width  = self.mImageWidth
            image_height = self.mImageHeight

            coord_x = random.randint( 0, self.mDefaultWidth - image_width )
            coord_y = random.randint( 0, self.mDefaultHeight - image_height )            
            
            self.mMapImage2Coords.append( Image(self.mFilenames[x], coord_x, coord_y, image_width, image_height, opacity = opacity ) )
            
        self.mDataWidth = self.mDefaultWidth
        self.mDataHeight = self.mDefaultHeight
    

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: png2svg.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-i", "--title", dest="title", type="string",
                      help="page title.")
    parser.add_option("-f", "--footer", dest="footer", type="string",
                      help="page footer.")
    parser.add_option("-l", "--layout", dest="layout", type="choice",
                      choices=("grid", "random", "random-layers",),
                      help="layout to choose.")
    parser.add_option( "--num-layers", dest="num_layers", type="int",
                      help="number of layers.")
    parser.add_option( "--canvas-width", dest="canvas_width", type="int",
                      help="canvas width [default: 1000].")
    parser.add_option( "--canvas-height", dest="canvas_height", type="int",
                      help="canvas height [default: 1000].")
    parser.add_option( "--image-width", dest="image_width", type="int",
                      help="image width.")
    parser.add_option( "--image-height", dest="image_height", type="int",
                      help="image height.")
    parser.add_option( "--force-square", dest="force_square", action="store_true",
                      help="force square layout.")

    
    parser.set_defaults(
        titles = "",
        title = "",
        footer = "",
        layout = "grid",
        force_square = False,
        num_layers = 2,
        image_width = 100,
        image_height = 100,
        canvas_width = 1000,
        canvas_height = 1000,
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    if len(args) > 0:
        ## read filenames from the command line
        filenames = args
    else:
        ## read filenames from stdin
        filenames = map(lambda x: x[:-1], sys.stdin.readlines())

    if options.loglevel >= 1:
        options.stdlog.write("# arranging %i images.\n" % len(filenames))

    if options.layout == "grid":
        plot = Grid( filenames, force_square = options.force_square )
    elif options.layout == "random":
        plot = RandomLayers( filenames,
                             num_layers = 1,
                             image_width = options.image_width,
                             image_height = options.image_height)
    elif options.layout == "random-layers":
        plot = RandomLayers( filenames,
                             num_layers = options.num_layers,
                             image_width = options.image_width,
                             image_height = options.image_height)
    else:
        raise "unknown layout %s" % options.layout

    plot.setCanvasSize( options.canvas_width, options.canvas_height )
    
    plot.initializePlot()

    plot.writeToFile(sys.stdout)
    
    E.Stop()

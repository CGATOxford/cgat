##########################################################################
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
##########################################################################
'''
optic/plot_duplications.py -
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

   python optic/plot_duplications.py --help

Type::

   python optic/plot_duplications.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import math
import tempfile
import bisect
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.SVGdraw as SVGdraw
import CGAT.TreeTools as TreeTools

USAGE = """python optic/plot_duplications.py [OPTIONS] < stdin > stdout

plots a wheel plot of duplications in SVG format.

The script requirse the following input files:

--contig-sizes: a filename with contig sizes to output. This is a tab-separated
    file with contig<tab>size on each line.

In standard input the duplications are listed as a tab-separated file
containing the following three fields:

1. cluster_id: a numerical identifier for the cluster

2. locations: a ";" separated list of duplicated genes with their
   locations. Each entry is of the following format:
   gene_id:chr:strand:start:end

3. tree: the tree of duplicated genes. This is used to color arcs by
   tree height.

The trees are in newick format.

Output: the svg formatted file is writtent to standard output. Comment lines
start with a "#". You can remove this by setting -v 0 or --verbose=0.
"""

# some definitions for the layout of the picture
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)
YELLOW = (255, 255, 0)
CYAN = (0, 255, 255)
PURPLE = (255, 0, 255)
GREY = (128, 128, 128)

MAX_GREY = 240
GREY_COLORS = map(lambda x: (x, x, x), range(0, MAX_GREY))

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

DEFAULT_COLOR_SEQUENCE = BLACK
DEFAULT_COLOR_ANNOTATION = BLACK


class DuplicationPlot:

    def __init__(self,
                 contigs,
                 map_contig2size,
                 num_entries,
                 ):

        self.mElements = {}

        # save contig sizes
        self.contigs = contigs
        self.mMapContig2Size = map_contig2size

        # get row size
        self.mNumEntries = num_entries

        # colours
        self.startColour = RED
        self.mMiddleColour = YELLOW
        self.mStopColour = GREEN
        self.mCircleColour = GREY

        # a space
        self.mSeparator = 10

        # widht and hight of a row/column
        self.mRowHeight = 200
        self.mColWidth = 200

        # maximum size of a box
        self.mMaxBoxSize = 100
        self.mRevertSize = True

        # Height and width get set
        self.mHeaderFontSize = self.mRowHeight / 2
        self.mHeaderFont = "Verdana"

        # write a grid line every five rows
        self.mRowTicks = 5
        self.mColTicks = 5

        # footer
        self.mFooterFrom = 10
        self.mFooterFontSize = self.mRowHeight / 2
        self.mFooterFont = "Verdana"
        self.mFooter = None

        # Title
        self.mTitleFontSize = self.mRowHeight
        self.mTitleFont = "Verdana"
        self.mTitle = None

        self.mRadius = 3000
        self.mRadiusIncrement = 40
        self.mRadiusStart = 3000

        self.mMinValue = 0.0
        self.mMaxValue = 1.0

        # planes
        self.mPlaneJoins = 4
        self.mPlaneGrid = 5

    #####################################################################
    def addElement(self, element, plane=0):
        """add element to list in plane.

        The list is later sorted by priority.
        """
        if plane not in self.mElements:
            self.mElements[plane] = []
        self.mElements[plane].append(element)

    #####################################################################
    def initializePlot(self):
        """set various coordinates in the plot."""

        # height of col header
        self.mHeaderHeight = self.mHeaderFontSize / 2

        if self.mTitle:
            self.mHeaderHeight += self.mSeparator + self.mTitleFontSize

        # width of row header
        self.mHeaderWidth = self.mHeaderFontSize / 2

        # height of footer:
        self.mFooterHeight = 2 * \
            max(self.mFooterFontSize, self.mMaxBoxSize) + self.mSeparator

        if self.mFooter:
            self.mFooterHeight += self.mSeparator + self.mFooterFontSize

        self.mTotalSize = reduce(
            lambda x, y: x + y, self.mMapContig2Size.values())
        # estimated height and width of data section
        self.mDataWidth = self.mNumEntries / 2 * \
            self.mRadiusIncrement + self.mRadiusStart
        self.mDataHeight = self.mDataWidth

        # build coordinates for rows and columns
        self.buildMap2Position()

        # build colour map
        self.buildColourMap()

        self.mPageWidth = self.mHeaderWidth + \
            self.mDataWidth + 1 * self.mSeparator
        self.mPageHeight = self.mHeaderHeight + self.mDataHeight + \
            self.mFooterHeight + 2 * self.mSeparator

        # smallest radius
        self.mRadiusStart = self.mRadius

        # radius to fall back on
        self.mRadiusFallBack = self.mRadius

        self.mDataMiddleX, self.mDataMiddleY = self.mDataWidth / \
            2, self.mDataHeight / 2

        # maximum coordinate of longest range for fallback
        self.mLastMax = 0

        # maximum used radius in slice
        self.mRadiusMax = self.mRadius

        # maximum coordinate of previous range
        self.mPreviousMax = 0

        self.mLastChr = None

        # size of marker
        self.mMarkerWidth = self.mRadiusIncrement / 4

        # keeps track of written separators, so that they do not get written
        # multiply
        self.mSeparators = {}

        # minimum distance for resolving dots
        self.mMinDistance = self.mTotalSize / 360.0

        # whether or not last entry was cis
        self.mPreviousCis = False

        # number format for legend
        self.mFormatNumberLegend = "%5.2f"

    #####################################################################
    def setTitle(self, title):
        """set title."""
        self.mTitle = title

    #####################################################################
    def setFooter(self, footer):
        """set footer."""
        self.mFooter = footer

    #####################################################################
    def writeTitle(self):
        """write title into plot."""

        if self.mTitle:
            e = SVGdraw.text(self.mPageWidth / 2,
                             self.mTitleFontSize,
                             self.mTitle,
                             self.mTitleFontSize,
                             self.mTitleFont,
                             stroke="rgb(%i,%i,%i)" % BLACK,
                             text_anchor="middle")

            self.addElement(e)

    #####################################################################
    def buildColourMap(self):
        """build map of heights to colours.

        This is two gradients:
        1: first half: start to middle
        2: second half: middle to end
        """

        num_steps = 100
        self.mColourIncrement = float(
            (self.mMaxValue - self.mMinValue)) / num_steps / 2

        self.mColours = []
        self.mColourThresholds = []

        value = self.mMinValue
        d = map(lambda x, y: (x - y) / num_steps,
                self.mMiddleColour, self.startColour)
        for x in range(num_steps):
            self.mColours.append((self.startColour[0] + x * d[0],
                                  self.startColour[1] + x * d[1],
                                  self.startColour[2] + x * d[2]))
            value += self.mColourIncrement
            self.mColourThresholds.append(value)

        self.mColours.append(self.mMiddleColour)
        value += self.mColourIncrement
        self.mColourThresholds.append(value)

        d = map(lambda x, y: (x - y) / num_steps,
                self.mStopColour, self.mMiddleColour)
        for x in range(num_steps):
            self.mColours.append((self.mMiddleColour[0] + x * d[0],
                                  self.mMiddleColour[1] + x * d[1],
                                  self.mMiddleColour[2] + x * d[2]))

            value += self.mColourIncrement
            self.mColourThresholds.append(value)

        self.mColours.append(self.mStopColour)

    #####################################################################
    def buildMap2Position(self):
        """build map of contigs to start."""
        self.contig2Position = {}

        l = 0
        for x in self.contigs:
            self.contig2Position[x] = l
            l += self.mMapContig2Size[x]

    #####################################################################
    def addSeparator(self):
        """add separator on circles."""
        if self.mRadius not in self.mSeparators:
            e = SVGdraw.circle(self.mDataMiddleX, self.mDataMiddleY, self.mRadius, fill="none",
                               stroke="rgb(%i,%i,%i)" % self.mCircleColour,
                               stroke_width=1)
            self.addElement(e, self.mPlaneGrid)
        self.mSeparators[self.mRadius] = 1

    #####################################################################
    def getPosition(self, chr, strand, residue):
        """return position on arc where a residue lies."""
        return self.contig2Position[chr] + residue

    #####################################################################
    def getAngle(self, pos):
        """return angle on arc for position pos."""

        return 360.0 * pos / self.mTotalSize

    #####################################################################
    def getPosOnArc(self, angle, radius):
        """return position on arc centered on x, y at angle and radius."""

        degPerRad = math.pi / 180.0

        dx = radius * math.cos(angle * degPerRad)
        dy = radius * math.sin(angle * degPerRad)

        return dx + self.mDataMiddleX, dy + self.mDataMiddleY

    #####################################################################
    def getPosOnArcForRadius(self, x, radius):
        """get pos on arc for coordinate x."""
        return self.getPosOnArc(self.getAngle(x), radius)

    #####################################################################
    def addDuplication(self, entries, map_gene2pos, height,
                       url=None,
                       with_separator=True,
                       link_to_previous=False,
                       quality2symbol={},
                       quality2mask={}):
        """add a dot in row/col."""

        mi, ma = None, 0

        pos = bisect.bisect(self.mColourThresholds, height)
        master_colour = self.mColours[pos]

        chrs = {}
        points = []
        for species, transcript, gene, quality in entries:

            chr, strand, first_res, last_res = map_gene2pos[gene]
            chrs[chr] = 1
            pos1 = self.getPosition(chr, strand, first_res)
            pos2 = self.getPosition(chr, strand, last_res)

            a = min(pos1, pos2)
            b = max(pos1, pos2)

            if mi is None:
                mi = a
            else:
                mi = min(a, mi)
            ma = max(b, ma)

            points.append((pos1, pos2, gene, quality))

        # decide whether we need to increment the radius
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
                    # overlap due to close proximitiy
                    self.mRadius += self.mRadiusIncrement
                    if with_separator:
                        self.addSeparator()

                    # true overlap
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
                if with_separator:
                    self.addSeparator()

        if cis and self.mPreviousCis and (link_to_previous or is_overlap):

            # print old_x1, old_y1, old_x2, old_y2, new_x1, new_y2, new_x2,
            # new_y2
            r1 = old_radius
            r2 = self.mRadius
            old_x1, old_y1 = self.getPosOnArcForRadius(self.mPreviousMin, r1)
            old_x2, old_y2 = self.getPosOnArcForRadius(self.mPreviousMax, r1)
            new_x1, new_y1 = self.getPosOnArcForRadius(mi, r2)
            new_x2, new_y2 = self.getPosOnArcForRadius(ma, r2)

            if link_to_previous:

                # print cis, entries, chrs

                d = SVGdraw.pathdata(old_x1, old_y1)
                d.ellarc(r1, r1, 0, 0, 1, old_x2, old_y2)
                d.line(new_x2, new_y2)
                d.ellarc(r2, r2, 0, 0, 1, new_x1, new_y1)
                d.line(old_x1, old_y1)

                e = SVGdraw.path(d,
                                 fill="rgb(%i,%i,%i)" % GREY,
                                 stroke="rgb(%i,%i,%i)" % GREY,
                                 stroke_width=1)

            else:

                # get points of center
                old_x1, old_y1 = self.getPosOnArcForRadius(
                    (self.mPreviousMax + self.mPreviousMin) / 2, r1)
                new_x1, new_y1 = self.getPosOnArcForRadius((ma + mi) / 2, r2)

                # lines for interleaved spans: skip
                if not ((self.mPreviousMin < mi and self.mPreviousMax > ma) or
                        (self.mPreviousMin > mi and self.mPreviousMax < ma)):

                    e = SVGdraw.line(old_x1, old_y1,
                                     new_x1, new_y1,
                                     fill="none",
                                     stroke="rgb(%i,%i,%i)" % GREY,
                                     stroke_width=self.mMarkerWidth / 2)
                else:
                    e = None

            if e:
                self.addElement(e, self.mPlaneJoins)

        self.mPreviousMin = mi
        self.mPreviousMax = ma
        self.mPreviousCis = cis

        self.mRadiusMax = max(self.mRadius, self.mRadiusMax)

        # draw points
        link_colour = master_colour
        link_width = 10
        for pos1, pos2, gene, quality in points:

            angle = self.getAngle((pos1 + pos2) / 2)

            x, y = self.getPosOnArc(angle, self.mRadius)

            try:
                symbol = quality2symbol[quality]
            except KeyError:
                symbol = "rect"

            if quality in quality2mask:
                colour = GREY
                link_colour = GREY
                link_width = 1
            else:
                colour = master_colour

            if symbol == "circle":
                ee = SVGdraw.circle(x, y, self.mMarkerWidth,
                                    fill="rgb(%i,%i,%i)" % colour,
                                    stroke="black",
                                    stroke_width=1)
            elif symbol == "rect":
                ee = SVGdraw.rect(x - self.mMarkerWidth / 2, y - self.mMarkerWidth / 2,
                                  self.mMarkerWidth, self.mMarkerWidth,
                                  fill="rgb(%i,%i,%i)" % colour,
                                  stroke="black",
                                  stroke_width=1)

            if url:
                e = SVGdraw.link(url % gene)
                e.addElement(ee)
            else:
                e = ee
            self.addElement(e)

        angle1 = self.getAngle(mi)
        angle2 = self.getAngle(ma)

        x1, y1 = self.getPosOnArc(angle1, self.mRadius)
        x2, y2 = self.getPosOnArc(angle2, self.mRadius)

        d = SVGdraw.pathdata(x1, y1)
        if cis:
            d.ellarc(self.mRadius, self.mRadius, 0, 0, 1, x2, y2)
        else:
            d.ellarc(self.mRadius * 2, self.mRadius * 2, 0, 0, 0, x2, y2)

        e = SVGdraw.path(d,
                         fill="none",
                         stroke="rgb(%i,%i,%i)" % link_colour,
                         stroke_width=link_width)

        self.addElement(e)

    #####################################################################
    def pushRadius(self):
        """push radius as fallback."""

        self.mRadiusFallBack = self.mRadiusMax + self.mRadiusIncrement
        self.mRadius = self.mRadiusFallBack
        # self.addSeparator()
        self.mLastMax = 0
        self.mPreviousMax = 0

    #####################################################################
    def writeColHeaders(self):
        """write row headers."""

        current_x = self.mColWidth / 2
        current_y = self.mHeaderHeight
        i = 0
        for name in self.mColNames:
            e = SVGdraw.text(current_x,
                             current_y,
                             name,
                             self.mHeaderFontSize,
                             self.mHeaderFont,
                             stroke="rgb(%i,%i,%i)" % BLACK,
                             text_anchor="start",
                             transform="rotate(-45,%i,%i)" % (current_x, current_y))

            self.addElement(e)

            current_x += self.mColWidth

    #####################################################################
    def writeGrid(self):
        """add grid lines."""

        # print last circle
        self.mRadius = self.mRadiusMax + self.mRadiusIncrement
        self.addSeparator()

        # print separators
        for contig, pos in self.contig2Position.items():

            pos = self.getPosition(contig, "+", 0)

            angle = self.getAngle(pos)

            x, y = self.getPosOnArc(angle, self.mRadius)

            e = SVGdraw.line(self.mDataMiddleX,
                             self.mDataMiddleY,
                             x,
                             y,
                             stroke="rgb(%i,%i,%i)" % BLACK)

            self.addElement(e, self.mPlaneGrid)

            x, y = self.getPosOnArc(
                angle, (self.mRadius + self.mRadiusStart) / 2)

            e = SVGdraw.text(x,
                             y,
                             contig,
                             self.mHeaderFontSize,
                             self.mHeaderFont,
                             stroke="rgb(%i,%i,%i)" % BLACK,
                             text_anchor="start",
                             transform="rotate(%i,%i,%i)" % (angle, x, y))

            self.addElement(e)

    #####################################################################
    def writeRowHeaders(self):
        """write row headers."""

        current_x = self.mDataWidth + self.mSeparator

        current_y = self.mHeaderHeight + self.mSeparator + self.mHeaderFontSize

        for name in self.mRowNames:
            e = SVGdraw.text(current_x,
                             current_y,
                             name,
                             self.mHeaderFontSize,
                             self.mHeaderFont,
                             stroke="rgb(%i,%i,%i)" % BLACK,
                             text_anchor="start")

            self.addElement(e)

            current_y += self.mRowHeight

        self.mHeaderWidth = max(
            map(len, self.mRowNames)) * self.mHeaderFontSize / 2

    #####################################################################
    def writeFooter(self):
        """write footer.

        The footer contains the legend.
        """
        current_x = self.mFooterFrom
        current_y = self.mHeaderHeight + self.mDataHeight + 2 * self.mSeparator

        self.mFooterBoxSize = 30
        self.mNumTicks = 20
        for x in range(len(self.mColourThresholds)):

            e = SVGdraw.rect(current_x,
                             current_y,
                             self.mFooterBoxSize,
                             self.mFooterBoxSize,
                             fill="rgb(%i,%i,%i)" % self.mColours[x],
                             stroke="rgb(%i,%i,%i)" % self.mColours[x])

            self.addElement(e)

            if x % self.mNumTicks == 0:

                e = SVGdraw.line(current_x,
                                 current_y,
                                 current_x,
                                 current_y + self.mFooterBoxSize,
                                 stroke="rgb(%i,%i,%i)" % BLACK,
                                 stroke_width=5)
                self.addElement(e)

                e = SVGdraw.text(current_x,
                                 current_y - self.mFooterBoxSize,
                                 self.mFormatNumberLegend % self.mColourThresholds[
                                     x],
                                 self.mFooterFontSize,
                                 self.mFooterFont,
                                 stroke="rgb(%i,%i,%i)" % BLACK,
                                 text_anchor="start")
                self.addElement(e)

            current_x += self.mFooterBoxSize

        ###########################################################
        if self.mFooter:

            current_y += max(self.mFooterFontSize,
                             self.mMaxBoxSize) + self.mSeparator

            e = SVGdraw.text(self.mPageWidth / 2,
                             current_y + self.mFooterFontSize,
                             self.mFooter,
                             self.mFooterFontSize,
                             self.mFooterFont,
                             stroke="rgb(%i,%i,%i)" % BLACK,
                             text_anchor="middle")

            self.addElement(e)

    #####################################################################
    def finalizePlot(self):
        """write remaining parts of the plot."""

        # corrected height and width of data section
        self.mDataWidth = 2 * self.mRadiusMax
        self.mDataHeight = self.mDataWidth

        self.writeGrid()
        self.writeTitle()
        self.writeFooter()

    #####################################################################
    def writeToFile(self, outfile):
        """write svg image to file.
        """
        self.finalizePlot()

        self.mRoot = SVGdraw.drawing()
        self.mDraw = SVGdraw.svg(
            (0, 0, self.mPageWidth, self.mPageHeight), "100%", "100%")

        kk = self.mElements.keys()
        kk.sort()
        kk.reverse()
        for k in kk:
            for e in self.mElements[k]:
                self.mDraw.addElement(e)

        self.mRoot.setSVG(self.mDraw)

        tfile = tempfile.mktemp()

        self.mRoot.toXml(tfile)

        lines = open(tfile, "r").readlines()

        outfile.write(string.join(lines, ""))
        outfile.write("\n")

        os.remove(tfile)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/plot_duplications.py 2782 2009-09-10 11:40:29Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-e", "--header-names", dest="headers", action="store_true",
                      help="first row is a header [ignored].")
    parser.add_option("-t", "--title", dest="title", type="string",
                      help="page title.")
    parser.add_option("-f", "--footer", dest="footer", type="string",
                      help="page footer.")
    parser.add_option("-c", "--contigs-tsv-file", dest="filename_contig_sizes", type="string",
                      help="filname with contig sizes.")
    parser.add_option("-r", "--territory-extension", dest="radius", type="int",
                      help="radius.")
    parser.add_option("-i", "--flank-increment-size", dest="radius_increment", type="int",
                      help="radius increment.")
    parser.add_option("-u", "--url", dest="url", type="string",
                      help="string to build url for annotation.")
    parser.add_option("--min-contig", dest="min_contig_size", type="string",
                      help="minimum contig size to delineate.")

    parser.add_option("--min-value", dest="min_value", type="float",
                      help="minimum branch length.")

    parser.add_option("--max-value", dest="max_value", type="float",
                      help="maximum branch length.")

    parser.set_defaults(
        filename_contig_sizes=None,
        headers=False,
        titles="",
        pattern_filename=None,
        title="",
        footer="",
        radius=3000,
        min_value=0.0,
        max_value=0.2,
        url=None,
        radius_increment=40,
        min_contig_size=10000,
        remove_empty_contigs=True,
        separator="|",
        quality2symbol={'CG': "circle", 'PG': "circle", 'SG': "circle"},
        quality2mask=(
            "RG", "CP", "PP", "SP", "RP", "CF", "PF", "SF", "UG", "UP", "UF", "BF", "UK"),
        sort_by_size=True,
        input_format="pairwise",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if options.filename_contig_sizes:
        map_contig2size = IOTools.ReadMap(open(options.filename_contig_sizes, "r"),
                                          map_functions=(str, int))

    # read data and get contigs that are used (i.e.: remove empty contigs)
    chrs = {}
    lines = sys.stdin.readlines()

    if options.remove_empty_contigs:
        for line in lines:
            if line[0] == "#":
                continue

            d = line[:-1].split("\t")

            cluster_id, in_locations, in_tree = d[:3]

            for l in in_locations.split(";"):
                gene_id, chr, strand, sbjct_from, sbjct_to = l.split(":")
                if chr not in map_contig2size:
                    continue
                chrs[chr] = 1
        for k in map_contig2size.keys():
            if k not in chrs:
                del map_contig2size[k]

    k = map_contig2size.keys()

    if len(k) == 0:
        E.Stop()
        sys.exit(0)

    k.sort()

    if options.sort_by_size:
        k.sort(lambda x, y: cmp(map_contig2size[x], map_contig2size[y]))

    plot = DuplicationPlot(k,
                           map_contig2size,
                           num_entries=0)

    plot.mRadiusIncrement = options.radius_increment
    plot.mRadius = options.radius
    plot.mMaxValue = options.max_value
    plot.mMinValue = options.min_value

    if options.title:
        plot.setTitle(options.title)
    if options.footer:
        plot.setFooter(options.footer)

    plot.initializePlot()

    data = []

    if options.input_format == "pairwise":

        # read data from pairwise analysis
        # format is: cluster_id, locations of duplications, tree of
        # duplications

        for line in lines:
            if line[0] == "#":
                continue

            d = line[:-1].split("\t")

            cluster_id, in_locations, in_tree = d[:3]

            mi, ma = 0, 0
            found = False
            n = 0
            chrs = {}
            for l in in_locations.split(";"):
                gene_id, chr, strand, sbjct_from, sbjct_to = l.split(":")
                if chr not in map_contig2size:
                    continue
                chrs[chr] = 1
                sbjct_from, sbjct_to = int(sbjct_from), int(sbjct_to)

                xi = plot.getPosition(chr, strand, sbjct_from)
                xa = plot.getPosition(chr, strand, sbjct_to)

                if not mi:
                    mi = xi
                else:
                    mi = min(mi, xi)

                n += 1
                ma = max(ma, xa)
                found = True

            if not found:
                continue
            cis = len(chrs) == 1
            if options.loglevel >= 2:
                options.stdlog.write("# adding duplications in cluster %s: %s with tree %s\n" % (
                    cluster_id, in_locations, in_tree))
            data.append((cis, n, mi, ma, cluster_id, in_locations, in_tree))

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

        if not map_gene2location:
            continue

        tree = TreeTools.Newick2Tree(in_tree)

        # the last subset is all nodes again.
        s = TreeTools.GetSubsets(tree)

        is_first = True
        for children, height, branchlength in s[:-1]:
            if len(children) == 1:
                continue
            c = map(lambda x: x.split(options.separator), children)
            plot.addDuplication(c, map_gene2location, height,
                                url=options.url,
                                with_separator=is_first,
                                link_to_previous=not is_first,
                                quality2symbol=options.quality2symbol,
                                quality2mask=options.quality2mask)
            is_first = False

    plot.writeToFile(sys.stdout)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

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
optic/plot_multiple_synteny.py -
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

   python optic/plot_multiple_synteny.py --help

Type::

   python optic/plot_multiple_synteny.py --help

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
import CGAT.Experiment as E
import CGAT.SVGdraw as SVGdraw
import CGAT.Synteny as Synteny
import pgdb

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

###################################################################

###################################################################


class PlotMultipleSynteny:

    def __init__(self,
                 orthologs,
                 ):

        self.mElements = {}

        # plot layout options

        # plot links between all genomes
        self.mPlotLinks = "all"

        # orient by first
        self.mOrientation = "local"

        # save orthologs
        self.mOrthologs = orthologs

        # font form factor. Use this to estimate text length
        self.mFontFactor = 2.0 / 3.0

        # footer
        self.mFooterFrom = 10
        self.mFooterFontSize = 12
        self.mFooterFont = "Verdana"
        self.mFooter = None

        # Title
        self.mTitleFontSize = 12
        self.mTitleFont = "Verdana"
        self.mTitle = None

        # Legend
        self.mLegendFontSize = 12
        self.mLegendFont = "Verdana"
        self.mLegendHorizontalSeparator = 10
        self.mLegendVerticalSeparator = 10

        # Header
        self.mHeaderFontSize = 24
        self.mHeaderFont = "Verdana"

        # Contig
        self.contigFontSize = 12
        self.contigFont = "Verdana"

        # Row/column labels
        self.mLabelFontSize = 20
        self.mLabelFont = "Verdana"

        self.mUseGeneOrder = True

        self.mHeaderHeight = 100
        self.mHeaderWidth = 100

        self.mFooterHeight = 100
        self.mFooterWidth = 100

        self.mBlockSize = 10
        self.mHorizontalSeparator = 100
        self.mVerticalSeparator = 10

        # smallest contig to plot should have at least 0.1 % entries
        # of maximum contig
        self.mMinRelativeContigSize = 0.1

        # maximum distance between synteny blocks on same contig
        self.mMaxSyntenyDifference = 10000000

        # species sort order
        self.mSpeciesList = None

        # url for genes
        self.mGeneUrl = "%s"

        # planes to use
        self.mPlaneLinks = 10

        self.mSymbols = ("filled-circle", "filled-box",
                         "open-circle", "open-box",
                         "open-lefttriangle", "filled-lefttriangle",
                         "open-righttriangle", "filled-righttriangle",
                         "open-uptriangle", "filled-uptriangle",
                         "open-downtriangle", "filled-downtriangle")

    #####################################################################
    def setSpeciesList(self, species):
        """set sort order for species."""
        self.mSpeciesList = species

    #####################################################################
    def setOrientation(self, orientation):
        """set sort order for species."""
        self.mOrientation = orientation

    #####################################################################
    def addElement(self, element, plane=0):
        """add element to list in plane.

        The list is later sorted by priority.
        """
        if plane not in self.mElements:
            self.mElements[plane] = []

        self.mElements[plane].append(element)

    #####################################################################
    def getOrthologGroups(self):
        """return list of ortholog groups."""
        groups = []
        for g in self.mOrthologs.values():
            groups += list(g.keys())

        return list(set(groups))

    #####################################################################
    def getGenes(self, species):
        """return list of genes for species."""
        genes = []
        if species not in self.mOrthologs:
            return []

        for ortho in self.mOrthologs[species].values():
            for gene in ortho:
                genes.append(gene)
        return genes

    #####################################################################
    def buildColourWheel(self):

        self.mColours = []

        nPalSize = 5
        nColorShift = 3

        pi = 3.1415926535

        # Generate color palette.
        # Colors represent the outer edge of a 24 bit (maximum) color palette.
        # Colors are equally spaced around the edge of the color circle.
        # Use CArray object to use functions and dynamic sizing.

        # Arrays to translate palette size and color shift
        # from parameter array integers. Initial integers
        # correspond to radio button values in parameter
        # dialog box.
        PalSizeArray = (24, 48, 96, 192, 384, 768)
        # Color shift blue, green, red, yellow, magenta, cyan
        ColorShiftArray = (0, 511, 255, 383, 127, 639)

        nPalSize = PalSizeArray[nPalSize]
        nColorShift = ColorShiftArray[nColorShift]

        # First range check
        if nPalSize > 768:
            nPalSize = 768
        if nPalSize < 24:
            nPalSize = 24

        # Now set the array size. May be one too big after modifications.
        # At least not harmful.

        # Then set the palette
        for n in range(0, nPalSize):

            # nFactor calculate position of n in maximized palette ( 768 total colors )
            # Assure that it runs 0 to 768
            nFactor = 768 * n / nPalSize
            nFactor = nFactor + nColorShift
            # Adjust for color shift
            if nFactor >= 768:
                nFactor = nFactor - 768

            # Color calculations include sine/cosine functions to create
            # colors at constant radius from center of color circle. This
            # creates a constant brightness level. Math not optimized, left in
            # readable state.
            if nFactor <= 256:
                # Red increasing, no green, blue decreasing
                nRed = int(255.0 * (math.sin((nFactor / 256.0) * pi / 2)))
                nGreen = 0
                nBlue = int(255.0 * (math.cos((nFactor / 256.0) * pi / 2)))

            elif ((nFactor > 256.0) and(nFactor <= 512)):
                # Red decreasing, green increasing, blue 0.
                nRed = int(
                    255.0 * (math.cos(((nFactor - 256.0) / 256.0) * pi / 2)))
                nGreen = int(
                    255.0 * (math.sin(((nFactor - 256.0) / 256.0) * pi / 2)))
                nBlue = 0

            elif nFactor > 512:
                # Red 0, green decreasing, blue increasing
                nRed = 0
                nGreen = int(
                    255.0 * (math.cos(((nFactor - 512) / 256.0) * pi / 2)))
                nBlue = int(
                    255.0 * (math.sin(((nFactor - 512) / 256.0) * pi / 2)))

            self.mColours.append((nRed, nGreen, nBlue))

    #####################################################################
    def buildColourLinear(self):

        groups = self.getOrthologGroups()

        num_steps = len(groups)

        self.startColour = RED
        self.mStopColour = YELLOW

        d = map(lambda x, y: (x - y) / num_steps,
                self.mStopColour, self.startColour)

        for x in range(num_steps):
            self.mColours.append((self.startColour[0] + x * d[0],
                                  self.startColour[1] + x * d[1],
                                  self.startColour[2] + x * d[2]))

    #####################################################################
    def buildColourMap(self):
        """build colour and symbol schema.
        """

        self.buildColourWheel()
        self.mMapGroup2Colour = {}

        groups = self.getOrthologGroups()

        increment = int(math.floor(len(self.mColours) / 16.0))
        offset = 0
        x = 0
        nsymbol = 0
        for group in groups:
            self.mMapGroup2Colour[group] = (
                self.mColours[x], self.mSymbols[nsymbol])
            x += increment
            if x >= len(self.mColours):
                x = offset
                nsymbol += 1
                if nsymbol >= len(self.mSymbols):
                    nsymbol = 0

    #####################################################################
    def initializePlot(self):
        """set various coordinates in the plot.

        Note:

        Width  = X = coordinate 1
        Height = Y = coordinate 2

        """

        if not self.mSpeciesList:
            self.mSpeciesList = list(self.mOrthologs.keys())

        # prepare data
        self.mNumSpecies = len(self.mSpeciesList)

        max_width = 0
        for species in self.mOrthologs.keys():
            max_width = max(
                max_width, self.getRowWidth(self.getGenes(species)))

        # calculate size of diagonal
        self.mDataWidth = max_width
        self.mDataHeight = self.mNumSpecies * \
            (self.mBlockSize + self.mHorizontalSeparator)

        # write row header - sets mHeaderWidth
        self.writeRowHeader()

        # set the page size
        self.mPageWidth = self.mHeaderWidth + \
            self.mDataWidth + self.mFooterWidth
        self.mPageHeight = self.mHeaderHeight + \
            self.mDataHeight + self.mFooterHeight

        # add plot elements
        self.buildColourMap()

        # build map of group to genes
        self.mMapGroup2Genes = {}

        for g in self.mOrthologs.values():
            for group_id, genes in g.items():
                if group_id not in self.mMapGroup2Genes:
                    self.mMapGroup2Genes[group_id] = []
                self.mMapGroup2Genes[group_id] += genes

        self.plotGenes()
        self.plotLinks()

    #####################################################################
    def getRowWidth(self, genes):
        """width width of a row."""
        return len(genes) * (self.mBlockSize + self.mVerticalSeparator)

    #####################################################################
    def addGroupStart(self, x, y, contig, start):
        """plot start of contig."""

        t = "%s:%i" % (contig, start)
        max_x = x + len(t) * self.contigFontSize * self.mFontFactor

        # plot line
        e = SVGdraw.line(x, y, max_x, y,
                         stroke="rgb(%i,%i,%i)" % BLACK,
                         )

        self.addElement(e)
        e = SVGdraw.text(x, y, t,
                         self.contigFontSize,
                         self.contigFont,
                         stroke="rgb(%i,%i,%i)" % BLACK,
                         text_anchor="left")

        self.addElement(e)

        return max_x

    def getBlocks(self, genes):
        """get blocks in genes."""
        blocks = []

        b = []
        last_contig = None
        last_end = None
        for gene in genes:
            if gene.contig != last_contig or abs(gene.mFrom - last_end) > self.mMaxSyntenyDifference:
                if b:
                    blocks.append(b)
                b = []

            last_contig = gene.contig
            last_end = gene.mTo
            b.append(gene)

        if b:
            blocks.append(b)

        return blocks

    #####################################################################
    def sortBlock(self, block, last_block):
        """sort blocks by reference"""

        # check if sort order should be reversed by counting the number of
        # consistent intervals
        groups1 = map(lambda x: x.mOrthologId, block)
        groups2 = map(lambda x: x.mOrthologId, last_block)

        def getDifferences(groups):
            last = 0
            d = []
            for x in groups:
                d.append(last - int(x))
                last = int(x)
            return d

        diffs1 = set(getDifferences(groups1))
        diffs2 = set(getDifferences(groups2))
        groups2.reverse()
        rdiffs2 = set(getDifferences(groups2))

        if len(diffs1.intersection(rdiffs2)) > len(diffs1.intersection(diffs2)):
            block.reverse()

    #####################################################################
    def getSortedGenes(self, genes, last_genes):
        """sort genes according to last_genes."""

        # split list of genes into blocks
        blocks = self.getBlocks(genes)
        last_blocks = self.getBlocks(last_genes)

        # get best matching blocks
        best_matches = []
        for x in range(len(blocks)):
            for y in range(len(last_blocks)):
                groups1 = set(map(lambda x: x.mOrthologId, blocks[x]))
                groups2 = set(map(lambda x: x.mOrthologId, last_blocks[y]))
                best_matches.append(
                    (y, -len(groups1.intersection(groups2)), x))

        best_matches.sort()

        new_genes = []
        assigned = set()
        last_y = None

        for y, l, x in best_matches:

            if l == 0:
                continue

            if x in assigned:
                continue
            if y == last_y:
                continue

            self.sortBlock(blocks[x], last_blocks[y])
            new_genes += blocks[x]
            assigned.add(x)
            last_y = y

        for x in range(len(blocks)):
            if x not in assigned:
                new_genes += blocks[x]

        return new_genes

    #####################################################################
    def getSortedGenesByBlockSize(self, genes):
        """sort genes according to last_genes."""

        # split list of genes into blocks
        blocks = self.getBlocks(genes)

        blocks.sort(lambda x, y: cmp(len(x), len(y)))
        blocks.reverse()

        new_genes = []
        for b in blocks:
            new_genes += b

        return new_genes

    #####################################################################
    def plotGene(self, x, y, gene=None, group_id=None):
        """plot a gene at x,y."""

        if gene:
            group_id = gene.mOrthologId

        colour, format = self.mMapGroup2Colour[group_id]

        filled, shape = format.split("-")

        if filled == "filled":
            fill = "rgb(%i,%i,%i)" % colour
            stroke = "rgb(%i,%i,%i)" % BLACK
            stroke_width = 1

        elif filled == "open":
            fill = "rgb(%i,%i,%i)" % WHITE
            stroke = "rgb(%i,%i,%i)" % colour
            stroke_width = self.mBlockSize / 2

        if shape == "circle":
            ee = SVGdraw.circle(x, y,
                                self.mBlockSize,
                                fill=fill,
                                stroke=stroke,
                                stroke_width=stroke_width)
        elif shape == "box":
            ee = SVGdraw.rect(x - self.mBlockSize, y - self.mBlockSize,
                              2 * self.mBlockSize, 2 * self.mBlockSize,
                              fill=fill,
                              stroke=stroke,
                              stroke_width=stroke_width)
        elif shape == "lefttriangle":
            ee = SVGdraw.polygon(((x - self.mBlockSize, y),
                                  (x + self.mBlockSize, y - self.mBlockSize),
                                  (x + self.mBlockSize, y + self.mBlockSize)),
                                 fill=fill,
                                 stroke=stroke,
                                 stroke_width=stroke_width)
        elif shape == "righttriangle":
            ee = SVGdraw.polygon(((x - self.mBlockSize, y - self.mBlockSize),
                                  (x + self.mBlockSize, y),
                                  (x - self.mBlockSize, y + self.mBlockSize)),
                                 fill=fill,
                                 stroke=stroke,
                                 stroke_width=stroke_width)
        elif shape == "uptriangle":
            ee = SVGdraw.polygon(((x, y - self.mBlockSize),
                                  (x + self.mBlockSize, y + self.mBlockSize),
                                  (x - self.mBlockSize, y + self.mBlockSize)),
                                 fill="rgb(%i,%i,%i)" % WHITE,
                                 stroke="rgb(%i,%i,%i)" % colour,
                                 stroke_width=self.mBlockSize / 2)
        elif shape == "downtriangle":
            ee = SVGdraw.polygon(((x, y + self.mBlockSize),
                                  (x + self.mBlockSize, y - self.mBlockSize),
                                  (x - self.mBlockSize, y - self.mBlockSize)),
                                 fill="rgb(%i,%i,%i)" % WHITE,
                                 stroke="rgb(%i,%i,%i)" % colour,
                                 stroke_width=self.mBlockSize / 2)

        return ee

    #####################################################################
    def plotGenes(self):
        """plot orthologs orthologwise."""

        # sort all entries by coordinates
        # and build coordinates
        map_gene2coord = {}

        y = self.mHeaderHeight

        self.mMapGene2Coord = {}

        reference_genes = None

        for species in self.mSpeciesList:
            genes = self.getGenes(species)
            genes.sort(
                lambda x, y: cmp((x.contig, x.mFrom), (y.contig, y.mFrom)))

            if reference_genes:
                genes = self.getSortedGenes(genes, reference_genes)
            else:
                genes = self.getSortedGenesByBlockSize(genes)

            x = self.getRowWidth(genes)

            x = self.mHeaderWidth + (self.mDataWidth - x) / 2

            last_contig, last_end = None, 0

            # build list of genes to skip (note: this can be improved
            # by working in blocks from the start.)
            blocks = self.getBlocks(genes)

            for block in blocks:

                nlinks = 0
                for gene in block:
                    nlinks += len(self.mMapGroup2Genes[gene.mOrthologId])
                if nlinks == len(block):
                    continue

                x = self.addGroupStart(x, y, gene.contig, gene.mFrom)

                # count number of links to be made
                for gene in block:

                    self.mMapGene2Coord[gene] = (x, y)

                    ee = self.plotGene(x, y, gene)

                    if self.mGeneUrl:
                        e = SVGdraw.link(self.mGeneUrl % gene)
                        e.addElement(ee)
                    else:
                        e = ee

                    self.addElement(e)
                    x += self.mBlockSize + self.mVerticalSeparator

                    last_contig, last_end = gene.contig, gene.mTo

            y += self.mBlockSize + self.mHorizontalSeparator

            # set orientation
            if self.mOrientation == "local":
                reference_genes = genes
            elif self.mOrientation == "first":
                if not reference_genes:
                    reference_genes = genes

    #####################################################################
    def plotLinks(self):

        if self.mPlotLinks == "all":
            for x in range(len(self.mSpeciesList) - 1):
                for y in range(x + 1, len(self.mSpeciesList)):
                    if abs(x - y) == 1:
                        self.plotLinksForPair(
                            self.mSpeciesList[x], self.mSpeciesList[y], BLACK, self.mPlaneLinks)
                    else:
                        self.plotLinksForPair(
                            self.mSpeciesList[x], self.mSpeciesList[y], GREEN, self.mPlaneLinks + 1)

        elif self.mPlotLinks == "pairs":
            for x in range(len(self.mSpeciesList) - 1):
                self.plotLinksForPair(
                    self.mSpeciesList[x], self.mSpeciesList[x + 1])

    #####################################################################
    def plotLinksForPair(self, species1, species2, colour=BLACK, plane=0):
        """plot links."""

        if species1 not in self.mOrthologs or species2 not in self.mOrthologs:
            return

        for group1 in self.mOrthologs[species1]:
            if group1 in self.mOrthologs[species2]:
                for g1 in self.mOrthologs[species1][group1]:
                    for g2 in self.mOrthologs[species2][group1]:
                        x1, y1 = self.mMapGene2Coord[g1]
                        x2, y2 = self.mMapGene2Coord[g2]

                        e = SVGdraw.line(x1, y1, x2, y2,
                                         stroke="rgb(%i,%i,%i)" % colour,
                                         )

                        self.addElement(e, plane=plane)

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
    def writeRowHeader(self):
        """write row header

        The row header contains the species names
        """
        x = 0
        y = self.mHeaderHeight

        max_l = 0

        for species in self.mSpeciesList:
            max_l = max(
                max_l, len(species) * self.mHeaderFontSize) * self.mFontFactor
            e = SVGdraw.text(x, y, species,
                             self.mHeaderFontSize,
                             self.mHeaderFont,
                             stroke="rgb(%i,%i,%i)" % BLACK,
                             text_anchor="left")

            self.addElement(e)

            y += self.mBlockSize + self.mHorizontalSeparator

        self.mHeaderWidth = max_l

    #####################################################################
    def writeFooter(self):
        """write footer.

        The footer contains the legend.
        """

        x = self.mHeaderWidth + self.mBlockSize

        group_ids = self.mMapGroup2Colour.keys()
        group_ids.sort()
        ncolumns = 10
        l = int(math.ceil(len(group_ids) / float(ncolumns)))

        max_d = 0

        for a in range(0, ncolumns):

            y = self.mHeaderHeight + self.mDataHeight
            for b in range(a * l, min(len(group_ids), (a + 1) * l)):

                group_id = group_ids[b]

                e = self.plotGene(x, y, group_id=group_id)
                self.addElement(e)

                t = str(group_id)

                max_d = max(
                    len(t) * self.mLegendFontSize * self.mFontFactor, max_d)

                e = SVGdraw.text(x + self.mBlockSize + self.mLegendVerticalSeparator,
                                 y,
                                 str(group_id),
                                 self.mLegendFontSize,
                                 self.mLegendFont,
                                 stroke="rgb(%i,%i,%i)" % BLACK,
                                 text_anchor="left")

                self.addElement(e)

                y += self.mLegendFontSize + self.mLegendHorizontalSeparator

            x += max_d + 2 * self.mBlockSize + \
                2 * self.mLegendVerticalSeparator

    #####################################################################
    def finalizePlot(self):
        """write remaining parts of the plot."""

        # corrected height and width of data section
# self.writeGrid()
# self.writeRowLabels()
# self.writeColLabels()

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

###################################################################


def GetMap(values):
    map_a2b = {}
    map_b2a = []

    for v in values:
        map_a2b[v] = len(map_b2a)
        map_b2a.append(v)

    return map_a2b, map_b2a


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/plot_multiple_synteny.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-t", "--title", dest="title", type="string",
                      help="page title.")
    parser.add_option("-f", "--footer", dest="footer", type="string",
                      help="page footer.")
    parser.add_option("-g", "--groups", dest="groups", type="int", action="append",
                      help="groups to extract from database.")
    parser.add_option("-r", "--synteny-radius", dest="synteny_radius", type="int",
                      help="number of genes to include of either side of a group in synteny plot.")
    parser.add_option("-n", "--table-groups", dest="table_name_groups", type="string",
                      help="table name with group information.")
    parser.add_option("--table-synteny", dest="table_name_synteny", type="string",
                      help="table name with synteny information.")
    parser.add_option("-s", "--schema", dest="schema", type="string",
                      help="schema name in psql database.")
    parser.add_option("-u", "--url", dest="url", type="string",
                      help="string to build url for annotation.")
    parser.add_option("--species", dest="species", type="string",
                      help="comma separted list of species order.")
    parser.add_option("--orientation", dest="orientation", type="choice",
                      choices=("local", "first"),
                      help="how to orient synteny blocks.")

    parser.set_defaults(
        filename_contig_sizes=None,
        headers=False,
        titles="",
        pattern_filename=None,
        title="",
        footer="",
        url=None,
        min_contig_size=10000,
        remove_empty_contigs=True,
        separator="|",
        species=None,
        quality2symbol={'CG': "circle", 'PG': "circle", 'SG': "circle"},
        quality2mask=(
            "RG", "CP", "PP", "SP", "RP", "CF", "PF", "SF", "UG", "UP", "UF", "BF", "UK"),
        sort_by_size=True,
        groups=[],
        table_name_synteny="synteny",
        schema=None,
        synteny_radius=5,
        orientation="local",
    )

    (options, args) = E.Start(
        parser, add_pipe_options=True, add_database_options=True)

    if options.species:
        options.species = options.species.split(",")

    if options.schema and options.table_name_synteny:
        dbhandle = pgdb.connect(options.psql_connection)
        statement = """SELECT
        group_id,
        schema,
        gene_id,
        sbjct_token AS contig,
        sbjct_strand AS strand,export_sbjct_genome_from AS start,
        export_sbjct_genome_to AS end
        FROM %s.%s
        WHERE synteny_id IN
        (SELECT * FROM %s.get_near_synteny(%%i,%i ))
        ORDER BY schema, sbjct_token, export_sbjct_genome_from""" % (options.schema,
                                                                     options.table_name_synteny,
                                                                     options.schema,
                                                                     options.synteny_radius)

        lines = []
        for group_id in options.groups:
            cc = dbhandle.cursor()
            cc.execute(statement % group_id)
            lines += cc.fetchall()
            cc.close()

        group_ids = list(set(map(lambda x: x[0], lines)))

        statement = """SELECT group_id,
        schema,
        gene_id,
        sbjct_token AS contig,
        sbjct_strand AS strand,export_sbjct_genome_from AS start,
        export_sbjct_genome_to AS end
        FROM %s.%s WHERE group_id IN ('%s')""" % (options.schema, options.table_name_synteny,
                                                  "','".join(map(str, group_ids)))

        cc = dbhandle.cursor()
        cc.execute(statement)
        lines = cc.fetchall()
        cc.close()

        orthologs = Synteny.ReadOrthologsPerSpecies(lines)
    else:
        orthologs = Synteny.ReadOrthologsPerSpecies(sys.stdin)

    plot = PlotMultipleSynteny(orthologs)

    plot.setSpeciesList(options.species)
    plot.setOrientation(options.orientation)

    plot.initializePlot()

    plot.writeToFile(sys.stdout)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

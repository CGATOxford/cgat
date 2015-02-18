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
optic/plot_synteny.py -
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

   python optic/plot_synteny.py --help

Type::

   python optic/plot_synteny.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import tempfile
import CGAT.Experiment as E
import CGAT.SVGdraw as SVGdraw
import CGAT.Synteny as Synteny
import CGAT.IndexedFasta as IndexedFasta

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

###################################################################


class PlotSynteny:

    """syntey plot with gene order."""

    def __init__(self,
                 orthologs1, orthologs2,
                 sorted_contigs1, sorted_contigs2,
                 contigs1, contigs2
                 ):

        self.mElements = {}

        # save orthologs
        self.mOrthologs1 = orthologs1
        self.mOrthologs2 = orthologs2

        # save contigs
        self.mSortedContigs1 = sorted_contigs1
        self.mSortedContigs2 = sorted_contigs2

        self.mContentsContig1 = contigs1
        self.mContentsContig2 = contigs2

        # footer
        self.mFooterFrom = 10
        self.mFooterFontSize = 12
        self.mFooterFont = "Verdana"
        self.mFooter = None

        # Title
        self.mTitleFontSize = 12
        self.mTitleFont = "Verdana"
        self.mTitle = None

        # Row/column labels
        self.mLabelFontSize = 20
        self.mLabelFont = "Verdana"

        self.mHeaderHeight = 100
        self.mHeaderWidth = 100

        self.mFooterHeight = 100
        self.mFooterWidth = 100

        # size of one dot
        self.mDotRadius = 4

        # smallest contig to plot should have at least 0.1 % entries
        # of maximum contig
        self.mMinRelativeContigSize = 0.1

    #####################################################################
    def addElement(self, element, plane=0):
        """add element to list in plane.

        The list is later sorted by priority.
        """
        if plane not in self.mElements:
            self.mElements[plane] = []

        self.mElements[plane].append(element)

    #####################################################################
    def getDiagonal(self, axis, contig_contents, sorted_contigs):
        """get length of diagnal.

        a. sum number of genes per contig.
        b. sum of residues per contig.

        do not add contigs, which are less than mMinContigSize% of largest contigs size
        """
        l = 0
        map_contig2start = {}
        map_contig2size = {}
        max_n = 0

        for contig in sorted_contigs:
            n = len(contig_contents[contig])
            if max_n:
                if max_n * self.mMinRelativeContigSize > n:
                    continue
            else:
                max_n = n
            map_contig2start[contig] = l
            map_contig2size[contig] = n
            l += n

        return l, map_contig2start, map_contig2size

    #####################################################################
    def getCoords1(self, ortholog):
        """get coordinates for an ortholog on first axis."""
        return self.mHeaderWidth + self.mMapContig2Start1[ortholog.contig] + ortholog.mRank

    #####################################################################
    def getCoords2(self, ortholog):
        """get coordinates for an ortholog on second axis."""
        return self.mHeaderHeight + self.mMapContig2Start2[ortholog.contig] + ortholog.mRank

    #####################################################################
    def assignCoordinates(self):
        """assign coordinates to orthologs."""

        # set rank for orthologs within a contig
        Synteny.SetRankToPosition(self.mOrthologs1, self.mContentsContig1)
        Synteny.SetRankToPosition(self.mOrthologs2, self.mContentsContig2)

    #####################################################################
    def initializePlot(self):
        """set various coordinates in the plot.

        Note:

        Width  = X = coordinate 1
        Height = Y = coordinate 2

        """

        # assign coordinates to orthologs
        self.assignCoordinates()

        # calculate size of diagonal
        self.mDataWidth,  self.mMapContig2Start1, self.mMapContig2Size1 = self.getDiagonal(
            0, self.mContentsContig1, self.mSortedContigs1)
        self.mDataHeight, self.mMapContig2Start2, self.mMapContig2Size2 = self.getDiagonal(
            1, self.mContentsContig2, self.mSortedContigs2)

        # set the page size
        self.mPageWidth = self.mHeaderWidth + \
            self.mDataWidth + self.mFooterWidth
        self.mPageHeight = self.mHeaderHeight + \
            self.mDataHeight + self.mFooterHeight

        # set the font size factor
        self.mFontSizeFactor = int(
            max(1.0, float(min(self.mPageWidth, self.mPageHeight)) / 1000))

        self.plotOrthologs()

    #####################################################################

    def plotOrthologs(self):
        """plot orthologs orthologwise."""

        # Do the plot
        for ortholog_id in self.mOrthologs1.keys():

            for o1 in self.mOrthologs1[ortholog_id]:

                try:
                    x = self.getCoords1(o1)
                except KeyError:
                    continue

                for o2 in self.mOrthologs2[ortholog_id]:

                    try:
                        y = self.getCoords2(o2)
                    except KeyError:
                        continue

                    e = SVGdraw.circle(x, y,
                                       self.mDotRadius,
                                       fill="rgb(%i,%i,%i)" % BLACK,
                                       stroke="rgb(%i,%i,%i)" % BLACK,
                                       stroke_width=1)

                    self.addElement(e)

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
    def writeColLabels(self):
        """write column headers."""

        y = self.mHeaderHeight

        for contig in self.mSortedContigs1:

            if contig not in self.mMapContig2Start1:
                continue

            x = self.mHeaderWidth + self.mMapContig2Start1[contig]
            e = SVGdraw.text(x,
                             y,
                             contig,
                             self.mLabelFontSize * self.mFontSizeFactor,
                             self.mLabelFont,
                             stroke="rgb(%i,%i,%i)" % BLACK,
                             text_anchor="start")

            self.addElement(e)

    #####################################################################
    def writeRowLabels(self):
        """write column headers."""

        x = self.mHeaderWidth

        for contig in self.mSortedContigs2:

            if contig not in self.mMapContig2Start2:
                continue

            y = self.mHeaderHeight + \
                self.mMapContig2Start2[contig] + self.mMapContig2Size2[contig]
            e = SVGdraw.text(x,
                             y,
                             contig,
                             self.mLabelFontSize * self.mFontSizeFactor,
                             self.mLabelFont,
                             stroke="rgb(%i,%i,%i)" % BLACK,
                             text_anchor="start",
                             transform="rotate(-90,%i,%i)" % (x, y))

            self.addElement(e)

    #####################################################################
    def writeGrid(self):
        """add grid lines."""

        min_x = self.mHeaderWidth
        max_x = min_x + self.mDataWidth
        min_y = self.mHeaderHeight
        max_y = min_y + self.mDataHeight

        for contig in self.mSortedContigs1:
            if contig not in self.mMapContig2Start1:
                continue
            x = self.mMapContig2Start1[contig]

            e = SVGdraw.line(min_x + x, min_y, min_x + x, max_y,
                             stroke="rgb(%i,%i,%i)" % GREEN,
                             )
            self.addElement(e)

        self.addElement(SVGdraw.line(max_x, min_y, max_x, max_y,
                                     stroke="rgb(%i,%i,%i)" % GREEN))

        for contig in self.mSortedContigs2:
            if contig not in self.mMapContig2Start2:
                continue
            y = self.mMapContig2Start2[contig]
            e = SVGdraw.line(min_x, min_y + y, max_x, min_y + y,
                             stroke="rgb(%i,%i,%i)" % GREEN,
                             )
            self.addElement(e)

        self.addElement(SVGdraw.line(min_x, max_y, max_x, max_y,
                                     stroke="rgb(%i,%i,%i)" % GREEN))

    #####################################################################
    def writeFooter(self):
        """write footer.

        The footer contains the legend.
        """
        return

    #####################################################################
    def finalizePlot(self):
        """write remaining parts of the plot."""

        # corrected height and width of data section
        self.writeGrid()
        self.writeRowLabels()
        self.writeColLabels()

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


class PlotSyntenyGenomic(PlotSynteny):

    """synteny plot with genomic coordinates."""

    def __init__(self, contig_sizes1, contig_sizes2, *args, **kwargs):

        self.contigSizes1 = contig_sizes1
        self.contigSizes2 = contig_sizes2

        PlotSynteny.__init__(self, *args)

        self.mDataWidth = 1000
        self.mDataHeight = 1000

    #####################################################################
    def getDiagonal(self, axis, contig_contents, sorted_contigs):
        """get length of diagnal.

        sum of residues per contig.

        do not add contigs, which are less than mMinContigSize% of largest contigs size
        """
        total = 0
        map_contig2start = {}
        map_contig2size = {}

        if axis == 0:
            contig_sizes = self.contigSizes1
            sf = self.mScaleFactor1
        else:
            contig_sizes = self.contigSizes2
            sf = self.mScaleFactor2

        for contig in sorted_contigs:
            length = contig_sizes[contig] * sf
            map_contig2start[contig] = total
            map_contig2size[contig] = length
            total += length

        return total, map_contig2start, map_contig2size

    #####################################################################
    def getCoords1(self, ortholog):
        """get coordinates for an ortholog on first axis."""
        return self.mHeaderWidth + self.mMapContig2Start1[ortholog.contig] + ortholog.mRank

    #####################################################################
    def getCoords2(self, ortholog):
        """get coordinates for an ortholog on second axis."""
        return self.mHeaderHeight + self.mMapContig2Start2[ortholog.contig] + ortholog.mRank

    #####################################################################
    def assignCoordinates(self):
        """assign coordinates to orthologs with in a contig.

        use mean genomic location of an ortholog.
        """

        def getTotalContigLength(contig_sizes, contents_contig):
            """also updates contents_contig."""

            min_length = max(contig_sizes.values()) * \
                self.mMinRelativeContigSize
            total = 0
            for contig, length in contig_sizes.items():
                if length < min_length:
                    if contig in contents_contig:
                        del contents_contig[contig]
                    continue
                total += length

            return total

        self.mScaleFactor1 = float(
            self.mDataWidth) / getTotalContigLength(self.contigSizes1, self.mContentsContig1) * self.mDotRadius
        self.mScaleFactor2 = float(
            self.mDataHeight) / getTotalContigLength(self.contigSizes2, self.mContentsContig2) * self.mDotRadius

        for contig, oo in self.mContentsContig1.items():
            for o in oo:
                o.mRank = int((o.mFrom + o.mTo) / 2.0) * self.mScaleFactor1

        for contig, oo in self.mContentsContig2.items():
            for o in oo:
                o.mRank = int((o.mFrom + o.mTo) / 2.0) * self.mScaleFactor2


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/plot_synteny.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-e", "--header-names", dest="headers", action="store_true",
                      help="first row is a header [ignored].")
    parser.add_option("-t", "--title", dest="title", type="string",
                      help="page title.")
    parser.add_option("-f", "--footer", dest="footer", type="string",
                      help="page footer.")
    parser.add_option("--genome-file1", dest="genome_file1", type="string",
                      help="filename with genome information for first set of orthologs.")
    parser.add_option("--genome-file2", dest="genome_file2", type="string",
                      help="filename with genome information for second set of orthologs.")
    parser.add_option("-c", "--coordinates", dest="coordinates", type="choice",
                      choices=("genes", "genomic"),
                      help="coordinate system to use.")
    parser.add_option("-u", "--url", dest="url", type="string",
                      help="string to build url for annotation.")
    parser.add_option("--min-contig", dest="min_contig_size", type="string",
                      help="minimum contig size to delineate.")

    parser.set_defaults(
        genome_file1=None,
        genome_file2=None,
        use_coordinates=False,
        headers=False,
        titles="",
        pattern_filename=None,
        title="",
        footer="",
        url=None,
        min_contig_size=10000,
        remove_empty_contigs=True,
        separator="|",
        quality2symbol={'CG': "circle", 'PG': "circle", 'SG': "circle"},
        quality2mask=(
            "RG", "CP", "PP", "SP", "RP", "CF", "PF", "SF", "UG", "UP", "UF", "BF", "UK"),
        sort_by_size=True,
        coordinates="genes",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if len(args) != 2:
        raise "please supply two files with synteny information."

    options.filename1, options.filename2 = args

    orthologs1 = Synteny.ReadOrthologs(open(options.filename1, "r"))
    orthologs2 = Synteny.ReadOrthologs(open(options.filename2, "r"))

    # make sure, that there are no orthologs without a corresponding
    # entry in the other set.
    Synteny.CleanOrthologs(orthologs1, orthologs2)

    # get sorted list of contigs and list of orthologs per contig
    sorted_contigs1, sorted_contigs2, contigs1, contigs2 = Synteny.SortAssignOrthologs(
        orthologs1, orthologs2)

    if options.coordinates == "genes":
        plot = PlotSynteny(orthologs1,
                           orthologs2,
                           sorted_contigs1,
                           sorted_contigs2,
                           contigs1,
                           contigs2)
    elif options.coordinates == "genomic":

        if not options.genome_file1 or not options.genome_file2:
            raise "please supply two genomes for contig sizes (option --genome-file1 and --genome-file2)"

        contig_sizes1 = IndexedFasta.IndexedFasta(
            options.genome_file1).getContigSizes()
        contig_sizes2 = IndexedFasta.IndexedFasta(
            options.genome_file2).getContigSizes()

        plot = PlotSyntenyGenomic(contig_sizes1,
                                  contig_sizes2,
                                  orthologs1,
                                  orthologs2,
                                  sorted_contigs1,
                                  sorted_contigs2,
                                  contigs1,
                                  contigs2)

    plot.initializePlot()

    plot.writeToFile(sys.stdout)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

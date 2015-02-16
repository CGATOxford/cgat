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
optic/analyze_orthology_multiple.py -
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

   python optic/analyze_orthology_multiple.py --help

Type::

   python optic/analyze_orthology_multiple.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import re
import time
import cStringIO

from types import *

import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools
import CGAT.Tree as Tree
import CGAT.SVGTree as SVGTree
import CGAT.SVGdraw as SVGdraw
import CGAT.IOTools as IOTools
import CGAT.Stats as Stats
import scipy
import scipy.stats
import numpy
import CGAT.Histogram as Histogram

import CGAT.SVGDuplicationsWheel as SVGDuplicationsWheel

# Do full import - sort out names later
from CGAT.TreeReconciliation import *

SVG_MAP_TYPE2COLOUR = {"Speciation": SVGTree.GREEN,
                       "SpeciationDeletion": SVGTree.OLIVE,
                       "Transcripts": SVGTree.BLUE,
                       "DuplicationLineage": SVGTree.CYAN,
                       "Duplication": SVGTree.PLUM,
                       "DuplicationDeletion": SVGTree.ORANGE,
                       "DuplicationInconsistency": SVGTree.BROWN,
                       "Outparalogs": SVGTree.YELLOW,
                       "InconsistentTranscripts": SVGTree.RED,
                       "Masked": SVGTree.GREY,
                       "Inconsistency": SVGTree.PURPLE}

###################################################################
###################################################################
###################################################################
# map of categories to quality indices
###################################################################
MAP_FILTER2QUALITY = {'all': set(("CG", "SG", "PG", "RG", "UG", "CP", "SP", "PP", "RP", "UP", "SF", "CF", "PF", "UF", "BF")),
                      'genes': set(("CG", "SG", "PG", "UG", "SF", "CF", "PF", "UF")),
                      'pseudogenes': set(("RG", "CP", "SP", "PP", "RP", "UP", "BF"))}

###################################################################
###################################################################
###################################################################
# contigs with name are considered junk. Compare in lower case.
###################################################################
MAP_CONTIG2JUNK = {'chr3l_random': 1,
                   'chr3r_random': 1,
                   'chr3h_random': 1,
                   'chr2r_random': 1,
                   'chr2l_random': 1,
                   'chr2h_random': 1,
                   'chr2h': 1,
                   'chr3h': 1,
                   'chr4h': 1,
                   'chrxh': 1,
                   'chryh': 1,
                   'chr4_random': 1,
                   'chru_random': 1,
                   'chru': 1,
                   'chruh': 1,
                   'chrm': 1,
                   'chrx_random': 1,
                   'chrxh_random': 1,
                   'chry_random': 1,
                   'chryh_random': 1,
                   }

# to do: add dpse: unknown_group_...

###################################################################


class Annotation:

    def __init__(self, identifier, short, description):
        self.mIdentifier = identifier
        self.mShort = short
        self.mDescription = description


def readAnnotationInterpro(infile):
    map_id2annotation = {}

    for line in infile:
        if line[0] == "#":
            continue
        data = line[:-1].split("\t")
        id, identifier, short, description = data[1], data[3], data[4], data[5]
        if short == "NULL":
            continue

        if id not in map_id2annotation:
            map_id2annotation[id] = []

        map_id2annotation[id].append(
            Annotation(identifier, short, description))

    return map_id2annotation

###################################################################


class Location:

    def __init__(self):
        pass

    def __str__(self):
        return "%s\t%s\t%i\t%i" % (self.contig, self.strand, self.mFrom, self.mTo)

###################################################################


class Duplication:

    """container for duplication information."""
    format_branch_length = "%6.4f"

    def __init__(self):

        self.mGeneTree = None
        self.mGeneNode = None
        self.mSpeciesNode = None
        self.mSpeciesTree = None
        self.mDistance2Root = None
        self.mDistance2Min = None
        self.mDistance2Max = None
        self.mRelativeHeight = None
        self.mTotalDistance = None
        self.mRelativeDistance = None
        self.mDistance2Parent = None
        self.mDistances2Children = None
        self.mNGenes = None
        self.mNPseudogene = None
        self.mNOthers = None
        self.mNMaxContigs = None
        self.mNSpecies = 0
        self.mNChildren = 0

    def __str__(self):

        return "\t".join((self.mGeneTree,
                          str(self.mGeneNode),
                          str(self.mSpeciesTree),
                          str(self.mSpeciesNode),
                          str(self.mNSpecies),
                          str(self.mNChildren),
                          self.format_branch_length % self.mDistance2Root,
                          self.format_branch_length % self.mDistance2Min,
                          self.format_branch_length % self.mDistance2Max,
                          self.format_branch_length % self.mRelativeHeight,
                          self.format_branch_length % self.mTotalDistance,
                          self.format_branch_length % self.mRelativeDistance,
                          self.format_branch_length % self.mDistance2Parent,
                          ",".join(
                              [self.format_branch_length % x for x in self.mDistance2Children]),
                          str(self.mNGenes),
                          str(self.mNPseudogenes),
                          str(self.mNOthers),
                          str(self.mNMaxContigs)))

    def readFromString(self, line):
        """fill data from string."""
        data = line.split("\t")

        gene_tree, gene_node, species_tree, species_node, \
            nspecies, nchildren, \
            d2root, d2min, d2max, \
            relative_height, total_distance, rel_distance, \
            distance_to_parent, distances_to_children, \
            ngenes, npseudogenes, nothers, max_ncontigs = data[:18]

        self.mGeneTree = gene_tree

        (self.mGeneNode, self.mSpeciesNode, self.mSpeciesTree,
         self.mNSpecies, self.mNChildren,
         self.mNGenes, self.mNPseudogenes, self.mNOthers, self.mNMaxContigs) = \
            map(int, (gene_node, species_node, species_tree,
                      nspecies, nchildren,
                      ngenes, npseudogenes, nothers, max_ncontigs))

        self.mDistance2Root, self.mDistance2Min, self.mDistance2Max, self.mRelativeHeight, self.mTotalDistance, self.mRelativeDistance, self.mDistance2Parent = \
            map(float, (d2root, d2min, d2max, relative_height,
                total_distance, rel_distance, distance_to_parent))

        self.mDistance2Children = map(float, distances_to_children.split(","))

    def getHeader(self):
        return "\t".join((
            "gene_tree", "gene_node", "species_tree", "species_node",
            "nspecies", "nchildren",
            "d2root", "d2min", "d2max",
            "height", "distance",
            "rel_distance", "parent_distance",
            "children_distances",
            "ngenes", "npseudogenes", "nothers", "max_contigs"))

###################################################################


class SVGRuler:

    mFontSize = 12
    mFont = "Verdana"
    mSeparator = 5
    mFormat = "%5.2f"
    mColour = SVGTree.BLACK
    mTickSize = 5

    def __init__(self):
        pass

    def getHeight(self):
        return self.mTickSize * 3

    def getElements(self, x1, y, x2,
                    min_value, max_value,
                    value_increment=None,
                    ruler_elements=("ruler", "right-ticks", "scale")):
        """add ruler at positions x1, y, x2.

        [y is top of ruler]
        """
        ee = []

        value_range = max_value - min_value

        if value_increment is None:
            # quick fix for kaks

            if value_range > 1:
                value_increment = 0.5
            elif value_range > 0.1:
                value_increment = 0.05
            else:
                value_increment = 0.01

        increment = float(value_increment) / value_range * (x2 - x1)

        if "ruler" in ruler_elements:
            # full length ruler with tick marks and labels
            e = SVGdraw.line(x1, y + self.mTickSize + 1, x2, y + self.mTickSize + 1,
                             stroke="rgb(%i,%i,%i)" % self.mColour,
                             stroke_width=1)
            ee.append(e)

        if "right-ticks" in ruler_elements:

            x = x2
            while x >= x1:
                e = SVGdraw.line(x, y, x, y + 2 * self.mTickSize + 1,
                                 stroke="rgb(%i,%i,%i)" % self.mColour,
                                 stroke_width=1)
                ee.append(e)
                x -= increment

            y += 2 * self.mTickSize + 1 + self.mSeparator

        if "left-ticks" in ruler_elements:

            x = x1
            while x <= x2:
                e = SVGdraw.line(x, y, x, y + 2 * self.mTickSize + 1,
                                 stroke="rgb(%i,%i,%i)" % self.mColour,
                                 stroke_width=1)
                ee.append(e)
                x += increment

            y += 2 * self.mTickSize + 1 + self.mSeparator

        if "scale" in ruler_elements:

            e = SVGdraw.line(x2,
                             y + self.mTickSize + 1,
                             x2 - increment,
                             y + self.mTickSize + 1,
                             stroke="rgb(%i,%i,%i)" % self.mColour,
                             stroke_width=1)

            ee.append(e)

            e = SVGdraw.line(x2,
                             y,
                             x2,
                             y + 2 * self.mTickSize + 1,
                             stroke="rgb(%i,%i,%i)" % self.mColour,
                             stroke_width=1)

            ee.append(e)

            e = SVGdraw.line(x2 - increment,
                             y,
                             x2 - increment,
                             y + 2 * self.mTickSize + 1,
                             stroke="rgb(%i,%i,%i)" % self.mColour,
                             stroke_width=1)

            ee.append(e)

            e = SVGdraw.text(x2 - increment / 2,
                             y + 2 * self.mTickSize + 1 + self.mFontSize,
                             self.mFormat % value_increment,
                             self.mFontSize,
                             self.mFont,
                             stroke="rgb(%i,%i,%i)" % self.mColour,
                             text_anchor="middle")

            ee.append(e)

        return ee


###################################################################
class ColourScale:

    """build and serve colour scale."""

    def __init__(self, start_colour, end_colour, steps=100):
        self.mPalette = []

        rstart, gstart, bstart = start_colour
        rend, gend, bend = end_colour

        rinc = float((rend - rstart)) / (steps - 1)
        ginc = float((gend - gstart)) / (steps - 1)
        binc = float((bend - bstart)) / (steps - 1)

        r, g, b = rstart, gstart, bstart
        for x in range(steps):

            self.mPalette.append((int(r), int(g), int(b)))
            r += rinc
            g += ginc
            b += binc

        self.mPalette.append((rend, gend, bend))

    def __getitem__(self, x):
        return self.mPalette[x]

###################################################################
# Branch decorator for species trees


class BranchDecoratorHorizontalDistributions(SVGTree.BranchDecoratorHorizontal):

    mBoxHeight = 20

    def __init__(self, tree, map_node2counts, *args, **kwargs):
        SVGTree.BranchDecoratorHorizontal.__init__(self, tree, *args, **kwargs)

        self.mMapNode2Counts = map_node2counts

        self.mMapType2Colour = SVG_MAP_TYPE2COLOUR

        # contents of top/bottom bar for branches to internal nodes
        self.mInternalContents = (
            ('Duplication',),
            ('DuplicationDeletion',),
            ('Inconsistency', 'InconsistentTranscripts'),
        )

        # contents of top/bottom bar for branches to external nodes
        self.mExternalContents = (
            ('DuplicationLineage', ),
            ('InconsistentTranscripts',),
        )

        self.mSpeciationContents = ('Speciation', 'SpeciationDeletion')

        # maximum distance
        self.mMaxDistance = None

    def getValues(self, node_id, contents):

        if node_id not in self.mMapNode2Counts:
            return []
        vals = []
        for c in contents:
            if c in self.mMapNode2Counts[node_id]:
                vals += self.mMapNode2Counts[node_id][c]

        return vals

    def getElements(self, node_id, x1, x2, y):

        e = []

        node = self.mTree.node(node_id)

        if node_id not in self.mMapNode2Counts:
            # no info at node
            return SVGTree.BranchDecoratorHorizontal.getElements(self, node_id, x1, x2, y)

        elif node.succ:
            # Internal nodes
            contents = self.mInternalContents
        else:
            # Terminal nodes
            contents = self.mExternalContents

        values = []
        colours = []
        # get speciation time point
        parent = node.prev
        v = self.getValues(parent, self.mSpeciationContents)
        if not v:
            return SVGTree.BranchDecoratorHorizontal.getElements(self, node_id, x1, x2, y)

        values.append(v)
        colours.append(self.mMapType2Colour[self.mSpeciationContents[0]])

        # get all values and colour scales
        for cc in contents:
            v = self.getValues(node_id, cc)
            if not v:
                continue
            values.append(v)
            colours.append(self.mMapType2Colour[cc[0]])

        min_value = min(map(lambda x: min(x), values))
        if self.mMaxDistance is None:
            max_value = max(map(lambda x: numpy.median(x), values)) * 4
        elif self.mMaxDistance == 0.0:
            max_value = min(
                max(map(lambda x: max(x), values)), self.mMaxDistance)
        else:
            max_value = min(max(map(lambda x: max(x), values)))

        yy = y - ((len(values) / 2) + 1) * self.mBoxHeight

        for x in range(len(values)):

            e += DistributionBox().getElements(x1, yy, x2, yy + self.mBoxHeight,
                                               values=values[x],
                                               min_value=min_value,
                                               max_value=max_value,
                                               colour_distn=colours[x])
            yy += self.mBoxHeight

        e += SVGRuler().getElements(x1, yy, x2,
                                    min_value=min_value,
                                    max_value=max_value)

        return e

    def getHeight(self, node_id):
        node = self.mTree.node(node_id)
        if node.succ:
            return self.mBoxHeight * (1 + len(self.mInternalContents)) * 3
        else:
            return self.mBoxHeight * (1 + len(self.mExternalContents)) * 3

###################################################################


class DistributionBox:

    mSeparator = 5

    mFontSize = 12
    mFont = "Verdana"
    mFontColour = SVGTree.BLACK

    def __init__(self):
        pass

    def getElements(self, x1, y1, x2, y2, values,
                    min_value=None,
                    max_value=None,
                    num_bins=100,
                    num_colours=100,
                    add_counts=True,
                    colour_distn=SVGTree.BLUE,
                    colour_median=SVGTree.RED):
        """plot a box with coordinates x,y,x2,y2.

        Colour box with gradient according to values.
        """

        self.mValues = values

        if min_value is None:
            min_value = min(values)
        if max_value is None:
            max_value = max(values)

        value_range = max_value - min_value
        width = x2 - x1
        height = y2 - y1

        e = []

        # plot distribution
        stats = Stats.DistributionalParameters(values)

        increment_values = float(value_range) / num_bins
        increment_coords = float(width) / num_bins
        bins = [min_value + increment_values *
                x for x in range(0, num_bins + 1)]
        hist = scipy.stats.histogram2(values, bins)

        max_height = max(hist)

        scale = ColourScale(SVGTree.WHITE, colour_distn, num_colours + 1)

        x = x2
        for b in range(num_bins):
            x -= increment_coords
            c = int((float(hist[b]) / max_height) * num_colours)
            e.append(SVGdraw.rect(int(x), y1,
                                  increment_coords, height,
                                  stroke=None,
                                  fill="rgb(%i,%i,%i)" % scale[c]))

        # plot median bar for speciation
        x = x2 - width * ((stats['median'] - min_value) / value_range)

        e.append(SVGdraw.line(x, y1, x, y2,
                              stroke="rgb(%i,%i,%i)" % SVGTree.RED,
                              stroke_width=1))

        # plot box
        ee = SVGdraw.rect(x1, y1, width, height,
                          stroke="rgb(%i,%i,%i)" % SVGTree.BLACK,
                          fill="rgb(%i,%i,%i)" % SVGTree.WHITE)

        ee.attributes['fill-opacity'] = 0
        e.append(ee)

        # plot counts
        if add_counts:
            ee = SVGdraw.text(x1 + width / 4, y1 + height,
                              "%i" % (len(values)),
                              self.mFontSize,
                              self.mFont,
                              stroke="rgb(%i,%i,%i)" % self.mFontColour,
                              fill="rgb(%i,%i,%i)" % self.mFontColour,
                              text_anchor="center")

            e.append(ee)

        return e

###################################################################


class NodeDecoratorDuplications(SVGTree.NodeDecorator):

    """class for decorating external nodes. Colour by duplication type.
    """

    mRadius = 6
    mFontSize = 12
    mFont = "Verdana"
    mSeparator = 5

    def __init__(self, tree, node_types, add_ids=False, add_legend=True, *args, **kwargs):
        SVGTree.NodeDecorator.__init__(self, tree, *args, **kwargs)
        self.mMapType2Colour = SVG_MAP_TYPE2COLOUR
        self.mNodeTypes = node_types
        self.mAddIds = add_ids
        self.mAddLegend = add_legend

    def getElements(self, node_id, x, y):
        t = self.mNodeTypes[node_id].mType
        if t not in self.mMapType2Colour:
            self.mMapType2Colour[t] = SVGTree.COLOURS[
                len(self.mMapType2Colour) % len(SVGTree.COLOURS)]

        e = []
        e.append(SVGdraw.circle(x, y,
                                self.mRadius,
                                stroke="rgb(%i,%i,%i)" % self.mMapType2Colour[
                                    t],
                                fill="rgb(%i,%i,%i)" % self.mMapType2Colour[t]))

        if self.mAddIds:
            e.append(SVGdraw.text(x + self.mRadius, y + self.mRadius,
                                  str(node_id),
                                  self.mFontSize,
                                  self.mFont,
                                  stroke="rgb(%i,%i,%i)" % (0, 0, 0)))

        return e

    def getHeight(self, node_id):
        return 0

    def getWidth(self, node_id):
        return 0

    def getLegend(self, x, y):
        """add legend to plot starting at x and y.
        """
        e = []
        if self.mAddLegend:
            keys = self.mMapType2Colour.keys()
            keys.sort()

            for k in keys:
                v = self.mMapType2Colour[k]

                e.append(SVGdraw.circle(x, y,
                                        self.mRadius,
                                        stroke="rgb(%i,%i,%i)" % v,
                                        fill="rgb(%i,%i,%i)" % v))

                e.append(SVGdraw.text(x + self.mRadius + self.mSeparator,
                                      y,
                                      k,
                                      self.mFontSize,
                                      self.mFont,
                                      stroke="rgb(%i,%i,%i)" % SVGTree.BLACK,
                                      text_anchor="left"))

                y += self.mFontSize + self.mSeparator

        return e, x, y

    def getFooterHeight(self):
        if self.mAddLegend:
            return len(self.mMapType2Colour) * (max(self.mFontSize, self.mRadius) + self.mSeparator)
        else:
            return 0

###################################################################


class NodeDecoratorCounts(SVGTree.NodeDecorator):

    """class for decorating internal/external nodes.

    Add counts to branch
    """

    mFontSize = 12
    mFont = "Verdana"
    mFontColour = SVGTree.BLACK

    mRadius = 6
    mSeparator = 5

    # maxinum size of bar in pixels
    mBarWidth = 100
    mBarHeight = 15

    def __init__(self, tree, map_node2counts,
                 mode="bars",
                 bar_scale_factor=0.0,
                 *args, **kwargs):
        """colours a ','-separated tuples of RGB values."""

        SVGTree.NodeDecorator.__init__(self, tree, *args, **kwargs)
        self.mMapNode2Counts = map_node2counts

        self.mMapType2Colour = {"Speciation": SVGTree.GREEN,
                                "SpeciationDeletion": SVGTree.OLIVE,
                                "Transcripts": SVGTree.BLUE,
                                "DuplicationLineage": SVGTree.GREY,
                                "Duplication": SVGTree.PLUM,
                                "DuplicationDeletion": SVGTree.ORANGE,
                                "DuplicationInconsistency": SVGTree.BROWN,
                                "Outparalogs": SVGTree.YELLOW,
                                "InconsistentTranscripts": SVGTree.RED,
                                "Inconsistency": SVGTree.PURPLE}

        # node types to write
        self.mNodeTypes = ("DuplicationLineage",
                           "Duplication",
                           "DuplicationDeletion",
                           "DuplicationInconsistency",
                           "InconsistentTranscripts",
                           "Inconsistency")

        if mode not in ("text", "bars"):
            raise "unknown mode %s" % mode

        self.mMode = mode

        self.mBarScaleFactor = bar_scale_factor

    def getElements(self, node_id, x, y):
        """write a list of counts for each node.
        """
        e = []

        if node_id in self.mMapNode2Counts:

            node_counts = self.mMapNode2Counts[node_id]

            y -= self.mBarHeight * len(self.mNodeTypes) / 2

            if self.mBarScaleFactor == 0.0:
                m = 0.0
                for node_type in self.mNodeTypes:
                    if node_type in node_counts:
                        m = max(m, len(node_counts[node_type]))
                # if maximum is zero, no values need to be output
                if m == 0:
                    return e
                scale_factor = float(self.mBarWidth) / m
            else:
                scale_factor = self.mBarScaleFactor

            for node_type in self.mNodeTypes:

                if node_type in node_counts:

                    vals = node_counts[node_type]

                    colour = self.mMapType2Colour[node_type]

                    if self.mMode == "text":

                        ee = SVGdraw.text(x, y,
                                          "%s-%i" % (node_type, len(vals)),
                                          self.mFontSize,
                                          self.mFont,
                                          stroke="rgb(%i,%i,%i)" % colour,
                                          fill="rgb(%i,%i,%i)" % colour,
                                          text_anchor="left")
                        e.append(ee)

                    elif self.mMode == "bars":

                        bar_width = int(len(vals) * scale_factor)

                        ee = SVGdraw.rect(x, y, bar_width, self.mBarHeight,
                                          stroke="rgb(%i,%i,%i)" % colour,
                                          fill="rgb(%i,%i,%i)" % colour)

                        e.append(ee)

                        ee = SVGdraw.text(x + bar_width + 1, y + self.mBarHeight,
                                          "%i" % (len(vals)),
                                          self.mFontSize,
                                          self.mFont,
                                          stroke="rgb(%i,%i,%i)" % self.mFontColour,
                                          fill="rgb(%i,%i,%i)" % self.mFontColour,
                                          text_anchor="left")

                        e.append(ee)

                y += self.mBarHeight

        return e

    def getHeight(self, node_id):
        return self.mFontSize

    def getWidth(self, node_id):
        return len(self.mTree.node(node_id).data.taxon) * self.mFontSize

    def getLegend(self, x, y):
        """add legend to plot starting at x and y."""

        e = []
        for k, v in self.mMapType2Colour.items():

            e.append(SVGdraw.circle(x, y,
                                    self.mRadius,
                                    stroke="rgb(%i,%i,%i)" % v,
                                    fill="rgb(%i,%i,%i)" % v))

            e.append(SVGdraw.text(x + self.mRadius + self.mSeparator,
                                  y,
                                  k,
                                  self.mFontSize,
                                  self.mFont,
                                  stroke="rgb(%i,%i,%i)" % SVGTree.BLACK,
                                  text_anchor="left"))

            y += self.mFontSize + self.mSeparator

        return e, x, y

    def getFooterHeight(self):
        return len(self.mMapType2Colour) * (max(self.mFontSize, self.mRadius) + self.mSeparator)

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def getInfile(filename):
    if filename == "-":
        return sys.stdin
    else:
        return open(filename, "r")

opened_files = {}

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def getFile(options, section):
    """returns tuple with file and flag of whether file is new.
    """

    global opened_files

    if options.output_pattern:

        if section in opened_files:
            return opened_files[section], False
        else:
            if "%s" in options.output_pattern:
                fn = options.output_pattern % section
            else:
                fn = options.output_pattern

            if options.loglevel >= 1:
                options.stdout.write(
                    "# output for section %s goes to %s\n" % (section, fn))

            if os.path.exists(fn):
                outfile, rc = open(fn, "a"), False
            else:
                outfile, rc = open(fn, "w"), True

            opened_files[section] = outfile
            return outfile, rc
    else:
        if section in opened_files:
            return options.stdout, False
        else:
            opened_files[section] = 1
            return options.stdout, True


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def printSummary(counts_categories, counts_trees, counts_nodes, duplications,
                 reference_height,
                 gene_tree, gene_tree_species_list,
                 species_nexus, species_list,
                 options, prefix_header, prefix_row):
    """print summary information.


    * nspecies
    * species distribution
    * Number of out-paralogs
    * Number of species with lineage specific distributions
    * Number of inconsistencies
    * Number of deletions
    * 1:1 orthologs?
    """
    # counts per tree and type
    outfile, is_new = getFile(options, "summary")

    if is_new:
        outfile.write("%s%s\n" % (prefix_header,
                                  "\t".join(("notus",
                                             "nspecies",
                                             "has_outgroup",
                                             "is_full",
                                             "ninconsistencies",
                                             "refheight",
                                             "ndup_total",
                                             "ndup_internal",
                                             "\t".join(map(lambda x: "ndup_%s" % x, species_list))))))

    # count lineage and non-lineage specific duplication events
    lineage_specific_duplications = {}
    for s in species_list:
        lineage_specific_duplications[s] = 0
    ninternal_duplications = 0
    nduplications = 0
    for duplication in duplications:
        species_tree = species_nexus.trees[duplication.mSpeciesTree]
        node = species_tree.node(duplication.mSpeciesNode)
        if node.succ == []:
            lineage_specific_duplications[node.data.taxon] += 1
        else:
            ninternal_duplications += 1
        nduplications += 1

    nspecies = len(gene_tree_species_list)

    gt_set = set(gene_tree_species_list)
    st_set = set(species_list)
    if options.outgroup_species:
        ou_set = set(options.outgroup_species)
    else:
        ou_set = set()
    st_set_no_ou = st_set.difference(ou_set)

    def yes_no(x):
        if x:
            return "1"
        else:
            return "0"

    has_outgroup = len(gt_set.intersection(ou_set)) > 0
    is_full = len(st_set_no_ou.difference(gt_set)) == 0

    notus = len(gene_tree.get_terminals())

    nproblems = sum([len(counts_categories[x])
                    for x in options.nodetypes_inconsistency])

    outfile.write("%s%s\n" % (prefix_row,
                              "\t".join((
                                  str(notus),
                                  str(nspecies),
                                  yes_no(has_outgroup),
                                  yes_no(is_full),
                                  str(nproblems),
                                  options.format_branch_length % reference_height,
                                  str(nduplications),
                                  str(ninternal_duplications),
                                  "\t".join(map(str, [lineage_specific_duplications[s] for s in species_list]))))))

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def printResolution(resolution, options, prefix_header, prefix_row):
    """print resolution data."""

    outfile, is_new = getFile(options, "resolution")
    outfile.write("prefix\tgene_node\ttype\ttree\td_np\td_nl\td_nr\n")
    for cluster, gene_node_id, best_type, best_tree, d_np, d_nl, d_nr in resolution:
        outfile.write("\t".join((cluster,
                                 str(gene_node_id),
                                 best_type,
                                 str(best_tree),
                                 options.format_branch_length % d_np,
                                 options.format_branch_length % d_nl,
                                 options.format_branch_length % d_nr)) + "\n")
    if outfile != options.stdout:
        outfile.close()

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def printCounts(counts_categories, counts_trees, counts_nodes,
                species_nexus,
                options, prefix_header, prefix_row):
    """print counts for each section."""

    # counts per tree and type
    outfile, is_new = getFile(options, "trees")
    if is_new:
        outfile.write("%stree\ttotal\t%s\n" %
                      (prefix_header, "\t".join(sorted(options.priority.keys()))))

    for k, counts in enumerate(counts_trees):
        values = [counts[x] for x in sorted(options.priority.keys())]
        outfile.write("%s%s\t%i\t%s\n" %
                      (prefix_row, k, sum(values), "\t".join(map(str, values))))

    # counts per category
    outfile, is_new = getFile(options, "categories")
    if is_new:
        outfile.write("%stype\tcounts\theights\n" % prefix_header)
    for k, c in counts_categories.items():
        outfile.write("%s%s\t%i\t%s\n" % (prefix_row, k, len(c), ",".join(
            map(lambda x: options.format_branch_length % x, c))))

    # counts per node and tree
    outfile, is_new = getFile(options, "nodes")
    if is_new:
        outfile.write("%stree\tnode\ttype\t%s\theights\ttaxa\n" % (
            prefix_header, "\t".join(Stats.DistributionalParameters().getHeaders())))

    for tree in range(len(species_nexus.trees)):
        for node, types in counts_nodes[tree].items():
            for type, heights in types.items():
                s = Stats.DistributionalParameters(heights)
                outfile.write("%s%i\t%i\t%s\t%s\t%s\t%s\n" % (prefix_row,
                                                              tree,
                                                              node,
                                                              type,
                                                              str(s),
                                                              ",".join(
                                                                  map(lambda x: options.format_branch_length % x, heights)),
                                                              ",".join(species_nexus.trees[tree].get_taxa(node))))


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
def readCountsNodes(infile, nspecies_trees):
    """print counts for each section."""

    counts_nodes = [{} for x in range(nspecies_trees)]

    take = None
    for line in infile:
        if line[0] == "#":
            continue

        data = line[:-1].split("\t")

        if not take:
            if data[0] == "prefix":
                if data[1] == "total":
                    take = range(2, len(data))
                else:
                    take = range(1, len(data))
            else:
                raise "can't parse header %s" % data[0]

        else:
            tree = int(data[take[0]])
            node = int(data[take[1]])
            node_type = data[take[2]]
            heights = data[take[-2]]
            taxa = data[take[-1]]

            a = counts_nodes[tree]
            if node not in a:
                a[node] = {}
            a = a[node]
            a[node_type] = map(float, heights.split(","))

    return counts_nodes
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def printBestNodes(best_nodes, options, prefix_header, prefix_row):
    """print section with best nodes.
    """
    outfile, is_new = getFile(options, "details")
    if is_new:
        outfile.write(
            "%sgenenode\tspeciesnode\ttree\ttype\theight\tspecies\totus\n" % prefix_header)

    for best in best_nodes:

        if best is None:
            continue
        if best.mType == "Leaf":
            continue

        best_tree = best.mSpeciesTree

        # Output best choice
        species_tree = species_nexus.trees[best_tree]

        outfile.write("%s%i\t%i\t%i\t%s\t%s\t" % (
            prefix_row,
            best.mGeneNode, best.mSpeciesNode, best_tree, best.mType,
            options.format_branch_length % ((min_branch_lengths[best.mGeneNode] + max_branch_lengths[best.mGeneNode]) / 2)))

        # output species joined (node: species tree differs)
        # for s in species_tree.node(best.mSpeciesNode).succ:
        outfile.write(
            "\t" + ",".join(species_nexus.trees[best_tree].get_taxa(best.mSpeciesNode)))

        # output species for the genes joined
        # for s in .succ:
        outfile.write("\t" + ",".join(gene_tree.get_taxa(best.mGeneNode)))
        outfile.write("\n")

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def printClusters(clusters, options, prefix_prefix, otus):
    """print clusters from subtree calculation.
    """
    if options.subtrees_trees:
        outfile_subtrees, is_new = getFile(options, "subtrees")
    else:
        outfile_subtrees = None

    if options.subtrees_identifiers:
        outfile_subids, is_new = getFile(options, "subids")
    else:
        outfile_subids = None

    nclustered_taxa = sum(map(lambda x: len(x), clusters))
    written_subtrees = 0
    nsubtrees = 0

    for cluster in clusters:
        nsubtrees += 1
        if outfile_subids:
            for otu in cluster:
                outfile_subids.write("%s\t%s%i\n" %
                                     (otu, prefix_prefix, nsubtrees))
        if outfile_subtrees:
            subtree = TreeTools.GetSubtree(gene_tree, gene_tree.root)
            TreeTools.PruneTree(subtree, cluster)
            outfile_subtrees.write(">%s%i\n%s\n" %
                                   (prefix_prefix, nsubtrees, TreeTools.Tree2Newick(subtree)))
            written_subtrees += 1

    if options.loglevel >= 1:
        options.stdlog.write("# written: %i/%i subtrees (%5.2f%%) %i/%i (%5.2f%%) taxa\n" % (written_subtrees, len(clusters),
                                                                                             100.0 * written_subtrees /
                                                                                             len(clusters),
                                                                                             nclustered_taxa, len(
                                                                                                 otus),
                                                                                             100.0 * nclustered_taxa / len(otus)))

    if outfile_subtrees:
        outfile_subtrees.close()
    if outfile_subids:
        outfile_subids.close()

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def printSVGSpeciesTrees(trees, all_counts_nodes, options):
    """svg output for a species tree.

    The tree is annotated with node information.
    Branch lengths are set to median speciation events, if available, otherwise
    they are set to the median of all values.
    """
    for ntree in range(len(species_nexus.trees)):

        tree = trees[ntree]

        counts_nodes = total_counts_nodes[ntree]

        # set branch lengths in species tree
        if options.svg_branch_lengths == "contemporary":
            # set branch lengths such that all tips are contemporaneous
            pass
        elif options.svg_branch_lengths == "uniform":
            # set branch lengths to same value
            for id, node in tree.chain.items():
                node.data.branchlength = 1.0
        elif options.svg_branch_lengths == "median":
            # set branch lengths to median branch length of speciation events
            for id, node in tree.chain.items():
                if id in counts_nodes and "Speciation" in counts_nodes[id]:
                    v = numpy.median(counts_nodes[id]["Speciation"])
                    for s in node.succ:
                        tree.node(s).data.branchlength = v

        plot = SVGTree.SVGTree(tree)

        d = NodeDecoratorCounts(tree, counts_nodes)
        plot.setDecoratorInternalNodes(d)
        plot.setDecoratorHorizontalBranches(
            BranchDecoratorHorizontalDistributions(tree, counts_nodes))

        plot.initializePlot()

        if options.output_pattern_svg:
            if "%s" in options.output_pattern_svg:
                fn = re.sub("%s", tree.name, options.output_pattern_svg)
            else:
                fn = options.output_pattern_svg

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# svg output for tree %s goes to %s.\n" % (tree.name, fn))

            outfile = open(fn, "w")
        else:
            outfile = options.stdout

        plot.writeToFile(outfile)

        if outfile != options.stdout:
            outfile.close()


# ------------------------------------------------------------------------
# save previously loaded data
map_species2url = None
map_species2colour = None
map_taxon2url = None
map_species2name = None

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def printSVGGeneTree(tree, options, extract_species):
    """print svg gene tree.
    """
    global map_species2url
    global map_species2colour
    global map_taxon2url
    global map_species2name

    if options.filename_species2url:
        if not map_species2url:
            map_species2url = IOTools.ReadMap(
                open(options.filename_species2url, "r"))

            class MapperTaxon:

                def __init__(self, map_species2url, separator):
                    self.mMapSpecies2URL = map_species2url
                    self.mSeparator = separator

                def __call__(self, taxon):
                    data = taxon.split(self.mSeparator)
                    if len(data) > 3:
                        species, prediction_id, gene_id = data[:3]
                    elif len(data) == 2:
                        species, gene_id = data
                    else:
                        raise "can't understand taxon id %s" % taxon

                    if species in self.mMapSpecies2URL:
                        params = {'species': species,
                                  'gene': gene_id}
                        return self.mMapSpecies2URL[species] % params
                    else:
                        return None

            map_taxon2url = MapperTaxon(map_species2url, options.separator)

    else:
        map_species2url = None
        map_taxon2url = None

    if options.filename_species2colour:
        if not map_species2colour:
            map_species2colour = IOTools.ReadMap(
                open(options.filename_species2colour, "r"))
    else:
        map_species2colour = None

    if options.filename_species2name:
        if not map_species2name:
            map_species2name = IOTools.ReadMap(
                open(options.filename_species2name, "r"))
    else:
        map_species2name = None

    plot = SVGTree.SVGTree(tree)

    node_types = {}
    for best in best_nodes:

        if best is None:
            continue
        if best.mType == "Leaf":
            continue

        node_types[best.mGeneNode] = best

    plot.setDecoratorExternalNodes(SVGTree.NodeDecoratorBySpecies(tree, extract_species=extract_species,
                                                                  map_species2colour=map_species2colour,
                                                                  map_species2name=map_species2name,
                                                                  map_taxon2url=map_taxon2url))

    plot.setDecoratorInternalNodes(NodeDecoratorDuplications(tree, node_types,
                                                             add_ids=options.svg_add_ids,
                                                             add_legend=options.svg_add_legend))

    if clusters:
        map_id2cluster = {}
        for c in range(len(clusters)):
            for id in clusters[c]:
                map_id2cluster[id] = c
        plot.addBorderDecorator(SVGTree.BorderDecoratorClusters(tree,
                                                                map_id2cluster))

    if options.filename_locations:
        plot.addBorderDecorator(SVGTree.BorderDecoratorLocations(tree,
                                                                 map_id2location,
                                                                 max_separation=options.max_separation,
                                                                 extract_species=extract_species,
                                                                 map_species2colour=map_species2colour,
                                                                 print_location=options.svg_print_location))

    plot.initializePlot()

    if options.output_format == "table":

        # write output in tabular format with tree name, tree and svg
        svg_file = cStringIO.StringIO()
        plot.writeToFile(svg_file)

        options.stdout.write("%s\t%s\t%s\n" % (tree.name,
                                               TreeTools.Tree2Newick(tree),
                                               re.sub("[\t\n]", " ", svg_file.getvalue())))

    else:
        if options.output_pattern_svg:
            if "%s" in options.output_pattern_svg:
                fn = re.sub("%s", tree.name, options.output_pattern_svg)
            else:
                fn = options.output_pattern_svg

            if options.loglevel >= 1:
                options.stdlog.write(
                    "# svg output for tree %s goes to %s.\n" % (tree.name, fn))

            outfile = open(fn, "w")
        else:
            outfile = options.stdout

        plot.writeToFile(outfile)

        if outfile != options.stdout:
            outfile.close()


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
def getAnalysisSets(gene_tree, extract_quality, options):
    """return sets of otus in gene tree according to filtering critera.

    returns a tuple of analysis_set, gene_set, pseudogene_set, other_set
    """
    # get list of nodes to analyze, classify by quality
    genes = []
    pseudogenes = []
    othergenes = []
    for node_id in gene_tree.get_terminals():
        node = gene_tree.node(node_id)
        try:
            quality = extract_quality(node.data.taxon)
        except IndexError:
            # default quality is gene
            quality = "CG"

        if quality in MAP_FILTER2QUALITY['genes']:
            genes.append(node_id)
        elif quality in MAP_FILTER2QUALITY['pseudogenes']:
            pseudogenes.append(node_id)
        else:
            othergenes.append(node_id)

    # define analysis set
    if options.filter_quality == "genes":
        analysis_set = genes
    elif options.filter_quality == "pseudogenes":
        analysis_set = pseudogenes
    elif options.filter_quality == "all":
        analysis_set = genes + pseudogenes + othergenes
    else:
        raise "unknown filter: %s" % options.filter_quality

    return analysis_set, genes, pseudogenes, othergenes


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def analyzeDuplicationData(duplications, species_nexus, options, extract_gene=None, extract_species=None,
                           gene_nexus=None, map_id2location=None):
    """analyze duplication data from infile and perform analysis.
    """

    def getDuplicationCounts(duplications):
        counts_per_node = {}

        mode = options.filter_quality

        for duplication in duplications:

            # ignore non-default tree duplications
            if duplication.mSpeciesTree != 0:
                continue

            if mode == "genes":
                if duplication.mNPseudogenes or duplication.mNOthers:
                    continue
            elif mode == "pseudogenes":
                if duplication.mGenes or duplication.mNOthers:
                    continue

            try:
                counts_per_node[duplication.mSpeciesNode].append(
                    duplication.mRelativeDistance)
            except KeyError:
                counts_per_node[duplication.mSpeciesNode] = [
                    duplication.mRelativeDistance]

        return counts_per_node

    # ------------------------------------------------------------------------
    def printDuplicationHistograms(outfile, counts, species_tree, options, cumulate=False, normalize=False):

        histograms = []

        headers = []

        values_all, values_lineage, values_internal = [], [], []

        for t in sorted(counts.keys()):

            h = Histogram.Calculate(counts[t],
                                    increment=options.hist_dups_bin_size,
                                    min_value=options.hist_dups_min_value,
                                    max_value=options.hist_dups_max_value,
                                    no_empty_bins=True)

            if cumulate:
                Histogram.cumulate(h)
            if normalize:
                Histogram.normalize(h)

            histograms.append(h)

            taxa = species_tree.get_taxa(t)
            headers.append("-".join(map(lambda x: x[:x.find("_")], taxa)))

            if len(taxa) == 1:
                values_lineage += counts[t]
            else:
                values_internal += counts[t]

            values_all += counts[t]

        for title, values in (("all", values_all), ("lineage", values_lineage), ("internal", values_internal)):
            h = Histogram.Calculate(values,
                                    increment=options.hist_dups_bin_size,
                                    min_value=options.hist_dups_min_value,
                                    max_value=options.hist_dups_max_value,
                                    no_empty_bins=True)
            if cumulate:
                Histogram.cumulate(h)
            if normalize:
                Histogram.normalize(h)

            histograms.append(h)
            headers.append(title)

        hh = Histogram.Combine(histograms, missing_value="na")

        outfile.write("%s\t%s\n" % ("bin", "\t".join(headers)))
        Histogram.Write(outfile, hh, format_value=options.format_branch_length)

    # ------------------------------------------------------------------------
    def printDuplicationStats(outfile, counts, species_tree, options):
        """print summary statistics for duplication counts."""

        headers = []

        outfile.write("node\t%s\n" %
                      ("\t".join(Stats.DistributionalParameters().getHeaders())))

        for t in sorted(counts.keys()):

            s = Stats.DistributionalParameters(counts[t])
            s.setFormat(options.format_branch_length)

            outfile.write("%s\t%s\n" %
                          ("-".join(species_tree.get_taxa(t)), str(s)))

    # ------------------------------------------------------------------------
    def printDuplicationsMembers(outfile, duplications, species_tree, options, extract_gene, extract_species, gene_nexus, map_id2location):
        """print list of Bmembers involved in duplications.
        """

        map_name2tree = {}
        for t in gene_nexus.trees:
            map_name2tree[t.name] = t

        nskipped = 0
        mode = options.filter_quality

        # group duplications trees and subgroups
        for duplication in duplications:

            # ignore non-default tree duplications
            if duplication.mSpeciesTree != 0:
                continue

            if mode == "genes":
                if duplication.mNPseudogenes or duplication.mNOthers:
                    continue
            elif mode == "pseudogenes":
                if duplication.mGenes or duplication.mNOthers:
                    continue

            try:
                gene_tree = map_name2tree[duplication.mGeneTree]
            except KeyError:
                nskipped += 1
                continue

            taxa = species_tree.get_taxa(duplication.mSpeciesNode)
            otus = gene_tree.get_taxa(duplication.mGeneNode)

            outfile.write(">%s\n" % duplication.mGeneTree)
            for otu in otus:
                outfile.write("%s\t%s\t%s\t%s\n" % (extract_species(otu),
                                                    extract_gene(otu),
                                                    otu,
                                                    map_id2location[otu]))

    # ------------------------------------------------------------------------
    def printDuplicationsLocal(outfile, duplications, species_tree, options, extract_gene, extract_species, gene_nexus, map_id2location):
        """print local duplications. These duplications are lineage-specific and
        not adjacent without any intermittent genes.
        """

        map_name2tree = {}
        for t in gene_nexus.trees:
            map_name2tree[t.name] = t

        nskipped = 0
        mode = options.filter_quality

        # container for children that should not be output
        visited = {}

        # group duplications trees and subgroups
        for duplication in duplications:

            # ignore non-default tree duplications
            if duplication.mSpeciesTree != 0:
                continue

            if mode == "genes":
                if duplication.mNPseudogenes or duplication.mNOthers:
                    continue
            elif mode == "pseudogenes":
                if duplication.mGenes or duplication.mNOthers:
                    continue

            try:
                gene_tree = map_name2tree[duplication.mGeneTree]
            except KeyError:
                nskipped += 1
                continue

            # only take lineage specific duplications
            if duplication.mNSpecies > 1:
                continue

            # only print out nodes from highest part of tree
            if duplication.mGeneTree in visited and duplication.mGeneNode in visited[duplication.mGeneTree]:
                continue

            # print duplications only from the highest node,
            # Note: this only works, if duplications are sorted
            # from high up in the tree to down in the tree
            if duplication.mGeneTree not in visited:
                visited[duplication.mGeneTree] = set()

            visited[duplication.mGeneTree].add(duplication.mGeneNode)
            for c in gene_tree._walk(duplication.mGeneNode):
                assert(c not in visited[duplication.mGeneTree])
                visited[duplication.mGeneTree].add(c)

            taxa = species_tree.get_taxa(duplication.mSpeciesNode)
            otus = gene_tree.get_taxa(duplication.mGeneNode)

            otus.sort(
                lambda x, y: cmp(map_id2location[x].mId, map_id2location[y].mId))

            def printOTU(last_otu):
                try:
                    outfile.write("%s\t%s\t%s\t%i\t%s\n" % (extract_species(last_otu),
                                                            extract_gene(
                                                                last_otu),
                                                            last_otu,
                                                            map_id2location[
                                                                last_otu].mId,
                                                            map_id2location[last_otu]))
                except AttributeError:
                    options.stderr.write("# unable to parse: %s\n" % last_otu)
                    options.stderr.flush()

            last_otu = otus[0]
            last_loc = map_id2location[last_otu]
            first = True

            for otu in otus[1:]:

                loc = map_id2location[otu]

                is_local = abs(loc.mId - last_loc.mId) <= 1

                if is_local:
                    if first:
                        outfile.write(">%s\n" % duplication.mGeneTree)
                        printOTU(last_otu)
                    printOTU(otu)
                    first = False
                else:
                    first = True

                last_otu = otu
                last_loc = loc

    # ------------------------------------------------------------------------
    def printDuplicationsSVG(duplications, species_tree, options, extract_gene, extract_species, gene_nexus, map_id2location):
        """print wheel plots for duplications.
        """

        # setup all species
        species_list = species_tree.get_taxa()

        if not options.filename_contig_sizes:
            raise "need location to filename with contig sizes."

        mode = options.filter_quality

        if options.filename_species2url:
            map_species2url = IOTools.ReadMap(
                open(options.filename_species2url, "r"))
        else:
            map_species2url = None

        map_name2tree = {}
        for t in gene_nexus.trees:
            map_name2tree[t.name] = t

        ###################################################################
        ###################################################################
        ###################################################################
        # go through duplications and collect contigs for
        # each species and subtrees
        ###################################################################
        data = {}
        contigs = {}
        for x in species_list:
            data[x] = []
            contigs[x] = set()

        ninput, nskipped, nnested, noutput = 0, 0, 0, 0

        # sort by otus, so that largest trees come first for a particular gene tree
        # smaller subtrees can then be ignored.
        duplications.sort(
            lambda x, y: cmp((x.mGeneTree, x.mNChildren), (y.mGeneTree, y.mNChildren)))
        duplications.reverse()

        last_tree = None
        otu_set = set()

        for duplication in duplications:

            ninput += 1

            if options.loglevel >= 3:
                options.stdlog.write(
                    "# processing duplication %s\n" % (str(duplication)))

            # ignore non-default tree duplications
            if duplication.mSpeciesTree != 0:
                nskipped += 1
                continue

            if mode == "genes":
                if duplication.mNPseudogenes or duplication.mNOthers:
                    continue
            elif mode == "pseudogenes":
                if duplication.mGenes or duplication.mNOthers:
                    continue

            try:
                gene_tree = map_name2tree[duplication.mGeneTree]
            except KeyError:
                nskipped += 1
                continue

            # write only the largest tree for a set of nested duplications
            if last_tree != duplication.mGeneTree:
                otu_set = set()
            last_tree = duplication.mGeneTree

            otus = gene_tree.get_taxa(duplication.mGeneNode)

            if len(otu_set.intersection(set(otus))) > 0:
                nnested += 1
                continue

            otu_set = otu_set.union(set(otus))

            is_lineage = duplication.mNSpecies == 1

            if not is_lineage:
                continue

            for species in set(map(extract_species, otus)):

                tree = TreeTools.GetSubtree(gene_tree, duplication.mGeneNode)
                this_otus = filter(
                    lambda x: extract_species(x) == species, otus)

                c = set()
                min_pos, max_pos = 0, 0
                gene_set = set()

                new_otus = []
                for otu in this_otus:
                    location = map_id2location[otu]
                    gene = extract_gene(otu)

                    if gene in gene_set:
                        continue
                    gene_set.add(gene)
                    new_otus.append(otu)
                    contigs[species].add(location.contig)

                    if max_pos == 0:
                        min_mos, max_pos = location.mFrom, location.mTo
                    else:
                        min_pos = min(location.mFrom, min_pos)
                        max_pos = max(location.mTo, max_pos)

                    c.add(location.contig)

                TreeTools.PruneTree(tree, new_otus)

                data[species].append(
                    (len(c), len(new_otus), min_pos, max_pos, new_otus, tree))

            noutput += 1

        if options.loglevel >= 1:
            options.stdout.write("# svg wheel plots: ninput=%i, nskipped=%i, nnested=%i, noutput=%i\n" % (
                ninput, nskipped, nnested, noutput))

        ###################################################################
        ###################################################################
        ###################################################################
        # create plots, add duplications and write them.
        ###################################################################
        for species in species_list:

            if len(data[species]) == 0:
                continue

            filename = re.sub("%s", species, options.filename_contig_sizes)
            map_contig2size = IOTools.ReadMap(open(filename, "r"),
                                              map_functions=(str, int))

            #########################################################
            # prepare contig list
            # for now: only those with duplications
            contig_list = list(contigs[species])

            # create plot
            plot = SVGDuplicationsWheel.DuplicationPlot(contig_list,
                                                        map_contig2size,
                                                        num_entries=0,
                                                        template=options.svg_template)

            plot.mRadiusIncrement = options.svg_wheel_plot_radius_increment
            plot.mRadius = options.svg_wheel_plot_radius
            plot.mMaxValue = options.svg_wheel_plot_max_value
            plot.mMinValue = options.svg_wheel_plot_min_value

            plot.initializePlot()

            #########################################################
            if map_species2url:
                url = map_species2url[species]
            else:
                url = None

            #########################################################
            # add duplications
            last_ndups = 0

            data[species].sort()

            for cis, ndups, min_pos, max_pos, otus, tree in data[species]:

                if ndups != last_ndups:
                    plot.pushRadius()
                    plot.addSeparator()

                last_ndups = ndups

                map_gene2location = {}

                for otu in otus:
                    location = map_id2location[otu]
                    if location.contig not in map_contig2size:
                        continue
                    gene = extract_gene(otu)
                    map_gene2location[gene] = (location.contig,
                                               location.strand,
                                               location.mFrom,
                                               location.mTo)

                if not map_gene2location:
                    continue

                # the last subset is all nodes
                s = TreeTools.GetSubsets(tree)

                is_first = True
                for children, height, branchlength in s:
                    if len(children) == 1:
                        continue
                    if options.loglevel >= 5:
                        options.stdlog.write(
                            "# species %s - %s: adding duplication to plot:\n" % (species, tree.name))
                        for child in children:
                            options.stdlog.write(
                                "#    %s\t%s\n" % (child, str(map_id2location[child])))

                    c = map(lambda x: x.split(options.separator), children)
                    plot.addDuplication(c, map_gene2location, height,
                                        url=url,
                                        with_separator=is_first,
                                        link_to_previous=not is_first,
                                        quality2symbol=options.svg_wheel_plot_quality2symbol,
                                        quality2mask=options.svg_wheel_plot_quality2mask)
                    is_first = False

            if options.output_pattern_svg:
                if "%s" in options.output_pattern_svg:
                    fn = re.sub("%s", species, options.output_pattern_svg)
                else:
                    fn = options.output_pattern_svg

                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# svg output for wheel %s goes to %s.\n" % (tree.name, fn))
                    options.stdlog.flush()
                outfile = open(fn, "w")
            else:
                outfile = options.stdout

            plot.writeToFile(outfile)

            if outfile != options.stdout:
                outfile.close()

    # ------------------------------------------------------------------------
    def annotateDuplications(outfile, duplications, species_tree, options, extract_gene, gene_nexus, map_id2location):
        """annotate duplications.

        ## NB: specific to fly project - dmel specific duplications
        """
        map_id2annotations = readAnnotationInterpro(
            open(options.filename_map_annotation_interpro, "r"))

        mode = options.filter_quality

        map_name2tree = {}
        for t in gene_nexus.trees:
            map_name2tree[t.name] = t

        dups_per_tree = {}
        nskipped = 0

        # group duplications trees and subgroups
        for duplication in duplications:

            # ignore non-default tree duplications
            if duplication.mSpeciesTree != 0:
                continue

            if mode == "genes":
                if duplication.mNPseudogenes or duplication.mNOthers:
                    continue
            elif mode == "pseudogenes":
                if duplication.mGenes or duplication.mNOthers:
                    continue

            taxa = species_tree.get_taxa(duplication.mSpeciesNode)

            if len(taxa) > 1:
                continue
            if taxa[0] != "dmel_vs_dmel4":
                continue

            try:
                gene_tree = map_name2tree[duplication.mGeneTree]
            except KeyError:
                nskipped += 1
                continue

            try:
                dups_per_tree[gene_tree.name].append(duplication)
            except KeyError:
                dups_per_tree[gene_tree.name] = [duplication]

        if options.loglevel >= 1:
            options.stdlog.write(
                "# analysing duplications in %i clusters - %i skipped.\n" % (len(dups_per_tree), nskipped))

        options.stdout.write("\t".join(("gene_tree",
                                        "dist2min",
                                        "dist2max",
                                        "notus",
                                        "otus",
                                        "locations",
                                        "nannotations",
                                        "identifiers",
                                        "descriptions")) + "\n")

        for gene_tree_name, dups in dups_per_tree.items():

            if options.loglevel >= 2:
                options.stdlog.write(
                    "# %s: %i duplications\n" % (gene_tree.name, len(dups)))

            gene_tree = map_name2tree[gene_tree_name]

            # find node with maximum number of children
            dups.sort(lambda x, y: cmp(x.mNChildren, y.mNChildren))

            duplication = dups[-1]

            taxa = gene_tree.get_taxa(duplication.mGeneNode)

            annotations = {}
            locations = []
            found_genes = set()
            genes = []
            for t in taxa:
                gene = extract_gene(t)

                # only one transcript per gene
                if gene in found_genes:
                    continue

                genes.append(gene)
                if map_id2location:
                    l = map_id2location[t]
                    locations.append(
                        ":".join(map(str, (l.contig, l.strand, l.mFrom, l.mTo))))

                found_genes.add(gene)
                try:
                    aa = map_id2annotations[gene]
                except KeyError:
                    continue

                for a in aa:
                    annotations[a.mIdentifier] = a

            # remove redundant identifiers
            identifiers = annotations.keys()
            identifiers.sort()
            descriptions = []
            for id in identifiers:
                descriptions.append(annotations[id].mDescription)

            options.stdout.write("\t".join((gene_tree.name,
                                            options.format_branch_length % duplication.mDistance2Min,
                                            options.format_branch_length % duplication.mDistance2Max,
                                            str(len(genes)),
                                            ";".join(genes),
                                            ";".join(locations),
                                            str(len(identifiers)),
                                            ";".join(identifiers),
                                            ";".join(descriptions))) + "\n")

    # get distribution of minimum/maximum separation per node type
    counts_per_node = getDuplicationCounts(duplications)

    for method in options.analyze_duplication_data:

        section = "duplications_%s" % (method)

        outfile, is_new = getFile(options, section)

        if method == "cumul-norm-hist":
            printDuplicationHistograms(outfile, counts_per_node, species_nexus.trees[0],
                                       options,
                                       cumulate=True, normalize=True)

        elif method == "cumul-hist":
            printDuplicationHistograms(outfile, counts_per_node, species_nexus.trees[0],
                                       options,
                                       cumulate=True)

        elif method == "stats":
            printDuplicationStats(outfile, counts_per_node, species_nexus.trees[0],
                                  options)

        elif method == "annotate":
            annotateDuplications(outfile, duplications, species_nexus.trees[
                                 0], options, extract_gene, gene_nexus, map_id2location)

        elif method == "members":
            printDuplicationsMembers(outfile, duplications, species_nexus.trees[
                                     0], options, extract_gene, extract_species, gene_nexus, map_id2location)

        elif method == "svg-wheel-plot":
            printDuplicationsSVG(duplications, species_nexus.trees[
                                 0], options, extract_gene, extract_species, gene_nexus, map_id2location)

        elif method == "local-dups":
            printDuplicationsLocal(outfile, duplications, species_nexus.trees[0], options, extract_gene, extract_species, gene_nexus,
                                   map_id2location)

        if outfile != options.stdout:
            outfile.close()

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def analyzeResolutionData(infile, options):
    """read resolution data from infile and perform analysis (counts/histograms).
    """

    def getCounts(data, field, f):
        counts = {}

        for cluster, gene_node_id, type, tree, d_np, d_nl, d_nr in data:
            if field == 0:
                a = type
            else:
                a = tree

            try:
                counts[a].append(f(d_np, d_nl, d_nr))
            except KeyError:
                counts[a] = [f(d_np, d_nl, d_nr)]

        return counts

    # ------------------------------------------------------------------------
    def printResolutionCounts(outfile, data, field, f, options, prefix_header="", prefix_row=""):

        counts = getCounts(data, field, f)

        outfile.write("%s%s\t%s\n" % (prefix_header,  "type", "\t".join(
            Stats.DistributionalParameters().getHeaders())))
        for t in sorted(counts.keys()):

            s = Stats.DistributionalParameters(counts[t])
            s.setFormat(options.format_branch_length)
            outfile.write("%s%s\t%s\n" % (prefix_row, t, str(s)))

    # ------------------------------------------------------------------------
    def printResolutionHistograms(outfile, data, field, f, options):

        counts = getCounts(data, field, f)

        histograms = []

        for t in sorted(counts.keys()):
            h = Histogram.Calculate(counts[t],
                                    increment=options.hist_dist_bin_size,
                                    min_value=options.hist_dist_min_value,
                                    max_value=options.hist_dist_max_value,
                                    no_empty_bins=True)
            histograms.append(h)

        combined_histogram = Histogram.Combine(histograms)
        outfile.write("%s\t%s\n" % ("bin", "\t".join(sorted(counts.keys()))))
        Histogram.Write(outfile, combined_histogram)

    # ------------------------------------------------------------------------
    def middle(*args):
        return sorted(args)[1]

    ###########################################################################
    ###########################################################################
    ###########################################################################
    # main function body
    ###########################################################################
    total_resolution = []
    for line in infile:
        if line[0] == "#":
            continue
        if line[:4] == "type":
            continue

        cluster, gene_node_id, type, tree, d_np, d_nl, d_nr = line[
            :-1].split("\t")[:7]

        d_np, d_nl, d_nr = map(float, (d_np, d_nl, d_nr))

        total_resolution.append(
            (cluster, gene_node_id, type, tree, d_np, d_nl, d_nr))

    # get distribution of minimum/maximum separation per node type
    for method in options.analyze_resolution_data:

        section = "resolution_%s" % (method)

        outfile, is_new = getFile(options, section)

        if method == "stats":
            printResolutionCounts(sys.stdout, total_resolution, 2, min, options,
                                  prefix_header="prefix\t", prefix_row="min\t")

        elif method == "histograms":
            for field, name in ((2, "category"), (3, "tree")):
                for f, method in ((min, "min"), (middle, "middle")):
                    section = "hist_%s_%s" % (name, method)
                    outfile, is_new = getFile(options, section)
                    printResolutionHistograms(
                        outfile, total_resolution, field, f, options)
                    if outfile != options.stdout:
                        outfile.close()

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def setupSpeciesTrees(species_tree_lines, options):
    """setup species trees.

    Species trees are pruned before reconciliation to include only species that are present
    in the gene tree. Thus they have to be recreated before each new tree.
    """
    species_nexus = TreeTools.Newick2Nexus(species_tree_lines)
    Tree.updateNexus(species_nexus)
    nspecies_trees = len(species_nexus.trees)

    if options.svg_otus:
        for species_tree in species_nexus.trees:
            TreeTools.PruneTree(species_tree, options.svg_otus)

    return species_nexus


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
def annotateDuplicationEvents(best, gene_tree, species_nexus, best_nodes,
                              counts_duplications, duplications,
                              min_distance, max_distance, distance_to_root,
                              extract_species, extract_quality,
                              map_id2location={}):
    """annotate relative location of a duplication event.

    updates counts and returns True, if a value has been added.

    A duplication event is bracketed by speciation events.
    Locate the relative position of the duplication event on the
    tree compared to the surrounding speciation events.

           /-------------------------------- Species X, Y, Z
          /            /-------------- Species A
         /      /-----S2
        /      /       \-------------- Species B
       /      /
    --S1-----D1        /--------------- Species A
              \------S3
                     \---------------  Species B

    From this arrangement, the relative duplication time of D1 is then |S1-D1| / f( |S1-S2|, |S1-S3|),
    where f is min, max or mean. This then allows to date duplications on a certain
    branch in the species tree.

    Computational details
       * There can be several duplication events between S1 and S2/S3
          * these would be interesting cases, ancient but non-recent expansions
       * S2 and S3 can also be speciations with deletions if there are more
         than three child species.
       * The calculation for terminal duplications is trivial.

    The pseudogene status of a duplication is decided as follows:
       * all children are pseudogenes: pseudogene
       * no  children are pseudogenes: gene
       * some children are pseudogenes: mixed

    """

    def stop_function(node_id):
        if best_nodes[node_id].mType in options.nodetypes_speciation:
            return True

    # find parent speciation
    parent, distance_to_parent = TreeTools.GetParentNodeWhereTrue(
        best.mGeneNode, gene_tree, stop_function)

    # find child speciation
    children = TreeTools.GetChildNodesWhereTrue(
        best.mGeneNode, gene_tree, stop_function)

    gene_node = best.mGeneNode

    best_tree = best.mSpeciesTree

    # assign a duplication to first child node that is the same as the
    # the duplication node. If it is a lineage specific duplication,
    # the tree does not matter.
    if best.mType != "DuplicationLineage":
        for child, distance in children:
            if best_nodes[child].mSpeciesTree == best_tree:
                child_species_node = best_nodes[child].mSpeciesNode
                break
        else:
            if options.loglevel >= 1:
                options.stdlog.write("# tree %s: %s at %i found no speciation from same tree.\n" % (
                    gene_tree.name, best.mType, gene_node))
            return False
    else:
        child_species_node = best_nodes[children[0][0]].mSpeciesNode

    # used to perform various checks whether trees were conguent
    # not necessary
    if False:
        # check if we have the same species tree between parent and children
        if best.mType != "DuplicationLineage":
            if best_nodes[parent].mSpeciesTree != best_tree:
                if options.loglevel >= 1:
                    options.stdlog.write("# tree %s: %s at %i between speciation nodes from different trees.\n" % (
                        gene_tree.name, best.mType, gene_node))
                return False

            for child, distance in children:
                if best_nodes[child].mSpeciesTree != best_tree:
                    if options.loglevel >= 1:
                        options.stdlog.write("# tree %s: %s at %i between speciation nodes from different trees.\n" % (
                            gene_tree.name, best.mType, gene_node))
                    return False

            # check, if the two speciation nodes are adjacent and the
            # the child nodes correspond
            species_tree = species_nexus.trees[best_tree]
            parent_species_node = best_nodes[parent].mSpeciesNode
            child_species_node = best_nodes[children[0][0]].mSpeciesNode

            # The following tests catch all duplications with deletions (speciation
            # nodes of children will differ). Turned off, might add option
            # later.
            if False:
                if species_tree.node(child_species_node).prev != parent_species_node:
                    species_tree.display()
                    if options.loglevel >= 1:
                        options.stdlog.write("# tree %s: %s at %i: speciation events in tree at %i and %i not adjacent.\n" %
                                             (gene_tree.name, best.mType, best.mGeneNode, parent_species_node, child_species_node))
                    return False

                # check if all the children are the same
                for child, distance in children:

                    species_tree.display()
                    if child_species_node != best_nodes[child].mSpeciesNode:
                        if options.loglevel >= 1:
                            options.stdlog.write("# tree %s: %s at %i: incongruent speciation events at %i and %i.\n" %
                                                 (gene_tree.name, best.mType, best.mGeneNode, child_species_node, best_nodes[child].mSpeciesNode))

                        return False

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # compute relative position of duplication event
    # Values are recorded for the parent node, because
    # child nodes might be different.
    max_distance_to_child = max(map(lambda x: x[1], children))
    total_distance = distance_to_parent + max_distance_to_child

    if total_distance == 0.0:
        if options.loglevel >= 3:
            options.stdlog.write("# tree %s: %s at %i: zero distance.\n" %
                                 (gene_tree.name, best.mType, best.mGeneNode))
        return False

    rel_distance = 100.0 * (max_distance_to_child / total_distance)
    if child_species_node not in counts_duplications[best_tree]:
        counts_duplications[best_tree][child_species_node] = []

    counts_duplications[best_tree][child_species_node].append(rel_distance)

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # Compute various properties based on children
    taxa = gene_tree.get_taxa(gene_node)

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # compute number of genes/pseudogenes among the children
    taxa_codes = map(extract_quality, taxa)

    ngenes, npseudogenes, nothers = 0, 0, 0
    for child in taxa_codes:
        if child in MAP_FILTER2QUALITY['genes']:
            ngenes += 1
        elif child in MAP_FILTER2QUALITY['pseudogenes']:
            npseudogenes += 1
        else:
            nothers += 1

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # obtain information about locations of duplicated genes
    locations = {}
    for t in taxa:

        location = map_id2location[t]
        species = extract_species(t)

        try:
            locations[species].append((t, location))
        except KeyError:
            locations[species] = [(t, location)]

    nspecies = len(locations.keys())

    max_contigs = 0
    for species, l in locations.items():
        contigs = set()
        for t, ll in l:
            contigs.add(ll.contig)

        max_contigs = max(max_contigs, len(contigs))

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # save results
    d = Duplication()
    d.mGeneTree = gene_tree.name
    d.mGeneNode = best.mGeneNode
    d.mSpeciesNode = best.mSpeciesNode
    d.mSpeciesTree = best.mSpeciesTree
    d.mNSpecies = nspecies
    d.mNChildren = len(taxa)
    d.mDistance2Root = distance_to_root[gene_node]
    d.mDistance2Min = min_distance[gene_node]
    d.mDistance2Max = max_distance[gene_node]
    d.mRelativeHeight = max_distance_to_child
    d.mTotalDistance = total_distance
    d.mRelativeDistance = rel_distance
    d.mDistance2Parent = distance_to_parent
    d.mDistance2Children = tuple([x[1] for x in children])
    d.mNGenes = ngenes
    d.mNPseudogenes = npseudogenes
    d.mNOthers = nothers
    d.mNMaxContigs = max_contigs

    duplications.append(d)

    return True


def readSpeciesTrees(options):
    """read species tree/trees."""

    if options.filename_species_tree:
        species_tree_lines = open(
            options.filename_species_tree, "r").readlines()
    elif options.species_tree:
        species_tree_lines = options.species_tree
    else:
        raise "please supply a species tree."

    species_nexus = TreeTools.Newick2Nexus(species_tree_lines)

    Tree.updateNexus(species_nexus)
    nspecies_trees = len(species_nexus.trees)

    if options.loglevel >= 1:
        options.stdlog.write("# read %i species trees.\n" % nspecies_trees)

        if options.loglevel >= 10:
            for x in range(nspecies_trees):
                options.stdlog.write(
                    "# %i: %s\n" % (x, TreeTools.Tree2Newick(species_nexus.trees[x])))

        options.stdlog.flush()

    return species_nexus, species_tree_lines


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
def processGeneTrees(gene_nexus, options):
    """preprocess gene trees before analysis

    - root them if necessary.
    - remove unwanted taxa.
    - remove unwanted trees.
    """
    ninput = 0
    nskipped_filter = 0
    nfiltered = 0

    new_trees = []

    for gene_tree in gene_nexus.trees:

        ninput += 1

        xname = re.sub("_tree.*", "", gene_tree.name)
        xname = re.sub("subtree_", "", xname)

        if filter_positives and xname not in filter_positives:
            nskipped_filter += 1
            continue

        #######################################################################
        #######################################################################
        #######################################################################
        # apply filters to gene tree
        #######################################################################
        filterTree(gene_tree, options, map_id2location)

        otus = TreeTools.GetTaxa(gene_tree)

        if len(otus) <= 1:
            nfiltered += 1
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# tree %s: less than two species (%i) after filtering - skipped.\n" % (gene_tree.name, len(otus)))
                options.stdlog.flush()
            continue

        # check, if only outgroups
        if options.outgroup_species and not set(map(extract_species, otus)).difference(options.outgroup_species):
            nfiltered += 1
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# tree %s: only outgroups after filtering - skipped.\n" % gene_tree.name)
            continue

        #######################################################################
        #######################################################################
        #######################################################################
        # reroot gene tree, if outgroups have been given.
        #######################################################################

        if options.outgroup_species:
            rerootTree(gene_tree, extract_species, options)

        new_trees.append(gene_tree)

    gene_nexus.trees = new_trees

    if options.loglevel >= 1:
        options.stdlog.write("# after processing trees: ninput=%i, nskipped_filter=%i, nfiltered=%i, noutput=%i\n" %
                             (ninput, nskipped_filter, nfiltered, len(gene_nexus.trees)))

    return

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/analyze_orthology_multiple.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("--species-regex", dest="species_regex", type="string",
                      help="regular expression to extract species from identifier.")

    parser.add_option("--gene-regex", dest="gene_regex", type="string",
                      help="regular expression to extract gene from identifier.")

    parser.add_option("-s", "--filename-species-tree", dest="filename_species_tree", type="string",
                      help="filename with species tree.")

    parser.add_option("-g", "--filename-gene-trees", dest="filename_gene_trees", type="string",
                      help="filename with gene trees [default: %default].")

    parser.add_option("--filename-filter-positives", dest="filename_filter_positives", type="string",
                      help="filename with positive list of trees to analyze.")

    parser.add_option("--filename-species2colour", dest="filename_species2colour", type="string",
                      help="filename with map of species to colours. If not given, random colours are assigned to species.")

    parser.add_option("--filename-species2name", dest="filename_species2name", type="string",
                      help="filename for mapping species names.")

    parser.add_option("-t", "--species-tree", dest="species_tree", type="string",
                      help="species tree.")

    parser.add_option("-e", "--locations-tsv-file", dest="filename_locations", type="string",
                      help="filename with map of transcript information to location information.")

    parser.add_option("--no-create", dest="create", action="store_false",
                      help="do not create files, but append to them.")

    parser.add_option("--max-separation", dest="max_separation", type="int",
                      help="maximum allowable separation between syntenic segments for border plot (set to 0, if syntey is enough).")

    parser.add_option("--filename-species2url", dest="filename_species2url", type="string",
                      help="filename with mapping information of species to URL.")

    parser.add_option("--column-prefix", dest="prefix", type="string",
                      help="prefix to add as first column.")

    parser.add_option("--outgroup-species", dest="outgroup_species", type="string",
                      help="species to used as outgroups. Separate multiple species by ','.")

    parser.add_option("--keep-outgroup-clusters", dest="keep_outgroup_clusters", type="string",
                      help="Keep clusters that are exclusively outgroups.")

    parser.add_option("--subtrees-trees", dest="subtrees_trees", action="store_true",
                      help="write trees for subtrees.")

    parser.add_option("--subtrees-identifiers", dest="subtrees_identifiers", action="store_true",
                      help="write identifiers of subtrees.")

    parser.add_option("--svg-add-ids", dest="svg_add_ids", action="store_true",
                      help="add node ids to svg plot.")

    parser.add_option("--svg-print-location", dest="svg_print_location", action="store_true",
                      help="add location information in gene trees.")

    parser.add_option("--svg-no-legend", dest="svg_add_legend", action="store_false",
                      help="do not add legend to plots.")

    parser.add_option("--svg-otus", dest="svg_otus", type="string",
                      help="otus to output in svg species tree.")

    parser.add_option("--svg-branch-lengths", dest="svg_branch_lengths", type="choice",
                      choices=("contemporary", "uniform", "median"),
                      help="branch lengths in species tree.")

    parser.add_option("--print-totals", dest="print_totals", action="store_true",
                      help="output totals sections.")

    parser.add_option("--print-subtotals", dest="print_subtotals", action="store_true",
                      help="output subtotals sections.")

    parser.add_option("--print-summaries", dest="print_summaries", action="store_true",
                      help="output summary sections.")

    parser.add_option("--print-best", dest="print_best", action="store_true",
                      help="output best node assignment for each node in gene tree.")

    parser.add_option("--print-svg", dest="print_svg", action="store_true",
                      help="output svg files.")

    parser.add_option("--print-species-svg", dest="print_species_svg",
                      action="store_true",
                      help="output species svg files.")

    parser.add_option(
        "--output-filename-pattern", dest="output_pattern", type="string",
        help="""output pattern for separate output of sections [default: %default].
        Set to None, if output to stdout. Can contain one %s to be substituted with section.""")

    parser.add_option("--output-pattern-svg", dest="output_pattern_svg", type="string",
                      help="filename for svg output. If it contains %s, this is replaced by gene_tree name.")

    parser.add_option("--filename-node-types", dest="filename_node_types", type="string",
                      help="filename with node type information from a previous run.")

    parser.add_option("--filename-duplications", dest="filename_duplications", type="string",
                      help="filename with duplication information from previous run.")

    parser.add_option("--filename-species-list", dest="filename_species_list", type="string",
                      help="filename with species to use in the analysis.")

    parser.add_option("--filename-contig-sizes", dest="filename_contig_sizes", type="string",
                      help="filename with contig sizes for each species (used for wheel plot - should contain %s).")

    parser.add_option("--filename-output-gene-trees", dest="filename_output_gene_trees", type="string",
                      help="filename to which pre-processed gene trees are output.")

    parser.add_option("--skip-preprocessing", dest="skip_preprocessing", action="store_true",
                      help="skip preprocessing of gene trees.")

    parser.add_option("--analyze-resolution-data", dest="analyze_resolution_data", type="choice", action="append",
                      choices=("stats", "histograms"),
                      help="stdin is resolution data.")

    parser.add_option("--analyze-duplication-data", dest="analyze_duplication_data", type="choice", action="append",
                      choices=("stats", "cumul-hist", "cumul-norm-hist", "annotate", "members", "svg-wheel-plot",
                               "local-dups"),
                      help="""stdin is duplication data. Perform the following analyses:
members:    print duplicated genes for each cluster.
local-dups: print locally duplicated genes in each cluster. Local duplications are lineage specific and are not
                separated by any other gene.
""")

    parser.add_option("--filter-quality", dest="filter_quality", type="choice",
                      choices=("all", "genes", "pseudogenes"),
                      help="filter predictions by gene type.")

    parser.add_option("--filter-location", dest="filter_location", type="choice",
                      choices=("all", "local", "non-local", "cis", "unplaced"),
                      help="filter predictions by location.")

    parser.add_option("--remove-unplaced", dest="remove_unplaced", action="store_true",
                      help="remove predictions on unplaced contigs.")

    parser.add_option("--svg-template", dest="svg_template", type="choice",
                      choices=("screen", "publication", "poster"),
                      help="template to use for figures. The template changes colours and line styles.")

    parser.add_option("--output-format", dest="output_format", type="choice",
                      choices=("plain", "table"),
                      help="set output format. Same analysis can be put in tabular format.")

    parser.set_defaults(
        filter_quality="all",
        filter_location="all",
        remove_unplaced=False,
        species_regex="^([^|]+)",
        gene_regex="^[^|]+[|]([^|]+)",
        filename_species_tree=None,
        priority={"Leaf": 0,
                   "Speciation": 1,
                  "SpeciationDeletion": 2,
                   "Transcripts": 3,
                  "DuplicationLineage": 4,
                  "Duplication": 5,
                  "DuplicationDeletion": 6,
                  "DuplicationInconsistency": 7,
                  "Outparalogs": 8,
                  "InconsistentTranscripts": 9,
                  "Inconsistency": 10,
                  "Masked": 11},
        species_tree=None,
        filename_species2colour=None,
        filename_species2name=None,
        filename_locations=None,
        max_separation=0,
        filename_species2url=None,
        separator="|",
        prefix=None,
        output_pattern=None,
        output_pattern_svg=None,
        output_format="plain",
        outgroup_species=None,
        keep_outgroup_clusters=False,
        svg_add_ids=False,
        svg_branch_lengths="median",
        svg_otus=None,
        svg_print_location=False,
        subtrees=False,
        print_svg=False,
        print_subtotals=False,
        print_totals=False,
        print_best=False,
        print_summaries=False,
        subtrees_identifiers=False,
        create=True,
        min_branch_length=0.00,
        filename_node_types=None,
        filename_duplications=None,
        format_branch_length="%6.4f",
        nodetypes_inconsistency=("InconsistentTranscripts", "Inconsistency"),
        nodetypes_duplication=(
            "Duplication", "DuplicationLineage", "DuplicationDeletion"),
        nodetypes_speciation=("Speciation", "SpeciationDeletion"),
        analyze_resolution_data=None,
        analyze_duplication_data=None,
        hist_dist_bin_size=0.01,
        hist_dist_min_value=0,
        hist_dist_max_value=None,
        hist_dups_bin_size=1.0,
        hist_dups_min_value=0,
        hist_dups_max_value=None,
        filename_filter_positives=None,
        filename_species_list=None,
        filename_map_annotation_interpro="/net/cpp-data/backup/andreas/projects/flies/data_1v5/interpro.list",
        filename_gene_trees="-",
        filename_output_gene_trees=None,
        filename_contig_size=None,
        skip_preprocessing=False,
        svg_wheel_plot_radius=1500,
        svg_wheel_plot_min_value=0.0,
        svg_wheel_plot_max_value=0.2,
        svg_wheel_plot_radius_increment=40,
        svg_wheel_plot_quality2symbol={
            'CG': "circle", 'PG': "circle", 'SG': "circle"},
        svg_wheel_plot_quality2mask=(
            "RG", "CP", "PP", "SP", "RP", "CF", "PF", "SF", "UG", "UP", "UF", "BF", "UK"),
        svg_template="screen",
        svg_add_legend=True,
    )

    (options, args) = E.Start(parser, add_csv_options=True)

    if options.analyze_resolution_data:
        analyzeResolutionData(sys.stdin, options)
        E.Stop()
        sys.exit(0)

    if options.outgroup_species:
        options.outgroup_species = set(options.outgroup_species.split(","))

    if options.svg_otus:
        options.svg_otus = set(options.svg_otus.split(","))

    #########################################################################
    #########################################################################
    #########################################################################
    # read species trees
    #########################################################################
    species_nexus, species_tree_lines = readSpeciesTrees(options)

    nspecies_trees = len(species_nexus.trees)

    rx_species = re.compile(options.species_regex)
    extract_species = lambda x: rx_species.match(x).groups()[0]

    if options.gene_regex:
        rx_gene = re.compile(options.gene_regex)
        extract_gene = lambda x: rx_gene.match(x).groups()[0]
    else:
        extract_gene = None

    def extract_quality(x):
        try:
            q = x.split(options.separator)[3]
        except IndexError:
            q = "CG"
        return q

    if options.filename_species_list:
        species_list, nerrors = IOTools.ReadList(
            open(options.filename_species_list, "r"))
        # prune species trees
        for tree in species_nexus.trees:
            TreeTools.PruneTree(tree, species_list)
        # write them to lines, so that they can be recreated.
        species_tree_lines = TreeTools.Nexus2Newick(species_nexus).split("\n")
        # re-read them so that node numbering is consistent
        species_nexus = setupSpeciesTrees(species_tree_lines, options)
    else:
        species_list = set()
        for tree in species_nexus.trees:
            species_list = species_list.union(set(tree.get_taxa()))
        species_list = list(species_list)

    #########################################################################
    #########################################################################
    #########################################################################
    # read positive list of malis
    #########################################################################
    if options.filename_filter_positives:
        filter_positives, nerrors = IOTools.ReadList(
            open(options.filename_filter_positives, "r"))
        filter_positives = set(filter_positives)
    else:
        filter_positives = None

    #########################################################################
    #########################################################################
    #########################################################################
    # read location info
    #########################################################################
    if (options.remove_unplaced or options.filter_location != "all") and not options.filename_locations:
        raise "please supply a file with location information."

    if options.filename_locations:

        t0 = time.time()

        if options.loglevel >= 1:
            options.stdlog.write(
                "# reading locations from %s\n" % options.filename_locations)
            options.stdlog.flush()

        map_id2location = TreeReconciliation.readLocations(
            open(options.filename_locations, "r"), extract_species)

        if options.loglevel >= 1:
            options.stdlog.write("# read %i locations in %i seconds\n" % (
                len(map_id2location), time.time() - t0))
            options.stdlog.flush()
    else:
        map_id2location = {}

        if options.loglevel >= 1:
            options.stdlog.write("# reading gene trees.\n")
            options.stdlog.flush()

    #########################################################################
    #########################################################################
    #########################################################################
    # read gene trees
    #########################################################################
    if options.filename_gene_trees:

        t0 = time.time()

        infile = getInfile(options.filename_gene_trees)

        gene_nexus = TreeTools.Newick2Nexus(infile)

        if infile != sys.stdin:
            infile.close()

        if options.loglevel >= 1:
            options.stdlog.write("# read %i gene trees from %s in %i seconds.\n" % (
                len(gene_nexus.trees), options.filename_gene_trees, time.time() - t0))
            options.stdlog.flush()
    else:
        gene_nexus = None

    #########################################################################
    #########################################################################
    #########################################################################
    # preprocess gene trees (filtering/rerooting)
    #########################################################################
    if not options.skip_preprocessing and gene_nexus is not None and \
       len(gene_nexus.trees) > 0:

        t0 = time.time()

        processGeneTrees(gene_nexus, options)

        if options.loglevel >= 1:
            options.stdlog.write("# processed %i gene trees in %i seconds.\n" % (
                len(gene_nexus.trees), time.time() - t0))
            options.stdlog.flush()

    if options.filename_output_gene_trees:
        if options.loglevel >= 1:
            options.stdlog.write("# writing %i gene trees to %s.\n" % (
                len(gene_nexus.trees), options.filename_output_gene_trees))
            options.stdlog.flush()

        outfile = open(options.filename_output_gene_trees, "w")
        outfile.write(TreeTools.Nexus2Newick(gene_nexus, with_names=True))
        outfile.close()

    #########################################################################
    #########################################################################
    #########################################################################
    # perform analysis on previous runs
    #########################################################################
    if options.analyze_duplication_data and options.filename_duplications:

        infile = getInfile(options.filename_duplications)
        duplications = readDuplications(infile)
        if infile == sys.stdin:
            infile.close()

        analyzeDuplicationData(duplications,
                               species_nexus,
                               options,
                               extract_gene,
                               extract_species,
                               gene_nexus,
                               map_id2location)
        E.Stop()
        sys.exit(0)

    #########################################################################
    #########################################################################
    #########################################################################
    # print table headers
    #########################################################################
    if options.output_format == "table":
        if options.print_svg:
            options.stdout.write("name\ttree\tsvg\n")

    #########################################################################
    #########################################################################
    #########################################################################
    # setup total counts data structures
    #########################################################################
    total_counts_trees = []

    for x in range(nspecies_trees):
        a = {}
        for y in options.priority.keys():
            a[y] = 0
        total_counts_trees.append(a)

    # heights per tree and node and type
    total_counts_nodes = [{} for x in range(nspecies_trees)]

    # relative heights per tree and node for duplications
    total_counts_duplications = [{} for x in range(nspecies_trees)]

    # list of duplications tuple is: best, height, total
    total_duplications = []

    total_counts_categories = {}
    for x in options.priority.keys():
        total_counts_categories[x] = []

    # resolution data
    total_resolution = []

    if not options.filename_node_types:

        #######################################################################
        #######################################################################
        #######################################################################
        # delete output files
        #######################################################################
        if options.create and options.output_pattern:
            for section in ("details", "subtrees", "subids", "details", "trees", "nodes", "categories"):
                fn = options.output_pattern % section
                if os.path.exists(fn):
                    if options.loglevel >= 1:
                        options.stdlog.write("# deleting file %s.\n" % fn)
                    os.remove(fn)

        #######################################################################
        #######################################################################
        #######################################################################
        # main loop over gene trees
        #######################################################################
        ninput, nfiltered, noutput = 0, 0, 0
        nduplications, nskipped_duplications = 0, 0
        nskipped_filter = 0
        nskipped_refheight = 0
        nskipped_analysis_set = 0

        for gene_tree in gene_nexus.trees:

            ninput += 1

            ###################################################################
            ###################################################################
            ###################################################################
            # get identifier for this tree and update prefixes accordingly
            ###################################################################
            if options.prefix:
                if len(gene_nexus.trees) > 0:
                    prefix_header = "prefix1\tprefix2\t"
                    prefix_row = options.prefix + "\t" + gene_tree.name + "\t"
                    prefix_prefix = options.prefix + "_" + gene_tree.name + "_"
                    prefix_name = options.prefix + "_" + gene_tree.name
                else:
                    prefix_header = "prefix\t"
                    prefix_row = options.prefix + "\t"
                    prefix_prefix = options.prefix + "_"
                    prefix_name = options.prefix
            else:
                if len(gene_nexus.trees) > 0:
                    prefix_header = "prefix\t"
                    prefix_row = gene_tree.name + "\t"
                    prefix_prefix = gene_tree.name + "\t"
                    prefix_name = gene_tree.name
                else:
                    prefix_header, prefix_row, prefix_prefix, prefix_name = "", "", "", ""

            ###################################################################
            ###################################################################
            ###################################################################
            # get analysis set
            ###################################################################
            analysis_set, gene_set, pseudogene_set, other_set = getAnalysisSets(
                gene_tree, extract_quality, options)

            if len(analysis_set) == 0:
                nskipped_analysis_set += 1
                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# tree %s: empty analysis set - skipped.\n" % gene_tree.name)
                continue

            ###################################################################
            ###################################################################
            ###################################################################
            # compute various distance parameters for each node
            ###################################################################
            min_branch_lengths, max_branch_lengths = TreeTools.GetBranchLengths(
                gene_tree)
            distance_to_root = TreeTools.GetDistanceToRoot(gene_tree)

            reference_height = getReferenceHeight(distance_to_root,
                                                  gene_tree,
                                                  gene_set,
                                                  options,
                                                  extract_species,
                                                  method="median")

            if reference_height is None:
                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# tree %s: reference height not computable or 0 - skipped.\n" % gene_tree.name)
                nskipped_refheight += 1
                continue

            otus = TreeTools.GetTaxa(gene_tree)

            ###################################################################
            ###################################################################
            ###################################################################
            # restrict analysis to analysis set
            ###################################################################
            if options.filter_quality != "all":

                TreeTools.PruneTree(
                    gene_tree, set([gene_tree.node(x).data.taxon for x in analysis_set]))

                nstart = len(otus)

                if len(otus) <= 1:
                    nfiltered += 1
                    if options.loglevel >= 1:
                        options.stdlog.write(
                            "# tree %s: empty after filtering by analysis set - skipped.\n" % gene_tree.name)
                        options.stdlog.flush()
                    continue

                # check, if only outgroups
                if options.outgroup_species and not set(map(extract_species, otus)).difference(options.outgroup_species):
                    if not options.keep_outgroup_clusters:
                        nfiltered += 1
                        if options.loglevel >= 1:
                            options.stdlog.write(
                                "# tree %s: only outgroups after filtering - skipped.\n" % gene_tree.name)
                        continue

                if options.loglevel >= 1:
                    options.stdlog.write("# tree %s: analysis_set filtering - before: genes=%i, pseudogenes=%i, other=%i, total=%i\n" % (
                        gene_tree.name,
                        len(gene_set), len(pseudogene_set), len(other_set),
                        nstart))

                    options.stdlog.write("# tree %s: analysis_set filtering - after: analysis_set=%i, total=%i, reference_height=%f\n" % (
                        gene_tree.name,
                        len(analysis_set), len(otus), reference_height))

            gene_tree_species_list = list(set(map(extract_species, otus)))

            ###################################################################
            ###################################################################
            ###################################################################
            # re-read species trees
            ###################################################################
            species_nexus = setupSpeciesTrees(species_tree_lines, options)

            ###################################################################
            ###################################################################
            ###################################################################
            # apply RIO tree reconciliation
            ###################################################################
            node_types = [[] for x in range(TreeTools.GetSize(gene_tree))]

            for s in range(nspecies_trees):

                species_tree = species_nexus.trees[s]

                if options.loglevel >= 2:
                    options.stdlog.write("# analysis with tree %i\n" % s)

                TreeTools.PruneTree(species_tree, gene_tree_species_list)

                t = TreeTools.ReconciliateByRio(gene_tree, species_tree, extract_species, extract_gene,
                                                outgroup_species=options.outgroup_species,
                                                min_branch_length=options.min_branch_length)

                for n in range(len(node_types)):
                    node_types[n].append(t[n])

            ###################################################################
            # perform various analysis

            ###################################################################
            ###################################################################
            ###################################################################
            # find best node for node in species tree
            ###################################################################
            best_nodes = []

            for n in range(len(node_types)):

                # node_type for reference (first tree) to get gene tree info.
                # This is the same for all node_types
                node_type = node_types[n][0]
                if not node_type:
                    best_nodes.append(None)
                    continue

                if options.loglevel >= 2:
                    options.stdout.write("# multiple assignments for node %i\t%s\n" %
                                         (n,
                                          options.format_branch_length % ((min_branch_lengths[node_type.mGeneNode] + max_branch_lengths[node_type.mGeneNode]) / 2)))
                    options.stdout.write("# %s\t%s\n" % ("\t" * 4,
                                                         "\t".join(map(lambda x: x.mType, node_types[n]))))
                    options.stdlog.flush()

                    # output species for the genes joined
                    for s in gene_tree.node(node_type.mGeneNode).succ:
                        options.stdout.write(
                            "#" + "\t" * 5 + ",".join(tuple(set(map(extract_species, gene_tree.get_taxa(s))))) + "\n")

                # get the best node
                # the best node is the one with highest priority
                ##
                best_tree = 0
                best = node_types[n][0]
                for x in range(1, nspecies_trees):

                    this = node_types[n][x]
                    if options.priority[best.mType] > options.priority[this.mType]:
                        best_tree = x
                        best = this

                best.mSpeciesTree = best_tree
                best_nodes.append(best)

            if options.print_best:
                printBestNodes(best_nodes, options, prefix_header, prefix_row)

            ###################################################################
            ###################################################################
            ###################################################################
            # build subtrees
            ###################################################################
            if options.subtrees_trees or options.subtrees_identifiers:
                clusters = extractSubtrees(gene_tree, extract_species, options)
                printClusters(clusters, options, prefix_prefix, otus)
            else:
                clusters = None

            ###################################################################
            ###################################################################
            ###################################################################
            # compute counts over best nodes for categories, trees, etc.
            ###################################################################
            counts_trees = []
            for x in range(nspecies_trees):
                a = {}
                for y in options.priority.keys():
                    a[y] = 0
                counts_trees.append(a)

            # heights per tree and node and type
            counts_nodes = [{} for x in range(nspecies_trees)]

            # relative heights per tree and node for duplications
            counts_duplications = [{} for x in range(nspecies_trees)]

            duplications = []

            counts_categories = {}
            for x in options.priority.keys():
                counts_categories[x] = []

            for best in best_nodes:

                if best is None:
                    continue
                if best.mType == "Leaf":
                    continue

                best_tree = best.mSpeciesTree

                species_tree = species_nexus.trees[best_tree]

                height = max_branch_lengths[best.mGeneNode]

                counts_categories[best.mType].append(height)

                counts_trees[best_tree][best.mType] += 1

                if best.mSpeciesNode not in counts_nodes[best_tree]:
                    counts_nodes[best_tree][best.mSpeciesNode] = {}

                if best.mType not in counts_nodes[best_tree][best.mSpeciesNode]:
                    counts_nodes[best_tree][best.mSpeciesNode][best.mType] = []

                counts_nodes[best_tree][best.mSpeciesNode][
                    best.mType].append(height)

                # record branch length of node. Do not include terminal nodes
                node = gene_tree.node(best.mGeneNode)
                if node.succ:
                    total_resolution.append((gene_tree.name, best.mGeneNode,
                                             best.mType, best_tree, node.data.branchlength,
                                             gene_tree.node(
                                                 node.succ[0]).data.branchlength,
                                             gene_tree.node(node.succ[1]).data.branchlength))

                ###############################################################
                ###############################################################
                ###############################################################
                # annotate duplication nodes with relative tree height
                ###############################################################
                if best.mType in options.nodetypes_duplication:
                    retval = annotateDuplicationEvents(best,
                                                       gene_tree, species_nexus,
                                                       best_nodes,
                                                       counts_duplications,
                                                       duplications,
                                                       min_branch_lengths,
                                                       max_branch_lengths,
                                                       distance_to_root,
                                                       extract_species,
                                                       extract_quality,
                                                       map_id2location=map_id2location)

                    if retval:
                        nduplications += 1
                    else:
                        nskipped_duplications += 1

            if options.print_subtotals:
                printCounts(counts_categories, counts_trees, counts_nodes,
                            species_nexus, options, prefix_header, prefix_row)

            # update total counts
            appendCounts(total_counts_categories, counts_categories)
            appendCounts(total_counts_nodes, counts_nodes)
            appendCounts(total_counts_duplications, counts_duplications)
            addCounts(total_counts_trees, counts_trees)
            total_duplications += duplications

            ###################################################################
            ###################################################################
            ###################################################################
            # print summary
            ###################################################################
            if options.print_summaries:
                printSummary(counts_categories, counts_trees, counts_nodes, duplications,
                             reference_height,
                             gene_tree, gene_tree_species_list,
                             species_nexus, species_list,
                             options, prefix_header, prefix_row)

            ###################################################################
            ###################################################################
            ###################################################################
            # output annotated tree to svg file
            ###################################################################
            if options.print_svg:
                printSVGGeneTree(gene_tree, options, extract_species)

            noutput += 1

        #######################################################################
        #######################################################################
        #######################################################################
        # print totals
        #######################################################################
        if len(gene_nexus.trees) > 0:

            species_nexus = setupSpeciesTrees(species_tree_lines, options)

            if options.prefix:
                prefix_header = "prefix1\tprefix2\t"
                prefix_row = options.prefix + "\t" + "total" + "\t"
                prefix_prefix = options.prefix + "_" + "total" + "_"
                prefix_name = options.prefix + "_" + "total"
            else:
                prefix_header = "prefix\t"
                prefix_row = "total" + "\t"
                prefix_prefix = "total" + "_"
                prefix_name = "total"

            if options.print_totals:
                printCounts(total_counts_categories, total_counts_trees, total_counts_nodes,
                            species_nexus,
                            options, prefix_header, prefix_row)

            if total_resolution and options.print_totals:
                printResolution(
                    total_resolution, options, prefix_header, prefix_row)

            if total_duplications and options.print_totals:
                printDuplications(total_duplications, species_nexus, options)

        if options.loglevel >= 1:
            options.stdlog.write("# ninput=%i, nfiltered=%i, nskipped_filter=%i, noutput=%i, nskipped_duplications=%i, nskipped_refheight=%i, nskipped_analysis_set=%i, nduplications=%i\n" %
                                 (ninput, nfiltered, nskipped_filter, noutput, nskipped_duplications, nskipped_refheight, nskipped_analysis_set, nduplications))
            options.stdlog.flush()
    else:
        if options.filename_node_types:
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# reading previous counts from %s.\n" % options.filename_node_types)

            total_counts_nodes = readCountsNodes(
                open(options.filename_node_types, "r"), nspecies_trees)

        elif options.filename_duplications:

            total_duplications = readDuplications(
                open(options.filename_duplications, "r"))

    #########################################################################
    #########################################################################
    #########################################################################
    # output annotated species trees
    #########################################################################
    if total_counts_nodes:

        species_nexus = setupSpeciesTrees(species_tree_lines, options)

        #######################################################################
        #######################################################################
        #######################################################################
        # print svg annotated species trees
        #######################################################################
        if options.print_species_svg:
            printSVGSpeciesTrees(
                species_nexus.trees, total_counts_nodes, options)

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

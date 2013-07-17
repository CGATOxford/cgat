################################################################################
#   Gene prediction pipeline 
#
#   $Id: plot_tree.py 2782 2009-09-10 11:40:29Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
tree2svg.py - plot a tree in svg format
========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

plot a tree in New-Hampshire format. The output is svg.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Examples
--------

The following plots a tree with branch length from :file:`ks.tree` adding
branch length information from :file:`kaks.tree`, `ks2.tree` and `kaks2.tree`::

   cat ks.tree kaks.tree ks2.tree kaks2.tree |\
   python plot_tree.py --verbose=0 --show-branchlengths --annotation=kaks --annotation=master --annotation=kaks --font-style-tips=italic > tree.svg

Adding tables to trees
----------------------

Branches can be decorated with tables, for example::

   cat test.tree tables.tree | python plot_tree.py -v 0 --filename-tables=test.table --annotation=tables

Here, the file :file:`test.tree` is::

   (((A:1,B:1):1,C:2),D:3);

The file :file:`tables.tree` links branches to tables. The table number is given as the branch length::

   (((A:5,B:1):2,C:3),D:4);

Note that table ids need to start with 1, but need not be consecutive numbers. Finally, the tables are
in the file :file:`test.table` given as argument to the options ``--filename-tables``.

   table=3
   header1 header2 header3
   0row1   col1    col2
   0row2   col1    col2
   table=1
   header1 header2 header3
   1row1    col1    col2
   1row2    col1    col2
   table=2
   header1 header2 header3
   2row1    col1    col2
   2row2    col1    col2

A table is started by the line ``table=#`` where ``#`` is an integer number. The actual table follows as
a tab-separated table with the first line being interpreted as a header. Lines starting with ``#`` are
ignored.

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

from types import *

import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools
import CGAT.IOTools as IOTools
import CGAT.SVGTree as SVGTree
import CGAT.SVGdraw as SVGdraw


###################################################################
###################################################################
###################################################################
##
###################################################################
class NodeDecoratorSupportPieChart(SVGTree.NodeDecorator):
    """class for decorating internal nodes using a
    pie-chart.

    Uses the support information in the tree.
    """
    
    mRadius = 40
    mFontSize = 12
    mFont = "Verdana"
    mSeparator = 5
    mStrokeWidth = 1
    
    def __init__(self, tree, *args, **kwargs):
        
        SVGTree.NodeDecorator.__init__( self, tree, *args, **kwargs)

    def getElements( self, node_id, x, y):

        node = self.mTree.node(node_id)
        e = []

        p = node.data.support

        if p == 1.0:
            e.append( SVGdraw.circle( x, y, self.mRadius,
                                      stroke = "rgb(%i,%i,%i)" % SVGTree.BLACK,
                                      fill = "rgb(%i,%i,%i)" % SVGTree.RED))
        elif p > 0.0:
            e.append( SVGdraw.circle( x, y, self.mRadius,
                                      stroke = "rgb(%i,%i,%i)" % SVGTree.BLACK,
                                      fill = "rgb(%i,%i,%i)" % SVGTree.WHITE))

            d = SVGdraw.pathdata( x, y )

            d.line( x + self.mRadius, y )

            angle = 360.0 * p
            dx = self.mRadius * math.cos( math.radians(angle) ) + x
            dy = self.mRadius * math.sin( math.radians(angle) ) + y

            if p <= 0.5:
                d.ellarc( self.mRadius, self.mRadius, 0, 0, 1, dx, dy )
            else:
                d.ellarc( self.mRadius, self.mRadius, 0, 1, 1, dx, dy )            

            e.append( SVGdraw.path( d,
                                    stroke = "rgb(%i,%i,%i)" % SVGTree.RED,
                                    fill = "rgb(%i,%i,%i)" % SVGTree.RED,
                                    stroke_width = self.mStrokeWidth ) )

        else:
            pass

        return e

    def getHeight( self, node_id ):
        return 0

    def getWidth( self, node_id ):
        return 0

###################################################################
###################################################################
###################################################################
##
###################################################################
class BranchDecoratorTable( SVGTree.BranchDecoratorHorizontal ):
    """branch decorator - add a table onto a branch length.
    
    The table will be plotted below each branch.

    This decorator requires labeled branches within a tree.
    """

    mBranchLengthFormat = "%5.2f"
    mFontSize = 10
    
    def __init__(self, tree, filename, *args, **kwargs):
        SVGTree.BranchDecoratorHorizontal.__init__( self, tree, *args, **kwargs)
        self.mWritten = 0

        infile = open(filename, "r")
        self.mTables = {}
        table_id = None        
        for line in infile:
            if line.startswith("#"): continue
            if line.startswith("table="):
                if table_id: self.mTables[table_id] = table
                table_id = re.search("table=(\d+)", line).groups()[0]
                if int(table_id) == 0: raise ValueError( "table id 0 is invalid" )
                table = []
                continue
            table.append( line[:-1].split("\t") )
        if table_id: self.mTables[table_id] = table

        self.mColumnWidths = {}
        
        for id, table in self.mTables.iteritems():
            if len(table) == 0:
                raise ValueError("empty table %s" % id)
            
            column_widths = [0] * len( table[0] )
            
            for row in table:
                if len(column_widths) != len(row):
                    raise ValueError("table %s has unequal row lengths" % id )

                for x, col in enumerate(row):
                    column_widths[x] = max( column_widths[x], len(col) )
            
            self.mColumnWidths[id] = column_widths

    def getElements( self, node_id, x1, x2, y ):

        e = SVGTree.BranchDecoratorHorizontal.getElements( self, node_id, x1, x2, y )

        table_id = str(int(self.mTree.node(node_id).data.branchlength))

        if table_id not in self.mTables:
            return e

        startx = x1 + self.mFontSize
        y = y + self.mFontSize

        table, column_widths = self.mTables[table_id], self.mColumnWidths[table_id]

        font_weight = "bold"
        
        for r, row in enumerate(table):
            x = startx
            
            for c, col in enumerate(row):
                e.append( SVGdraw.text( x, y,
                                        col,
                                        20, 
                                        self.mFont,
                                        stroke = "rgb(%i,%i,%i)" % self.mFontColour,
                                        font_weight = font_weight,
                                        text_anchor = "left" ))
                x += column_widths[c] * self.mFontSize // 2
                
            y +=  self.mFontSize
            font_weight = "normal"
            
        return e
    
    def getHeight( self, node_id ):
        
        table_id = str(int(self.mTree.node(node_id).data.branchlength))

        if table_id in self.mTables:
            table_height = len( self.mTables[table_id]) * self.mFontSize
        else:
            table_height = 0
        
        return 5 + table_height
        
    def getWidth( self, node_id ):
        table_id = str(int(self.mTree.node(node_id).data.branchlength))
        if table_id in self.mTables:
            return 5 + sum(self.mColumnWidths[table_id])

        return 5
        
###------------------------------------------------------------------------------------------------
def main():

    parser = E.OptionParser( version = "%prog version: $Id: plot_tree.py 2782 2009-09-10 11:40:29Z andreas $")

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
    parser.add_option( "--support-style", dest="support_style", type="choice",
                       choices=("pie", "number"),
                      help="style for support information."  )
    parser.add_option( "--error-style", dest="error_style", type="choice",
                       choices=("pie", "number"),
                      help="style for error information."  )
    parser.add_option( "--branch-scale", dest="branch_scale", type="float",
                      help="branch length scale factor."  )
    parser.add_option( "--height-scale", dest="height_scale", type="float",
                      help="height scale factor."  )
    parser.add_option( "-a", "--annotations", dest="annotations", type="choice", action="append",
                       choices=("support", "error", "kaks", "master", "value", "tables" ),
                       help="annotations given by further trees."  )
    parser.add_option("--filename-tables", dest="filename_tables", type="string",
                       help="add tables from file (need also set options -a tables) [%default]")
    parser.add_option("--show-branchlengths", dest="show_branchlengths", action="store_true",
                      help="show branch lengths." )
    parser.add_option("--leaf-symbol", dest="plot_leaf_symbol", type="choice",
                      choices=("square", "circle" ),
                      help="Symbol for leaves." )
    parser.add_option("--font-size-branches", dest="font_size_branches", type="int",
                      help="set font size for branches." )
    parser.add_option("--font-size-tips", dest="font_size_tips", type="int",
                      help="set font size for tips." )
    parser.add_option("--font-style-tips", dest="font_style_tips", type="choice",
                      choices=("normal", "italic",),
                      help="set font style for tips." )
    parser.add_option("--filename-map", dest="filename_map", type="string",
                      help="filename with a name translation table." )
    parser.add_option("--filename-map-species2colour", dest="filename_colour_map", type="string",
                      help="filename with a map of species to colour." )
    parser.add_option("--no-leaf-labels", dest="plot_leaf_labels", action="store_false",
                      help="do not show labels at leafs." )
    parser.add_option("--no-ruler", dest="plot_ruler", action="store_false",
                      help="do not plot ruler." )
    
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
        support_style = None,
        error_style = "number",
        kaks_style = "number",
        annotations = None,
        show_branchlengths = False,
        branch_length_format = "%5.2f",
        font_size_tips = None,
        font_size_branches = None,
        font_style_tips = None,
        filename_map = None,
        filename_colour_map = None,
        plot_leaf_labels = True,
        plot_leaf_symbol = None,
        plot_ruler = True,
        filename_tables = None,
        )

    (options, args) = E.Start( parser, add_pipe_options = True )
    
    if options.filename_tree:
        tree_lines = open(options.filename_tree, "r").readlines()
    elif options.tree:
        tree_lines = options.tree
    else:
        tree_lines = sys.stdin.readlines()

    nexus = TreeTools.Newick2Nexus( tree_lines )
    master_tree = nexus.trees[0]

    if options.filename_map:
        map_names = IOTools.ReadMap( open(options.filename_map, "r") )

        for id, node in master_tree.chain.items():
            if node.data.taxon in map_names:
                node.data.taxon = map_names[node.data.taxon]

    if options.loglevel >= 2:
        master_tree.display()
    
    plot = SVGTree.SVGTree( master_tree )

    if options.branch_scale:
        plot.setBranchScale( options.branch_scale )
        
    if options.height_scale != None:
        plot.setHeightScale( options.height_scale )    

    if options.font_size_tips != None:
        plot.setFontSize( options.font_size_tips )

    if options.plot_ruler == False:
        plot.setRulerElements( [] )

    if options.show_branchlengths:
        b = SVGTree.BranchDecoratorHorizontalBranchLength( master_tree ) 
        if options.font_size_branches:
            b.setFontSize( options.font_size_branches )
        plot.setDecoratorHorizontalBranches( b )

    if options.colour_by_species:
        if options.filename_colour_map:
            map_species2colour = IOTools.ReadMap( open( options.filename_colour_map, "r"))
        else:
            map_species2colour = None
            
        rx = re.compile(options.species_regex)
        extract_species = lambda x: rx.search( x ).groups()[0]
        plot.setDecoratorExternalNodes( SVGTree.NodeDecoratorBySpecies( master_tree,
                                                                        plot_symbol = options.plot_leaf_symbol,
                                                                        plot_label = options.plot_leaf_labels,
                                                                        map_species2colour = map_species2colour,
                                                                        extract_species= extract_species ) )

    if options.font_style_tips:
        plot.getDecoratorExternalNodes().setFontStyle( options.font_style_tips )

    plot.getDecoratorExternalNodes().setPlotLabel( options.plot_leaf_labels )

    current_tree = 1
    
    ## add annotations by further trees given on the command line
    branch_length_annotations = []

    current_reference_tree = master_tree
    
    if options.annotations:
        for annotation in options.annotations:
            
            tree = nexus.trees[current_tree]
            
            if annotation == "support":

                tree.branchlength2support()
                for id, node in tree.chain.items():
                    node.data.branchlength = 1.0

                if options.support_style == "pie":
                    plot.setDecoratorInternalNodes( NodeDecoratorSupportPieChart( nexus.trees[current_tree] ) )
                    
            elif annotation == "error":

                if options.error_style == "number":
                    b = SVGTree.BranchDecoratorHorizontalBranchLengthError( current_reference_tree, tree )
                    if options.font_size_branches:
                        b.setFontSize( options.font_size_branches )
                    branch_length_annotations.append( b )

            elif annotation == "kaks":

                if options.kaks_style == "number":
                    b = SVGTree.BranchDecoratorHorizontalBranchLengthWithKaks( current_reference_tree, tree )
                    if options.font_size_branches:
                        b.setFontSize( options.font_size_branches )
                    branch_length_annotations.append( b )                        

            elif annotation == "value":

                b = SVGTree.BranchDecoratorHorizontalBranchLength( tree )
                if options.font_size_branches:
                    b.setFontSize( options.font_size_branches )
                branch_length_annotations.append( b )                        

            elif annotation == "master":
                current_reference_tree = tree

            elif annotation == "tables":
                b = BranchDecoratorTable( tree, filename= options.filename_tables )
                plot.setDecoratorHorizontalBranches( b )
                
            current_tree += 1

        if len(branch_length_annotations) == 1:
            b = branch_length_annotations[0]
        elif len(branch_length_annotations) == 2:
            b1, b2 = branch_length_annotations
            b1.setFontColour( SVGTree.BLUE )
            b2.setFontColour( SVGTree.RED )
            b = SVGTree.BranchDecoratorHorizontalAboveBelow( master_tree, b1, b2 )
        elif len(branch_length_annotations) > 2:
            raise "obtained more than three branch length annotations. Layout not implemented"

        plot.setDecoratorHorizontalBranches( b )
    
    plot.initializePlot()
    
    plot.writeToFile(sys.stdout)
    
    E.Stop()


if __name__ == "__main__":
    sys.exit(main())
    

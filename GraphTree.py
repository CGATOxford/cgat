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
GraphTree.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import sys

from Graph import *

import Phylolib

FORMAT = """
"""

def SortAByB( a, b):

    c = map( lambda x,y: (y, x), a, b)
    c.sort()
    return map(lambda x: x[1], c)

##-------------------------------------------------------------------------------
class GraphTree( Graph ):

    def __init__(self, tree, order = None, info = None, colours = None ):

        Graph.__init__(self,FORMAT)

        num_leaves = tree.getNumLeaves()
        
        if not order: order = tree.getNodesDepthFirstFinish()
        
        last_height = 0.0
        vertical_order = 1

        # sort arrays according to heights
        order = list(order)
        if colours: colours = list(colours)
        if info: info = list(info)

        heights = []
        for index in order:
            heights.append(tree.getHeight( index ))
            
        order = SortAByB(order, heights)

        for index in order:
            
            if colours:
                colour = colours[index]
            else:
                colour = None
                
            if last_height < tree.getHeight(index):
                vertical_order += 2
                last_height = tree.getHeight(index)

            if info:
                info1 = info[index]
            else:
                info1 = None
                
            n = Node( str(index),
                      str(index),
                      vertical_order = vertical_order,
                      info1 = info1,
                      colour = colour
                      )
                
            self.AddNode( n )
            
            if index >= num_leaves:
                e = Edge( str(tree.getLeftChild(index)), str(index) )
                self.AddEdge(e)
                e = Edge( str(tree.getRightChild(index)), str(index) )
                self.AddEdge(e)
                
    
##-------------------------------------------------------------------------------
class GraphTreeExplicit( Graph ):

    def __init__(self,
                 tree,
                 info = None,
                 colours = None ):

        Graph.__init__(self,FORMAT)

        num_leaves = tree.getNumLeaves()
        
        order = tree.getNodesDepthFirstFinish()
        
        coords = [0] * num_leaves * 2

        self.mStepSizeX = 100 # spacing of nodes on baseline
        self.mStepSizeY = 100 # horizontal spacing of nodes
        self.mBaseLineY = 0
        
        current_x = 0
        current_y = self.mBaseLineY + self.mStepSizeY

        for node in order:

            # if it is leave:
            if node < num_leaves:
                coords[node] = (current_x, self.mBaseLineY) 
                current_x += self.mStepSizeX
            else:
                left_child  = tree.getLeftChild(node)
                right_child = tree.getRightChild(node)
                ## print node, coords, left_child, right_child
                coords[node] = ((coords[left_child][0] + coords[right_child][0]) / 2, current_y)
                current_y += self.mStepSizeY

        for index in order:
            
            if colours:
                colour = colours[index]
            else:
                colour = None
                
            if info:
                info1 = info[index]
            else:
                info1 = None
                
            n = Node( str(index),
                      str(index),
                      location = coords[index],
                      info1 = info1,
                      colour = colour
                      )
                
            self.AddNode( n )
            
            if index >= num_leaves:
                e = Edge( str(tree.getLeftChild(index)), str(index) )
                self.AddEdge(e)
                e = Edge( str(tree.getRightChild(index)), str(index) )
                self.AddEdge(e)
                
        
if __name__ == "__main__":

    t = Phylolib.makeTree(5)

    t.joinNodes( 0, 1, 1.0, 1.0)
    t.joinNodes( 2, 3, 1.0, 1.0)
    t.joinNodes( 4, 5, 1.0, 1.0)
    t.joinNodes( 6, 7, 1.0, 1.0)
    
    g = GraphTree( t,
                   order = t.getNodesDepthFirstFinish(),
                   colours=["blue", "green", "red", "yellow", "black", "blue", "green", "red", "yellow"],
                   )
    
    g.Write( sys.stdout )

    










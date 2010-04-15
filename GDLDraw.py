####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: GDLDraw.py 2784 2009-09-10 11:41:14Z andreas $
##
##
####
####


##---------------------------------------------------------------------------
class Node:

    def __init__(self, title,
                 label = None,
                 colour = None,
                 info1 = None,
                 info2 = None,
                 vertical_order = None,
                 horizontal_order = None,
                 shape = None,
                 height = None,
                 width = None,
                 location = None):

        self.mLabel = label
        self.mTitle = title
        self.mInfo1 = info1
        self.mInfo2 = info2
        self.mVerticalOrder = vertical_order
        self.mHorizontalOrder = horizontal_order
        self.mColour = colour
        self.mShape = shape
        self.mHeight = height
        self.mWidth = width
        self.mLocation   = location

    def Write( self, outfile ):
        """write node to outputfile."""

        statement = ' title: "%s" ' % self.mTitle
        
        if self.mLabel != None:
            statement += ' label: "%s" ' % self.mLabel

        if self.mShape: statement += ' shape: %s ' % self.mShape
        if self.mInfo1: statement += ' info1: "%s" ' % self.mInfo1
        if self.mInfo2: statement += ' info2: "%s" ' % self.mInfo2        
        if self.mVerticalOrder: statement += ' vertical_order: %i ' % self.mVerticalOrder
        if self.mHorizontalOrder: statement += ' horizontal_order: %i ' % self.mHorizontalOrder        
        if self.mColour: statement += ' color: %s ' % self.mColour    
        if self.mHeight: statement += ' height: %i ' % self.mHeight
        if self.mWidth: statement += ' width: %i ' % self.mWidth        
        if self.mLocation: statement += ' loc: { x: %i y: %i } ' % self.mLocation
        
        outfile.write( '\tnode: { ' + statement + '}\n')

##---------------------------------------------------------------------------
class Edge:
    
    def __init__( self, source_name, target_name,
                  colour = None,
                  line_style = None,
                  type = 'edge',
                  priority = None,
                  edge_class = None):
        
        self.sourceName = source_name
        self.mTargetName = target_name
        self.mType       = type
        self.mColour     = colour
        self.mLineStyle  = line_style
        self.mPriority   = priority
        self.mEdgeClass  = edge_class
                  
    def Write( self, outfile ):
        """write node to outputfile."""

        statement = ' sourcename: "%s" targetname: "%s" ' % (self.sourceName, self.mTargetName)
        if self.mColour: statement += ' color: %s' % self.mColour
        if self.mLineStyle: statement += ' linestyle: %s' % self.mLineStyle
        if self.mPriority: statement += ' priority: %i' % self.mPriority
        if self.mEdgeClass: statement += ' class: %i' % self.mEdgeClass                
        
        outfile.write( '\t%s: { ' % self.mType + statement + '}\n')
        
##-------------------------------------------------------------------------------
class Graph:

    COLOURS = ('black', 'blue', 'red', 'green', 'yellow', 'cyan', 'magenta',
               'darkblue','lightblue','aquamarine','darkred','lightred','khaki',
               'darkgreen','lightgreen','purple','darkyellow','lightyellow','yellowgreen',
               'darkmagenta','lightmagenta','pink','darkcyan','lightcyan','orange',
               'gold','lilac','orchid','darkgrey','lightgrey','turquoise','white')

    def __init__(self, format = None):
        self.mFormat = format
        self.mEdges = []
        self.mNodes = []

    def Write( self, outfile ):
        """write graph to outputfile."""
        outfile.write( 'graph: {\n')
        if self.mFormat:
            outfile.write( self.mFormat + "\n" )

        for node in self.mNodes:
            node.Write( outfile )
        for edge in self.mEdges:
            edge.Write( outfile )

        outfile.write( '}\n' )

    def AddNode ( self, node ):
        self.mNodes.append(node)

    def AddEdge(self,edge):
        self.mEdges.append(edge)
        
    def Read( self, infile ):
        """read graph from inputfile."""

        self.mNodes = []
        self.mEdges = []
        
        for line in infile.readlines():
            pass
        
    








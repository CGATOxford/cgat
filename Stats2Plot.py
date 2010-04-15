####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Stats2Plot.py 2784 2009-09-10 11:41:14Z andreas $
##
##
####
####


# translate an ascii-table into an html-table

from Pairsdb import *

import sys
import os
import string
import re

from RSDB2PLOT import RSDB2PLOT
import GraphXY
import Plot

class RSDB2PLOT_test( RSDB2PLOT):

    def __init__( self, grap ):
        RSDB2PLOT.__init__( self, graph )

    def ParseColumns( self, columns ):
        return columns
    
    
class RSDB2PLOT_clustering ( RSDB2PLOT ):

    def __init__( self, graph, columns ):
        RSDB2PLOT.__init__( self, graph )
        self.mColumns = columns

    def ParseHeader( self, line ):
        columns = string.split( line[:-1], "\t")
        result = []
        for col in self.mColumns:
            result.append( columns[col] )
        
        self.mGraph.SetLegend( result )

    def ParseColumns( self, columns ):
        result = []
        for col in self.mColumns:
            result.append( columns[col] )
        return result
    
    
    
            
        

    

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
Stats2Plot.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import sys
import os
import string
import re

from Rsdb2Plot import Rsdb2Plot

class Rsdb2Plot_test( Rsdb2Plot):

    def __init__( self, grap ):
        Rsdb2Plot.__init__( self, graph )

    def ParseColumns( self, columns ):
        return columns
    
    
class Rsdb2Plot_clustering ( Rsdb2Plot ):

    def __init__( self, graph, columns ):
        Rsdb2Plot.__init__( self, graph )
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
    
    
    
            
        

    

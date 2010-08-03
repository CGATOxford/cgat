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
Rsdb2Plot.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import sys, os, string

class Rsdb2Plot:

    def __init__ (self, graph):
        self.mGraph = graph
        self.mUseIndex = 1

    def ParseColumns( self, columns ):
        return list(columns)

    def Write( self ):
        self.WritePreamble()            # skip

        self.WriteBody()

        self.WritePostamble()           # skip

    def WritePostamble( self ):
        """write postamble. Stop at '//', ignore '#' at first position."""
        while 1:
            line = sys.stdin.readline()
            
            if not line: break
            if line[0] == "#": continue
            if line[0:2] == "//": break
            
            print line
            
    def WritePreamble( self ):
        """write preamble. Stop at '//', ignore '#' at first position."""
        
        while 1:
            line = sys.stdin.readline()
            
            if not line: break
            if line[0] == "#": continue
            if line[0:2] == "//": break
            
    def ParseHeader(self, line ):
        legend = string.split( line[:-1], "\t")
        self.mGraph.SetLegend( legend )

    def WriteBody( self ):
        """assumes, first line is header, parse until '//', ignore '#'."""

        self.ParseHeader( sys.stdin.readline() )

        num_lines = 0
        total_lines = 0
        sys.stderr.write("parsing..")
        
        while 1:
            line = sys.stdin.readline()
            if not line: break
            if line[0] == "#": continue
            if line[0:2] == "//": break
            
            total_lines = total_lines + 1

            (columns) = string.split( line[:-1], "\t" )

            if not columns:
                break

            columns = self.ParseColumns( columns )

            if self.mUseIndex:
                columns.insert(0, str(total_lines))
            
            self.mGraph.AddPoint( columns )

        sys.stderr.write("done\n")

        









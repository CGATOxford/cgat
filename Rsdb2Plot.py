####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Rsdb2Plot.py 2784 2009-09-10 11:41:14Z andreas $
##
##
####
####

# feed an ascii table to grace

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

        









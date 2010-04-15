####
####
##
## Project PythonTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Rsdb2Html.py 2784 2009-09-10 11:41:14Z andreas $
##
##
####
####

"""
translate an ascii-table into an html-table
"""

import sys
import os
import string
import HTMLgen

class Rsdb2Html:

    def __init__ (self):
        self.mTitle   = ""
        self.mResource = "/homes/heger/Picasso/bin/HTMLgen.rc"
        self.mCaption = ""
        self.mHeading = ("",)
        self.mUseIndex = 1              # number rows
        self.mFullWidth = 0             # use full width
        
    def CreateTable( self ):
        
        self.table = HTMLgen.Table(self.mCaption)
        
        if self.mUseIndex:
            self.table.heading = (('No',) + tuple(self.mHeading))
        else:
            self.table.heading = (self.mHeading)

        if not self.mFullWidth:
            self.table.width = ''
            
        self.table.body = []

    def ParseColumns( self, columns ):
        return list(columns)

    def Write( self ):

        self.doc = HTMLgen.SeriesDocument(self.mResource)
        self.doc.title = self.mTitle
        
        # doc.subtitle = "Matches of Prosite-patterns in syster-clusters"
        print self.doc.html_head() 
        print self.doc.html_body_tag()

        self.WritePreamble()

        self.WriteBody()

        self.WritePostamble()

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
            
            print line
            
    def ParseHeader(self, line ):
        self.mHeading = string.split( line[:-1], "\t")
                
    def WriteBody( self ):
        """assumes, first line is header, parse until '//', ignore '#'."""

        self.ParseHeader( sys.stdin.readline() )

        self.CreateTable()

        num_lines = 0

        total_lines = 0
        sys.stderr.write("parsing..")
        
        while 1:
            line = sys.stdin.readline()
            
            if not line: break
            if line[0] == "#": continue
            if line[0:2] == "//": break
            
            if not (total_lines % 100) and total_lines > 0:
                sys.stderr.write( ".%i.." % total_lines)
                print self.table 
                self.table.body = []

            total_lines = total_lines + 1

            (columns) = string.split( line[:-1], "\t" )
            
            if not columns:
                break

            if self.mUseIndex:
                col_string = [str(total_lines),] + self.ParseColumns( columns )
            else:
                col_string = self.ParseColumns( columns )
                
            self.table.body.append( col_string ) 
        
        print self.table
        
        sys.stderr.write("done\n")

    

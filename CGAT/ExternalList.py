################################################################################
#
#   Tools
#
#   $Id: ExternalList.py 2784 2009-09-10 11:41:14Z andreas $
#
#   Copyright (C) 2007 Andreas Heger
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
"""
ExternalList.py - large disk-based lists
========================================

A list class that serves as a stand-in for python lists but
stores data on the file system.

Implements not all functionality of lists yet. Lists are stored
as tab-separated values for unix sort functionality, so all
elements in the list should have the same type.
"""
import os, sys, string, tempfile

class ListIterator:
    
    def __init__( self, file ):
        self.mFile = file

    def __iter__( self ):
        return self
    
    def next( self ):

        while 1:
            line = self.mFile.readline()
            if not line: raise StopIteration
            if line[0] != "#": break
            
        entry = PredictionParser.PredictionParserEntry( expand = 0 )
        entry.Read(line)
        return entry
    
class ExternalList:

    mMapKey2Field = {}

    def __init__( self ):
        self.mIsOpened = 0
    
    def __del__( self ):
        pass
    
    def open( self, filename = None, mode = "r"):
        if filename: self.mFilename = filename
        self.mFile = open( self.mFilename, mode )
        self.mIsOpened = 1

    def getFileName( self ):
        return self.mFilename
    
    def close( self ):
        self.mIsOpened = 0
        self.mFile.close()
        
    def append( self, entry ):
        self.mFile.write( str(entry) + "\n" )
        
    def sort( self, fields, offset = 0 ):
        """sort file according to criteria."""

        if self.mIsOpened:
            reopen = 1
            self.close()
        else:
            reopen = 0
        
        sort_criteria = []
        for field in fields:
            id, modifier = self.mMapKey2Field[ field ]
            field_id = 1 + id + offset
            sort_criteria.append( "-k%i,%i%s" % (field_id, field_id, modifier) )            
            
        outfile, filename_temp = tempfile.mkstemp()
        os.close(outfile)

        ## empty fields are a problem with sort, which by default uses whitespace to non-whitespace
        ## transitions as field separate. 
        ## make \t explicit field separator.
        statement = "sort -t'\t' %s %s > %s" % (string.join( sort_criteria, " " ), self.mFilename, filename_temp)
        exit_code = os.system(statement)
        if exit_code:
            raise "error while sorting, statement =%s" % statement
        os.system( "mv %s %s" % (filename_temp, self.mFilename))

        if reopen: self.open()
        
    def __iter__( self ):
        self.mFile.seek(0)
        return FileIterator( self.mFile )

if __name__ == "__main__":

    file = PredictionFile()

    file.open( "tmp.predictions" )

    file.sort( ("mQueryToken", "mSbjctToken"), offset = -1 )
    
    for entry in file:
        print str(entry)

        
    

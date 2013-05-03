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
GraphTools.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, time, tempfile, os, subprocess

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

    def __str__(self):
        return str(self.message)

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class RuntimeError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class ExternalGraph:

    mExecutableComponents = "ga_components"

    def __init__( self, filename = None):

        if filename:
           self.mInputFilename = filename
           self.mKeepInput = True
           self.mIsOpen = False
        else:
            x, self.mInputFilename = tempfile.mkstemp()
            os.close(x)
            self.mInputFile = open(self.mInputFilename, "w" )
            self.mIsOpen = True
            self.mKeepInput = False

    def __del__(self):
        
        if not self.mKeepInput:
            if os.path.exists( self.mInputFilename):
                os.remove( self.mInputFilename)
                
    def add_edge( self, node1, node2, weight = None ):
        """add edge to graph.
        """

        if self.mIsOpen:
            if weight != None:
                self.mInputFile.write( "%s\t%s\t%s\n" % (map(str,(node1,node2, weight))))
            else:
                self.mInputFile.write( "%s\t%s\n" % (str(node1),str(node2)))
        else:
            raise InputError( "Graph has already been finalized.")

    def finalize( self ):
        """
        """
        if self.mIsOpen:
            self.mInputFile.close()
            self.mIsOpen = False
            
    def connected_components( self ):
        """calculate connected components for graph."""

        self.finalize()        

        statement = "%s %s " % (self.mExecutableComponents, self.mInputFilename)
        
        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              close_fds = True)                              

        (out, err) = s.communicate()


        if s.returncode != 0:
            raise RuntimeError("Error %s while executing statement\n%s" % (err, statement))

        keep = False
        
        map_cluster2members = {}

        rx2 = re.compile("# total number of components")
        rx1 = re.compile("# token\tvertex")
        for line in out.split("\n"):
            if rx1.match(line):
                keep = 1
            elif rx2.match(line):
                break
            elif keep:
                id, cluster = line.split("\t")[:2]
                if cluster not in map_cluster2members:
                    map_cluster2members[cluster] = []
                map_cluster2members[cluster].append(id)

        return map_cluster2members.values()
                
        
        












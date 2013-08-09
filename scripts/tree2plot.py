################################################################################
#   Gene prediction pipeline 
#
#   $Id: tree2plot.py 2782 2009-09-10 11:40:29Z andreas $
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
"""
tree2plot.py - plot a tree as postscript (deprecated)
=====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Plot a tree in postscript format. This script has been replaced
by :mod:`tree2svg`.

Usage
-----

Example::

   python <script_name>.py --help

Type::

   python <script_name>.py --help

for command line help.

Documentation
-------------

Code
----

""" 

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

TEMPLATE = """
\\begindef

    \\name{tree0} % imported from NEXUS file
    \\version{1.0} % version of treegraph file format

    %
    % tree geometry
    %
    \\width{150}  % width in mm
    \\height{250}  % height in mm
    %     \\autolength  % ignore len{} in nodes
    %     \\variable
    
    %
    % layout
    %
    \\roundness{0.2}  % between 0.0 and 1.0
    \\thickness{0.3}  % in mm
    \\margin{0}{0}{0}{0} % left top right bottom, in mm
    \\style{default}{plain}{9.6}
    \\style{r}{plain}{12}
    
    %     \\style{u1}{plain}{9.6}
    %     \\style{d1}{italic}{9.6}
    %     \\style{u2}{bold}{8.4}
    %     \\style{br}{plain}{14.4}
    %     \\rulepos{5}{20}
    % position of rule
    %     \\separator{:}
    %     \\proof
    
    %
    % brackets
    %
    %     \\bracket{0}{n1}{n2}{example text}
    %     \\brace*{1.5}{root}{root}{example text}
    
    \\enddef
"""

class TreeGraphError(Exception):
    pass

def to_string(tree,
              branchlengths="d1",
              support="d2",
              length_format="%6.4f",
              support_format="%1.5f"):
    """Return a paup compatible tree line.

    branchlengths go into u2
    support goes into d2
    """
    def make_info_string(data,terminal=False):
        s = '\\len{%f}\n' % data.branchlength
        
        if branchlengths: # write only branchlengths, ignore support
            s += '\\%s{%s}\n' % (branchlengths, length_format % data.branchlength)
        if support:
            ## no support for terminal nodes
            if not terminal:
                s += '\\%s{%s}\n' % (support, length_format % data.support)            
        return s
    
    def newickize(node):
        """Convert a node tree to a newick tree recursively."""

        if not tree.node(node).succ:    #terminal
            return '\\r{%s}\n%s' % (tree.node(node).data.taxon, make_info_string(tree.node(node).data,terminal=True))
        else:
            return '%s(\n%s\n)\n' % (make_info_string(tree.node(node).data), ',\n'.join(map(newickize,tree.node(node).succ)))
        return subtree

    treeline='%s\n\label{root}\n(\n%s\n)\n' % (
        TEMPLATE,
        ','.join(map(newickize,tree.node(tree.root).succ)))
    
    return treeline 

class TreeGraph:

    mExecutable = "/net/cpp-group/bin/tgf"

    def __init__( self,
                  options="",
                  branchlengths = "u1",
                  support = "d2",
                  loglevel = 0):
        
        self.mOptions = options
        self.mBranchLengths = branchlengths
        self.mSupport = support
        self.mLogLevel = loglevel
        
    def CreateTemporaryFiles( self ):
        """create temporary files."""
        self.mTempDirectory = tempfile.mkdtemp()
        
        self.mFilenameTempInput = self.mTempDirectory + "/input"
        self.mFilenameTempOutput = self.mTempDirectory + "/input.svg"

    def DeleteTemporaryFiles( self ):
        """clean up."""
        os.remove( self.mFilenameTempInput )
        os.remove( self.mFilenameTempOutput )
        os.rmdir( self.mTempDirectory )

    def SetStderr( self, file = None):
        """set file for dumping stderr."""
        self.mStderr= file
        
    def WriteOutput( self, lines, filename_output = None):
        """write output to file.

        If file is not given, lines are written to stdout.
        """
        
        if filename_output:
            outfile = open(filename_output, "w")
        else:
            outfile = sys.stdout
            
        outfile.write( string.join( lines, "") )

        if filename_output:
            outfile.close()

    def ParseResult( self, trace_file = None, information_file = None ):

        result.Read( trace_file, information_file )
        
        return result
    
    def Run( self, tree ):

        self.CreateTemporaryFiles()

        tempfile = open(self.mFilenameTempInput, "w")        
        tempfile.write( to_string(input_tree,
                                  branchlengths = self.mBranchLengths,
                                  support = self.mSupport) )
        tempfile.close()

        if self.mLogLevel >= 2:
            os.system("cat %s" % self.mFilenameTempInput)
        
        statement = string.join( ( self.mExecutable,
                                   "-v",
                                   self.mFilenameTempInput ),
                                 " ")


        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = self.mTempDirectory,
                              close_fds = True)                              

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise TreeGraphError, "Error in calculating svg file\n%s" % err

        d = open(self.mFilenameTempOutput).readlines()

        self.DeleteTemporaryFiles()
        
        return "".join(d)


if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: tree2plot.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.set_defaults(
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    lines = filter( lambda x: x[0] != "#", sys.stdin.readlines())

    nexus = TreeTools.Newick2Nexus( lines )
    
    input_tree = nexus.trees[0]
    
    treegraph = TreeGraph( support = None, loglevel = options.loglevel)

    print treegraph.Run( input_tree )

    E.Stop()

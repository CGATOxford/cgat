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
cgat_make2help.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python cgat_make2help.py --help

Type::

   python cgat_make2help.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import optparse
import time

"""extract help test from a makefile
"""

import CGAT.Experiment as E
from CGAT.Local import getMakefiles, getScripts, getModules

class Parameter:
    def __init__(self, name, comment, default_value = "na"):
        self.mName = name
        self.mComment = comment
        self.mDefaultValue = value
        self.mValue = None        

    def printPretty( self ):
        lines = []

        if not self.mValue or self.mDefaultValue == self.mValue:
            value = "%s [default]" % self.mDefaultValue
        else:
            value = "%s [default=%s]" % (self.mValue, self.mDefaultValue )

        lines.append( "%s = %s" % (self.mName, value) )
        
        if self.mComment:
            lines.append( "" )
            lines.append( "  " + "\n  ".join( self.mComment ) )
        return "\n".join(lines)

class Target:
    def __init__(self, name, dependencies, comment, ):
        self.mName = name
        self.mComment = comment
        self.mDependencies = dependencies

    def printPretty( self ):
        lines = []
        
        lines.append( "%s" % self.mName )
        if self.mDependencies:
            lines.append( "" )
            lines.append( "  depends on: %s" % self.mDependencies )
        if self.mComment:
            lines.append( "" )
            lines.append( "  " + "\n  ".join( self.mComment ) )
        return "\n".join(lines)
            

if __name__  == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: cgat_make2help.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option( "-m", "--method", dest="method", type="choice",
                       help="method to use [t-test=t-test]",
                       choices=( "t-test", ) )
    parser.set_defaults(
        filename="Makefile",
        )
    
    (options, args) = E.Start( parser,
                                        quiet = True,
                                        add_pipe_options = True )

    src_dir = os.path.dirname( os.path.abspath( __file__  ) )

    makefiles = getMakefiles( (options.filename,), source_directory = src_dir, ignore_missing = True )

    ## variables starting with "PARAM"
    parameters = {}
    ## main targets: no patterns or no dependencies
    targets = []
    ## contain patterns and dependencies
    rules = []

    options.stdout.write( "Help for makefile %s\n\n" % os.path.abspath( options.filename ) )

    options.stdout.write( """This file is organized in three sections.

1. Primary Targets

   Primary targets are targets to be run by the users. Usually includes "all".

2. Parameters

   Parameters defined in the makefile.

3. Secondary targets and rules

   Rules defined in the makefile to build the primary targets.

""")

    for makefile in makefiles:
        
        infile = open(makefile, "r" )
        comment = []
        last_line = None
        
        first = True
        for line in infile:

            if last_line: 
                line = last_line + line
                last_line = None

            if not re.sub( "\s+", "", line): 
                first = False
                continue

            if line[-2] == "\\":
                last_line = line[:-2]
                continue
            
            if line[0] == "\t": continue

            if line[0] == "#":
                l = re.sub( "^#+\s*", "", line[:-1] )
                if l and not first:
                    comment.append( l )
                continue

            x = re.match("(PARAM_\S+[^?:])\s*\?*\s*=(.+)", line )
            if x:
                param,value = x.groups()
                if param not in parameters:
                    p = Parameter( param, comment, value )
                    parameters[param] = p 
                else:
                    p = parameters[param]
                    p.mComment += comment
                    p.mDefaultValue = value

            x = re.match("(PARAM_\S+[^?])\s*=(.+)", line )
            if x:
                param,value = x.groups()
                if param not in parameters:
                    p = Parameter( param, comment )
                    parameters[param] = p 
                else:
                    p = parameters[param]
                    p.mComment += comment
                    
                p.mValue = value

            x = re.match( "(\S+)\s*:\s*(.*)", line )
            if x:
                target, dependencies = x.groups()
                if "=" in target or "=" in dependencies: continue
                if target[0] in ".$": continue
                if "hook" in target: continue
                dependencies = re.sub( "\s+", " ", dependencies )
                if "." in target:
                    rules.append( Target( target, dependencies, comment ) )
                else:
                    targets.append( Target( target, dependencies, comment ) )

            comment = []
            first = False

        infile.close()

    options.stdout.write("\n\n Section 1: Primary targets\n" )
    options.stdout.write(" --------------------------\n\n" )
    options.stdout.write(""" Primary targets are targets to be run by the users. Usually includes "all".\n\n""")

    for x in targets:
        options.stdout.write( x.printPretty() + "\n\n" )

    options.stdout.write("\n\n Section 2: Parameters\n" )
    options.stdout.write(" --------------------------\n\n" )
    options.stdout.write(""" A list of parameters defined in the makefile.\n\n""")

    for x, y in sorted(parameters.items()):
        options.stdout.write( y.printPretty() + "\n\n\n" )        

    for x in targets:
        options.stdout.write( x.printPretty() + "\n\n" )

    options.stdout.write("\n\n Section 3: Secondary targets and rules\n" )
    options.stdout.write(" --------------------------\n\n" )
    options.stdout.write(""" Secondary targets are run by the pipeline.\n\n""")

    for x in rules:
        options.stdout.write( x.printPretty() + "\n\n" )
        
    options.stdout.write("\n\n Section 4: Contents\n" )
    options.stdout.write(" --------------------------\n\n" )
    options.stdout.write(" This help was created by scanning the contents of the following makefiles.\n" )

    for makefile in makefiles:
        options.stdout.write("    %s\n" % makefile )
        
    E.Stop()

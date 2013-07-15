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
analyze_logfiles.py - create summary from logfiles
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts are collection of logfiles and collates
summary information.

This script uses the ``# job finished`` tag that
is added by scripts using the module :mod:`Experiment`.

Usage
-----

Example::

   python analyze_logfiles.py --help

Type::

   python analyze_logfiles.py --help

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
import tempfile
import subprocess
import optparse
import gzip
import glob

from types import *

import CGAT.Experiment as Experiment
import CGAT.IOTools as IOTools

class LogFileData:
    mRegex = re.compile("# job finished in (\d+) seconds at (.*) --\s+([.\d]+)\s+([.\d]+)\s+([.\d]+)\s+([.\d]+)")
    mFormat = "%6.2f"
    mDivider = 1.0
    
    def __init__(self):

        self.mWall = 0
        self.mUser = 0
        self.mSys = 0
        self.mChildUser = 0
        self.mChildSys = 0
        self.mNChunks = 0

    def add( self, line ):

        if not self.mRegex.match(line): return
        t_wall, date, t_user, t_sys, t_child_user, t_child_sys = self.mRegex.match(line).groups()

        self.mWall += int(t_wall)
        self.mUser += float(t_user)
        self.mSys += float(t_sys)
        self.mChildUser += float(t_child_user)
        self.mChildSys += float(t_child_sys)
        self.mNChunks += 1

    def __getitem__( self, key ):
        
        if key == "wall":
            return self.mWall
        elif key == "user":
            return self.mUser
        elif key == "sys":
            return self.mSys
        elif key == "cuser":
            return self.mChildUser
        elif key == "csys":
            return self.mChildSys
        elif key == "nchunks":
            return self.mNChunks
        else:
            raise ValueError("key %s not found" % key)

    def __add__(self, other ):

        self.mWall += other.mWall
        self.mUser += other.mUser
        self.mSys += other.mSys
        self.mChildUser += other.mChildUser
        self.mChildSys += other.mChildSys
        self.mNChunks += other.mNChunks

        return self
    
    def __str__(self):
        
        return "%i\t%s" % (
            self.mNChunks,
            "\t".join( map( lambda x: self.mFormat % (float(x)/self.mDivider), \
                                (self.mWall, self.mUser, self.mSys,
                                 self.mChildUser, self.mChildSys ) ) ) )

    def getHeader(self):
        return "\t".join( ("chunks", "wall", "user", "sys", "cuser", "csys") )


class LogFileDataLines(LogFileData):
    """record lines."""
    def __init__(self):
        LogFileData.__init__(self)
        self.mNLines = 0

    def add( self, line ):
        if line[0] != "#":
            self.mNLines += 1
        else:
            return LogFileData.add( self, line )
    def __getitem__( self, key ):
        if key == "lines":
            return self.mNLines
        else:
            return LogFileData.__getitem__( self, key )
    def __add__(self, other ):
        self.mNLines += other.mNLines
        return LogFileData.__add__( self, other )

    def __str__(self):
        return "%s\t%i" % (LogFileData.__str__(self), self.mNLines )
                          
    def getHeader(self):
        return "%s\t%s" % (LogFileData.getHeader(self), "lines" )
    
if __name__ == "__main__":
    
    parser = E.OptionParser( version = "%prog version: $Id: analyze_logfiles.py 2781 2009-09-10 11:33:14Z andreas $" )

    parser.add_option( "-g", "--glob", dest="glob_pattern", type="string" ,
                       help="glob pattern to use for collecting files [%default].")

    parser.add_option( "-f", "--file-pattern", dest="file_pattern", type="string",
                       help="only check files matching this pattern [%default]." )

    parser.add_option( "-m", "--mode", dest="mode", type="choice",
                       choices = ("file", "node" ),
                       help="analysis mode [%default]." )

    parser.add_option( "-r", "--recursive", action="store_true",
                       help="recursively look for logfiles from current directory [%default]." )

    parser.set_defaults(
        truncate_sites_list = 0,
        glob_pattern = "*.log",
        mode = "file",
        recursive = False,
        )

    (options, args) = Experiment.Start( parser )

    if args:
        filenames = args
    elif options.glob_pattern:
        filenames = glob.glob( options.glob_pattern )

    if len(filenames) == 0:
        raise "no files to analyse"

    if options.mode == "file":
        totals = LogFileData()

        options.stdout.write( "file\t%s\n" % totals.getHeader() )

        for filename in filenames:
            if filename == "-":
                infile = sys.stdin
            elif filename[-3:] == ".gz":
                infile = gzip.open( filename, "r" )
            else:
                infile = open(filename, "r" )

            subtotals = LogFileData()
            for line in infile:
                subtotals.add( line )

            infile.close()

            options.stdout.write( "%s\t%s\n" % (filename, str(subtotals) ) )
            totals += subtotals

        options.stdout.write( "%s\t%s\n" % ("total", str(totals) ) )

    elif options.mode == "node":
        
        chunks_per_node = {}

        rx_node = re.compile( "# job started at .* \d+ on (\S+)" )

        for filename in filenames:
            if filename == "-":
                infile = sys.stdin
            elif filename[-3:] == ".gz":
                infile = gzip.open( filename, "r" )
            else:
                infile = open(filename, "r" )
                
            data = LogFileDataLines()
               
            for line in infile:
                
                if rx_node.match(line):
                    node_id = rx_node.match(line).groups()[0]
                    data = LogFileDataLines()
                    if node_id not in chunks_per_node:
                        chunks_per_node[node_id] = []
                    chunks_per_node[node_id].append( data )
                    continue

                data.add( line )

        options.stdout.write( "node\t%s\n" % data.getHeader() )
        total = LogFileDataLines()
                   
        for node, data in sorted(chunks_per_node.items()):
            subtotal = LogFileDataLines()
            for d in data:
                # options.stdout.write( "%s\t%s\n" % (node, str(d) ) )  
                subtotal += d

            options.stdout.write( "%s\t%s\n" % (node, str(subtotal) ) )  

            total += subtotal
            
        options.stdout.write( "%s\t%s\n" % ("total", str(total) ) )  

    Experiment.Stop()
    

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
matrix2tree.py - build tree from a distance matrix
==================================================

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

   python matrix2tree.py --help

Type::

   python matrix2tree.py --help

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
import getopt
import time
import optparse
import math
import tempfile
import subprocess
import random

from types import *

import CGAT.Experiment as E
import CGAT.TreeTools as TreeTools
import CGAT.Tree as Tree
import CGAT.IOTools as IOTools
import CGAT.WrapperPhylip as WrapperPhylip

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: matrix2tree.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-i", "--invert-map", dest="invert_map", action="store_true",
                      help="""invert map.""")

    parser.add_option( "--input-format", dest="input_format", type="choice",
                       choices=("phylip", "full") ,
                       help="""input format.""")

    parser.add_option("-t", "--filename-tree", dest="filename_tree", type="string",
                      help="""filename with tree to fit.""")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("nj", "kitsch", "fitch"),
                      help="""algorithm to run.""")

    parser.add_option("-e", "--replicates", dest="replicates", action="store_true",
                      help="replicates."  )

    parser.add_option("-r", "--root", dest="root", action="store_true",
                      help="midpoint root (if it is not rooted)."  )

    parser.add_option("-u", "--unroot", dest="unroot", action="store_true",
                      help="unroot tree (if it is rooted)."  )

    parser.add_option( "--skip-separators", dest="write_separators", action="store_false",
                      help="do not echo separators (starting with >)")
    
#    parser.add_option("-i", "--iterations", dest="iterations", type="int",
#                      help="number of iterations." )

    parser.add_option("-p", "--power", dest="power", type="float",
                      help="power." )

    parser.add_option( "--prune-tree", dest="prune_tree", action="store_true",
                       help="prune tree such to include only taxa which are part of the input matrix.")

    parser.add_option( "--add-random", dest="add_random", action="store_true",
                       help="add small random value to off-diagonal zero elements in matrix.")

    parser.add_option( "--pseudo-replicates", dest="pseudo_replicates", action="store_true",
                       help="add small random value to off-diagonal zero elements in matrix, even if they have no replicates.")

    parser.add_option( "--debug", dest="debug", action="store_true",
                       help="dump debug information." )


    parser.set_defaults(
        value = 0,
        method = "nj",
        input_format = "phylip",
        filename_tree = None,
        outgroup = None,
        replicates = False,
        root = False,
        unroot = False,
        power = 0,
        write_separators = True,
        prune_tree = False,
        add_random = False,
        debug = False,
        )

    (options, args) = E.Start( parser, add_pipe_options = True )

    phylip = WrapperPhylip.Phylip()

    if options.debug:
        phylip.setLogLevel( options.loglevel )
    
    phylip.setPruneTree( options.prune_tree )

    lines = filter( lambda x: x[0] != "#", sys.stdin.readlines())

    chunks = filter( lambda x: lines[x][0] == ">", range(len(lines)) )

    if not chunks:
        options.write_separators = False
        chunks = [-1]
        
    chunks.append( len(lines) )

    for x in range(len(chunks) -1 ):
        
        matrix = lines[chunks[x]+1:chunks[x+1]]

        ## parse phylip matrix
        if options.add_random:
            mm = []
            ids = []
            for l in range(1, len(matrix)):
                values = re.split("\s+", matrix[l][:-1])
                ids.append( values[0] )
                mm.append( map( lambda x: x.strip(), values[1:] ) )
                
            d = len(mm)
            if options.replicates:
                for row in range(d-1):
                    for col in range(row + 1, d):
                        cc = col * 2
                        rr = row * 2
                        if mm[row][cc] == "0" and mm[row][cc+1] != "0":
                            mm[row][cc+1] = "1"
                            mm[col][rr+1] = "1"
                            v = str(random.random() / 10000.0 )
                            mm[row][cc] = v
                            mm[col][rr] = v

            else:
                for row in range(d-1):
                    for col in range(row + 1, d):
                        if mm[row][col] == "0":
                            v = str(random.random() / 10000.0 )
                            mm[row][col] = v
                            mm[col][row] = v
                    
            matrix = ["%i\n" % d]
            for row in range(d):
                matrix.append( ids[row] + "    " + "    ".join(mm[row]) + "\n" )

        ## parse phylip matrix
        if options.pseudo_replicates:
            mm = []
            ids = []
            for l in range(1, len(matrix)):
                values = re.split("\s+", matrix[l][:-1])
                ids.append( values[0] )
                mm.append( map( lambda x: x.strip(), values[1:] ) )
                
            d = len(mm)
            if options.replicates:
                for row in range(d-1):
                    for col in range(row + 1, d):
                        cc = col * 2
                        rr = row * 2
                        if mm[row][cc+1] == "0":
                            mm[row][cc+1] = "1"
                            mm[col][rr+1] = "1"
                            v = str(random.random() / 10000.0 )
                            mm[row][cc] = v
                            mm[col][rr] = v
                        else:
                            mm[row][cc+1] = "100"
                            mm[col][rr+1] = "100"
            else:
                for row in range(d-1):
                    for col in range(row + 1, d):
                        if mm[row][col] == "0":
                            v = str(random.random() / 10000.0 )
                            mm[row][col] = v
                            mm[col][row] = v
                    
            matrix = ["%i\n" % d]
            for row in range(d):
                matrix.append( ids[row] + "    " + "    ".join(mm[row]) + "\n" )

        phylip.setMatrix( matrix )
        
        phylip_options = []

        if options.filename_tree:
            nexus = TreeTools.Newick2Nexus( open(options.filename_tree, "r") )
            ref_tree = nexus.trees[0]
            phylip.setTree( ref_tree )
            phylip_options.append( "U" )
        else:
            ref_tree = None
            
        if options.method == "nj":
            phylip.setProgram( "neighbor")
            
        elif options.method == "fitch":
            phylip.setProgram( "fitch")
            
        elif options.method == "kitsch":        
            phylip.setProgram( "kitsch")

        if options.replicates:
            phylip_options.append("S")

        if options.power > 0:
            phylip_options.append("P")
            phylip_options.append("%f" % options.power)        

        phylip_options.append("Y")        

        phylip.setOptions( phylip_options )

        result = phylip.run()

        ## root with outgroup
        if options.root:        
            if options.outgroup:
                pass
            # midpoint root
            else:
                for tree in result.mNexus.trees:
                    tree.root_midpoint()
                    
        ## explicitely unroot
        elif options.unroot:
            phylip.setOptions( ("Y", "W", "U", "Q" ) )
            phylip.setProgram( "retree" )
            for x in range( len(result.mNexus.trees)):
                phylip.setTree( result.mNexus.trees[x])
                xresult = phylip.run()
                result.mNexus.trees[x] = xresult.mNexus.trees[0]

        if options.write_separators:
            options.stdout.write( lines[chunks[x]] )

        if result.mNexus:
            options.stdout.write( TreeTools.Nexus2Newick( result.mNexus ) + "\n" )

        if options.loglevel >= 1:
            if ref_tree:
                nref = len(ref_tree.get_terminals() )
            else:
                nref = 0
            for tree in result.mNexus.trees:
                options.stdlog.write("# ninput=%i, nreference=%i, noutput=%i\n" % ( len(matrix)-1, nref, len(tree.get_terminals()) ) )

    E.Stop()

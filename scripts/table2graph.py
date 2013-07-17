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
"""
convert adjacency lists to edge lists
=====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

convert a graph in table format into a graph 
with links.

A table in graph format stores adjacency lists for each node, for example::

   node1   node2;node3;node4   value1;value2;value3

This will be converted to::

   node1   node2   value1
   node1   node3   value2
   node1   node4   value3

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os
import sys
import re
import optparse

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", usage = globals()["__doc__"] )
    
    parser.add_option("-e", "--headers", dest="headers", type="string",
                      help="',' separated list of node headers [default=%default]."  )

    parser.set_defaults(
        headers="node1,node2",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ## do sth
    ninput, nskipped, noutput = 0, 0, 0

    first = True
    for line in options.stdin:
        if line.startswith( "#" ): continue
        data = line[:-1].split("\t")
        if first:
            headers = options.headers.split(",")
            if len(data) >= 2: extra = "\t%s" % ("\t".join(data[2:]))
            else: extra = ""
            options.stdout.write("%s%s\n" % ("\t".join(headers), extra))
            first = False
            continue

        ninput += 1
        if len(data) < 2: continue
        values = [ x.split(";") for x in data[1:] if x != ""]
        if len(values) == 0: 
            nskipped += 1
            continue

        l = min([len(x) for x in values])
        assert l == max([len(x) for x in values]), "unequal number of fields in line '%s'" % line[:-1]
        node1 = [[data[0]]*l]
        for n in zip( *(node1 + values) ):
            options.stdout.write("\t".join(n) + "\n" )
            noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

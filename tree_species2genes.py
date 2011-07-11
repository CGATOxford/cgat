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
tree_species2genes.py - collapse leaves with the same species.
==============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

-p, --pattern-species=          regex pattern to extract species from identifier
-g, --genes=                    filename with list of genes per species

Usage
-----

Example::

   python tree_species2genes.py --help

Type::

   python tree_species2genes.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os, sys, string, re, getopt, tempfile, time, popen2

from Bio.Nexus import Nexus
from Bio.Nexus.Nodes import Node

import Experiment
import BlastAlignments
import alignlib
import Genomics

param_loglevel = 1

param_long_options=["verbose=", "help",
                    "pattern-species=",
                    "genes="]

param_short_options="v:hp:g:"

param_pattern_species = "^([^@:]+)[@:]"

param_filename_genes = None

if __name__ == "__main__":

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)
        elif o in ("-p", "--pattern-species"):
            param_pattern_species = a
        elif o in ("-g", "--genes"):
            param_filename_genes = a

    if not param_filename_genes:
        print "please specify a filename with genes"
        print USAGE
        sys.exit(1)
    
    print Experiment.GetHeader()
    print Experiment.GetParams()

    rx = re.compile( param_pattern_species )
    infile = open(param_filename_genes, "r")
    map_species2genes = {}
    for line in infile:
        if line[0] == "#": continue
        gene = line[:-1].split("\t")[0]
        species = rx.search( gene ).groups()[0]
        if species not in map_species2genes:
            map_species2genes[species] = []
        map_species2genes[species].append( gene )
    
    nexus = TreeTools.Newick2Nexus( sys.stdin )

    TreeTools.Species2Genes( nexus, map_species2genes )

    print TreeTools.Nexus2Newick( nexus )
    
    print Experiment.GetFooter()

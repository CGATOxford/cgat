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
analyze_codonbias_orthology.py - 
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

   python analyze_codonbias_orthology.py --help

Type::

   python analyze_codonbias_orthology.py --help

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
import math

"""analyze codonbias of orthologs.

The input looks like:

schema1 gene1 cai1 schema2, gene2, cai2
"""

import CGAT.Experiment as E

def GetThreshold( dbhandle, schema ):
    """get dominant set threshold for schema."""

    statement = "SELECT cai FROM %s.codonbias WHERE is_selected = TRUE ORDER BY cai DESC";
    
    ## get thresholds
    cc = dbhandle.cursor()
    cc.execute(statement)
    r = cc.fetchall()
    cc.close()

    return threshold

##--------------------------------------------------------------------------------------------
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: analyze_codonbias_orthology.py 2781 2009-09-10 11:33:14Z andreas $")
    
    parser.add_option( "--filters", dest="filters", type="string",
                      help="Filters to use for filtering sequences [all|codon1|codon2|codon3|44]." )
    parser.add_option( "--fields", dest="fields", type="string",
                      help="Fields to output [aligned|nunaligned1|nunaligned2|identical|transitions|transversions|jc69|t92]." )

    parser.set_defaults(
        filename_map = None,
        gap_char = "-",
        )
    
    (options, args) = E.Start( parser, add_pipe_options = True,
                                        add_psql_options = True)

    data = map(lambda x: x[:-1].split("\t"), filter( lambda x: x[0] != "#", sys.stdin.readlines()))

    schema1 = data[0][0]
    schema2 = data[0][3]

    dbhandle = pgdb.connect( options.psql_connection )
    
    threshold1 = GetThreshold( dhhandle, schema1)
    threshold1 = GetThreshold( dhhandle, schema2)
    

    
    

    
    E.Stop()

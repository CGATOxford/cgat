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
analyze_go.py - 
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

   python analyze_go.py --help

Type::

   python analyze_go.py --help

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

import CGAT.Experiment as E
import CGAT.Database as Database


def WriteBackground( go_type, options, suffix ):

    statement = """SELECT DISTINCTROW
    gene_stable_id, glook_%s_id, description, olook_evidence_code
    FROM %s.%s_go_%s__go_%s__main
    WHERE glook_%s_id IS NOT NULL
    GROUP BY gene_stable_id, glook_%s_id, description
    ORDER BY gene_stable_id
    """ % (go_type,
           options.database, options.species, go_type, go_type,
           go_type, go_type)
    
    result = dbhandle.Execute(statement).fetchall()

    outfile = open( "%s%s" % (options.prefix, suffix), "w" )
    
    for r in result:
        outfile.write( " : ".join( map(str, r) ) + "\n" )
        
    outfile.close()
    

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: analyze_go.py 309 2005-12-01 15:50:26Z andreas $")

    dbhandle = Database.Database()
    
    parser.add_option("-s", "--species", dest="species", type="string",
                      help="species to use." )

    parser.add_option("-p", "--prefix", dest="prefix", type="string",
                      help="prefix to use for temporary files." )

    parser.set_defaults( species = "dmelanogaster")
    parser.set_defaults( database = "ensembl_mart_31")
    parser.set_defaults( prefix = "dm_go_")    
    
    (options, args) = E.Start( parser, add_mysql_options = True )

    dbhandle.Connect( options )

    WriteBackground( "biol_process", options, "bp" )
    WriteBackground( "cell_location", options, "lm" )
    WriteBackground( "mol_function", options, "fm" )        


        

    

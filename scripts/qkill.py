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
qkill.py - kill jobs in the queue
==================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

kill jobs from the sun grid engine according to certain criteria.

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
import subprocess
import CGAT.Experiment as E
import xml.etree.ElementTree
import cStringIO as StringIO

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-p", "--pattern", dest="pattern", type="string",
                      help="jobs matching `pattern` in their job description will be killed [default=%default]."  )

    parser.add_option("-n", "--dry-run", dest="dry_run", action="store_true",
                      help="do dry run, do not kill [default=%default]."  )

    parser.set_defaults(
        pattern = None,
        dry_run = False,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    output = StringIO.StringIO(subprocess.Popen(["qstat", "-xml"], stdout=subprocess.PIPE).communicate()[0])
    
    tree = xml.etree.ElementTree.ElementTree( file = output )

    ntested = 0
    to_kill = set()

    if options.pattern:
        pattern = re.compile(options.pattern)
    else:
        pattern = None
    
    for x in tree.getiterator( "job_list"):
        ntested += 1
        id = x.find( "JB_job_number" ).text
        name = x.find( "JB_name" ).text
        if pattern and pattern.search( name ):
            to_kill.add( id )
            
    nkilled = len(to_kill)
    if not options.dry_run:
        p = subprocess.Popen(["qdel", ",".join( to_kill)], stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()
    
    E.info( "ntested=%i, nkilled=%i" % (ntested, nkilled) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

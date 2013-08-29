################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
bed.plot.py - plot 
=============================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Create genomic plots in a set of intervals.

Methods implemented

igv 
   plot through an active IGV instance. 

Usage
-----

Example::

   python bed2plot.py < in.bed 

Type::

   python script_template.py --help

for command line help.

Documentation
-------------

Command line options
--------------------

'''

import os
import sys
import re
import glob
import sqlite3

import CGAT.IGV as IGV
import CGAT.Bed as Bed

import CGAT.Experiment as E

def main( argv = sys.argv ):

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id$", usage = globals()["__doc__"] )

    parser.add_option( "-m", "--method", dest="method", type="choice",
                       choices = ("igv",),
                       help = "method to create plots with [%default]" )

    parser.add_option( "-d", "--snapshot-dir", dest="snapshotdir", type = "string",
                       help = "directory to save snapshots in [%default]" )

    parser.add_option( "-f", "--format", dest="format", type="choice",
                       choices = ("png", "eps", "svg"),
                       help = "output file format [%default]" )

    parser.add_option( "-o", "--host", dest="host", type = "string",
                       help = "host that IGV is running on [%default]" )

    parser.add_option( "-p", "--port", dest="port", type = "int",
                       help = "port that IGV listens at [%default]" )

    parser.add_option( "-e", "--extend", dest="extend", type = "int",
                       help = "extend each interval by a number of bases [%default]" )

    parser.add_option( "-x", "--expand", dest="expand", type = "float",
                       help = "expand each region by a certain factor [%default]" )

    parser.set_defaults(
        method = "igv",
        host = '127.0.0.1',
        port = 61111,
        snapshotdir = os.getcwd(),
        extend = 0,
        format = "png",
        expand = 1.0,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    E.info( "connection to session on %s:%s" % (options.host, options.port) )

    E.info( "saving images in %s" % options.snapshotdir)
    igv = IGV.IGV( host = options.host,
                   port = options.port,
                   snapshot_dir = os.path.abspath( options.snapshotdir ) )

    c = E.Counter()
    for bed in Bed.iterator(options.stdin):

        c.input += 1

        # IGV can not deal with white-space in filenames
        name = re.sub("\s", "_", bed.name)
        
        E.info( "going to %s:%i-%i for %s" % (bed.contig, bed.start, bed.end, name))

        start, end = bed.start, bed.end
        extend = options.extend
        if options.expand:
            d = end - start
            extend = max( extend, (options.expand * d - d) // 2 )

        start -= extend
        end += extend

        igv.go( "%s:%i-%i" % (bed.contig, start, end ))

        fn = "%s.%s" % (name, options.format)
        E.info( "writing snapshot to '%s'" % fn )
        igv.save( fn )
    
        c.snapshots +=1 

    E.info( c )
    E.Stop()

if __name__ == "__main__":
    sys.exit( main() )

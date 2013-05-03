################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: adda2coverage.py 2781 2009-09-10 11:33:14Z andreas $
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

:Author: Andreas Heger
:Release: $Id: adda2coverage.py 2781 2009-09-10 11:33:14Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os, sys, re, optparse

import Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: adda2coverage.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-f", "--filename-lengths", dest="filename_lengths", type="string",
                      help="filename with length information [default=%default]."  )

    parser.set_defaults(
        filename_lengths = "test",
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    map_id2length = {}
    
    for line in open(options.filename_lengths, "r"):
        if line.startswith("#"): continue
        if line.startswith("id\t"): continue
        id, length = line[:-1].split()[:2]
        map_id2length[bytes(id)] = int(length)

    E.info("read sequence length information for %i sequences" % len(map_id2length) )

    ## do sth
    ninput, nskipped, noutput = 0, 0, 0

    def iterator_domains( infile ):
        last = None
        for line in infile:
            if line.startswith("#"): continue
            if line.startswith("id\t"): continue
            id, start, end, family = line[:-1].split()
            if id != last:
                if last: yield domains
                domains = []
                last = id
            domains.append( (bytes(id), int(start), int(end), bytes(family) ) )
        yield domains

    options.stdout.write( "id\tcoverage\n" )
        
    for domains in iterator_domains( options.stdin ):
        ninput += 1
        id = domains[0][0]
        if id not in map_id2length:
            nskipped += 1
            E.warn( "length for sequence %s not known" % id )
            continue

        t = sum( [ x[2] - x[1] for x in domains ] )
        options.stdout.write( "%s\t%5.2f\n" % (id, 100.0 * t / map_id2length[id] ) )
            
        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

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
liftover.py - 
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

   python liftover.py --help

Type::

   python liftover.py --help

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

"""liftover coordinates using a liftover formatted file from the ucsc
"""

import CGAT.Experiment as E
import numpy

parser = E.OptionParser( version = "%prog version: $Id: liftover.py 2404 2009-01-12 11:22:36Z andreas $")

def readLiftOver( infile, chromosome,
                  chromosome_size = 250000000,
                  report_step = 1000000 ):
    """read a matrix. There probably is a routine for this in Numpy, which
    I haven't found yet.
    """

    if options.loglevel >= 2:
        print "## started reading mapping information"
        sys.stdout.flush()
        
    map_position   = numpy.zeros( (chromosome_size,), numpy.int) 
    
    ## signed character for chromosme, negative values for negative strands
    map_chromosome = numpy.zeros( (chromosome_size,), numpy.int8 )

    map_id2chromosome = ["",]
    map_chromosome2id = {}
    n = 0
    
    for line in infile:

        n += 1

        if not (n % report_step):
            if options.loglevel >= 2:
                print "# iteration %i" % n
                sys.stdout.flush()
                
        if line[:5] == "chain":
            ( chr_x, size_x, strand_x, first_x, last_x,
              chr_y, size_y, strand_y, first_y, last_y,
              dontknow) = line[:-1].split(" ")[2:]

            if strand_x == "-":
                raise "what shall I do with negative strands?"

            x = int(first_x)
            
            ## revert coordinates for negative strands (it seems that
            ## the mapping file uses reverse coordinates, while liftover
            ## output doesn't)
            ## add 1 to coordinates, because 0 is flag for unmappable.
            if strand_y == "-":
                invert = True
                ## no +1, because already one past current residue (due to open
                ## bracket)
                y = int(size_y) - int(first_y) 
            else:
                invert = False
                y = int(first_y) + 1
            
            if chr_x != chromosome:
                keep = False
            else:
                keep = True
                if options.loglevel >= 3:
                    print "# adding alignment", line[:-1]
                    
            continue
        
        elif line.strip() == "":
            keep = False
            continue
            
        elif keep:
            
            data = map(int, line[:-1].split("\t"))
            
            if len(data) == 3:
                size, increment_x, increment_y = data
            else:
                size, increment_x, increment_y = data[0], 0, 0

            # add position
            if invert:
                map_position[ x:x+size ] = numpy.arrayrange( y, y-size, -1 )
            else:
                map_position[ x:x+size ] = numpy.arrayrange( y, y+size, 1 )

            if chr_y not in map_id2chromosome:
                map_chromosome2id[chr_y] = len(map_id2chromosome)                 
                map_id2chromosome.append( chr_y )
                
            id = map_chromosome2id[chr_y]
            if strand_y == "-": id = -id

            # add chromsome
            map_chromosome[x:x+size] = id

            x += increment_x + size
            if invert:
                y -= increment_y + size
            else:
                y += increment_y + size

            if y < 0:
                raise "illegal mapping: %i -> %i for %s %s:%s-%s(%s) to %s %s: %s-%s(%s)" % ( x, y, 
                                                                                              chr_x, strand_x, first_x, last_x, size_x,
                                                                                              chr_y, strand_y, first_y, last_y, size_y )
            
    return map_position, map_chromosome, map_chromosome2id, map_id2chromosome

if __name__ == "__main__":

    parser.add_option("-c", "--chromosome", dest="chromosome", type="string",
                      help="chromosome to take." )

    parser.add_option("-m", "--map", dest="filename_map", type="string",
                      help="filename with mapping info.",
                      metavar = "FILE")

    parser.set_defaults( \
        filename_map = "",
        chromosome = None,
        )
    
    (options, args) = E.Start( parser )

    if options.filename_map == "":
        raise "please specify the file with the liftover mapping information."

    if not options.chromosome:
        raise "please give a chromosome."        
        
    map_position, map_chromosome, map_chromosome2id, map_id2chromosome = readLiftOver( open(options.filename_map, "r"), options.chromosome )

    l = 0
    
    for line in sys.stdin:

        if line[0] == "#": continue

        data = line[:-1].split("\t")

        chromosome = data[0]
        range_from, range_to = int(data[1]), int(data[2])
        l += 1

        if chromosome == options.chromosome:

            if options.loglevel >= 1:
                print "#", l, ":", line[:-1]
            
            for x in range(range_from, range_to):
                if map_position[x]:
                    id = map_chromosome[x]
                    if id > 0:
                        c = "+"
                    else:
                        c = "-"
                        id = -id
                        
                    print "%s\t%i\t%s\t%s\t%i" % (chromosome, x, map_id2chromosome[ id ], c, map_position[x] - 1 )
                else:
                    pass
                    # print "%s\t%i\tna" % (chromosome, x )                    

    E.Stop()

    

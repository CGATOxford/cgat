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
linezip.py - 
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

   python linezip.py --help

Type::

   python linezip.py --help

for command line help.

Documentation
-------------

Code
----

'''
import zlib
import sys
import optparse
import struct

import CGAT.Experiment as Experiment

if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: linezip.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-i", "--filename_input", dest="filename_input", type="string",
                      help="filename for compressed input." )

    parser.add_option("-o", "--filename_output", dest="filename_output", type="string",
                      help="filename for compressed output." )

    parser.add_option("-x", "--filename-index", dest="filename_index", type="string",
                      help="filename for indexing information (sorted key)." )

    parser.add_option("-c", "--column-index", dest="column_index", type="string",
                      help="column of sorted key to index on." )

    parser.add_option("-l", "--level", dest="compression_level", type="int",
                      help="compression level to choose [1-9]." )
    
    parser.add_option("-u", "--uncompress", dest="uncompress", action="store_true",
                      help="uncompress file and write to stdout." )

    parser.set_defaults(
        filename_index = None,
        filename_output = None,
        filename_input = None,
        column_index = None,
        level = 6,
        uncompress = False,
        report_step = 100 )

    header_bytes = 4
    max_line_length = 256**header_bytes

    (options, args) = Experiment.Start( parser )

    iteration = 0
    uncompressed, compressed = 0, 0 
        
    if options.uncompress:
        
        f = file( options.filename_input, "rb") 
        while 1:
            
            iteration += 1
            
            size = f.read(4)
            if not size: break
            size = struct.unpack( "l", size )[0]


            c = f.read( size )
            line = zlib.decompress( c ) + "\n"
            options.stdout.write( line )

            uncompressed += len(line)
            compressed += len(c) + header_bytes

            if options.loglevel >= 1 and iteration % options.report_step == 0:
                options.stdlog.write("# iteration %i: compression=%5.2f: %i/%i\n" % (iteration, 100.0 * compressed / uncompressed, compressed, uncompressed ) )

        if options.loglevel >= 1:
            options.stdlog.write("# final: compression=%5.2f: %i/%i\n" % (100.0 * compressed / uncompressed, compressed, uncompressed ) )
            
    else:
        if len(args) == 1:
            options.filename_output= args[0]

        outfile = open(options.filename_output, "w")

        for line in sys.stdin:

            if line[0] == "#": continue

            iteration += 1

            uncompressed += len(line)

            if len(line) >= max_line_length:
                raise "maximum line length exceeded for line starting with: %s" % line[:30]
            
            c = zlib.compress( line[:-1], options.level )
            compressed += len(c) + header_bytes

            outfile.write( struct.pack( "l",len(c) ) + c )

            if options.loglevel >= 1 and iteration % options.report_step == 0:
                options.stdlog.write("# iteration %i: compression=%5.2f: %i/%i\n" % (iteration, 100.0 * compressed / uncompressed, compressed, uncompressed ) )
                
        if options.loglevel >= 1:
            options.stdlog.write("# final: compression=%5.2f: %i/%i\n" % (100.0 * compressed / uncompressed, compressed, uncompressed ) )

        outfile.close()
    
    
        
        
        
    
    

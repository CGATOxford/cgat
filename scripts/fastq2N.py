################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
fastq2N.py -
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   This script takes as input a fastq file and converts a position (base) to an N
   call. This is to keep all reads in line should there be an overrepresentation
   of N calls in any samples.

Usage
-----

Example::

   python fastq2N.py --help

Type::

   python fastq2N.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse
import CGAT.Experiment as E
import CGAT.Fastq as Fastq
import gzip

# define functions to be used in main
def replace(fastqfile, baseToReplace):
    '''replaces the specified base with N'''

    # use gzip as default to open the fastq file
    outf = gzip.open("replaced_" + fastqfile, "w")
    fastq = gzip.open(fastqfile)
    iterator = Fastq.iterate(fastq)
    for record in iterator:
        x = list(record.seq)
        x[int(baseToReplace)] = "N"
        record.seq = "".join(x)
        outf.write("@" + record.identifier + "\n" + record.seq + "\n" + "+" + record.identifier + "\n" + record.quals + "\n")
        
        

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                             usage = globals()["__doc__"] )

    parser.add_option("-i", "--infile", dest="infile", type="string",
                      help="Input filename")
    parser.add_option("-b", "--base-position", dest = "base_position",
                      help = "base position in which to replace with an N")

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if not options.infile:
        raise ValueError ("need to specify input fastq file")

    if not options.base_position:
        raise ValueError ("need to specify which base to convert")
    
    # main function
    replace(options.infile, options.base_position)
    
    ## write footer and output benchmark information.
    E.Stop()
    

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


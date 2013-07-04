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
mali2bootstrap.py - create bootstrap samples from a mali
========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Read a multiple alignment in write out new alignments
of the same lengths with columns sampled from the source
multiple alignment.

The default sampling is with-replacement.

Usage
-----

Example::

   python mali2bootstrap.py --help

Type::

   python mali2bootstrap.py --help

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
import random
import optparse

import CGAT.Experiment as E
import CGAT.Mali as Mali

def getBootstrappedMali( mali, block_size = 1):
    """return a multiple alignment of same
    length with randomly sampled columns from mali

    block_size of 3 samples codons.
    """

    # build list of columns to extract
    l = mali.getWidth() / block_size

    columns = []
    for x in range(l):
        columns.append( random.randint(0,l-1) )
        
    return getMali( mali, columns, block_size)
    
def getSampledMali( mali, sample_size, block_size = 1 ):
    """return a multiple alignment with columns sampled
    from mali without replacement.

    block_size of 3 samples codons.
    """

    # build list of columns to sample from
    l = mali.getWidth() / block_size

    if sample_size >= l:
        raise ValueError( "sample size (%i) larger than number of columns (%i) to sample from" % (sample_size, l) )

    columns = list(range(0,l))
    random.shuffle( columns )

    return getMali( mali, columns[:sample_size], block_size)

def getMali( mali, columns, block_size = 1 ):
        
    new_mali = Mali.Mali()

    for id, val in mali.items():
        sequence = val.mString
        chars = []
        for c in columns:
            chars.append( sequence[c*block_size: c*block_size+block_size] )

        new_sequence = "".join( chars )
        new_mali.addSequence( id, 0, mali.countCharacters( new_sequence), new_sequence )

    return new_mali

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: mali2bootstrap.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal", "stockholm", "phylip" ),
                      help="input format of multiple alignment"  )
    
    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("plain", "fasta", "stockholm", "phylip"),
                      help="output format of multiple alignment"  )

    parser.add_option( "-p", "--output-filename-pattern", dest="output_filename_pattern", type="string",
                       help="pattern for output filenames. Should contain a %(id)i. If not given, the output is to stdout with --separator [default=%default]."  )

    parser.add_option("-n", "--samples", dest="samples", type="int",
                      help="number of samples."  )

    parser.add_option("-r", "--no-replacement", dest="no_replacement", type="int",
                      help="sample without replacement. The parameter gives the size of the multiple alignment [default=%default]."  )
    
    parser.add_option("-b", "--block-size", dest="block_size", type="int",
                      help="block size. Use 3 for sampling from codons."  )

    parser.set_defaults(
        input_format="fasta",
        output_format="fasta",
        samples = 10,
        block_size = 1,
        output_filename_pattern = None,
        no_replacement = None,
        separator = "//")

    (options, args) = E.Start( parser )

    mali = Mali.Mali()

    mali.readFromFile( sys.stdin, format = options.input_format )

    for x in range(options.samples):
        
        if options.no_replacement != None:
            new_mali = getSampledMali( mali, options.no_replacement, options.block_size )
        else:
            new_mali = getBootstrappedMali( mali, options.block_size )

        if options.output_filename_pattern:

            filename = options.output_filename_pattern % { "id": x+1}
            target_directory = os.path.dirname( filename )
            if not os.path.exists( target_directory ):
                os.makedirs( target_directory )
            outfile = open( filename, "w" )
            E.info( "creating mali %s" % filename )
        else:
            outfile = options.stdout

        new_mali.writeToFile( outfile, format = options.output_format )

        if outfile == sys.stdout:
            if options.separator and x < options.samples - 1:
                options.stdout.write( options.separator + "\n" )
        else:
            outfile.close()

    E.Stop()
        

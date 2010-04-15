################################################################################
#   Gene prediction pipeline 
#
#   $Id: concatenate_sequences.py 2782 2009-09-10 11:40:29Z andreas $
#
#   Copyright (C) 2006 Tyler ???? and Andreas Heger 
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
import os, sys, string, re, optparse, math, time, tempfile, subprocess

USAGE="""python %s [OPTIONS] in1 in2 [...]

concatenate sequences in files in1 and in2.

""" % sys.argv[0]

import Experiment
import IOTools
import Genomics
import FastaIterator
import Masker

##------------------------------------------------------------
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: concatenate_sequences.py 2782 2009-09-10 11:40:29Z andreas $", usage = USAGE)

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("in-order",),
                      help="method on how to decide which sequences to concatenate."  )
    
    parser.set_defaults(
        method = "in-order",
        )

    (options, args) = Experiment.Start( parser )

    if len(args) < 2:
        raise "please supply at least two filenames to concatenate."

    iterators = []
    for a in args:
        iterators.append(FastaIterator.FastaIterator( open(a,"r") ) )

    ninput, noutput, nerrors = 0, 0, 0
    
    while 1:

        sequences = []
        ids = []

        for iterator in iterators:
            try:
                cur_record = iterator.next()
            except StopIteration:
                break

            sequences.append( re.sub( " ", "", cur_record.sequence) )
            ids.append( cur_record.title )

        if not sequences: break
        ninput += 1

        if len(sequences) != len(iterators):
            raise "unequal number of sequences in files."

        noutput += 1

        options.stdout.write( ">%s\n%s\n" % (ids[0],
                                             "".join( sequences )))

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nerrors=%i\n" % (ninput, noutput, nerrors) )
        
    Experiment.Stop()
    
    

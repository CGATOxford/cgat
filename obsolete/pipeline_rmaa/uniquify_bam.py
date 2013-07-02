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
uniquify_bam.py.py - only keep unique reads
===========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------


Usage
-----

Example::

   python uniquify_bam.py.py --help

Type::

   python uniquify_bam.py.py --help

for command line help.

Documentation
-------------

Code
----
'''

import os, sys, re, optparse

import Experiment as E
import pysam

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id$",
                                    usage = globals()["__doc__"] )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) != 2 :
        raise ValueError( "please supply two BAM files.")

    samfile = pysam.Samfile( args[0], "rb" )

    readone=set()
    readtwo=set()
    removeone=set()
    removetwo=set()

    for read in samfile.fetch():
        if read.is_read1:
            if read.qname in readone:
                removeone.add(read.qname)
            readone.add(read.qname)
        else:
            if read.qname in readtwo:
                removetwo.add(read.qname)
            readtwo.add(read.qname)

    discarded=0
    samout=pysam.Samfile( args[1], mode='wb', template=samfile)
    for read in samfile.fetch():
        if (read.qname in removeone) and read.is_read1:
            discarded+=1
        elif (read.qname in removetwo) and read.is_read2:
            discarded+=1
        else:
            samout.write(read)
    samfile.close()
    samout.close()

    E.info( "%s of %s first reads removed; %s of %s second reads; %s of %s multi mapped reads which fell into %s positions (average of %s positions per read)" % \
                ( len(removeone), 
                  len(readone), 
                  len(removetwo), 
                  len(readtwo), 
                  len(removeone)+len(removetwo), 
                  len(readone)+len(readtwo), 
                  discarded, 
                  float(discarded)/(len(removeone)+len(removetwo))) )
    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


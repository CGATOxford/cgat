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
nr2table.py - convert description in nr fasta file to table
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This scripts converts headers in ncbi `nr` fasta files to
a tabular format.

Usage
-----

Example::

   python nr2table.py --help

Type::

   python nr2table.py --help

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
import math
import time
import tempfile
import subprocess


import CGAT.Experiment as E
import Bio
import CGAT.FastaIterator as FastaIterator

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: nr2table.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.set_defaults()

    (options, args) = E.Start( parser )
    
    iterator = FastaIterator.FastaIterator( sys.stdin )

    sequences = []
    
    ninput, noutput, nentries = 0, 0, 0

    options.stdout.write( "gid\tsrc\tacc\tprotid\tannotation\tspecies\n" )
    while 1:
        
        cur_record = iterator.next()

        if cur_record == None: break

        ninput += 1

        records = cur_record.title.split( chr(1) )
        for record in records:

            a, anno = re.search("(\S+)\s+(.+)", record).groups()

            vals = a.split("|")
            
            try:
                gi, gid, src, acc, protid = vals
            except ValueError:
                raise "parsing error for record: %s" % record
            
            try:
                annotation, species = re.search( "(.+)\s+\[(.*)\]", anno).groups()
            except AttributeError:
                annotation = anno.strip()
                species = ""

            annotation = re.sub( "\t", " ", annotation)

            options.stdout.write( "\t".join( (gid, src, acc, protid,
                                              annotation, species) ) + "\n" )

            nentries += 1
            
        noutput += 1
        
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nentries=%i\n" % (ninput, noutput, nentries ))
        
    E.Stop()
    
    

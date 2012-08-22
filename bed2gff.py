####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: bed2gff.py 2861 2010-02-23 17:36:32Z andreas $
##
##
####
####
"""
bed2gff.py - convert bed to gff
===============================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script converts a bed-formatted file to
a gff-formatted file.

Usage
-----

Example::

   python <script_name>.py --help

Type::

   python <script_name>.py --help

for command line help.

Documentation
-------------

Code
----

""" 
import sys, re, string, optparse, time, os, itertools, tempfile, subprocess, shutil

import Experiment as E
import Stats
import GFF, GTF, Bed
import IndexedFasta, IOTools

def main( argv = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id: bed2gff.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-a", "--as-gtf", dest="as_gtf", action="store_true",
                      help="output as gtf."  )

    parser.set_defaults( as_gtf = False,
                         id_format = "%08i",
                         test = None )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    as_gtf = options.as_gtf
    id_format = options.id_format

    if as_gtf:
        gff = GTF.Entry()
    else:
        gff = GFF.Entry()

    gff.source = "bed"
    gff.feature = "exon"

    ninput, noutput, nskipped = 0, 0, 0

    id = 0
    for bed in Bed.iterator( options.stdin ):

        ninput += 1

        gff.contig = bed.contig
        gff.start = bed.start 
        gff.end = bed.end
        if bed.mFields and len(bed.mFields) >= 3:
            gff.strand = bed.mFields[2]
        else: 
            gff.strand = "."

        if bed.mFields and len(bed.mFields) >= 2:
            gff.score = bed.mFields[1]
        
        
        if as_gtf:
            if bed.mFields:
                gff.gene_id = bed.mFields[0]
                gff.transcript_id = bed.mFields[0]
            else:
                id += 1
                gff.gene_id = id_format % id
                gff.transcript_id = id_format % id            
        else:
            if bed.mFields:
                gff.source = bed.mFields[0]
            
        options.stdout.write( str(gff) + "\n" )

        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )

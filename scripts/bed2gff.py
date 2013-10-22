"""
bed2gff.py - convert bed to gff/gtf
===================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Intervals BED GFF Conversion

Purpose
-------

This script converts a :term:`bed`-formatted file to a :term:`gff` or :term:`gtf`-formatted file.

It aims to populate the appropriate fields in the :term:`gff` file with columns in the :term:`bed` file.

If ``--is-gtf`` is set and a name column in the :term:`bed` file is present, its contents will be set 
as ``gene_id`` and ``transcript_id``. Otherwise, a numeric ``gene_id`` or ``transcript_id`` will be set 
according to ``--id-format``.

Usage
-----

Example::

   python bed2gff.py < in.bed > out.gff

Type::

   python bed2gff.py --help

for command line help.

Command line options
--------------------

""" 
import sys
import re
import string
import optparse
import time
import os
import itertools
import tempfile
import subprocess
import shutil

import CGAT.Experiment as E
import CGAT.Stats as Stats
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IOTools as IOTools

def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id: bed2gff.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-a", "--as-gtf", dest="as_gtf", action="store_true",
                      help="output as gtf."  )

    parser.add_option("-f", "--id-format", dest="id_format", type="string",
                      help="format for numeric identifier if --as-gtf is set and no name in bed file [%default]."  )

    parser.set_defaults( as_gtf = False,
                         id_format = "%08i",
                         test = None )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    as_gtf = options.as_gtf
    id_format = options.id_format

    if as_gtf:
        gff = GTF.Entry()
    else:
        gff = GTF.Entry()

    gff.source = "bed"
    gff.feature = "exon"

    ninput, noutput, nskipped = 0, 0, 0

    id = 0
    for bed in Bed.iterator( options.stdin ):

        ninput += 1

        gff.contig = bed.contig
        gff.start = bed.start 
        gff.end = bed.end
        if bed.fields and len(bed.fields) >= 3:
            gff.strand = bed.fields[2]
        else: 
            gff.strand = "."

        if bed.fields and len(bed.fields) >= 2:
            gff.score = bed.fields[1]
        
        if as_gtf:
            if bed.fields:
                gff.gene_id = bed.fields[0]
                gff.transcript_id = bed.fields[0]
            else:
                id += 1
                gff.gene_id = id_format % id
                gff.transcript_id = id_format % id            
        else:
            if bed.fields:
                gff.source = bed.fields[0]
            
        options.stdout.write( str(gff) + "\n" )

        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput, nskipped) )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )

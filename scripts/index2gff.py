################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: index2gff.py 2880 2010-04-07 08:44:13Z andreas $
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
index2gff.py - convert indexed fasta file to gff file
=====================================================

:Author: Andreas Heger
:Release: $Id: index2gff.py 2880 2010-04-07 08:44:13Z andreas $
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

""" 

import os
import sys
import re
import optparse

import CGAT.IndexedFasta as IndexedFasta
import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: index2gff.py 2880 2010-04-07 08:44:13Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    ## do sth
    ninput, nskipped, noutput = 0, 0, 0

    fasta = IndexedFasta.IndexedFasta( options.genome_file )

    entry = GTF.Entry()
    entry.start = 0
    entry.feature = "contig"
    entry.source = "genome"

    for contig, size in fasta.getContigSizes( with_synonyms = False ).iteritems():
        ninput += 1
        entry.contig = contig
        entry.end = int(size)
        options.stdout.write( "%s\n" % str(entry) )
        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


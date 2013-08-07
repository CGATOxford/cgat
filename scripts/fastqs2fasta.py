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
fastqs2fasta.py
=============================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script is used to interleave two fastq files into a single fasta file. Read1 is followed by 
read2 in the resultant file. Fastq files must be sorted by read identifier.
Writes the output to stdout.

Usage
-----

Example::

   python fastqs2fasta.py --help

Type::

   python fastqs2fasta.py --help

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
import itertools
import CGAT.IOTools as IOTools
import CGAT.Fastq as Fastq

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = optparse.OptionParser( version = "%prog version: $Id: script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-a", "--fastq1", dest="fastq1", type="string",
                      help="supply read1 fastq file"  )
    parser.add_option("-b", "--fastq2", dest="fastq2", type="string",
                      help="supply read2 fastq file"  )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    fastq1 = IOTools.openFile(options.fastq1)
    fastq2 = IOTools.openFile(options.fastq2)

    E.info("reading fastq files")
    for f1, f2 in itertools.izip(Fastq.iterate(fastq1), Fastq.iterate(fastq2)):
        assert f1.identifier[:-1] == f2.identifier[:-1], "fastq entries are not sorted"
        options.stdout.write(">%s\n%s\n>%s\n%s\n" % (f1.identifier, f1.seq, f2.identifier, f2.seq))

## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


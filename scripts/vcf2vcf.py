'''
vcf2vcf.py - manipulate vcf files
=================================

Purpose
-------

Manipulate vcf-formatted files.


Usage
-----

Example::

   cat in.vcf | python vcf2vcf.py - --reorder alphabetical > sorted.vcf

This command generates a sorted vcf with the sample columns in in.vcf in alphabetical order.

Type::

   python vcf2vcf.py --help

for command line usage.

Methods
-------

This script provides the following methods:

re-order
   reorder sample columns in vcf formatted file according to a given sort order

Documentation
-------------

This is a tool for manipulating vcf-formatted files.  The following options are available:

+-----------+-------------------------+
|--reorder  |reorders sample columns  |
+-----------+-------------------------+

Reorder
^^^^^^^

To sort sample columns into alphabetical order::

   cat example.vcf | python vcf2vcf.py - --reorder alphabetical

This will sort the columns in the example.vcf into the order "#CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,
SAMPLE_A, SAMPLE_B, SAMPLE_C, ..."

To specify a non-alphabetical order::

   cat example.vcf | python vcf2vcf.py - --reorder SAMPLE_C,SAMPLE_A,SAMPLE_B,...

This will sort the columns in the example.vcf into the order "#CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,  |
SAMPLE_C, SAMPLE_A, SAMPLE_B, ..."

Command line options
--------------------
'''

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
import CGAT.IOTools as IOTools
import CGAT.VCF as VCF
    
def main( argv = sys.argv ):

    parser = E.OptionParser( version = "%prog version: $Id: bed2bed.py 2861 2010-02-23 17:36:32Z andreas $", 
                             usage = globals()["__doc__"] )

    parser.add_option( "--reorder", dest="reorder", type="string",
                       help="reorder columns. Give column names as comma-separated list or specify ``alphabetical`` [default=%default]"  )

    parser.set_defaults( reorder = None )
    
    (options, args) = E.Start( parser, add_pipe_options = True )

    noutput = 0
    
    infile = VCF.VCFFile( options.stdin )

    if options.reorder: 
        order = options.reorder.split(",")
        if "alphabetical" in order:
            order = sorted(infile.samples)
    else:
        order = False

    infile.writeHeader( options.stdout, order = order )
        
    for vcf in infile:
        if order: vcf.order = order 
        options.stdout.write( str(vcf) + "\n" )

    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv ) )

####
####
##
##
## Copyright (C) 2008 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: bed2bed.py 2861 2010-02-23 17:36:32Z andreas $
##
##
####
####
'''
vcf2vcf.py - manipulate vcf files
=================================

Purpose
-------

manipulate vcf-formatted files.


Usage
-----

Type::

   python vcf2vcf.py --help

for command line usage.

Methods
-------

This script provides several methods:

re-order
   reorder columns in vcf formatted file according to a given sort order.

Code
----
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
import CGAT.VCF as VCF
    
def main( argv = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id: bed2bed.py 2861 2010-02-23 17:36:32Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "--reorder", dest="reorder", type="string",
                       help="reorder columns. Give column names as comma-separated list or specify ``alphabetical`` [default=%default]"  )

    parser.set_defaults( reorder = None )
    
    (options, args) = E.Start( parser, add_pipe_options = True )


    noutput = 0

    infile = VCF.VCFFile( sys.stdin )

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

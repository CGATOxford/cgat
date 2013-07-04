################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2011 David Sims
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
====================================================
Add a number of random reads to an aligned bam file
====================================================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

This program takes two input file: a Samtools indexed genome file and a bam file to add reads to.
It also takes several input parameters: number of reads to be added, the read length, the insert size and the insert size standard deviation
Currently the program only adds paired end reads, with both reads the same length
The default insert size is 250, and insert size standard diviation is 20.
The default read length is 50 and the number of reads to be added is 10,000.


Usage
=====

add_random_reads_to_bam.py -g <genome.fa> -b <input.bam>-o <output.bam> -i insert size -s insert size sd -r number of reads -l read length


Code
====

"""

# load modules

import logging as L
import CGAT.Experiment as E
import CGAT.Pipeline as P
import sys
import os
import re
import shutil
import itertools
import math
import glob
import time
import gzip
import collections
import random
import optparse
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: snp2table.py 2861 2010-02-23 17:36:32Z davids $", usage = globals()["__doc__"] )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with Samtools indexed genome [default=%default]."  )
    parser.add_option("-b", "--bam-file", dest="bam_file", type="string",
                      help="filename of bam to add reads to [default=%default]."  )
    parser.add_option("-i", "--insert-size", dest="isize", type="string",
                      help="Insert size [default=%default]."  )
    parser.add_option("-s", "--insert-sd", dest="isd", type="string",
                      help="Insert size standard deviation [default=%default]."  )
    parser.add_option("-r", "--nreads", dest="nreads", type="string",
                      help="Number of random reads to add [default=%default]."  )
    parser.add_option("-l", "--readlength", dest="readlength", type="string",
                      help="length of reads to generate [default=%default]."  )
    parser.add_option("-o", "--output", dest="output_file", type="string",
                      help="output filename  [default=%default]."  )

    parser.set_defaults(
        genome_file = None,
        bam_file = None,
        isize=250,
        isd=20,
        nreads=10000,
        readlength=50,
        output_file = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    # Generate random reads and add to bam
    track = P.snip(os.path.basename(options.bam_file),".bam")
    readlen = options.readlength
    isize = options.isize
    isd = options.isd
    nreads = options.nreads
    genome = options.genome_file
    bam = options.bam_file
    out = options.output_file
    statement = '''java -jar -Xmx2048m /ifs/apps/bio/simseq-72ce499/SimSeq.jar 
                       -1 %(readlen)s -2 %(readlen)s \
                       --error  /ifs/apps/bio/simseq-72ce499/examples/hiseq_mito_default_bwa_mapping_mq10_1.txt \
                       --error2 /ifs/apps/bio/simseq-72ce499/examples/hiseq_mito_default_bwa_mapping_mq10_2.txt \
                       --insert_size %(isize)s \
                       --insert_stdev %(isd)s \
                       --read_number %(nreads)s \
                       --read_prefix simseq_ \
                       --reference %(genome)s \
                       --duplicate_probability 0.0 \
                       --out simseq.sam > simseq.log; ''' % locals()
    P.run()
    statement = '''samtools view -bS -t %(genome)s.fai -o simseq.bam simseq.sam; 
                   samtools sort simseq.bam simseq.srt;
                   samtools sort %(bam)s %(track)s.srt;
                   samtools merge %(out)s %(track)s.srt.bam simseq.srt.bam;
                   samtools index %(out)s;'''
    P.run()

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


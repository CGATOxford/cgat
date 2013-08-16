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
radar.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python radar.py --help

Type::

   python radar.py --help

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

USAGE="""python %s [OPTIONS] < mali > mali

reformat and edit a multiple alignment.

Methods are:

translate:      translate sequences.

mark-codons:    add space after each codon
pseudo-codons:  translate, but keep register with codons
interleaved-codons: mix amino acids and codons

Parameters are given to the option parameters in a comma-separated list in the order
that the edit operations are called upon.

""" % sys.argv[0]

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Genomics as Genomics
import CGAT.FastaIterator as FastaIterator

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: radar.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("translate", "translate-till-stop",
                               "truncate-at-stop",
                               "mark-codons", "build-map",
                               "pseudo-codons", "interleaved-codons",
                               "remove-gaps",
                               "mask-seg", "mask-bias", "mask-codons",
                               "upper", "lower" ),
                      help="method to apply to sequences."  )
    
    parser.add_option("-p", "--parameters", dest="parameters", type="string",
                      help="parameter stack for methods that require one."  )

    parser.add_option("-t", "--type", dest="type", type="choice",
                      choices = ("aa", "na"),
                      help="sequence type (aa or na)."  )

    parser.set_defaults(
        only_codons = False,
        methods = "",
        parameters = "",
        type = "na",
        aa_mask_chars = "xX",
        aa_mask_char = "x",
        na_mask_chars = "nN",
        na_mask_char = "n",
        gap_chars = "-.",
        )

    (options, args) = E.Start( parser )
    
    iterator = FastaIterator.FastaIterator( sys.stdin)

    nseq = 0

    step_size = 3

    while 1:
        cur_record = iterator.next()
        
        if cur_record is None: break
        nseq += 1

        sequence = re.sub( " ", "", cur_record.sequence)
        l = len(sequence)

        if options.loglevel >= 1:
            options.stdlog.write("# analysing %s : length = %i\n" % (cur_record.title, l))

        nidentities = [0] * l

        for x in range(0, l, step_size):
            print x
            for d in range(0, l-x, step_size):
                for s1 in sequence[x:l-d]:
                    for s2 in sequence[x+d:l]:
                        if s1 == s2:
                            nidentities[d] += 1

        for x in range(l):
            print x, distances[x]
                
    E.Stop()

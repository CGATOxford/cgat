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
graph_reweight_links.py - 
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

   python graph_reweight_links.py --help

Type::

   python graph_reweight_links.py --help

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
import getopt
import tempfile
import time
import popen2

USAGE="""python %s [OPTIONS] < graph.in > graph.out

Version: $Id: graph_reweight_links.py 2782 2009-09-10 11:40:29Z andreas $

Wrapper for rescoring blast alignments.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-m, --method=                   method to use for scoring [bitscore|normalize]
-f, --self-score=               filename with self-scores
--lambda=                       lambda for bitscore calculation
--k=                            K for bitscore calculation
--expected=                     expected score in matrix

-d, --distance=                 convert to distance (# = maximum distance)

new scores are saved in third column (overwriting the E-Value)

Methods:

kimura                          = set third column to kimura two parameter score from pid
bitscore                        = set third column to bitscore from score.
normalize-product               = normalize third column by dividing with self-scores
                                new = old * old / self1 / self2
normalize-max                   = normalize third column by max(old/self1), max(old/self2)
normalize-min                   = normalize third column by min(old/self1), min(old/self2)
normalize-avg                   = normalize third column by avg((old/self1 + old/self2) / 2)

Format:

blast                           = Blast graph
edges                           = edge list (weighted)

""" % sys.argv[0]

param_long_options = ["help", "verbose=", "method=", "lambda=", "k=", "self-scores=", "expected=", "distance=",
                      "version"]

param_short_options = "hv:m:f:o:d:"

param_loglevel = 1
param_method = "bitscore"
param_filename_self_scores = None
param_format = "blast"

import CGAT.Experiment as E
import CGAT.BlastAlignments as BlastAlignments
import math

param_lambda = 0.267
param_K = 0.0410
param_expected = -0.5209
param_gop = -11.0
param_gep = -1.0
param_max_distance = 0.0

##-------------------------------------------------------------------------------
if __name__ == "__main__":

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "--version", ):
            print "version="
            sys.exit(0)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)
        elif o in ("-m", "--method"):
            param_method = a
        elif o in ("-f", "--self-scores"):
            param_filename_self_scores = a
        elif o == "--lambda":
            param_lambda = float(a)
        elif o == "--k":
            param_K = float(a)
        elif o == "--expected":
            param_expected = float(a)
        elif o in ("-d", "--distance"):
            param_max_distance = float(a)
            
    print E.GetHeader()
    print E.GetParams()

    if param_filename_self_scores:
        self_scores = {}
        for line in open(param_filename_self_scores, "r"):
            if line[0] == "#": continue
            d = line[:-1].split("\t")[:2]
            if d[0] not in self_scores: self_scores[d[0]] = 0.0
            self_scores[d[0]] = max( self_scores[d[0]], float(d[1]))

    if param_method == "kimura":
        a = 0
        f = lambda x: x < 0.85 and 0.0000001 - math.log( 1.0 - x - 0.2 * x * x ) or 5.2030

    elif param_method == "bitscore":
        a = 1        
        lK = math.log(param_K)
        l2 = math.log(2)
        f = lambda x: (param_lambda * x - lK) / l2
        
    elif param_method in ("normalize-product", "normalize-max", "normalize-min", "normalize-avg"):
        a = 2
        if param_method == "normalize-product":
            f = lambda x,y,z: x * x / self_scores[y] / self_scores[z]
        elif param_method == "normalize-max":
            f = lambda x,y,z: max( x / self_scores[y], x / self_scores[z])
        elif param_method == "normalize-min":
            f = lambda x,y,z: min( x / self_scores[y], x / self_scores[z])
        elif param_method == "normalize-avg":
            f = lambda x,y,z: (x / self_scores[y] + x / self_scores[z]) / 2.0
    else:
        raise "unknown method %s" % param_method
    
    ninput, noutput, nskipped, nfailed = 0, 0, 0, 0

    for line in sys.stdin:

        if line[0] == "#": continue
        
        token1, token2, score = line.split("\t")[:3]
        try:
            score = float(score)
        except ValueError:
            # ignore headers
            continue
        ninput += 1
        
        if a == 0:
            score = f( (100.0 - score) / 100.0 )
        elif a == 1:
            score = f( score )
        elif a == 2:
            score = f( score, token1, token2 )

        
        if param_max_distance:
            score = param_max_distance - score

        print string.join(map( str, (token1, token2, score) ), "\t") 

        noutput += 1
            
    print "# ninput=%i, noutput=%i, nskipped=%i, failed=%i" % (\
        ninput, noutput, nskipped, nfailed )
    
    print E.GetFooter()



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
graph_blast2pairs.py - 
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

   python graph_blast2pairs.py --help

Type::

   python graph_blast2pairs.py --help

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

Version: $Id: graph_blast2pairs.py 2782 2009-09-10 11:40:29Z andreas $

Wrapper for rescoring blast alignments.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-m, --method=                   method to use for scoring [bitscore|normalize]
-f, --self-score=               filename with self-scores
--lambda=                       lambda for bitscore calculation
--k=                            K for bitscore calculation
--expected=                     expected score in matrix
-a, --append                    append new columns at the end
--evalue-to-log                 convert evalue to log
--effective-length=             calculate evalue from bitscore based on effective sbjct length
--min-evalue=                   minimum evalue (after which it is truncated)
--version                       ouptut version

new scores are saved in third column (overwriting the E-Value)

Methods:

kimura                          = set third column to kimura two parameter score from pid
bitscore                        = set third column to bitscore from score.
normalize-product               = normalize third column by dividing with self-scores
                                new = old * old / self1 / self2
normalize-max                   = normalize third column by max(old/self1), max(old/self2)
normalize-scoredist-avg         = normalize third column method by Sonnhammer & Hollich (2005)
scoredist-avg                   = set third column to score by Sonnhammer & Hollich (2005)
scoredist-min                   = set third column to modified score Sonnhammer & Hollich (2005) (minimum score for
                                        estimating the upper bound).
gapless-score                   = score of alignment without the gaps
reset-evalue                    = reset evalue based on bitscore
""" % sys.argv[0]

param_long_options = ["help", "verbose=", "method=", "lambda=", "k=", "self-scores=", 
                      "expected=", "append", "evalue-to-log",
                      "effective-length=", "version" ]

param_short_options = "hv:m:f:o:a"


param_loglevel = 1
param_method = "bitscore"
param_filename_self_scores = None

import CGAT.Experiment as E
import CGAT.BlastAlignments as BlastAlignments
import math

param_lambda = 0.267
param_K = 0.0410
param_expected = -0.5209
param_gop = -11.0
param_gep = -1.0

param_append = False
param_evalue_to_log = False

param_effective_length = None

param_min_evalue = 1e-200

def CalculateGapScore( ali, gop, gep ):

    s = 0.0
    data = re.split( "[+-]", ali[1:])
    for x in data[1:-1:2]:
        s += gop + int(x) * gep
    return s

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
        elif o in ("-a", "--append"):
            param_append = True
        elif o == "--evalue-to-log":
            param_evalue_to_log = True
        elif o == "--effective-length":
            param_effective_length = int(a)
        elif o == "--min-evalue":
            param_min_evalue = float(a)
            
    if param_loglevel >= 1:
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
        
    elif param_method in ("normalize-product", "normalize-product-distance",
                          "normalize-max", "normalize-max-distance"):
        a = 2
        if param_method == "normalize-product":
            f = lambda x,y,z: x * x / self_scores[y] / self_scores[z]
        elif param_method == "normalize-product-distance":            
            f = lambda x,y,z: max( 0.0, 1.0 - x * x / self_scores[y] / self_scores[z])
        elif param_method == "normalize-max":
            f = lambda x,y,z: max( x / self_scores[y], x / self_scores[z])
        elif param_method == "normalize-max-distance":            
            f = lambda x,y,z: max( 0.0, 1.0 - max( x / self_scores[y], x / self_scores[z]))
            
    elif param_method == "normalize-scoredist-avg":
        a = 3
        f = lambda x,y,z,l: max(0.0, -100.0 * math.log( (x - param_expected * l) / ( (self_scores[y] + self_scores[z]) * 0.5 - param_expected * l ) ))

    elif param_method == "scoredist-avg":
        a = 4
        f = lambda x,y,z,l: max(0.0, -100.0 * math.log( (x - param_expected * l) / ( (self_scores[y] + self_scores[z]) * 0.5 - param_expected * l ) ))

    elif param_method == "scoredist-min":
        a = 4
        f = lambda x,y,z,l: max(0.0, -100.0 * math.log( (x - param_expected * l) / ( min(self_scores[y], self_scores[z]) - param_expected * l ) ))
    elif param_method == "gapless-score":
        a = 5
        f = lambda x,a,b: x - CalculateGapScore( a, param_gop, param_gep ) - CalculateGapScore( b, param_gop, param_gep)
    elif param_method == "reset-evalue":
        a = 6
        if param_evalue_to_log:
            # this way is less likely to underflow (2^-s might be zero for large s)
            me = math.log(param_min_evalue)
            l2 = math.log(2)
            le = math.log(param_effective_length)
            f = lambda s,m,n: max(me, -s * l2 + math.log(m) + le)
        else:
            f = lambda s,m,n: max(param_min_evalue, math.pow(2,-s) * m * n)
            
        param_evalue_to_log = False
    else:
        raise "unknown method %s" % param_method

    ninput, noutput, nskipped, nfailed = 0, 0, 0, 0

    for line in sys.stdin:

        if line[0] == "#": continue

        link = BlastAlignments.Link()
        link.Read( line )
        ninput += 1

        try:
            if a == 0:
                new_val = f( (100.0 - link.mPercentIdentity) / 100.0 )
            elif a == 1:
                new_val = f( link.score )
            elif a == 2:
                ## note: used to be evalue
                new_val = f( link.mBitScore, link.mQueryToken, link.mSbjctToken )
            elif a == 3:
                new_val = f( link.mEvalue, link.mQueryToken, link.mSbjctToken, max(link.mQueryTo - link.mQueryFrom, link.mSbjctTo - link.mSbjctFrom) + 1 )
            elif a == 4:
                new_val = f( link.score, link.mQueryToken, link.mSbjctToken, max(link.mQueryTo - link.mQueryFrom, link.mSbjctTo - link.mSbjctFrom) + 1 )
            elif a == 5:
                new_val = int(f( link.score, link.mQueryAli, link.mSbjctAli))
            elif a == 6:
                new_val = f( link.mBitScore, link.mQueryLength, param_effective_length)
                
        except KeyError:
            if param_loglevel >= 2:
                print "# Key error in line", line[:-1]
                nfailed += 1
            continue

        if param_evalue_to_log:
            link.mEvalue = math.log( link.mEvalue )
        
        if param_append:
            print str(link) + "\t" + str(new_val)
        else:
            if a in (0, 1, 2, 3, 4, 6):
                link.mEvalue = new_val
            elif a in (5,):
                link.score = new_val
            print str(link)
            
        noutput += 1
            
    print "# ninput=%i, noutput=%i, nskipped=%i, failed=%i" % (\
        ninput, noutput, nskipped, nfailed )

    if param_loglevel >= 1:    
        print E.GetFooter()



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
graph2besthits.py - 
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

   python graph2besthits.py --help

Type::

   python graph2besthits.py --help

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
import optparse
import math

import CGAT.Experiment as E

USAGE="""python %s [OPTIONS] < graph.in > graph.out

Version: $Id: graph2besthits.py 2782 2009-09-10 11:40:29Z andreas $

Parse graph and only output best hits for each query per genome. Note that
the graph needs to be sorted by query_token and sbjct_token.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-p, --pattern-genome=           pattern to extract genome
-m, --method=                   method [blast|score|pid]
-s, --score-factor=             thresholding on score (multiplicative)
-i, --pide-factor=              thresholding on pide (additive)
""" % sys.argv[0]


param_score_threshold_factor = 1.0
param_pide_threshold_factor = 0.0

class Link:

    def __init__(self):
        self.mQueryToken, self.mSbjctToken, self.score = "","", 0
        
    def Read( self, line ):
        self.mQueryToken, self.mSbjctToken, self.score = line.split("\t")[:3]
        self.score = float(self.score)

    def __str__(self):
        return "\t".join( (self.mQueryToken, self.mSbjctToken, str(self.score)) )
    
##-------------------------------------------------------------------------------
def PrintMatches( matches,
                  options,
                  filter=None,
                  score_threshold_factor = 1.0 ):
    """print best match per organism (sorted by E-Value).

    If there are several matches with the same E-Value, all are printed.
    """

    orgs = matches.keys()
    orgs.sort()

    if options.method == "distance":
        f1 = lambda x,y: cmp(x.score, y.score)
        f2 = lambda x,y: x.score <= y
    elif options.method == "similarity":
        f1 = lambda x,y:  -cmp(x.score, y.score)
        f2 = lambda x,y: x.score >= y
        raise "unknown method"

    nbest = options.nbest
    
    for org in orgs:
        
        if filter and filter == org: continue

        m = matches[org]
        
        score_threshold = int(m[0].score * options.score_threshold_factor)

        if options.loglevel >= 2:
            print "# %i matches for %s, min_score=%i" % (len(m), org, score_threshold)

        if options.loglevel >= 3:
            for mm in m:
                print str(mm)
            
        m.sort( f1 )

        s = m[0].score
        x = 0

        if nbest:
            for x in range(min(nbest, len(m))):            
                print str(m[x])
        elif options.score_threshold_factor != 1.0:
            ## take best match and all within score threshold factor
            for x in range(len(m)):
                if f2( m[x], score_threshold):
                    print str(m[x])
        else:
            ## write all with same primary score
            while x < len(m) and s == m[x].score and f2( m[x], score_threshold) :            
                print str(m[x])
                x += 1

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: graph2besthits.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-p", "--pattern-genome", dest="pattern_genome", type="string",
                      help="genome_pattern."  )
    
    parser.add_option("-m", "--method", dest="method", type="string",
                      help="method to use."  )

    parser.add_option("-n", "--nbest", dest="nbest", type="int",
                      help="use nbest links."  )
    
    parser.add_option("-s", "--keep-self", dest="keep_self", action="store_true",
                      help="apply to links within same genome."  )

    parser.set_defaults(
        pattern_genome = "^([^|]+)|",
        method = "distance",
        keep_self = False,
        nbest = 0,
        score_threshold_factor = 1.0,
        )

    (options, args) = E.Start( parser )
    
    ninput, noutput, nskipped, nfailed = 0, 0, 0, 0

    last_query_token = None
    rx = re.compile(options.pattern_genome)

    matches = {}
    
    for line in sys.stdin:

        if line[0] == "#": continue
        
        link = Link()
        try:
            link.Read( line )
        except ValueError:
            sys.stderr.write( "parsing error in line %s\n" % line[:-1] )
            continue
        
        ninput += 1
        if link.mQueryToken != last_query_token :
            if last_query_token:
                if not options.keep_self:
                    filter = rx.search(last_query_token).groups()[0]
                else:
                    filter = None
                PrintMatches( matches,
                              options,
                              filter = filter )
                noutput += 1
            matches = {}
            last_query_token = link.mQueryToken

        org = rx.search(link.mSbjctToken).groups()[0]
        
        if org not in matches: matches[org] = []
        matches[org].append( link )

    if not options.keep_self:
        filter = rx.search(last_query_token).groups()[0]
    else:
        filter = None

    PrintMatches( matches,
                  options,
                  filter = filter )
    
    noutput += 1
    
    print "# ninput=%i, noutput=%i" % (\
        ninput, noutput )
    
    E.Stop()




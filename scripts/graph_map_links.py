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
graph_map_links.py - 
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

   python graph_map_links.py --help

Type::

   python graph_map_links.py --help

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
import optparse
import hashlib

USAGE="""python %s [OPTIONS] < graph.in > graph.out

Version: $Id: graph_map_links.py 2782 2009-09-10 11:40:29Z andreas $

Map blast links.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-q, --map-query=                map to use for query
-s, --map-sbjct=                map to use for sbjct
-m, --multiple                  write multiple matches
-k, --keep-unmapped             keep entries not mapped
-i, --identity-map              mapping is done by identities
""" % sys.argv[0]

import CGAT.Experiment as E
import CGAT.BlastAlignments as BlastAlignments

def ReadIdentityMap( infile):
    """read identity map.

    multiple entries can either be separated over several lines
    or concatenated by semicolon in the same line.
    """
    
    m = {}
    for line in infile:

        if line[0] == "#": continue
        olds, news = map(lambda x: x.split(";"), line[:-1].split("\t")[:2])

        for old in olds:
            if old not in m: m[old] = []
            for new in news:
                m[old].append( new )
    return m

##-------------------------------------------------------------------------------
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: graph_map_links.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-q", "--map-query", dest="filename_map_query", type="string",
                      help="filename with queries to map." )
    parser.add_option("-s", "--map-sbjct", dest="filename_map_sbjct", type="string",
                      help="filename with queries to map." )
    parser.add_option("-m", "--multiple", dest="multiple", action="store_true",
                      help="map multiple options." )
    parser.add_option("-k", "--keep-unmapped", dest="keep_unmapped", action="store_true",
                      help="keep unmapped entries." )
    parser.add_option("-i", "--map-identity", dest="map_identity", action="store_true",
                      help="use identity maps (only identifier change)." )
    parser.add_option("-n", "--non-redundant", dest="non_redundant", action="store_true",
                      help="write only unique links (requires a lot of memory for large graphs).")
    parser.set_defaults( \
        filename_map_query = None,
        filename_map_sbjct = None,
        multiple = False,
        keep_unmapped = False,
        map_identity = False,
        report_step = 1000000,
        non_redundant = False)
        
    (options, args) = E.Start( parser )

    if options.filename_map_query:
        infile = open(options.filename_map_query, "r")        
        if options.map_identity:
            map_query = ReadIdentityMap( infile )
        else:
            map_query = BlastAlignments.ReadMap( infile, options.multiple )
        infile.close()
        if options.loglevel >= 1:
            print "# read maps for %i queries." % len(map_query)
    else:
        map_query = None

    if options.filename_map_sbjct:
        if options.filename_map_sbjct == options.filename_map_query:
            map_sbjct = map_query
        else:
            infile = open(options.filename_map_sbjct, "r")                
            if options.map_identity:
                map_sbjct = ReadIdentityMap( infile )            
            else:
                map_sbjct = BlastAlignments.ReadMap( open(options.filename_map_sbjct, "r"), options.multiple)
            infile.close()
        if options.loglevel >= 1:
            print "# read maps for %i sbjcts." % len(map_sbjct)
    else:
        map_sbjct = None

    nfailed = 0
    ninput = 0
    nskipped = 0
    noutput = 0

    # number of identical/mapped links
    nsame, nmapped = 0, 0

    printed = {}
    
    map = BlastAlignments.Map()

    for line in sys.stdin:
        
        if line[0] == "#": continue

        data = line[:-1].split("\t")

        map.Read( line )
        skip = False
        ninput += 1

        if options.loglevel >= 2:
            print "#", str(map)

        if options.loglevel >= 2 and ninput % options.report_step == 0:
            sys.stderr.write( "# progress: ninput=%i, noutput=%i, nhash=%i\n" % (ninput, noutput, len(printed)) )
        
        if options.multiple:
            skip = False
            if map_query != None:
                if map.mQueryToken in map_query:
                    mq = map_query[map.mQueryToken]
                else:
                    skip = True
            else:
                mq = [ None ]                                

            if map_sbjct != None:                
                if map.mSbjctToken in map_sbjct:
                    ms = map_sbjct[map.mSbjctToken]
                else:
                    skip = True
            else:
                ms = [ None ]
                
            if skip:
                nskipped += 1
                continue
            
            if options.map_identity:

                ## only if non_redundant is set, do global comparison
                if not options.non_redundant: printed = {}
                            
                new_map = map.GetClone()
                do_redundant = len(mq) > 1 or len(ms) > 1                
                for q in mq:
                    for s in ms:
  
                        new_map.mQueryToken = q
                        new_map.mSbjctToken = s

                        ## check for non-redundant links for 1:many or many:many mappings
                        if do_redundant:
                            key = "%s-%i-%i-%s-%i-%i" % (new_map.mQueryToken, new_map.mQueryFrom, new_map.mQueryTo,
                                                         new_map.mSbjctToken, new_map.mSbjctFrom, new_map.mSbjctTo)
                        
                            # hash key to save space
                            hkey = hashlib.md5(key).digest()

                            if hkey in printed: continue
                            
                            printed[hkey] = 1

                        print string.join( [str(new_map)]+ data[9:], "\t")
                        noutput += 1
                        if new_map.mQueryToken == map.mQueryToken and \
                           new_map.mSbjctToken == map.mSbjctToken:
                            nsame += 1
                        else:
                            nmapped += 1
                        
            else:
                for q in mq:
                    for s in ms:
                        new_map = map.GetClone()

                        if options.loglevel >= 2:
                            print "#", str(q)
                            print "#", str(s)

                        is_ok = new_map.MapAlignment( q, s )

                        if not is_ok:
                            nfailed += 1
                        else:
                            print string.join( [str(new_map)]+ data[9:], "\t")
                            noutput += 1
        else:

            if map_query != None:
                if map.mQueryToken in map_query:
                    mq = map_query[map.mQueryToken ]
                else:
                    mq = None
                    skip = True
            else:
                mq = None

            if map_sbjct != None:
                if map.mSbjctToken in map_sbjct:
                    ms = map_sbjct[map.mSbjctToken ]
                else:
                    ms = None
                    skip = True
            else:
                ms = None
                
            if skip and not options.keep_unmapped:
                nskipped += 1
                continue

            if options.loglevel >= 2:
                print "#", str(mq)
                print "#", str(ms)

            if mq or ms:
                is_ok = map.MapAlignment( mq, ms )
            else:
                is_ok = True

            if not is_ok:
                nfailed += 1
            else:
                print string.join( [str(map)]+ data[9:], "\t")
                noutput += 1

    print "# ninput=%i, noutput=%i, nskipped=%i, nfailed=%i, nsame=%i, nmapped=%i" % (\
        ninput, noutput, nskipped, nfailed, nsame, nmapped )
    
    E.Stop()



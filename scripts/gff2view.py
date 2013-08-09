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
gff2view.py - view gff intervals in UCSC genome browser
=======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads intervals from a :term:`gff` formatted file and
opens the corresponding location in the UCSC genome browser. If there
are multiple intervals, multiple windows/tabs will be opened, so be
careful.

Usage
-----

Example::

   python gff2view.py --help

Type::

   python gff2view.py --help

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

import webbrowser

USAGE="""python %s [OPTIONS] < stdin

open a browser window at several locations. This script interfaces
with UCSC and genome browser.

TODO: add support for ensembl
""" % sys.argv[0]

import CGAT.GTF as GTF
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: gff2view.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-t", "--target", dest="target", type="choice",
                      choices = ("ucsc", "gbrowser"),
                      help="target location to open [%default]."  )

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome."  )

    parser.add_option( "--is-gtf", dest="is_gtf", action="store_true",
                      help="input is gtf."  )

    parser.add_option("-f", "--add-flank", dest="flank", type="int",
                      help="add # nucleotides for each region."  )

    parser.add_option("-z", "--zoom", dest="zoom", type="float",
                      help="zoom out (# > 1) or in (# < 1)."  )

    parser.add_option("-c", "--chunk-size", dest="chunk_size", type="int",
                      help="number of tabs to display in one go."  )

    parser.add_option( "--ucsc-assembly", dest="ucsc_assembly", type="string",
                      help="ucsc assembly."  )

    parser.add_option( "--ucsc-user-tracks", dest="ucsc_user_tracks", type="string",
                      help="ucsc user tracks."  )

    parser.add_option( "--gbrowser-assembly", dest="gbrowser_assembly", type="string",
                      help="gbrowser assembly."  )

    parser.add_option( "--randomize", dest="randomize", action="store_true",
                       help="randomize input [%default]" )
    
    parser.set_defaults(
        ucsc_assembly = "ponAbe2",
        ucsc_url = "http://genome.ucsc.edu/cgi-bin/hgTracks",
        gbrowser_assembly = "Songbird",
        gbrowser_url = "http://genserv.anat.ox.ac.uk/cgi-bin/devel/gbrowse",
        genome_file = None,
        ucsc_custom_annotation = "http://wwwfgu.anat.ox.ac.uk/~andreas/ucsc_tracks/%s",
        ucsc_user_tracks = None,
        flank = None,
        zoom = None,
        chunk_size = 50,
        is_gtf = False,
        target = "ucsc",
        randomize = False, 
        joined = False,
        )

    (options, args) = E.Start( parser )

    if len(args) != 1:
        print USAGE
        raise "please specify the gff file to open."
    
    if options.is_gtf:
        entry_iterator = GTF.iterator
        chunk_iterator = GTF.flat_gene_iterator
    else:
        entry_iterator = GTF.iterator
        if options.joined:
            chunk_iterator = GTF.joined_iterator
        else:
            chunk_iterator = GTF.chunk_iterator

    if len(args) == "0" or args[0] == "-":
        iterator = chunk_iterator( entry_iterator( sys.stdin) )
    else:
        iterator = chunk_iterator( entry_iterator( open(args[0], "r") ) )

    nopened = 0

    # b = webbrowser.get( "konqueror" )
    b = webbrowser.get( "firefox" )

    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta( options.genome_file )
    else:
        fasta = None

    if options.ucsc_user_tracks:
        annotations = "hgt.customText=%s" % (options.ucsc_custom_annotation % options.ucsc_user_tracks)
    else:
        annotations = None

    for chunk in iterator:
        start = min( [ x.start for x in chunk ] )
        end = max( [ x.end for x in chunk ] )

        if options.flank: 
            start -= options.flank
            end += options.flank
            
        if options.zoom:
            s = end - start
            d = options.zoom * s - s
            start -= d
            end += d
            
        start = max( 0, start )

        contig = chunk[0].contig

        if fasta:
            contig = fasta.getToken( contig )
            end = min( end, fasta.getLength( contig ) )
            
        if len(contig) < 3: contig = "chr%s" % contig

        if options.target == "ucsc":
            url_options = [ "db=%s" % options.ucsc_assembly,
                            "position=%s:%i-%i" % (contig, start, end) ]

            if annotations: url_options.append( annotations )

            url = "%s?%s" % (options.ucsc_url, 
                             "&".join( url_options) )
        elif options.target == "gbrowser":
            url = "%s/%s?name=%s:%i..%i" % (options.gbrowser_url, 
                                            options.gbrowser_assembly,
                                            contig, 
                                            start, 
                                            end )

        print "# opening browser window for:"
        print "#", url

        if nopened % options.chunk_size == 0: 
            if nopened != 0:
                x = raw_input('showing %i - hit return to continue:' % options.chunk_size )
            b.open_new(url)
            first = False
        else:
            b.open_new_tab( url )

        nopened += 1
        
    E.Stop()

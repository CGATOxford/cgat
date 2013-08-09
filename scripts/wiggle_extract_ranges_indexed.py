#!/usr/bin/env python

"""
wiggle_extract_ranges_indexed.py - 
===============================================================================

Reads a list of intervals and a wiggle file. Produces a list of 
values similar to wiggle_to_simple.py. 

The difference is that this file allows indexed access to a
set of wiggle files indexed by wiggle_build_index.py

It is assumed that each file `wiggle_fname` has a corresponding `wiggle_fname`.index
file.

WARNING: bz2/bz2t support and file cache support are new and not as well
         tested.

Note: this requires a patched bx/wiggle.py

usage: %prog wiggle_fname1 wiggle_fname2 ... [options] < interval_file
   -v, --version:    Output version
   -s, --src=s:      Use this src for all intervals
   -p, --prefix=p:   Prepend this to each src before lookup
   -S, --strand:     Strand is included as an additional column.
   -C, --usecache:   Use a cache that keeps blocks of the WIGGLE files in memory (requires ~20MB per WIGGLE)
"""

import psyco_full

from bx.cookbook import doc_optparse

import CGAT.Wiggle as Wiggle
import sys
import os.path

def main():
    # Parse Command Line
    options, args = doc_optparse.parse( __doc__ )
    if options.version: return

    try:
        wiggle_files = args
        if options.src: fixed_src = options.src
        else: fixed_src = None
        if options.prefix: prefix = options.prefix
        else: prefix = None
        do_strand = bool( options.strand )
        use_cache = bool( options.usecache )
    except:
        doc_optparse.exit()

    # Open indexed access to wiggles
    index = Wiggle.WiggleMultiIndexedAccess( wiggle_files,
                                             keep_open = True,
                                             use_cache=use_cache )

    for line in sys.stdin:
        strand = "+"
        fields = line.split()
        if fixed_src:
            src, start, end = fixed_src, int( fields[0] ), int( fields[1] )
            if do_strand: strand = fields[2]
        else:
            src, start, end = fields[0], int( fields[1] ), int( fields[2] )
            if do_strand: strand = fields[3]
        if prefix: src = prefix + src
        blocks = index.get( src, start, end )
        for x, values in blocks:
            for v in values:
                print "%s\t%i\t%f" % (src, x, v)
                x += 1
        
if __name__ == "__main__": main()

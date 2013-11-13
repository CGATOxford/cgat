#!/usr/bin/env python

"""
wiggle_build_index.py - build an index for a wiggle track
=========================================================

If index_file is not provided wiggle_file.index is used.

usage: %prog wiggle_file index_file
--version: version string

"""

import psyco_full

from bx.cookbook import doc_optparse
from bx import interval_index_file
import bx.wiggle

import sys
import os.path
from bx.misc.seekbzip2 import SeekableBzip2File

def main():

    # Parse command line

    options, args = doc_optparse.parse( __doc__ )
    if options.version: return

    try:
        wiggle_file = args[0]
        # If it appears to be a bz2 file, attempt to open with table
        if wiggle_file.endswith( ".bz2" ):
            table_file = wiggle_file + "t"
            if not os.path.exists( table_file ):
                doc_optparse.exit( "To index bz2 compressed files first "
                                   "create a bz2t file with bzip-table." )
            # Open with SeekableBzip2File so we have tell support
            wiggle_in = SeekableBzip2File( wiggle_file, table_file )
            # Strip .bz2 from the filename before adding ".index"
            wiggle_file = wiggle_file[:-4]
        elif wiggle_file.endswith( ".lzo" ):
            from bx.misc.seeklzop import SeekableLzopFile
            table_file = wiggle_file + "t"
            if not os.path.exists( table_file ):
                doc_optparse.exit( "To index lzo compressed files first "
                                   "create a lzot file with bzip-table." )
            # Open with SeekableBzip2File so we have tell support
            wiggle_in = SeekableLzopFile( wiggle_file, table_file )
            # Strip .lzo from the filename before adding ".index"
            wiggle_file = wiggle_file[:-4]
        else:
            wiggle_in = open( wiggle_file )
        # Determine the name of the index file
        if len( args ) > 1:
            index_file = args[1]
        else:
            index_file = wiggle_file + ".index"
    except:
        doc_optparse.exception()

    indexes = interval_index_file.Indexes()

    # Can't use the iterator, as there is no next() and thus
    # no way to access the positions. The following code is 
    # modified from wiggle.py
    last_chrom = None
    start = None
    end = None
    first_pos = None

    # always for wiggle data
    strand = '+'

    mode = "bed"

    while 1:
        pos = wiggle_in.tell()
        line = wiggle_in.readline()
        if not line: break

        if line.isspace() or line.startswith( "track" ) or line.startswith( "#" ) or line.startswith( "browser" ):
            continue
        elif line.startswith( "bed" ):
            indexes.add( fields[0], int( fields[1] ), int( fields[2] ), pos )
        elif line.startswith( "variableStep" ) or line.startswith( "fixedStep"):
            if first_pos != None:
                indexes.add( last_chrom, start, end, first_pos )
            first_pos = pos
            header = bx.wiggle.parse_header( line )
            last_chrom = header['chrom']
            start = int(header['start']) - 1
            end = start
            current_step = None
            if 'span' in header: 
                current_span = int( header['span'] )
            else: 
                current_span = 1
            if 'step' in header: 
                current_step = int( header['step'] )

            if line.startswith( "variableStep" ):
                mode = "variableStep"
            else:
                mode = "fixedStep"
        elif mode == "variableStep":
            fields = line.split()
            end = int( fields[0] ) - 1 + current_span
        elif mode == "fixedStep":
            end += current_step
        else:
            raise "Unexpected input line: %s" % line.strip()

    out = open( index_file, 'w' )
    indexes.write( out )
    out.close()

if __name__ == "__main__": main()



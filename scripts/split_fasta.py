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
split_fasta.py - 
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

   python split_fasta.py --help

Type::

   python split_fasta.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import os
import optparse
import CGAT.FastaIterator as FastaIterator
import CGAT.IOTools as IOTools
import CGAT.Experiment as E

class Files:

    mFiles = {}
    
    def __init__(self,
                 output_pattern   = None,
                 skip_identifiers = False):

        self.mOutputPattern = output_pattern
        self.mSkipIdentifiers = skip_identifiers
        self.mCounts = {}
        
    def __del__(self):
        """close all open files."""
        for file in self.mFiles.values():
            file.close()

    def GetFile( self, identifier ):
        return identifier

    def GetFilename( self, identifier ):
        """get filename for an identifier."""

        if self.mOutputPattern:
            return re.sub( "%s", str(identifier), self.mOutputPattern )
        else:
            return identifier

    def OpenFile( self, filename, mode = "w" ):
        """open file.
        
        If file is in a new directory, create directories.
        """
        if mode in ("w", "a"):
            dirname = os.path.dirname(filename)
            if dirname and not os.path.exists( dirname ):
                os.makedirs( dirname )
                
        return open( filename, mode)
        
    def Write( self, identifier, sequence ):

        filename = files.GetFilename( identifier )
        
        if filename not in self.mFiles:
            
            if len(self.mFiles) > 1000:
                for f in self.mFiles.values(): f.close()
                self.mFiles = {}
                
            self.mFiles[filename] = self.OpenFile( filename, "a" )

        if self.mSkipIdentifiers:
            self.mFiles[filename].write( "%s\n" % (sequence.sequence))
        else:
            self.mFiles[filename].write( ">%s\n%s\n" % (sequence.title, sequence.sequence))

        if filename not in self.mCounts:
            self.mCounts[filename] = 0
        self.mCounts[ filename ] += 1

    def DeleteFiles( self, min_size = 0 ):
        """delete all files below a minimum size."""

        ndeleted = 0
        for filename, counts in self.mCounts.items():
            if counts < min_size:
                os.remove( filename )
                ndeleted += 1
                
        return ndeleted

class FilesChunks(Files):

    def __init__(self,
                 chunk_size, **kwargs):

        Files.__init__(self, **kwargs )
        self.mChunkSize = chunk_size
        self.mFilename = 0
        
    def GetFilename(self, identifier):

        if not self.mFilename or self.mCounts[self.mFilename] % self.mChunkSize == 0:
            self.mFilename = re.sub( "%s", str(len(self.mCounts) + 1), self.mOutputPattern )
            
        return self.mFilename
        
            
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: split_fasta.py 1714 2007-12-11 16:51:12Z andreas $")
    
    parser.add_option("-f", "--file", dest="input_filename", type="string",
                      help="input filename. If not given, stdin is used.",
                      metavar="FILE" )

    parser.add_option ("-i", "--input-pattern", dest="input_pattern", type="string",
                       help="input pattern. Parses description line in order to extract id." )

    parser.add_option ("-o", "--output-pattern", dest="output_pattern", type="string",
                       help="output pattern. Gives filename for a given sequence." )

    parser.add_option ("-n", "--num-sequences", dest="num_sequences", type="int",
                       help="split by number of sequences (not implemented yet)." )

    parser.add_option ("-m", "--map", dest="map_filename", type="string",
                       help="map filename. Map identifiers to filenames",
                       metavar="FILE" )

    parser.add_option ("-s", "--skip-identifiers", dest="skip_identifiers", action="store_true",
                       help="do not write identifiers.",
                       metavar="FILE" )

    parser.add_option ("--min-size", dest="min_size", type="int",
                       help="minimum cluster size.")

    parser.set_defaults( \
        input_filename = None,
        map_filename = None,
        skip_identifiers = False,
        input_pattern = "^(\S+)",
        min_size = 0,
        num_sequences = None,
        output_pattern = "%s" )

    (options, args) = E.Start( parser ) 

    if options.input_filename:
        infile = IOTools.openFile( options.input_filename, "r")
    else:
        infile = sys.stdin

    if options.map_filename:
        map_id2filename = IOTools.ReadMap( open( options.map_filename, "r") )
    else:
        map_id2filename = {}
        
    if options.num_sequences:
        files = FilesChunks( chunk_size = options.num_sequences,
                             output_pattern = options.output_pattern,
                             skip_identifiers = options.skip_identifiers )
        
    else:
        files = Files( output_pattern = options.output_pattern,
                       skip_identifiers = options.skip_identifiers )

    if options.input_pattern:
        rx = re.compile( options.input_pattern )
    else:
        rx = None

    ninput = 0
    noutput = 0
    identifier = None
    chunk = 0
    
    for seq in FastaIterator.iterate( infile ):
        
        ninput += 1

        if rx:
            try:
                identifier = rx.search(seq.title).groups()[0]
            except AttributeError:
                print "# parsing error in description line %s" % (seq.title)
        else:
            identifier = seq.title

        if map_id2filename:
            if identifier in map_id2filename:
                identifier = map_id2filename[identifier]
            else:
                continue

        files.Write( identifier, seq )
        noutput += 1

    if options.input_filename:
        infile.close()

    ## delete all clusters below a minimum size
    ## Note: this has to be done at the end, because
    ## clusters sizes are only available once both the fasta
    ## file and the map has been parsed.
    if options.min_size:
        ndeleted = files.DeleteFiles( min_size = options.min_size )
    else:
        ndeleted = 0
        
    if options.loglevel >= 1:
        print "# input=%i, output=%i, ndeleted=%i" % (ninput, noutput, ndeleted)
        
    E.Stop()

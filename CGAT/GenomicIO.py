################################################################################
#   Gene prediction pipeline 
#
#   $Id: GenomicIO.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2007 Andreas Heger
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
GenomicIO.py - Subroutines for working on I/O of large genomic files
====================================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

I tried the Biopython parser, but it was too slow for large genomic chunks.
"""

import os, sys, array, string, re, types, tempfile
USAGE="""python IndexedFasta.py [options] name files

Index fasta formatted files to create a database called "name".
"""

import os, sys, array, string, re, types, optparse

global_last_sbjct_token = None
global_last_sbjct_from = None
global_last_sbjct_to = None

global_index = {}

##------------------------------------------------------------
class MySequence(array.array):

    def __init__(self, *args):
        self.mType = args[0]
        array.array.__init__(self, *args)

    def __getslice__( self, *args):
        """return slice as a string."""
        return array.array.__getslice__(self, *args).tostring()

    def __setslice__( self, start, end, sub):
        """set slice start:end from a string sub."""
        return array.array.__setslice__(self, start, end, array.array( self.mType, sub ))

    def __str__(self):
        return self.tostring()

##------------------------------------------------------------
def index_file( filenames, db_name ):
    """index file/files.

    Two new files are create - db_name.fasta and db_name.idx
    """

    if db_name + ".fasta" in filenames:
        raise ValueError( "database (%s.fasta) is part of input set." % db_name)
        
    outfile_index = open( db_name + ".idx", "w")
    outfile_fasta = open( db_name + ".fasta", "w" )

    if type(filenames) == types.StringType:
        filenames = [filenames]

    identifiers = {}
    first = True
    lsequence = 0
    for filename in filenames:

        infile = open( filename, "r")
        
        for line in infile:

            pos = outfile_fasta.tell()        

            if line[0] == "#":  continue
            
            if line[0] == ">" :
                identifier = re.split("\s", line[1:-1])[0]
                if identifier in identifiers:
                    raise ValueError, "%s occurs more than once in %s and %s" % (identifiers[identifier], filename)
                identifiers[identifier] = filename
                
                if not first:
                    outfile_fasta.write("\n")
                    outfile_index.write("%i\n" % lsequence)
                first = False
                    
                outfile_fasta.write( "%s" % line)
                outfile_index.write( "%s\t%i\t%i\t" % (identifier, pos, outfile_fasta.tell()) )
                lsequence = 0
            else:
                s = re.sub( "\s", "", line[:-1] )
                lsequence += len(s)
                outfile_fasta.write( s )
                
    outfile_index.write("%i\n" % lsequence)
    
##------------------------------------------------------------
def _load_index( filename ):
    """load index from file."""

    global global_index
    for line in open(filename, "r"):
        
        (identifier, offset1, offset2, lsequence) = line[:-1].split("\t")

        offset1 = int(offset1)
        offset2 = int(offset2)
        lsequence = int(lsequence)
        
        global_index[identifier] = (offset1, offset2, lsequence)

##------------------------------------------------------------
def index_exists( filename ):
    """check if a certain file has been indexed."""
    return os.path.exists( filename + ".idx" )

##------------------------------------------------------------
def getSequence( db_name,
                 sbjct_token, sbjct_strand,
                 sbjct_from, sbjct_to,
                 as_array = False,
                 forward_coordinates = False):
    """get genomic fragment.
    """

    if not IndexIsLoaded: _load_index( db_name )

    if sbjct_token not in global_index:
        raise KeyError, "%s not in index" % sbjct_token

    is_negative = sbjct_strand in ("-", 0, "0", "-1")
    pos1, pos2, lsequence = global_index[sbjct_token][1]
    infile = open( db_name + ".fasta" )

    if sbjct_to == 0: sbjct_to = lsquence
    if sbjct_from <= 0: sbjct_from = 0

    if is_negative and not forward_coordinates:
        first_pos = pos2 + max( lsequence - sbjct_to, 0 ) 
        last_pos  = pos2 + lsequence - sbjct_from
    else:
        first_pos = pos2 + max(0, sbjct_from)
        last_pos  = pos2 + min(sbjct_to, lsequence)

    infile.seek( first_pos )
    
    p = MySequence( "c" )

    p.fromfile( infile, last_pos - first_pos)

    if is_negative:
        p.reverse()

    if as_array:
        return p
    else:
        return p[:]

def splitFasta( infile, chunk_size, dir="/tmp", pattern=None ):
    """split a fasta file into a subset of files.

    If pattern is not given, random file names are chosen.
    """

    n = 0
    chunk = 0
    
    def __getFilename( chunk ):
        if pattern:
            outname = pattern % chunk
            outfile = os.open( outname, "w" )
        else:
            (outfile, outname) = tempfile.mkstemp( dir=dir )
        return (outfile, outname)

    outfile, outname = __getFilename(chunk)
    filenames = [ outname ]
    noutput = 0
    for line in infile:
        if line[0] == "#": continue
        if line[0] == ">":
            n += 1
            if n > chunk_size:
                os.close(outfile)
                n = 1               
                outfile, outname = __getFilename(chunk)
                filenames.append( outname )
            noutput += 1
            
        os.write( outfile, line )
        
    os.close(outfile)

    if noutput == 0:
        os.remove( outname )
        return []

    return filenames

class IndexedFasta:

    def __init__( self, dbname ):
        self.mDbname = dbname
        self.mIsLoaded = False

    def __getitem__(self, key ):
        return self.getSequence( key, "+", 0, 0, as_array = True )
        
    def loadIndex( self ):

        self.mIndex = {}
        for line in open(self.mDbname + ".idx", "r"):
            
            (identifier, offset1, offset2, lsequence) = line[:-1].split("\t")
            
            offset1 = int(offset1)
            offset2 = int(offset2)
            lsequence = int(lsequence)
            
            self.mIndex[identifier] = (offset1, offset2, lsequence)

        self.mIsLoaded = True

    def getLength( self, sbjct_token ):
        """return sequence length for sbjct_token."""
        if not self.mIsLoaded: self.loadIndex()
        
        return self.mIndex[sbjct_token][2]

    def getContigs( self ):
        """return hash with contig sizes."""
        if not self.mIsLoaded: self.loadIndex()
        return self.mIndex.keys()

    def getContigSizes( self ):
        """return hash with contig sizes."""
        if not self.mIsLoaded: self.loadIndex()
        contig_sizes = {}
        for key, val in self.mIndex.items():
            contig_sizes[key] = val[2]
        return contig_sizes

    def getSequence( self,
                     sbjct_token, sbjct_strand,
                     sbjct_from, sbjct_to,
                     converter = None,
                     as_array = False):
        """get genomic fragment.
        """

        if not self.mIsLoaded: self.loadIndex()

        if sbjct_token not in self.mIndex:
            raise KeyError, "%s not in index" % sbjct_token

        is_negative = sbjct_strand in ("-", 0, "-1", "0")
        pos1, pos2, lsequence = self.mIndex[sbjct_token]
        infile = open( self.mDbname + ".fasta", "r" )

        if sbjct_to == 0: sbjct_to = lsequence
        if sbjct_to > lsequence: sbjct_to = lsequence
        if sbjct_from <= 0: sbjct_from = 0

        if converter:
            first_pos, last_pos = converter( sbjct_from, sbjct_to,
                                             sbjct_strand in ("+", 1, "1"),
                                             lsequence )
        else:
            first_pos, last_pos = sbjct_from, sbjct_to

        if is_negative:
            first_pos, last_pos = lsequence - last_pos, lsequence - first_pos

        first_pos += pos2
        last_pos += pos2
        
        infile.seek( first_pos )

        p = MySequence( "c" )

        p.fromfile( infile, last_pos - first_pos)

        if is_negative:
            p.reverse()            
            p = MySequence("c", string.translate( p[:], string.maketrans("ACGTacgt", "TGCAtgca") ) )

        if as_array:
            return p
        else:
            return p[:]

## converter functions. There are a few of them, but
## I did not want to call functions within functions.
def __one_forward_closed(x, y, c, l):
    x -= 1
    if not c: x, y = l - y, l - x
    return x, y 
def __zero_forward_closed(x, y, c, l):
    y += 1
    if not c: x, y = l - y, l - x
    return x, y 
def __one_both_closed(x, y, c  = None, l = None):
    return x-1, y
def __zero_both_closed(x, y, c = None, l = None):
    return x, y + 1
def __one_forward_open(x, y, c, l):
    x -= 1
    y -= 1
    if not c: x, y = l - y, l - x
    return x, y 
def __zero_forward_open(x, y, c, l):
    if not c: x, y = l - y, l - x
    return x, y 
def __one_both_open(x, y, c  = None, l = None):
    return x-1, y-1
def __zero_both_open(x, y, c = None, l = None):
    return x, y

def getConverter( format ):
    """return a converter function for
    converting various coordinate schemes into
    0-based, both strand, closed-open ranges.

    converter functions have the parameters
    x, y, s, l: with x and y the coordinates of
    a sequence fragment, s the strand (True is positive)
    and l being the length of the contig.
    """

    data = set(format.split("-"))

    if "one" in data:
        if "forward" in data:
            if "closed" in data:
                return __one_forward_closed                
            else:
                return __one_forward_open
        else:
            if "closed" in data:
                return __one_both_closed
            else:
                return __one_both_open
    else:
        if "forward" in data:
            if "closed" in data:
                return __zero_forward_closed
            else:
                return __zero_forward_open
        else:
            if "closed" in data:
                return __zero_both_closed
            else:
                return __zero_both_open                

if __name__ == "__main__":

    import Experiment
    
    parser = E.OptionParser( version = "%prog version: $Id: GenomicIO.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    (options, args) = Experiment.Start( parser )

    if len(args) < 2:
        print USAGE
        sys.exit(1)

    if options.loglevel >= 1:
        options.stdlog.write("# indexing the following files: %s\n" % (",".join( args[1:] ) ))
        options.stdlog.write("# creating database %s\n" % args[0])
        options.stdlog.flush()
        
    index_file( args[1:], args[0] )
    
    Experiment.Stop()

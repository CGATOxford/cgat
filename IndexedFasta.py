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
IndexedFasta.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, array, string, re, types, optparse, time, struct, math, tarfile, logging
import platform, anydbm

USAGE="""python %s [options] name [ files ]

Index fasta formatted files to create a database called "name".

Example:

  python %s oa_ornAna1_softmasked /net/cpp-mirror/ucsc/ornAna1/bigZips/ornAna1.fa.gz > oa_ornAna1_softmasked.log

""" % (sys.argv[0], sys.argv[0] )

# import psyco
# psyco.full()

import math
import random
import zlib
import gzip
import cStringIO
import Experiment as E
from AString import AString

##------------------------------------------------------------
class Uncompressor:
    def __init__(self, filename, unmangler):
        self.mFile = open(filename, "rb" )
        self.mUnMangler = unmangler
        
    def read( self, block_size, indices, start, end ):
        """read an uncompressed block from start:end.
        
        The compressed chunk starts at first_pos.
        NOTE: This is poorly implemented - do better.
        """

        # skip over uncompressed blocks
        d = int(math.floor(float(start) / block_size) )
        r = start % block_size
        assert( d < len(indices) )
        self.mFile.seek( indices[d] )

        # read x bytes of compressed data, at least one full chunk.
        nchunks = int(math.ceil( float((r+end-start)) / block_size) )

        fragments = []
        for x in range(d, d+nchunks):
            s = self.mFile.read( indices[x+1] - indices[x] )
            fragments.append(self.mUnMangler( s ))
        u = "".join(fragments)

        assert len(u) >= end - start, "fragment smaller than requested size: %i > %i-%i=%i" % (len(u), end,start,end-start)
        
        return u[r:r+end-start]

##------------------------------------------------------------    
def writeFragments( outfile_fasta, 
                    outfile_index,
                    fragments, 
                    mangler, size,
                    write_all = False):
    """write mangled fragments to *outfile_fasta* in chunks of *size*
    updating *outfile_index*.

    returns part of last fragment that has not been written and is 
    less than *size* and the number of fragments output.

    If *write_all* is True, all of the fragments are written to
    the file and the last file position is added to *outfile_index*
    as well.
    """

    s = "".join(fragments)
    rest = len(s) % size
    if len(s) > size:
        for x in range(0, len(s)-rest, size):
            outfile_index.write( "\t%i" % outfile_fasta.tell() )
            outfile_fasta.write( mangler(s[x:x+size]) )
        
    if rest:
        if write_all:
            outfile_index.write("\t%i" % outfile_fasta.tell() )
            outfile_fasta.write( mangler( s[-rest:] ) )
            outfile_index.write("\t%i" % outfile_fasta.tell() )            
            return ""
        else:
            return s[-rest:]
    else:
        return ""

def gzip_mangler(s):

    xfile = cStringIO.StringIO()

    gzipfile = gzip.GzipFile( fileobj = xfile, mode="wb" )
    gzipfile.write( s )
    gzipfile.close()
    
    m = xfile.getvalue()
    xfile.close()
    return m

def gzip_demangler(s):

    gzipfile = gzip.GzipFile( fileobj = cStringIO.StringIO(s), mode="rb" )
    m = gzipfile.readline()
    return m

##------------------------------------------------------------
class Translator:
    """translate a sequence."""
    def __init__(self): 
        self.mRegEx = re.compile( " +" )

    def __call__(self, sequence):
        return "".join( self.mMapScore2Char[self.mMapScore2Score[int(x)]] for x in self.mRegEx.split( sequence.strip() ) ) 

    def translate( self, sequence ):
        raise NotImplementedError( "translate not implemented.")

class TranslatorPhred(Translator):
    """translate phred quality scores."""
    def __init__(self, *args, **kwargs):
        Translator.__init__(self,*args,**kwargs)
        self.mMapScore2Char = [ chr(33 + x) for x in range( 0, 93) ]
        self.mMapScore2Score = range( 0, 93 )

    def translate( self, sequence ):
        return array.array( "I", ( ord(x)-33 for x in sequence ) )

class TranslatorSolexa(Translator):
    """translate solexa quality scores."""
    def __init__(self, *args, **kwargs):
        Translator.__init__(self,*args,**kwargs)
        self.mMapScore2Char = [ chr(64 + x) for x in range( 0, 128) ]
        self.mMapScore2Score =[ int(10.0 * math.log( 1.0 + 10 ** (x / 10.0)) / math.log(10)+.499) for x in range(-64,65) ]

    def translate( self, sequence ):
        raise NotImplementedError( "translate not implemented.")
        return array.array( "i", ( ord(x) - 64 ) )

class TranslatorRange200(Translator):
    """translate pcap quality scores.

    For example for PCAP scores.

    These scores range from 0 to 100 and are the
    "a weighted sum of input base quality values
    (Huang and Madan 1999)

    The numerical values from 0 to 200 are stored
    as values form 33 to 233
    "
    """
    def __init__(self, *args, **kwargs):
        Translator.__init__(self,*args,**kwargs)
        self.mMapScore2Char = [ chr(33 + x) for x in range( 0, 200) ]

    def __call__(self, sequence):
        try:
            return "".join( self.mMapScore2Char[int(x)] for x in self.mRegEx.split( sequence.strip() ) ) 
        except ValueError, msg:
            raise ValueError( msg + " parsing error in fragment: %s" % sequence )

    def translate( self, sequence ):
        return array.array( "I", ( ord(x) - 33 for x in sequence ) )

class TranslatorBytes(Translator):
    """output binary values as bytes permitting values from 0 to 255

    Note the resulting file will not be iterable as newline is not
    a record-separator any more.
    """
    def __init__(self, *args, **kwargs):
        Translator.__init__(self,*args,**kwargs)

    def __call__(self, sequence):
        try:
            return "".join( chr(int(x)) for x in self.mRegEx.split( sequence.strip() ) ) 
        except ValueError, msg:
            print "parsing error in line: %s" % sequence
            print "message=%s" % str(msg)
            return ""
        
    def translate( self, sequence ):
        return array.array( "I", ( ord(x) for x in sequence ) )

##------------------------------------------------------------
class MultipleFastaIterator:

    def __init__(self, filenames, regex_identifier = None ):

        if type(filenames) == types.StringType:
            self.mFilenames = [filenames]
        else:
            self.mFilenames = filenames

        self.mRegexIdentifier = regex_identifier
        self.mIterator = self._iterate()
        
    def __iter__(self):
        return self

    def next(self):
        try:
            return self.mIterator.next( )
        except StopIteration:
            return None

    def _iterate( self ):
        """iterate over muliple files."""
        
        def _iter( infile ):

            identifier = None

            for line in infile:
                if line.startswith("#"):  continue
                if line.startswith(">"):

                    if self.mRegexIdentifier:
                        try:
                            identifier = re.search(regex_identifier, line[1:-1]).groups()[0]
                        except AttributeError:
                            raise ValueError("could not parse identifier from line %s - check the input" % line[1:-1])
                    else:
                        identifier = re.split("\s", line[1:-1])[0]

                else:
                    if not identifier:
                        raise ValueError("refusing to emit sequence without identifier - check the input")
                    yield identifier, line.strip()

        for filename in self.mFilenames:
            if filename.endswith( "tar.gz" ):
                tf = tarfile.open( filename, "r" )
                for f in tf:
                    b, ext = os.path.splitext( f.name )
                    if ext.lower() in ( ".fasta", ".fa" ):
                        E.info( "extracting %s" % f.name)
                        infile = tf.extractfile( f )
                        for x in _iter( infile ): yield x
                    else:
                        E.info( "skipping %s" % f.name )

                tf.close()
                continue
            elif filename.endswith(".gz"):
                infile = gzip.open( filename, "r" )
            elif filename == "-":
                infile = sys.stdin
            else:
                infile = open( filename, "r")

            for x in _iter( infile ): yield x
            if filename != "-": infile.close()

        raise StopIteration

##------------------------------------------------------------
def createDatabase( db, iterator,
                    force = False,
                    synonyms = None,
                    compression = None,
                    random_access_points = None,
                    regex_identifier = None,
                    clean_sequence = False,
                    ignore_duplicates = False,
                    allow_duplicates = False,
                    translator = None ):
    """index files in filenames to create database.

    Two new files are created - db.fasta and db_name.idx

    If compression is enabled, provide random access points
    every # bytes.

    Dictzip is treated as an uncompressed file.

    regex_identifier: pattern to extract identifier from description line.
    If None, the part until the first white-space character is used.

    translator: specify a translator
    """

    if db.endswith( ".fasta"): db = db[:-len(".fasta")]

    if compression:
        if compression == "lzo":
            import lzo
            def lzo_mangler( s ): return lzo.compress(s, 9)
            mangler = lzo_mangler
            db_name = db + ".lzo"
            write_chunks = True
        elif compression == "zlib":
            def zlib_mangler( s ): return zlib.compress( s, 9)
            mangler = zlib_mangler
            db_name = db + ".zlib"
            write_chunks = True     
        elif compression == "gzip":
            mangler = gzip_mangler
            db_name = db + ".gz"
            write_chunks = True     
        elif compression == "dictzip":
            import dictzip
            mangler = lambda x: x
            db_name = db + ".dz"
            write_chunks = False
        elif compression == "bzip2":
            import bz2
            def bzip_mangler( s ): return bz2.compress( s, 9)
            mangler = bzip_mangler
            db_name = db + ".bz2"
            write_chunks = True
        elif compression == "debug":
            mangler = lambda x: x
            db_name = db + ".debug"
            write_chunks = True
        elif compression == "rle":
            import RLE
            mangler = RLE.compress
            db_name = db + ".rle"
            write_chunks = True
        else:
            raise ValueError("unknown compression library: %s" % compression)

        index_name = db + ".cdx"

        if write_chunks and random_access_points == None or random_access_points <= 0:
            raise ValueError("specify chunksize in --random-access-points")
        
    else:
        mangler = lambda x: x
        db_name = db + ".fasta"
        write_chunks = False
        index_name = db + ".idx"
    
    if os.path.exists( db_name ) and not force:
        raise ValueError( "database %s already exists." % db_name )

    if os.path.exists( index_name ) and not force:
        raise ValueError( "database index %s already exists." % index_name )
    
    outfile_index = open( index_name, "w" )
    if compression == "dictzip":
        import dictzip
        if random_access_points == None or random_access_points <= 0:
            raise ValueError("specify dictzip chunksize in --random-access-points")
        outfile_fasta = dictzip.open( db_name, "wb", buffersize=1000000, chunksize=random_access_points )
        compression = None
    else:
        outfile_fasta = open( db_name, "wb" )

    identifiers = {}
    lsequence = 0
    identifier_pos, sequence_pos = 0, 0

    translation = string.maketrans("xX", "nN")

    fragments = []
    lfragment = 0
    last_identifier = None

    while 1:

        try:
            result = iterator.next()
        except StopIteration: 
            break

        if not result: break

        identifier, fragment = result

        if identifier != last_identifier:

            ## check for duplicate identifiers
            if identifier in identifiers:
                if ignore_duplicates:
                    raise ValueError, "ignore duplicates not implemented"
                elif allow_duplicates:
                    # the current implementation will fail if the same identifiers
                    # are directly succeeding each other
                    # better: add return to iterator that indicates a new identifier
                    out_identifier = identifier + "_%i" % (identifiers[identifier])
                    identifiers[identifier] += 1
                    identifiers[out_identifier] = 1
                else:
                    raise ValueError, "%s occurs more than once" %\
                        (identifier,)
            else:
                identifiers[identifier] = 1
                out_identifier = identifier

            if last_identifier:
                if write_chunks:
                    writeFragments( outfile_fasta, outfile_index, 
                                    fragments, mangler,
                                    size = random_access_points, 
                                    write_all = True )
                    
                    fragments = []
                    lfragment = 0
                else:
                    outfile_fasta.write( "\n" )
                        
                outfile_index.write("\t%i\n" % lsequence)

            # write identifier
            identifier_pos = outfile_fasta.tell()
            outfile_fasta.write( mangler(">%s\n" % out_identifier) )
            sequence_pos = outfile_fasta.tell()
                
            outfile_index.write( "%s\t%i" % (out_identifier,
                                             identifier_pos ) )
            if write_chunks:
                outfile_index.write( "\t%i" % random_access_points )
            else:
                outfile_index.write( "\t%i" % sequence_pos )

            fragments = []
            lsequence = 0
            last_identifier = identifier
                
        if translator:
            s = translator( fragment )
        else:
            s = re.sub( "\s", "", fragment.strip() )
            if clean_sequence:
                s = s.translate( translation )
                        
        lsequence += len(s)
                
        if write_chunks:
            fragments.append(s)
            lfragment += len(s)
            if lfragment > random_access_points:
                rest = writeFragments( outfile_fasta, 
                                       outfile_index,
                                       fragments, 
                                       mangler, 
                                       size = random_access_points,
                                       write_all = False)
                fragments = [rest]
                lfragment = len(rest)
        else:
            outfile_fasta.write( mangler(s) )
                    
    if write_chunks:
        writeFragments( outfile_fasta, outfile_index, fragments, mangler, 
                        size = random_access_points, write_all = True )
    else:
        outfile_fasta.write( "\n" )
            
    outfile_index.write("\t%i\n" % lsequence )

    # add synonyms for the table
    if synonyms:
        for key, vals in synonyms.items():
            for val in vals:
                outfile_index.write( "%s\t%s\n" % (key, val) )

# map of names
# order is suffix data, suffix index, noSeek
NAME_MAP={
    'uncompressed' : ('fasta', 'idx', False),
    'lzo'          : ('lzo',   'cdx', True ),    
    'dictzip'      : ('dz',    'idx', False ),
    'zlib'         : ('zlib',  'cdx', True ),
    'gzip'         : ('gzip',  'cdx', True ),
    'bzip2'        : ('bz2',   'cdx', True ),
    'debug'        : ('debug', 'cdx', True ),
    }

PREFERENCES=('uncompressed', 'lzo', 'dictzip', 'zlib', 'gzip', 'bzip2', 'debug')

class IndexedFasta:
    """an indexed fasta file."""

    def __init__( self, dbname ):

        if dbname.endswith( ".fasta"): dbname = dbname[:-len(".fasta")]

        for x in PREFERENCES:
            d =  "%s.%s" % (dbname, NAME_MAP[x][0] )
            i =  "%s.%s" % (dbname, NAME_MAP[x][1] )
            if os.path.exists( d ) and os.path.exists( i ):
                self.mMethod = x
                self.mDbname = d
                self.mNameIndex = i
                self.mNoSeek = NAME_MAP[x][2]
                break
        else:
            raise KeyError, "unknown database %s" % dbname 
        
        self.mIsLoaded = False
        self.mSynonyms = {} 
        self.mConverter = None
        self.mIndex = {}
        self.mTranslator = None

    def __len__(self):
        """return the number of sequences in fasta file."""
        if not self.mIsLoaded: self.__loadIndex()
        return len(self.mIndex)

    def __contains__(self, contig):
        if not self.mIsLoaded: self.__loadIndex()
        return contig in self.mIndex or contig in self.mSynonyms

    def __getitem__(self, key ):
        """return full length sequence."""
        return self.getSequence( key, "+", 0, 0, as_array = True )
        
    def __loadIndex( self, compress=False ):
        """load complete index into memory.

        if compress is set to true, the index will not be loaded,
        but a compressed index will be created instead.
        """

        if self.mMethod == "uncompressed":
            self.mDatabaseFile = open( self.mDbname, "r" )
        elif self.mMethod == "dictzip":
            import dictzip
            self.mDatabaseFile = dictzip.GzipFile( self.mDbname)
        elif self.mMethod == "lzo":
            import lzo
            self.mDatabaseFile = Uncompressor( self.mDbname, lzo.decompress )
        elif self.mMethod == "gzip":
            self.mDatabaseFile = Uncompressor( self.mDbname, gzip_demangler )
        elif self.mMethod == "zlib":
            self.mDatabaseFile = Uncompressor( self.mDbname, zlib.decompress )
        elif self.mMethod == "bzip2":
            import bz2
            self.mDatabaseFile = Uncompressor( self.mDbname, bz2.decompress )
        elif self.mMethod == "debug":
            self.mDatabaseFile = Uncompressor( self.mDbname + ".debug", lambda x: x )            

        filename_index = self.mNameIndex + ".dbm"
        
        if compress:
            if os.path.exists( filename_index ):
                raise OSError( "file %s already exists" % filename_index )
            self.mIndex = anydbm.open( filename_index, "n" )
        elif os.path.exists( filename_index ):
            self.mIndex = anydbm.open( filename_index, "r" )
            self.mIsLoaded = True
            return
        else:
            self.mIndex = {}
            
        for line in open(self.mNameIndex, "r"):

            data = line[:-1].split("\t")

            if len(data) == 2:
                # ignore synonyms of non-existent contigs
                identifier = data[1]
                if data[0] not in self.mIndex: continue
                self.mSynonyms[identifier] = data[0]
            else:
                ## index with random access points
                if len(data) > 4:
                    (identifier, pos_id, block_size, lsequence) = data[0], int(data[1]), int(data[2]), int(data[-1])
                    points = map(int, data[3:-1])
                    self.mIndex[identifier] = (pos_id, block_size, lsequence, points)
                else:
                    (identifier, pos_id, pos_seq, lsequence) = data[0], int(data[1]), int(data[2]), int(data[-1])
                    self.mIndex[identifier] = struct.pack( "QQi", pos_id, pos_seq, lsequence)
                    
        # Treat common cases of naming incompatibilites like chr1 = 1. 
        # Truncate or add known prefixes
        def __addSynonyms( src, target ):
            for p in ("chr", "contig", "scaffold", "Chr" ):
                if src.startswith(p):
                    k = src[len(p):]
                    if k not in self.mIndex: self.mSynonyms[k]= target
                    # add lower/upper-case version
                    k = src[0].upper() + src[1:]
                    if k not in self.mIndex: self.mSynonyms[k]= target
                    k = src[0].lower() + src[1:]
                    if k not in self.mIndex: self.mSynonyms[k]= target
                    break
            else:
                for p in ("chr", "contig", "scaffold" ):
                    k = "%s%s" % (p,src)
                    if k not in self.mIndex: 
                        self.mSynonyms[k]= target

        k = self.mSynonyms.items()

        # fix the ambiguity between chrMT and chrM between UCSC and ENSEMBL
        if "chrM" in self.mIndex and "chrMT" not in self.mIndex:
            self.mSynonyms[ "chrMT" ] = "chrM"
        elif "chrM" not in self.mSynonyms and "chrMT" in self.mSynonyms:
            self.mSynonyms[ "chrM" ] = "chrMT"

        for key, val in k: __addSynonyms( key, val )

        # add pointers to self
        for key in self.mIndex.keys(): __addSynonyms( key, key )

        self.mIsLoaded = True

    def setTranslator(self, translator = None ):
        """set the :class:`Translator` to use."""
        self.mTranslator = translator
            
    def getDatabaseName( self ):
        """returns the name of the database."""
        return self.mDbname

    def getToken( self, contig ):
        """check if token is in index."""
        if not self.mIsLoaded: self.__loadIndex()

        if contig in self.mSynonyms:
            contig = self.mSynonyms[contig]

        if contig not in self.mIndex:
            raise KeyError, "%s not in index" % contig

        return contig

    def getLength( self, contig ):
        """return sequence length for sbjct_token."""
        if not self.mIsLoaded: self.__loadIndex()
        return struct.unpack( "QQi", self.mIndex[self.getToken(contig)] )[2]

    def getLengths( self ):
        """return all sequence lengths."""
        if not self.mIsLoaded: self.__loadIndex()
        return [ struct.unpack( "QQi", x)[2] for x in self.mIndex.values() ]

    def compressIndex( self ):
        """compress index.
        Creates a database interface to an index.
        """
        self.__loadIndex( compress=True )

    def getContigs( self ):
        """return a list of contigs (no synonyms)."""
        if not self.mIsLoaded: self.__loadIndex()
        return self.mIndex.keys()

    def getContigSizes( self, with_synonyms = True ):
        """return hash with contig sizes including synonyms."""
        if not self.mIsLoaded: self.__loadIndex()

        contig_sizes = {}
        for key, val in self.mIndex.items():
            contig_sizes[key] = self.getLength( key )

        if with_synonyms:
            for key, val in self.mSynonyms.items():
                contig_sizes[key] = self.getLength( val )

        return contig_sizes

    def setConverter( self, converter ):
        """set converter from coordinate system to 0-based, both strand, open/closed
        coordinate system."""
        self.mConverter = None

    def getSequence( self,
                     contig, 
                     strand = "+", 
                     start = 0, 
                     end = 0,
                     converter = None,
                     as_array = False):
        """get a genomic fragment.

        A genomic fragment is identified by the coordinates
        contig, strand, start, end.

        The converter function supplied translated these coordinates
        into 0-based coordinates.

        If as_array is set to true, return the AString object. This might
        be beneficial for large sequence chunks. If as_array is set to False,
        return a python string.
        """

        contig = self.getToken( contig )

        data = self.mIndex[contig]
        # dummy is
        # -> pos_seq for seekable streams
        # -> block_size for unseekable streams
        try:
            pos_id, dummy, lsequence = struct.unpack( "QQi", data )
        except struct.error:
            pos_id, dummy, lsequence, points = data

        pos_seq = dummy
        block_size = dummy
        
        if end == 0: end = lsequence
        
        if end > lsequence:
            raise ValueError("3' coordinate on %s out of bounds: %i > %i" % (contig, end, lsequence))
        if start < 0:
            raise ValueError("5' coordinate on %s out of bounds: %i < 0" % (contig, start))

        if converter:
            first_pos, last_pos = converter( start, end,
                                             str(strand) in ("+", "1"),
                                             lsequence )
        elif self.mConverter:
            first_pos, last_pos = self.mConverter( start, end,
                                                   str(strand) in ("+", "1"),
                                                   lsequence )
        else:
            first_pos, last_pos = start, end
            if str(strand) in ("-", "0", "-1"):
                first_pos, last_pos = lsequence - last_pos, lsequence - first_pos
                
        if first_pos == last_pos: return ""

        assert first_pos < last_pos, "first position %i is larger than last position %i " % (first_pos, last_pos)
        
        p = AString()
        
        if self.mNoSeek:
            ## read directly from position
            p.fromstring( self.mDatabaseFile.read( block_size, data[3], first_pos, last_pos) )
        else:
            first_pos += pos_seq
            last_pos += pos_seq

            self.mDatabaseFile.seek( first_pos )
            p.fromstring( self.mDatabaseFile.read( last_pos - first_pos ) )

        if str(strand) in ("-", "0", "-1"):
            p.reverse()            
            p = AString( string.translate( p[:],
                                           string.maketrans("ACGTacgt", "TGCAtgca") ) )

        if self.mTranslator:
            return self.mTranslator.translate( p )
        elif as_array:
            return p
        else:
            # cast to string
            return p[:]

    def getRandomCoordinates( self, size ):
        """returns coordinates for a random fragment of size #.

        The coordinates are forward/reverse.

        Default sampling mode:

        Each residue has the same probability of being
        in a fragment. Thus, the fragment can be smaller than
        size due to contig boundaries.
        """
        if not self.mIsLoaded: self.__loadIndex()

        token = random.choice( self.mIndex.keys() )        
        strand = random.choice( ("+", "-") )
        data = self.mIndex[token]
        pos_id, pos_seq, lcontig = struct.unpack( "QQi", data )
        rpos = random.randint( 0, lcontig )
        if size >= lcontig:
            start = 0
            end = lcontig
        else:
            if random.choice( ("True", "False") ):
                start = rpos
                end = min(rpos + size, lcontig)
            else:
                start = max(0, rpos - size)
                end = rpos
            
        return token, strand, start, end

###############################################################################
###############################################################################
###############################################################################
## converter functions. Some code duplication could be avoided but 
## I preferred to keep the functions lean.
###############################################################################        
def __one_forward_closed(x, y, c, l):
    """convert coordinates to zero-based, both strand, open/closed coordinates.
    
     Parameters are from, to, is_positive_strand, length of contig.
    """
    x -= 1
    if not c: x, y = l - y, l - x
    return x, y 
def __zero_forward_closed(x, y, c, l):
    """convert coordinates to zero-based, both strand, open/closed coordinates.
    
     Parameters are from, to, is_positive_strand, length of contig.
    """
    y += 1
    if not c: x, y = l - y, l - x
    return x, y 
def __one_both_closed(x, y, c = None, l = None):
    """convert coordinates to zero-based, both strand, open/closed coordinates.
    
     Parameters are from, to, is_positive_strand, length of contig.
    """
    return x - 1, y
def __zero_both_closed(x, y, c = None, l = None):
    """convert coordinates to zero-based, both strand, open/closed coordinates.
    
     Parameters are from, to, is_positive_strand, length of contig.
    """
    return x, y + 1

def __one_forward_open(x, y, c, l):
    """convert coordinates to zero-based, both strand, open/closed coordinates.
    
     Parameters are from, to, is_positive_strand, length of contig.
    """
    x -= 1
    y -= 1
    if not c: x, y = l - y, l - x
    return x, y 
def __zero_forward_open(x, y, c, l):
    """convert coordinates to zero-based, both strand, open/closed coordinates.
    
     Parameters are from, to, is_positive_strand, length of contig.
    """
    if not c: x, y = l - y, l - x
    return x, y 
def __one_both_open(x, y, c = None, l = None):
    """convert coordinates to zero-based, both strand, open/closed coordinates.
    
     Parameters are from, to, is_positive_strand, length of contig.
    """
    return x - 1, y - 1
def __zero_both_open(x, y, c = None, l = None):
    """convert coordinates to zero-based, both strand, open/closed coordinates.
    
    Parameters are from, to, is_positive_strand, length of contig.
    """
    return x, y

def getConverter( format ):
    """return a converter function for converting various
    coordinate schemes into 0-based, both strand, closed-open ranges.

    converter functions have the parameters
    x, y, s, l: with x and y the coordinates of
    a sequence fragment, s the strand (True is positive)
    and l being the length of the contig.

    Format is a "-" separated combination of the keywords
    "one", "zero", "forward", "both", "open", "closed"
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

## Test function for benchmarking purposes
def benchmarkRandomFragment( fasta, size ):
    """returns a random fragment of size."""

    contig, strand, start, end = fasta.getRandomCoordinates( size )
    s = fasta.getSequence( contig, strand, start, end )
    return s

def verify( fasta1, fasta2, num_iterations, fragment_size,
            stdout = sys.stdout, quiet = False ):
    """verify two databases.

    Get segment from fasta1 and check for presence in fasta2.
    """
    if not quiet:
        options.stdout.write("verifying %s and %s using %i random segments of length %i\n" %\
                             (fasta1.getDatabaseName(),
                              fasta2.getDatabaseName(),
                              num_iterations,
                              fragment_size ))
        options.stdout.flush()
    nerrors = 0
    for x in range(num_iterations):
        contig, strand, start, end = fasta1.getRandomCoordinates( fragment_size )
        s1 = fasta1.getSequence(contig,strand,start,end)
        s2 = fasta2.getSequence(contig,strand,start,end)
        if s1 != s2:
            if not quiet:
                options.stdout.write("discordant segment: %s:%s:%i:%i\n%s\n%s\n" %\
                                     (contig, strand, start, end, s1, s2) )
            nerrors += 1
    return nerrors

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

def parseCoordinates( s ):
    '''parse a coordinate string.'''

    if ":" in s:
        d = s.strip().split(":")
        if len(d) == 4:
            contig, strand, start, end = d
        elif len(d) == 2:
            contig = d[0]
            strand = "+"
            if ".." in d[1]:
                start, end = d[1].split("..")
            else:
                start = d[1]
    else:
        contig = s
        strand = "+"
        start = "0" 
        # full sequence
        end = "0"

    start = int( re.sub( ",", "", start) )
    if end: end = int( re.sub( ",", "", end) )
    else: end = start + 1
    
    return contig, strand, start, end

def main():

    import Experiment as E

    parser = optparse.OptionParser( version = "%prog version: $Id: IndexedFasta.py 2801 2009-10-22 13:40:39Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "-e", "--extract", dest="extract", type="string",
                       help="""extract region for testing purposes. Format is contig:strand:from:to. The default coordinates are
0-based open/closed coordinates on both strands. For example, chr1:+:10:12 will return bases 11 to 12 on chr1.""" )

    parser.add_option( "-c", "--compression", dest="compression", type="choice",
                       choices=("lzo", "zlib", "gzip", "dictzip", "bzip2", "debug"),
                       help="compress database [default=%default]." )

    parser.add_option( "--random-access-points", dest="random_access_points", type="int",
                       help="save random access points every # number of nucleotides [default=%default]." )

    parser.add_option( "-i", "--input-format", dest="input_format", type="choice",
                       choices=("one-forward-open", "zero-both-open" ),
                       help="coordinate format of input [default=%default]." )

    parser.add_option( "-s", "--synonyms", dest="synonyms", type="string",
                       help="list of synonyms, comma separated with =, for example, chr1=chr1b [default=%default]" )

    parser.add_option( "-b", "--benchmark", dest="benchmark", action="store_true",
                       help="benchmark time for read access [default=%default]." )
    
    parser.add_option( "--benchmark-num-iterations", dest="benchmark_num_iterations", type="int",
                       help="number of iterations for benchmark [default=%default]." )

    parser.add_option( "--benchmark-fragment-size", dest="benchmark_fragment_size", type="int",
                       help="benchmark: fragment size [default=%default]." )

    parser.add_option( "--verify", dest="verify", type="string",
                       help="verify against other database [default=%default].")

    parser.add_option( "-a", "--clean-sequence", dest="clean_sequence", action="store_true",
                       help="remove X/x from DNA sequences - they cause errors in exonerate [default=%default]." )

    parser.add_option( "--allow-duplicates", dest="allow_duplicates", action="store_true",
                       help="allow duplicate identifiers. Further occurances of an identifier are suffixed by an '_%i' [default=%default]." )

    parser.add_option( "--regex-identifier", dest="regex_identifier", type="string",
                       help="regular expression for extracting the identifier from fasta description line [default=%default]." )

    parser.add_option( "--compress-index", dest="compress_index", action="store_true",
                       help="compress index [default=%default]." )

    parser.add_option( "--force", dest="force", action="store_true",
                       help="force overwriting of existing files [default=%default]." )

    parser.add_option( "-t", "--translator", dest="translator", type="choice",
                       choices=("solexa", "phred", "bytes", "range200" ),
                       help="translate numerical quality scores [default=%default]." )


    parser.set_defaults(
        extract = None,
        input_format = "zero-both-open",
        benchmark_fragment_size = 1000,
        benchmark_num_iterations = 1000000,
        benchmark = False,
        compression = None,
        random_access_points = 0,
        synonyms = None,
        verify = None,
        verify_num_iterations = 100000,
        verify_fragment_size = 100,
        clean_sequence = False,
        allow_duplicates = False,
        regex_identifier = None,
        compress_index = False,
        force = False,
        translator = None )
    
    (options, args) = E.Start( parser )

    if options.synonyms:
        synonyms = {}
        for x in options.synonyms.split(","):
            a,b = x.split("=")
            a = a.strip()
            b = b.strip()
            if a not in synonyms: synonyms[a] = []
            synonyms[a].append( b )
    else:
        synonyms = None

    if options.translator:
        if options.translator == "phred":
            options.translator = TranslatorPhred()
        elif options.translator == "solexa":
            options.translator = TranslatorSolexa()
        elif options.translator == "bytes":
            options.translator = TranslatorBytes()
        elif options.translator == "range200":
            options.translator = TranslatorRange200()
        else:
            raise ValueError("unknown translator %s" % options.translator)
        
    if options.extract:
        fasta = IndexedFasta( args[0] )
        fasta.setTranslator( options.translator )
        converter = getConverter( options.input_format )
        
        contig, strand, start, end = parseCoordinates( options.extract )
        sequence = fasta.getSequence( contig, strand,
                                      start, end,
                                      converter = converter )
        options.stdout.write( ">%s\n%s\n" % \
                              ( options.extract, sequence ) )
    elif options.benchmark:
        import timeit
        timer = timeit.Timer( stmt="benchmarkRandomFragment( fasta = fasta, size = %i)" % (options.benchmark_fragment_size),
                              setup="""from __main__ import benchmarkRandomFragment,IndexedFasta\nfasta=IndexedFasta( "%s" )""" % (args[0] ) )

        t = timer.timeit( number = options.benchmark_num_iterations )
        options.stdout.write("iter\tsize\ttime\n" )
        options.stdout.write("%i\t%i\t%i\n" % (options.benchmark_num_iterations, options.benchmark_fragment_size, t ) )
    elif options.verify:
        fasta1 = IndexedFasta( args[0] ) 
        fasta2 = IndexedFasta( options.verify )
        nerrors1 = verify( fasta1, fasta2,
                           options.verify_num_iterations,
                           options.verify_fragment_size,
                           stdout=options.stdout )
        options.stdout.write("errors=%i\n" % (nerrors1) )        
        nerrors2 = verify( fasta2, fasta1,
                           options.verify_num_iterations,
                           options.verify_fragment_size,
                           stdout=options.stdout )
        options.stdout.write("errors=%i\n" % (nerrors2) )        
    elif options.compress_index:
        fasta = IndexedFasta( args[0] ) 
        fasta.compressIndex()
    else:
        if options.loglevel >= 1:
            options.stdlog.write("# creating database %s\n" % args[0])            
            options.stdlog.write("# indexing the following files: \n# %s\n" %\
                                 (" \n# ".join( args[1:] ) ))
            options.stdlog.flush()

            if synonyms:
                options.stdlog.write("# Applying the following synonyms:\n" )
                for k,v in synonyms.items():
                    options.stdlog.write( "# %s=%s\n" % (k, ",".join(v) ) )
                options.stdlog.flush()
        if len(args) < 2:
            print USAGE
            sys.exit(1)
            
        iterator = MultipleFastaIterator( args[1:],     
                                          regex_identifier = options.regex_identifier )
        createDatabase( args[0], 
                        iterator, 
                        synonyms = synonyms,
                        random_access_points = options.random_access_points,
                        compression = options.compression,
                        clean_sequence = options.clean_sequence,
                        allow_duplicates = options.allow_duplicates,
                        translator = options.translator,
                        force = options.force )

    
    E.Stop()

if __name__ == "__main__": main()

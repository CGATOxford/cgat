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
Fastq.py - methods for dealing with fastq files
===============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''

from math import log
import CGAT.Experiment as E
# see http://en.wikipedia.org/wiki/FASTQ_format
# ranges are conservative - they are open-ended
RANGES = {
    'sanger' : (33, 75) ,
    'illumina-1.8' : (33, 76) ,
    'solexa' : (59, 106) ,
    'phred64' : (64, 106) ,
    }

class Record:

    def __init__( self, identifier, seq, quals, format = None ):
        self.identifier, self.seq, self.quals, format = identifier, seq, quals, format
        self.format = None

    def __str__(self ):
        return "@%s\n%s\n+\n%s" % (self.identifier, self.seq, self.quals )

    def guessFormat( self ):
        '''return quality score format - 
        might return several if ambiguous.'''
        
        c = [ord(x) for x in self.quals ]
        mi, ma = min(c), max(c)
        r = []
        for format, v in RANGES.iteritems():
            m1, m2 = v
            if mi >= m1 and ma < m2: r.append( format )
        return r

    def trim( self, trim3, trim5 = 0 ):
        self.seq = self.seq[trim5:-trim3]
        self.quals = self.quals[trim5:-trim3]

    def toPhred( self ):
        '''return qualities as a list of phred-scores.'''
        assert self.format != None, "format needs to be set for conversion"
        if self.format == "sanger":
            return [ ord(x) - 33 for x in self.quals ]
        elif self.format == "illumina-1.8":
            return [ ord(x) - 33 for x in self.quals ]
        elif self.format == "solexa":
            # from -5 to 40 (i.e., can be negative)
            log10x = log(10.0) + .499
            return [ int(10.0 * log( 1.0 + 10 ** (ord(x) / 10.0), 10) / log10x ) for x in self.quals ]
        elif self.format == "phred64":
            return [ ord(x) - 64 for x in self.quals ]

    def fromPhred( self, quals, format ):
        '''set qualities from a list of phred-scores.'''
        self.format = format
        # -1 for color space fastq file
        assert len(quals) == len(self.seq) or len(quals) == len(self.seq) - 1
        if self.format == "sanger":
            self.quals = "".join( [ chr(33 + x) for x in quals ] )
        elif self.format == "illumina-1.8":
            return [ chr(33 + x) for x in quals ]
        elif self.format == "solexa":
            log10x = log(10.0, 10) / 10.0
            q = [int(10.0 * (log( 10 ** (x * log10x) - 1.0, 10 ))) for x in quals ] 
            self.quals = "".join( [ chr(64 + x) for x in q ] )
        elif self.format == "phred64":
            self.quals = "".join( [ chr(64 + x) for x in quals ] )
        elif self.format == "integer":
            self.quals = " ".join( map( str, quals ) )

def iterate( infile ):
    '''iterate over contents of fastq file.'''

    while 1:
        line1 = infile.readline()
        if not line1: break
        if not line1.startswith('@'):
            raise ValueError( "parsing error: expected '@' in line %s" % line1 )
        line2 = infile.readline()
        line3 = infile.readline()
        if not line3.startswith('+'):
            raise ValueError( "parsing error: expected '+' in line %s" % line3 )
        line4 = infile.readline()
        # incomplete entry
        if not line4: raise ValueError( "incomplete entry for %s" % line1)
        
        yield Record( line1[1:-1], line2[:-1], line4[:-1])


def iterate_guess( infile, max_tries = 10000, guess = None):
    '''iterate over contents of fastq file.

    guess quality format.
    '''
    quals = set( RANGES.keys() )
    cache = []
    myiter = iterate(infile)
    lengths = []
    for c, record in enumerate(myiter):
        quals.intersection_update( set(record.guessFormat()) )
        if len(quals) == 0:
            raise ValueError( "could not guess format - ranges incompatible." )
        if len(quals) == 1:
            break
        cache.append( record )
        lengths.append( len(record.seq) )
        if c > max_tries: 
            break

    if len(quals) == 1:
        ref_format = list(quals)[0]
    elif guess in quals:
        E.warn("multiple input formats possible: %s. Continuing with %s" % (", ".join(quals), guess))
        ref_format = guess
    else:
        raise ValueError( "could not guess format - could be one of %s." % str(quals) )

    for r in cache:
        r.format = ref_format
        yield r

    for r in myiter:
        r.format = ref_format
        yield r

def iterate_convert( infile, format, max_tries = 10000, guess = None ):
    '''iterate over contents of fastq file.

    guess quality format and set it to new format.
    '''
    
    quals = set( RANGES.keys() )
    cache = []
    myiter = iterate(infile)
    lengths = []
    for c, record in enumerate(myiter):
        quals.intersection_update( set(record.guessFormat()) )

        if len(quals) == 0:
            raise ValueError( "could not guess format - ranges incompatible." )
        if len(quals) == 1:
            break
        cache.append( record )

        if c > max_tries: 
            break

    if len(quals) == 1:
        ref_format = list(quals)[0]
    elif guess in quals:
        E.warn("multiple input formats possible: %s. Continuing with %s" % (", ".join(quals), guess))
        ref_format = guess
    else:
        raise ValueError( "could not guess format - could be one of %s. If you know the format use the --format option" % str(quals) )
            
    for r in cache:
        r.format = ref_format
        r.fromPhred( r.toPhred(), format )
        yield r

    for r in myiter:
        r.format = ref_format
        r.fromPhred( r.toPhred(), format )
        yield r
         

def guessFormat( infile, max_lines = 10000, raises = True ):
    '''guess format of FASTQ File.

    If *raises* , a ValueError is raised if there is not a single format.
    
    '''
    
    quals = set( RANGES.keys() )
    myiter = iterate(infile)
    for c, record in enumerate(myiter):
        quals.intersection_update( set(record.guessFormat()) )
        if len(quals) == 0:
            raise ValueError( "could not guess format - ranges incompatible." )
        if len(quals) == 1:
            break
        if c > max_lines: 
            break

    if len(quals) == 1:
        return list(quals)[0]
    elif raises == False:
        return quals
    else:
        raise ValueError( "could not guess format - could be one of %s." % str(quals) )


def getOffset( format, raises = True ):
    '''returns the ASCII offset for a certain format.

    If *raises* , a ValueError is raised if there is not a single offset.

    Otherwise, a minimum offset is returned.
    '''
    if type(format) in (set, list, tuple):
        offsets = set([ RANGES[f][0] for f in format ])
        if len(offsets) == 1:
            return list(offsets)[0]
        elif raises == False:
            return min(offsets)
        else:
            raise ValueError("inconsistent offsets for multiple formats: %s" % offsets)

    return RANGES[format][0]



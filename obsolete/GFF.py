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
"""
GFF.py - Classes and methods for dealing with gff formatted files
=================================================================

.. note::
   This module is deprecated - use fastgtf instead.

The default GFF version is now 2. This module needs updating to

   * deal with multiple values
   * correct quoting of free text values
   
Currently GFF and GTF is handled separately. Should be merged.

"""

import string, sys, re, types, collections
import IndexedGenome

def toDot( v ):
    '''convert value to '.' if None'''
    if v == None: return "." 
    else: return str(v)

def quote( v ):
    '''return a quoted attribute.'''
    if type(v) in types.StringTypes:
        return '"%s"' % v
    else: 
        return str(v)

class Entry:
    """read/write gff formatted entry.

    The coordinates are kept internally in python coordinates (0-based, open-closed), but are
    output as inclusive 1-based coordinates according to

    http://www.sanger.ac.uk/Software/formats/GFF/
    """

    def __init__(self):
        self.contig = "."
        self.source = "."
        self.feature = "."
        self.frame = "."
        self.start = 0
        self.end = 0
        self.score = "."
        self.strand = "."
        self.mAttributes = {}

    def read( self, line ):
        """read gff entry from line.
        
        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
        """
        
        data = line[:-1].split("\t")

        try:
            (self.contig, self.source, self.feature,
             self.start, self.end, self.score, self.strand,
             self.frame ) = data[:8]
        except ValueError, msg:
            raise ValueError( "parsing error in line %s: %s" % (line[:-1], msg) )
        
        ## note: frame might be .
        (self.start, self.end) = map(int, (self.start, self.end))
        self.start -= 1

        self.__parseInfo( data[8], line )

    def __parseInfo( self, attributes, line ):
        """parse attributes.
        """
        # remove comments
        attributes = attributes.split( "#" )[0]
        # separate into fields
        fields = map( lambda x: x.strip(), attributes.split(";")[:-1])
        self.mAttributes = {}

        for f in fields:
            
            d = map( lambda x: x.strip(), f.split(" "))
            n,v = d[0], d[1]
            if len(d) > 2: v = d[1:]

            if v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            else:
                ## try to convert to a value
                try:
                    v = float( v )
                    v = int( v )
                except ValueError:
                    pass
                except TypeError:
                    pass

            self.mAttributes[n] = v

    def __str__(self):

        aa = []
        for k,v in self.mAttributes.items():
            if type(v) == types.StringType:
                aa.append( '%s "%s"' % (k,v) )
            else:
                aa.append( '%s %s' % (k,str(v)) )

        attributes = "; ".join( aa )
        if aa: attributes += ";"
        
        return "\t".join(map(str, (self.contig, self.source,
                                     self.feature, 
                                     self.start + 1, self.end,
                                     self.score, self.strand, self.frame,
                                     attributes))) 
        
        return cmp( (self.contig, self.strand, self.start),
                    (other.contig, other.strand, other.start))

    def invert( self, lcontig ):
        """invert genomic coordinates from
        forward to reverse coordinates and back.
        """

        if self.strand in ("-", "0", "-1"):
            x = min(self.start, self.end)
            y = max(self.start, self.end)
                
            self.start = lcontig - y
            self.end   = lcontig - x

    def copy( self, other ):
        """fill from other entry."""
        self.contig = other.contig
        self.source = other.source
        self.feature = other.feature
        self.start = other.start
        self.end = other.end
        self.score = other.score
        self.strand = other.strand
        self.frame = other.frame
        self.mAttributes = copy.copy(other.mAttributes)
        
    def clearAttributes( self ):
        self.mAttributes = {}

    def addAttribute( self, key, value = None):
        self.mAttributes[key] = value

    def __cmp__(self, other):
        ## note: does compare by strand as well!
        return cmp( (self.contig, self.strand, self.start),
                    (other.contig, other.strand, other.start))

    def __setitem__(self, key, value):
        self.addAttribute( key, value )
        
    def __getitem__(self, key):
        return self.mAttributes[key]

    def __contains__(self, key):
        return key in self.mAttributes
        
    def hasOverlap( self, other, min_overlap = 0 ):
        """returns true, if overlap with other entry.
        """
        return (self.contig == other.contig and self.strand == other.strand and \
                    min(self.end, other.end) - max(self.start, other.start) > min_overlap)

    def isIdentical( self, other, max_slippage = 0 ):
        """returns true, if self and other overlap completely.
        """
        return (self.contig == other.contig and \
                    self.strand == other.strand and \
                    abs(self.end - other.end) < max_slippage and \
                    abs(self.start - other.start < max_slippage) )

    def isHalfIdentical( self, other, max_slippage = 0 ):
        """returns true, if self and other overlap.
        """
        return (self.contig == other.contig and \
                    self.strand == other.strand and \
                    ( abs(self.end - other.end) < max_slippage or \
                          abs(self.start - other.start < max_slippage) ) )
        
def Overlap( entry1, entry2, min_overlap = 0 ):
    """returns true, if entry1 and entry2 overlap.
    """

    return (entry1.contig == entry2.contig and entry1.strand == entry2.strand and \
            min(entry1.end, entry2.end) - max(entry1.start, entry2.start) > min_overlap)

def Identity( entry1, entry2, max_slippage = 0 ):
    """returns true, if entry1 and entry2 overlap.
    """

    return (entry1.contig == entry2.contig and \
            entry1.strand == entry2.strand and \
            abs(entry1.end - entry2.end) < max_slippage and \
            abs(entry1.start - entry2.start < max_slippage) )

def HalfIdentity( entry1, entry2, max_slippage = 0 ):
    """returns true, if entry1 and entry2 overlap.
    """

    return (entry1.contig == entry2.contig and \
            entry1.strand == entry2.strand and \
            ( abs(entry1.end - entry2.end) < max_slippage or \
              abs(entry1.start - entry2.start < max_slippage) ) )

def readFromFile( infile ):
    """read gtf from file."""
    result = []
    for gff in iterator(infile):
        result.append( gff )
    return result

def CombineOverlaps( old_gff, method = "combine" ):
    """combine overlapping entries for a list of gffs.

    method can be any of combine|longest|shortest
    only the first letter is important.
    """

    old_gff.sort( lambda x,y: cmp( (x.contig, x.strand, x.start, x.end),
                                   (y.contig, y.strand, y.start, y.end) ) )
    
    new_gff = []

    last_e = old_gff[0]

    for e in old_gff[1:]:
        if not Overlap( last_e, e):
            new_gff.append( last_e )
            last_e = e
        else:
            if method[0] == "c":
                last_e.start = min( last_e.start, e.start )
                last_e.end = max( last_e.end, e.end )
                last_e.mInfo += " ; " + e.mInfo

    new_gff.append( last_e )
    
    return new_gff
    
def SortPerContig( gff ):
    """sort gff entries per contig and return a dictionary mapping a
    contig to the begin of the list."""
    map_contig2start = {}

    if len(gff) == 0: return map_contig2start
    
    gff.sort( lambda x,y: cmp( x.contig, y.contig) )

    last_contig = None
    start = 0
    for x in range( len(gff)):
        if last_contig != gff[x].contig:
            map_contig2start[last_contig] = (start, x)
            start = x
            last_contig = gff[x].contig

    map_contig2start[last_contig] = (start, x)

    return map_contig2start

def readAsIntervals( gff_iterator, with_values = False, with_records = False,
                     use_strand = False ):
    """read tuples of (start, end) from a GFF file.

    If with_values is True, a value is added to the tuples.
    If with_records is True, the record is added to the tuples.
    with_values and with_records are exclusive.
    
    If use_strand is True, intervals will be grouped by contig and strand.
    The default is to group by contig only.

    Returns a dictionary of intervals by contig.
    """

    assert not (with_values and with_records), "both with_values and with_records are true."
    intervals = collections.defaultdict(list)
    if use_strand:
        keyf = lambda x: (x.contig, x.strand)
    else:
        keyf = lambda x: x.contig

    if with_values:
        for gff in gff_iterator:
            intervals[keyf(gff)].append( (gff.start,gff.end,gff.score) )
    elif with_records:
        for gff in gff_iterator:
            intervals[keyf(gff)].append( (gff.start,gff.end,gff) )
    else:
        for gff in gff_iterator:
            intervals[keyf(gff)].append( (gff.start,gff.end) )
    return intervals

def iterator( infile ):
    """a simple iterator over all entries in a file."""
    ntracks = 0
    while 1:
        line = infile.readline()
        if not line: return
        if line[0] == "#": continue
        if len(line.strip()) == 0: continue
        if line.startswith("track"):
            ntracks += 1
            if ntracks > 1: raise ValueError( "more than one track in %s" % infile)
            continue
        gff = Entry()
        gff.read( line )
        yield gff

def chunk_iterator( gffs ):
    """iterate over the contents of a gff file.

    return entries as single element lists
    """
    for gff in gffs:
        yield [gff,]

def iterator_contigs( gffs ):
    """iterate over contigs.
    
    TODO: implement as coroutines
    """

    last_contig, data = None, []
    found = set()
    for gff in gffs:
        if last_contig != gff.contig:
            if last_contig: yield last_contig, data
            last_contig = gff.contig
            assert last_contig not in found, "input not sorted by contig."
            found.add(last_contig) 

        data.append( gff )

    if last_contig: yield last_contig, data

def joined_iterator(gffs, group_field = None):
    """iterate over the contents of a gff file.

    return a list of entries with the same group id.
    Note: the entries have to be consecutive in
    the file.
    """
    last_group_id = None
    matches = []

    for gff in gffs:
        
        if group_field:
            group_id = gff.mFields[group_field]
        else:
            group_id = gff.mAttributes
                
        if last_group_id != group_id:
            if last_group_id: yield matches
            matches = []
            last_group_id = group_id
                
        matches.append( gff )

    if last_group_id: yield matches

def iterator_filtered( gff_iterator, feature = None, source = None, contig = None, interval = None, strand = None):
    """iterate over the contents of a gff file.

    yield only entries for a given feature
    """
    if interval: start, end = interval
    if strand == ".": strand = None

    for gff in gff_iterator:
        if feature and gff.feature != feature: continue
        if source and gff.source != source: continue
        if contig and gff.contig != contig: continue
        if strand and gff.strand != strand: continue
        if interval and min(end,gff.end) - max(start,gff.start) < 0: continue
        yield gff

def iterator_overlaps( gff_iterator, min_overlap = 0 ):
    """iterate over gff file and return a list of features that
    are overlapping.

    The input should be sorted by contig,start
    """
    
    last = gff_iterator.next()
    matches = [last]
    end = last.end
    for this in gff_iterator:
        if last.contig != this.contig or \
                end - this.start <= min_overlap:
            yield matches
            matches = [this]
            end = this.end
            last = this
            continue

        assert last.start <= this.start, "input file needs to be sorted by contig, start:\n%s\n%s\n" % (str(last), str(this) )
        matches.append( this )
        last = this
        end = max( end, this.end )

    yield matches

def readAndIndex( iterator, with_value = True ):
    '''read from gtf stream and index.

    returns an :class:`IndexedGenome.IndexedGenome`
    '''

    if with_value:
        index = IndexedGenome.IndexedGenome()
        for gtf in iterator:
            index.add( gtf.contig, gtf.start, gtf.end, gtf )
    else:
        index = IndexedGenome.Simple()
        for gtf in iterator:
            index.add( gtf.contig, gtf.start, gtf.end )

    return index


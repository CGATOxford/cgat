"""Classes and methods for dealing with gtf formatted files.

The default GTF version is 2.2.
"""

import string, sys, re, copy, types, collections
import Intervals
import Genomics
import GFF
import IndexedGenome
import fastgtf

class Error(Exception):
    """Base class for exceptions in this module."""
    def __str__(self):
        return str(self.message)
    def _get_message(self, message): return self._message
    def _set_message(self, message): self._message = message
    message = property(_get_message, _set_message)

class ParsingError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class Entry:
    """read/write gtf formatted entry.

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
        self.gene_id = None
        self.transcript_id =  None
        self.attributes = {}

    def read( self, line ):
        """read gff entry from line.
        
        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
        """
        
        data = line[:-1].split("\t")

        try:
            (self.contig, self.source, self.feature,
             self.start, self.end, self.score, self.strand,
             self.frame ) = data[:8]
        except ValueError:
            raise ValueError( "parsing error in line `%s`" % line )

        ## note: frame might be .
        (self.start, self.end) = map(int, (self.start, self.end))
        self.start -= 1

        self.parseInfo( data[8], line )

    def parseInfo( self, attributes, line ):
        """parse attributes.
        """
        # remove comments
        attributes = attributes.split( "#" )[0]
        # separate into fields
        fields = map( lambda x: x.strip(), attributes.split(";")[:-1])
        self.attributes = {}
        
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
                
            if n == "gene_id": 
                self.gene_id = v     
            elif n == "transcript_id": 
                self.transcript_id = v
            else: 
                self.attributes[n] = v

        if not self.gene_id:
            raise ParsingError( "missing attribute 'gene_id' in line %s" % line)
        if not self.transcript_id:
            raise ParsingError( "missing attribute 'transcript_id' in line %s" % line)

    def getAttributeField( self, full = True ):
        aa = []
        for k,v in self.attributes.items():
            if type(v) == types.StringType:
                aa.append( '%s "%s"' % (k,v) )
            else:
                aa.append( '%s %s' % (k,str(v)) )

        if full:
            return"; ".join( ['gene_id "%s"' % self.gene_id, 
                                     'transcript_id "%s"' % self.transcript_id ] +\
                                        aa )
        else:
            return "; ".join( aa )

        return attributes

    def __cmp__(self, other):
        ## note: does compare by strand as well!
        return cmp( (self.contig, self.strand, self.start),
                    (other.contig, other.strand, other.start))

    def __str__(self):
        
        def _todot( val ):
            if val == None: return "."
            else: return str(val)

        attributes = self.getAttributeField( )
        return "\t".join(map(str, (self.contig, self.source,
                                   self.feature, 
                                   self.start + 1, self.end,
                                   _todot(self.score), 
                                   self.strand, 
                                   self.frame,
                                   attributes))) + ";" 

    def invert( self, lcontig ):
        """invert genomic coordinates from
        forward to reverse coordinates and back.
        """

        if self.strand in ("-", "0", "-1"):
            x = min(self.start, self.end)
            y = max(self.start, self.end)
                
            self.start = lcontig - y
            self.end   = lcontig - x

    def fromGFF( self, other, gene_id, transcript_id ):
        """fill from other entry."""
        self.contig = other.contig
        self.source = other.source
        self.feature = other.feature
        self.start = other.start
        self.end = other.end
        self.score = other.score
        self.strand = other.strand
        self.frame = other.frame
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        return self

    def copy( self, other ):
        """fill from other entry.
        works both if other is :class:`GTF.Entry` or 
        :class:`fastgtf.GTFProxy`
        """
        self.contig = other.contig
        self.source = other.source
        self.feature = other.feature
        self.start = other.start
        self.end = other.end
        self.score = other.score
        self.strand = other.strand
        self.frame = other.frame
        self.gene_id = other.gene_id
        self.transcript_id = other.transcript_id
        self.attributes = copy.copy(other.asDict())
        return self

    def asDict( self ):
        '''return attributes as a dictionary.'''
        return self.attributes

    def clearAttributes( self ):
        self.attributes = {}

    def addAttribute( self, key, value = None):
        self.attributes[key] = value

    def __setitem__(self, key, value):
        self.addAttribute( key, value )
        
    def __getitem__(self, key):
        return self.attributes[key]

    def __contains__(self, key):
        return key in self.attributes
        
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

def asRanges( gffs, feature = None ):
    """return ranges within a set of gffs.

    Overlapping intervals are merged.
    """

    if type(feature) == types.StringType:
        gg = filter( lambda x: x.feature == feature, gffs )
    elif feature:
        gg = filter( lambda x: x.feature in feature, gffs )
    else:
        gg = gffs[:]
    
    r = [ (g.start, g.end) for g in gg ]
    return Intervals.combine( r )

def readFromFile( infile ):
    """read gtf from file."""
    result = []
    for gff in fastgtf.iterator( infile ):
        result.append( gff )
    return result

def iterator( infile ):
    """return a simple iterator over all entries in a file."""
    return fastgtf.iterator(infile)

def chunk_iterator( gff_iterator ):
    """iterate over the contents of a gff file.

    return entries as single element lists
    """
    for gff in gffs_iterator:
        yield [gff,]

def transcript_iterator(gff_iterator, strict = True ):
    """iterate over the contents of a gtf file.

    return a list of entries with the same transcript id.

    Note: the entries for the same transcript have to be consecutive 
    in the file.
    """
    last = None
    matches = []
    found = set()

    for gff in gff_iterator:
        this = gff.transcript_id + gff.gene_id
        if last != this:
            if last: yield matches
            matches = []
            assert not strict or this not in found, "duplicate entry: %s" % this
            found.add( this )
            last = this
        matches.append( gff )

    if last: yield matches

def gene_iterator( gff_iterator, strict = True ):
    """iterate over the contents of a gtf file.

    return a list of transcripts with the same gene id. 

    Note: the entries have to be consecutive in the file, i.e,
    first sorted by transcript and then by gene id.
    """
    last = None
    matches = []
    found = set()
    for gffs in transcript_iterator( gff_iterator, strict ):
        
        if last != gffs[0].gene_id:
            if last: yield matches
            matches = []
            last = gffs[0].gene_id
            assert not strict or last not in found, "duplicate entry %s" % last 
            found.add(last)

        matches.append( gffs )

    if last: yield matches

def flat_gene_iterator( gff_iterator, strict = True ):
    """iterate over the contents of a gtf file.

    return a list of entries with the same gene id.

    Note: the entries have to be consecutive in the file, i.e,
    sorted by gene_id
    """

    last = None
    matches = []
    found = set()
    for gff in gff_iterator:
        
        if last != gff.gene_id:
            if last: yield matches
            matches = []
            last = gff.gene_id
            assert not strict or last not in found, "duplicate entry %s" % last 
            found.add(last)
        matches.append( gff )

    if last: yield matches

def merged_gene_iterator( gff_iterator ):
    """iterate over the contents of a gtf file.

    Each gene is merged into a single entry spanning the whole
    stretch that a gene covers.

    Note: the entries have to be consecutive in the file, i.e,
    sorted by gene_id
    """
    for m in flat_gene_iterator( gff_iterator ):
        gff = Entry()
        gff.copy( m[0] )
        gff.start = min( [x.start for x in m ] )
        gff.end = max( [x.end for x in m ] )
        yield gff

def iterator_filtered( gff_iterator, feature = None, source = None):
    """iterate over the contents of a gff file.

    yield only entries for a given feature
    """
    for gff in gff_iterator:
        if feature and gff.feature != feature: continue
        if source and gff.source != source: continue
        yield gff

def iterator_sorted_chunks( gff_iterator, sort_by = "location-start" ):
    """iterate over chunks in a sorted order

    sort_by can be "location-start"
    """
    
    ## get all chunks and annotate with sort order
    if sort_by == "location-start":
        chunks = ([ (x[0].contig, min( [y.start for y in x] ), x) for x in gff_iterator ] )
    else:
        raise ValueError( "unknown sort order %s" % sort_by )
    chunks.sort()
    for contig, start, chunk in chunks:
        yield chunk

def iterator_min_feature_length( gff_iterator, min_length, feature="exon" ):
    """select only those genes with a minimum length of a given feature."""
    for gffs in gff_iterator:
        intervals = [ (x.start, x.end) for x in gffs if x.feature == feature ]
        intervals = Intervals.combine( intervals )
        t = sum( ( x[1] - x[0] for x in intervals ) )
        if t >= min_length: yield gffs

def toIntronIntervals( chunk ):
    '''convert a set of gtf elements within a transcript to intron coordinates.

    Will raise an error if more than one transcript is submitted.
    
    Note that coordinates will still be forward strand coordinates
    '''
    if len(chunk) == 0: return []
    t = set([ x.transcript_id for x in chunk ])
    contig, strand, transcript_id = chunk[0].contig, chunk[0].strand, chunk[0].transcript_id
    for gff in chunk:
        assert gff.strand == strand, "features on different strands."
        assert gff.contig == contig, "features on different contigs."
        assert gff.transcript_id == transcript_id, "more than one transcript submitted"

    intervals = Intervals.combine( [ (x.start, x.end) for x in chunk ] )
    return Intervals.complement( intervals )
    

    
def toSequence( chunk, fasta ):
    """convert a list of gff attributes to a single sequence.
    
    This function ensures correct in-order concatenation on
    positive/negative strand. Overlapping regions are merged.
    """
    if len(chunk) == 0: return ""

    contig, strand = chunk[0].contig, chunk[0].strand 

    for gff in chunk:
        assert gff.strand == strand, "features on different strands."
        assert gff.contig == contig, "features on different contigs."
    
    intervals = Intervals.combine( [ (x.start, x.end) for x in chunk ] )
    lcontig = fasta.getLength(contig)
    positive = Genomics.IsPositiveStrand( strand )

    if not positive: 
        intervals = [ (lcontig - end, lcontig - start) for start,end in intervals ]
        intervals.reverse()

    s = [ fasta.getSequence( contig, strand, start, end ) for start, end in intervals ]

    return "".join(s)

##----------------------------------------------------------------------------------------
def iterator_overlapping_genes( gtf_iterator, min_overlap = 0 ):
    """return overlapping genes."""
    for gene in flat_gene_iterator(gtf_iterator):
        gene.sort( key=lambda x: x.start )
        genes.append( gene[0].contig, gene[0].start, gene )
        
    genes.sort()

    contig, last_end, ovl = None, 0, 0, []

    for this_contig, this_start, gene in genes:
        if this_contig != last_contig or last_end < this_start:
            if last_contig: yield ovl
            ovl, last_contig = [], this_contig
            last_end = 0

        last_end = max( last_end, gene[-1].end )
       
    if last_contig: yield ovl

##----------------------------------------------------------------------------------------
def iterator_transcripts2genes( gtf_iterator, min_overlap = 0 ):
    """cluster transcripts by exon overlap.

    The gene id is set to the first transcript encountered of a gene.
    If a gene stretches over several contigs, subsequent copies are
    appended a number.
    """

    map_transcript2gene = {}
    gene_ids = collections.defaultdict( list )

    for chunk in GFF.iterator_overlaps( gtf_iterator ):
        transcript_ids = list(set( [x.transcript_id for x in chunk ] ))
        contig = chunk[0].contig

        # have any of these already encountered?
        for x in transcript_ids:
            if x in map_transcript2gene:
                gene_id = map_transcript2gene[x]
                break
        else:
            # arbitrarily pick one
            gene_id = transcript_ids[0]

        if gene_id not in gene_ids:
            gene_ids[gene_id].append( contig )
            index = 0
        else:
            try:
                index = gene_ids[gene_id].index(contig)
            except ValueError:
                index = len(gene_ids[gene_id])
                gene_ids[gene_id].append( contig )

        for x in transcript_ids: map_transcript2gene[x] = gene_id
        
        if index: gene_id += ".%i" % index
        for x in chunk: x.gene_id = gene_id
        yield chunk

def readAsIntervals( gff_iterator, 
                     with_values = False, 
                     with_records = False,
                     merge_genes = False,
                     with_gene_id = False,
                     with_transcript_id = False ):
    """read tuples of (start, end) from a GTF file.

    If with_values is True, a value is added to the tuples.
    If with_records is True, the record is added to the tuples.
    If with_gene_id is True, the gene_id is added to the tuples
    with_values and with_records are exclusive.
    
    Ignores strand and everything but the coordinates.

    Returns a dictionary of intervals by contig.
    """

    assert not (with_values and with_records), "both with_values and with_records are true."
    intervals = {}

    if merge_genes:
        it = merged_gene_iterator( gff_iterator )
    else:
        it = gff_iterator
        
    if with_values:
        for gff in it:
            contig, start, end = gff.contig, gff.start, gff.end
            if contig not in intervals: intervals[contig] = []
            intervals[contig].append( (start,end,gff.score) )
    elif with_records:
        for gff in it:
            contig, start, end = gff.contig, gff.start, gff.end
            if contig not in intervals: intervals[contig] = []
            intervals[contig].append( (start,end,gff) )
    elif with_gene_id:
        for gff in it:
            contig, start, end = gff.contig, gff.start, gff.end
            if contig not in intervals: intervals[contig] = []
            intervals[contig].append( (start,end,gff.gene_id) )
    elif with_transcript_id:
        for gff in it:
            contig, start, end = gff.contig, gff.start, gff.end
            if contig not in intervals: intervals[contig] = []
            intervals[contig].append( (start,end,gff.transcript_id) )
    else:
        for gff in gff_iterator:
            contig, start, end = gff.contig, gff.start, gff.end
            if contig not in intervals: intervals[contig] = []
            intervals[contig].append( (start,end) )
        
    return intervals

def readAndIndex( iterator ):
    '''read from gtf stream and index.

    returns an :class:`IndexedGenome.IndexedGenome`
    '''

    index = IndexedGenome.IndexedGenome()
    for gtf in iterator:
        index.add( gtf.contig, gtf.start, gtf.end, gtf )

    return index

"""GTF.py - Classes and methods for dealing with GTF/GFF formatted files
=====================================================================

The coordinates are kept internally in python coordinates (0-based,
open-closed), but are output as inclusive 1-based coordinates
according to http://www.sanger.ac.uk/Software/formats/GFF/.

The default GTF version is 2.2.

This module uses pysam_ to provide the principal engine for iterating over
files (:func:`iterate`). As a consequence, the returned objects are of
type :func:`pysam.GTFProxy`.

The class defined in this model :class:`Entry` is useful for re-formatting
records.

Apart from basic iteration, this module provides the following utilities:

* Additional iterators for grouping/modifying :term:`GTF` formatted files:
  :func:`track_iterator`, :func:`chunk_iterator`, :func:`iterator_contigs`,
  :func:`transcript_iterator`, :func:`joined_iterator`, :func:`gene_iterator`,
  :func:`flat_gene_iterator`, :func:`merged_gene_iterator`,
  :func:`iterator_filtered`, :func:`iterator_sorted_chunks`,
  :func:`iterator_min_feature_length`, :func:`iterator_sorted`
  :func:`iterator_overlapping_genes`, :func:`iterator_transcripts2genes`
  :func:`iterator_overlaps`

* Compare intervals: :func:`Identity`, :func:`HalfIdentity`, :func:`Overlap`

* Read GTF formatted files and optionally index them: :func:`readFromFile`,
  :func:`readAsIntervals`, :func:`readAndIndex`

* Manipulate lists of GTF records: :func:`asRanges`, :func:`CombineOverlaps`,
  :func:`SortPerContig`, :func:`toIntronIntervals`, :func:`toSequence`

"""

import collections
from CGAT import Intervals as Intervals
from CGAT import Genomics as Genomics
from CGAT import IndexedGenome as IndexedGenome
import pysam
from CGAT import IOTools as IOTools


def iterator(infile):
    """return a simple iterator over all entries in a file."""
    return pysam.tabix_iterator(infile, pysam.asGTF())


def track_iterator(infile):
    """a simple iterator over all entries in a file."""
    # note: taken from GFF.py
    ntracks = 0
    while 1:
        line = infile.readline()
        if not line:
            return
        if line[0] == "#":
            continue
        if len(line.strip()) == 0:
            continue
        if line.startswith("track"):
            ntracks += 1
            if ntracks > 1:
                raise ValueError("more than one track in %s" % infile)
            continue
        gff = Entry()
        gff.read(line)
        yield gff


def chunk_iterator(gff_iterator):
    """iterate over the contents of a gff file.

    return entries as single element lists
    """
    for gff in gff_iterator:
        yield [gff, ]


def iterator_contigs(gffs):
    """iterate over contigs.

    TODO: implement as coroutines
    """

    last_contig, data = None, []
    found = set()
    for gff in gffs:
        if last_contig != gff.contig:
            if last_contig:
                yield last_contig, data
            last_contig = gff.contig
            assert last_contig not in found, "input not sorted by contig."
            found.add(last_contig)

        data.append(gff)

    if last_contig:
        yield last_contig, data


def transcript_iterator(gff_iterator, strict=True):
    """iterate over the contents of a gtf file.

    return a list of entries with the same transcript id.

    Any features without a transcript_id will be ignored.

    The entries for the same transcript have to be consecutive
    in the file. If *strict* is set an AssertionError will be
    raised if that is not true.
    """
    last = None
    matches = []
    found = set()

    for gff in gff_iterator:
        # ignore entries without transcript or gene id
        try:
            this = gff.transcript_id + gff.gene_id
        except KeyError:
            continue
        except AttributeError:
            continue

        if last != this:
            if last:
                yield matches
            matches = []
            assert not strict or this not in found, \
                "duplicate entry: %s" % this
            found.add(this)
            last = this
        matches.append(gff)
    if last:
        yield matches


def joined_iterator(gff_iterator, group_field=None):
    """iterate over the contents of a gff file.

    return a list of entries with the same group id.
    Note: the entries have to be consecutive in
    the file.
    """
    last_group_id = None
    matches = []

    if group_field is None:
        group_function = lambda x: x.attributes
    elif group_field == "gene_id":
        group_function = lambda x: x.gene_id
    elif group_field == "transcript_id":
        group_function = lambda x: x.transcript_id
    else:
        group_function = lambda x: x[group_field]

    for gff in gff_iterator:

        group_id = group_function(gff)

        if last_group_id != group_id:
            if last_group_id:
                yield matches
            matches = []
            last_group_id = group_id

        matches.append(gff)

    if last_group_id:
        yield matches


def gene_iterator(gff_iterator, strict=True):
    """iterate over the contents of a gtf file.

    return a list of transcripts with the same gene id.

    Note: the entries have to be consecutive in the file, i.e,
    first sorted by transcript and then by gene id.

    Genes with the same name on different contigs are resolved
    separately in *strict* = False.
    """
    last = Entry()
    matches = []
    found = set()
    for gffs in transcript_iterator(gff_iterator, strict):

        if last.gene_id != gffs[0].gene_id or last.contig != gffs[0].contig:
            if matches:
                yield matches
            matches = []
            last = gffs[0]
            assert not strict or last.gene_id not in found, \
                "duplicate entry %s" % last
            found.add(last.gene_id)

        matches.append(gffs)

    if matches:
        yield matches


def flat_gene_iterator(gff_iterator, strict=True):
    """iterate over the contents of a gtf file.

    return a list of entries with the same gene id.

    Note: the entries have to be consecutive in the file, i.e,
    sorted by gene_id

    Genes with the same name on different contigs are resolved
    separately in *strict* = False

    """

    last = Entry()
    matches = []
    found = set()
    for gff in gff_iterator:

        if last.gene_id != gff.gene_id or last.contig != gff.contig:
            if matches:
                yield matches
            matches = []
            last = gff
            if strict:
                assert last.gene_id not in found, \
                    "duplicate entry %s" % last.gene_id
                found.add(last.gene_id)
        matches.append(gff)

    if matches:
        yield matches


def merged_gene_iterator(gff_iterator):
    """iterate over the contents of a gtf file.

    Each gene is merged into a single entry spanning the whole
    stretch that a gene covers.

    Note: the entries have to be consecutive in the file, i.e,
    sorted by gene_id
    """
    for m in flat_gene_iterator(gff_iterator):
        gff = Entry()
        gff.copy(m[0])
        gff.start = min([x.start for x in m])
        gff.end = max([x.end for x in m])
        yield gff


def iterator_filtered(gff_iterator, feature=None,
                      source=None, contig=None,
                      interval=None, strand=None):
    """iterate over the contents of a gff file.

    yield only entries for a given feature
    """
    if interval:
        start, end = interval
    if strand == ".":
        strand = None

    for gff in gff_iterator:
        if feature and gff.feature != feature:
            continue
        if source and gff.source != source:
            continue
        if contig and gff.contig != contig:
            continue
        if strand and gff.strand != strand:
            continue
        if interval and min(end, gff.end) - max(start, gff.start) < 0:
            continue
        yield gff


def iterator_sorted_chunks(gff_iterator, sort_by="contig-start"):
    """iterate over chunks in a sorted order

    sort_by can be

    contig-start
       sort by position ignoring the strand
    contig-strand-start
       sort by position taking the strand into account
    contig-strand-start-end
       intervals with the same start position will be sorted by end position

    returns the chunks.
    """

    # get all chunks and annotate with sort order
    if sort_by == "contig-start":
        chunks = ([(x[0].contig, min([y.start for y in x]), x)
                   for x in gff_iterator])
        chunks.sort()
        for contig, start, chunk in chunks:
            chunk.sort(key=lambda x: (x.contig, x.start))
            yield chunk
    elif sort_by == "contig-strand-start":
        chunks = ([(x[0].contig, x[0].strand, min([y.start for y in x]), x)
                   for x in gff_iterator])
        chunks.sort()
        for contig, start, strand, chunk in chunks:
            chunk.sort(key=lambda x: (x.contig, x.strand, x.start))
            yield chunk
    elif sort_by == "contig-strand-start-end":
        chunks = ([(x[0].contig, x[0].strand, min([y.start for y in x]), x)
                   for x in gff_iterator])
        chunks.sort()
        for contig, start, strand, chunk in chunks:
            chunk.sort(key=lambda x: (x.contig, x.strand, x.start, x.end))
            yield chunk
    else:
        raise ValueError("unknown sort order %s" % sort_by)


def iterator_min_feature_length(gff_iterator, min_length, feature="exon"):
    """select only those genes with a minimum length of a given feature."""
    for gffs in gff_iterator:
        intervals = [(x.start, x.end) for x in gffs if x.feature == feature]
        intervals = Intervals.combine(intervals)
        t = sum((x[1] - x[0] for x in intervals))
        if t >= min_length:
            yield gffs


def iterator_sorted(gff_iterator, sort_order="gene"):
    '''sort input and yield sorted output.'''
    entries = list(gff_iterator)
    if sort_order == "gene":
        entries.sort(key=lambda x: (x.gene_id, x.contig, x.start))
    elif sort_order == "gene+transcript":
        entries.sort(
            key=lambda x: (x.gene_id, x.transcript_id, x.contig, x.start))
    elif sort_order == "contig+gene":
        entries.sort(
            key=lambda x: (x.contig, x.gene_id, x.transcript_id, x.start))
    elif sort_order == "transcript":
        entries.sort(key=lambda x: (x.transcript_id, x.contig, x.start))
    elif sort_order == "position":
        entries.sort(key=lambda x: (x.contig, x.start))
    elif sort_order == "position+gene":
        entries.sort(key=lambda x: (x.gene_id, x.start))
        genes = list(flat_gene_iterator(entries))
        genes.sort(key=lambda x: (x[0].contig, x[0].start))
        entries = IOTools.flatten(genes)
    elif sort_order == "gene+exon":
        entries.sort(key=lambda x: (x.gene_id, x.exon_number))

    for entry in entries:
        yield entry


def iterator_overlapping_genes(gtf_iterator, min_overlap=0):
    """return overlapping genes."""

    genes = []
    for gene in flat_gene_iterator(gtf_iterator):
        gene.sort(key=lambda x: x.start)
        genes.append(gene[0].contig, gene[0].start, gene)

    genes.sort()

    last_contig, last_end, ovl = None, 0, []

    for this_contig, this_start, gene in genes:
        if this_contig != last_contig or last_end < this_start:
            if last_contig:
                yield ovl
            ovl, last_contig = [], this_contig
            last_end = 0

        last_end = max(last_end, gene[-1].end)

    if last_contig:
        yield ovl


def iterator_transcripts2genes(gtf_iterator, min_overlap=0):
    """cluster transcripts by exon overlap.

    The gene id is set to the first transcript encountered of a gene.
    If a gene stretches over several contigs, subsequent copies are
    appended a number.
    """

    map_transcript2gene = {}
    gene_ids = collections.defaultdict(list)

    for chunk in iterator_overlaps(gtf_iterator):
        transcript_ids = list(set([x.transcript_id for x in chunk]))
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
            gene_ids[gene_id].append(contig)
            index = 0
        else:
            try:
                index = gene_ids[gene_id].index(contig)
            except ValueError:
                index = len(gene_ids[gene_id])
                gene_ids[gene_id].append(contig)

        for x in transcript_ids:
            map_transcript2gene[x] = gene_id

        if index:
            gene_id += ".%i" % index
        for x in chunk:
            x.gene_id = gene_id
        yield chunk


def iterator_overlaps(gff_iterator, min_overlap=0):
    """iterate over gff file and return a list of features that
    are overlapping.

    The input should be sorted by contig,start
    """

    last = next(gff_iterator)
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

        assert last.start <= this.start, \
            "input file needs to be sorted by contig, start:\n%s\n%s\n" % \
            (str(last), str(this))
        matches.append(this)
        last = this
        end = max(end, this.end)

    yield matches


def Overlap(entry1, entry2, min_overlap=0):
    """returns true, if entry1 and entry2 overlap by a minimum number of
    residues.

    """

    return (entry1.contig == entry2.contig and
            entry1.strand == entry2.strand and
            min(entry1.end, entry2.end) -
            max(entry1.start, entry2.start) > min_overlap)


def Identity(entry1, entry2, max_slippage=0):
    """returns true, if entry1 and entry2 are (almost) identical, allowing
    a small amount of slippage at either end.
    """

    return (entry1.contig == entry2.contig and
            entry1.strand == entry2.strand and
            abs(entry1.end - entry2.end) < max_slippage and
            abs(entry1.start - entry2.start < max_slippage))


def HalfIdentity(entry1, entry2, max_slippage=0):
    """returns true, if entry1 and entry2 overlap and at least one end is
    within `max_slippage` residues.
    """

    return (entry1.contig == entry2.contig and
            entry1.strand == entry2.strand and
            (abs(entry1.end - entry2.end) < max_slippage or
             abs(entry1.start - entry2.start < max_slippage)))


def asRanges(gffs, feature=None):
    """return ranges within a set of gffs.

    Overlapping intervals are merged.

    The returned intervals are sorted.
    """

    if isinstance(feature, str):
        gg = [x for x in gffs if x.feature == feature]
    elif feature:
        gg = [x for x in gffs if x.feature in feature]
    else:
        gg = gffs[:]

    r = [(g.start, g.end) for g in gg]
    return Intervals.combine(r)


def CombineOverlaps(old_gff, method="combine"):
    """combine overlapping entries for a list of gffs.

    method can be any of combine|longest|shortest
    only the first letter is important.
    """

    old_gff.sort(lambda x, y: cmp((x.contig, x.strand, x.start, x.end),
                                  (y.contig, y.strand, y.start, y.end)))

    new_gff = []

    last_e = old_gff[0]

    for e in old_gff[1:]:
        if not Overlap(last_e, e):
            new_gff.append(last_e)
            last_e = e
        else:
            if method[0] == "c":
                last_e.start = min(last_e.start, e.start)
                last_e.end = max(last_e.end, e.end)
                last_e.attributes += " ; " + e.attributes

    new_gff.append(last_e)

    return new_gff


def SortPerContig(gff):
    """sort gff entries per contig and return a dictionary mapping a
    contig to the begin of the list."""
    map_contig2start = {}

    if len(gff) == 0:
        return map_contig2start

    gff.sort(lambda x, y: cmp(x.contig, y.contig))

    last_contig = None
    start = 0
    for x in range(len(gff)):
        if last_contig != gff[x].contig:
            map_contig2start[last_contig] = (start, x)
            start = x
            last_contig = gff[x].contig

    map_contig2start[last_contig] = (start, x)

    return map_contig2start


def toIntronIntervals(chunk):
    '''convert a set of gtf elements within a transcript to intron coordinates.

    Will use first transcript_id found.

    Note that coordinates will still be forward strand coordinates
    '''
    if len(chunk) == 0:
        return []
    contig, strand, transcript_id = (chunk[0].contig,
                                     chunk[0].strand,
                                     chunk[0].transcript_id)
    for gff in chunk:
        assert gff.strand == strand, "features on different strands."
        assert gff.contig == contig, "features on different contigs."

    intervals = Intervals.combine([(x.start, x.end)
                                   for x in chunk if x.feature == "exon"])
    return Intervals.complement(intervals)


def toSequence(chunk, fasta):
    """convert a list of gff attributes to a single sequence.

    This function ensures correct in-order concatenation on
    positive/negative strand. Overlapping regions are merged.
    """
    if len(chunk) == 0:
        return ""

    contig, strand = chunk[0].contig, chunk[0].strand

    for gff in chunk:
        assert gff.strand == strand, "features on different strands."
        assert gff.contig == contig, "features on different contigs."

    intervals = Intervals.combine([(x.start, x.end) for x in chunk])
    lcontig = fasta.getLength(contig)
    positive = Genomics.IsPositiveStrand(strand)

    if not positive:
        intervals = [(lcontig - end, lcontig - start)
                     for start, end in intervals]
        intervals.reverse()

    s = [fasta.getSequence(contig, strand, start, end)
         for start, end in intervals]

    return "".join(s)


def readFromFile(infile):
    """read records from file and return as list."""
    result = []
    for gff in pysam.tabix_iterator(infile, pysam.asGTF()):
        result.append(gff)
    return result


def readAsIntervals(gff_iterator,
                    with_values=False,
                    with_records=False,
                    merge_genes=False,
                    with_gene_id=False,
                    with_transcript_id=False,
                    use_strand=False):
    """read tuples of (start, end) from a GTF file.

    This method ignores everything else but the coordinates.

    The `with_values`, `with_gene_id` and `with_records` options are
    exclusive.

    Arguments
    ----------
    gff_iterator : iterator
       Iterator yielding GTF records.
    with_values
       If True, the content of the score field is added to the tuples.
    with_records
       If True, the entire record is added to the tuples.
    merge_genes
       If true, the GTF records are passed through the :func:
       `merged_gene_iterator` iterator first.
    with_gene_id
       If True, the gene_id is added to the tuples.
    with_transcript_id
       If True, the transcript_ids are added to the tuples.
    use_strand
       If true, intervals will be grouped by contig and strand.
       The default is to group by contig only.

    Returns a dictionary of intervals by contig.

    """

    assert not (
        with_values and with_records),\
        "both with_values and with_records are true."
    intervals = collections.defaultdict(list)

    if merge_genes:
        it = merged_gene_iterator(gff_iterator)
    else:
        it = gff_iterator

    if use_strand:
        keyf = lambda x: (x.contig, x.strand)
    else:
        keyf = lambda x: x.contig

    if with_values:
        for gff in it:
            intervals[keyf(gff)].append((gff.start, gff.end, gff.score))
    elif with_records:
        for gff in it:
            intervals[keyf(gff)].append((gff.start, gff.end, gff))
    elif with_gene_id:
        for gff in it:
            intervals[keyf(gff)].append((gff.start, gff.end, gff.gene_id))
    elif with_transcript_id:
        for gff in it:
            intervals[keyf(gff)].append(
                (gff.start, gff.end, gff.transcript_id))
    else:
        for gff in gff_iterator:
            intervals[keyf(gff)].append((gff.start, gff.end))

    return intervals


def readAndIndex(iterator, with_value=True):
    '''read from gtf stream and index.

    Returns
    -------

    index :
        an object of type :class:`IndexedGenome.IndexedGenome`
    '''

    if with_value:
        index = IndexedGenome.IndexedGenome()
        for gtf in iterator:
            index.add(gtf.contig, gtf.start, gtf.end, gtf)
    else:
        index = IndexedGenome.Simple()
        for gtf in iterator:
            index.add(gtf.contig, gtf.start, gtf.end)

    return index


class Error(Exception):
    """Base class for exceptions in this module."""

    def __str__(self):
        return str(self.message)

    def _get_message(self, message):
        return self._message

    def _set_message(self, message):
        self._message = message
    message = property(_get_message, _set_message)


class ParsingError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


def toDot(v):
    '''convert value to '.' if None'''
    if v is None:
        return "."
    else:
        return str(v)


def quote(v):
    '''return a quoted attribute.'''
    if type(v) in (str,):
        return '"%s"' % v
    else:
        return str(v)


class Entry:
    """representation of a :term:`GTF` formatted entry.

    Attributes
    ----------
    contig : string
       Chromosome/contig
    source : string
       The GTF source field
    feature : string
       The GTF feature field
    frame : string
       The frame
    start : int
       Start coordinate in 0-based coordinates, half-open coordinates
    end : int
       End coordinate in 0-based coordinates, half-open coordinates
    score : float
       Score associated with feature
    strand : string
       Strand of feature
    gene_id : string
       Gene identifier of feature. Not present for :term:`GFF` formatted
       data.
    transcript_id : string
       Transcript identifier of feature. Not present for :term:`GFF` formatted
       data.
    attributes : dict
       Dictionary of additional attributes in the GFF/GTF record (last column)
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
        self.transcript_id = None
        self.attributes = collections.OrderedDict()

    def read(self, line):
        """read gff entry from line in GTF/GFF format.

        <seqname> <source> <feature> <start> <end> <score> \
              <strand> <frame> [attributes] [comments]
        """

        data = line[:-1].split("\t")

        try:
            (self.contig, self.source, self.feature,
             self.start, self.end, self.score, self.strand,
             self.frame) = data[:8]
        except ValueError:
            raise ValueError("parsing error in line `%s`" % line)

        # note: frame might be .
        (self.start, self.end) = list(map(int, (self.start, self.end)))
        self.start -= 1

        self.parseInfo(data[8], line)

    def parseInfo(self, attributes, line):
        """parse attributes.

        This method will set the gene_id and transcript_id attributes
        if present.
        """

        # remove comments
        attributes = attributes.split("#")[0]

        # separate into fields
        # Fields might contain a ";", for example in ENSEMBL GTF file
        # for mouse, v78:
        # ...; transcript_name "TXNRD2;-001"; ....
        # The current heuristic is to split on a semicolon followed by a
        # space, which seems to be part of the specification, see
        # http://mblab.wustl.edu/GTF22.html
        fields = [x.strip() for x in attributes[:attributes.rfind(";")].split("; ")]
        self.attributes = collections.OrderedDict()

        for f in fields:

            d = [x.strip() for x in f.split(" ")]

            n, v = d[0], " ".join(d[1:])
            if len(d) > 2:
                v = d[1:]

            if v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            else:
                # try to convert to a value
                try:
                    v = float(v)
                    v = int(v)
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
            raise ParsingError("missing attribute 'gene_id' in line %s" % line)
        if not self.transcript_id:
            raise ParsingError(
                "missing attribute 'transcript_id' in line %s" % line)

    def getAttributeField(self, full=True):
        aa = []
        for k, v in list(self.attributes.items()):
            if isinstance(v, str):
                aa.append('%s "%s"' % (k, v))
            elif isinstance(v, list) or isinstance(v, tuple):
                aa.append('%s "%s"' % (k, " ".join(v)))
            else:
                aa.append('%s %s' % (k, str(v)))

        if full:
            return"; ".join(['gene_id "%s"' % self.gene_id,
                             'transcript_id "%s"' % self.transcript_id] +
                            aa)
        else:
            return "; ".join(aa)

        return self.attributes

    def __str__(self):
        def _todot(val):
            if val is None:
                return "."
            else:
                return str(val)

        attributes = self.getAttributeField()
        return "\t".join(map(str, (self.contig, self.source,
                                   self.feature,
                                   self.start + 1, self.end,
                                   _todot(self.score),
                                   self.strand,
                                   self.frame,
                                   attributes))) + ";"

    def invert(self, lcontig):
        """invert genomic coordinates from forward to reverse coordinates and
        back.

        Arguments
        ---------
        lcontig : int
           Length of the chromosome that the feature resides on.
        """

        if self.strand in ("-", "0", "-1"):
            x = min(self.start, self.end)
            y = max(self.start, self.end)

            self.start = lcontig - y
            self.end = lcontig - x

    def fromGTF(self, other, gene_id=None, transcript_id=None):
        """fill record from other GFF/GTF entry.

        The optional attributes are not copied.
        """
        self.contig = other.contig
        self.source = other.source
        self.feature = other.feature
        self.start = other.start
        self.end = other.end
        self.score = other.score
        self.strand = other.strand
        self.frame = other.frame
        if gene_id is not None:
            self.gene_id = gene_id
        else:
            self.gene_id = other.gene_id
        if transcript_id is not None:
            self.transcript_id = transcript_id
        else:
            self.transcript_id = other.transcript_id

        # capture None values in frame and strand and
        # revert them to correct . format

        if self.frame is None:
            self.frame = "."
        if self.strand is None:
            self.strand = "."

        return self

    def fromBed(self, other, **kwargs):
        """fill record from a bed entry."""
        self.contig = other.contig
        self.source = kwargs.get("source", "bed")
        self.feature = kwargs.get("feature", "interval")
        self.start = other.start
        self.end = other.end
        self.score = other.score
        self.strand = other.strand
        self.frame = kwargs.get("frame", ".")
        self.gene_id = kwargs.get("gene_id", None)
        self.transcript_id = kwargs.get("transcript_id", None)
        return self

    def copy(self, other):
        """fill from other entry.

        This method works if other is :class:`GTF.Entry` or
        :class:`pysam.GTFProxy`.
        """
        self.contig = other.contig
        self.source = other.source
        self.feature = other.feature
        self.start = other.start
        self.end = other.end
        self.score = other.score
        self.strand = other.strand
        self.frame = other.frame

        # capture if frame and strand are None then
        # set to appropriate . format
        if self.frame is None:
            self.frame = "."
        if self.strand is None:
            self.strand = "."

        # gene_id and transcript_id can be optional
        try:
            self.gene_id = other.gene_id
        except AttributeError:
            pass
        try:
            self.transcript_id = other.transcript_id
        except AttributeError:
            pass

        self.attributes = collections.OrderedDict(other.asDict().items())
        # from gff - remove gene_id and transcript_id from attributes

        try:
            del self.attributes["gene_id"]
            del self.attributes["transcript_id"]
        except KeyError:
            pass

        return self

    def asDict(self):
        '''return attributes as a dictionary.'''
        return self.attributes

    def clearAttributes(self):
        self.attributes = collections.OrderedDict()

    def addAttribute(self, key, value=None):
        self.attributes[key] = value

    def __cmp__(self, other):
        # note: does compare by strand as well!
        return cmp((self.contig, self.strand, self.start),
                   (other.contig, other.strand, other.start))

    # python 3 compatibility
    def __lt__(self, other):
        # note: does compare by strand as well!
        return (self.contig, self.strand, self.start) < \
            (other.contig, other.strand, other.start)

    def __setitem__(self, key, value):
        self.addAttribute(key, value)

    def __getitem__(self, key):
        return self.attributes[key]

    def __contains__(self, key):
        return key in self.attributes

    def hasOverlap(self, other, min_overlap=0):
        """returns true, if overlap with other entry.
        """
        return (self.contig == other.contig and
                self.strand == other.strand and
                min(self.end, other.end) -
                max(self.start, other.start) > min_overlap)

    def isIdentical(self, other, max_slippage=0):
        """returns true, if self and other overlap completely.
        """
        return (self.contig == other.contig and
                self.strand == other.strand and
                abs(self.end - other.end) < max_slippage and
                abs(self.start - other.start < max_slippage))

    def isHalfIdentical(self, other, max_slippage=0):
        """returns true, if self and other overlap.
        """
        return (self.contig == other.contig and
                self.strand == other.strand and
                (abs(self.end - other.end) < max_slippage or
                 abs(self.start - other.start < max_slippage)))

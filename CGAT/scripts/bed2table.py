'''Bed2table.py - annotate intervals
=================================

:Tags: Genomics Intervals Summary

Purpose
-------

This script takes a bed-formatted file as input and annotates each interval.
Possible annotators are (see option '--counter'):

overlap

    compute overlap with intervals in other bed file. If the other bed
    file contains tracks, the overlap is computed per track.

peaks

    compute peak location in intervals. Requires one or more
    bam-files. This counter can also count within an secondary set of
    bam-files (--control-bam-file) and add this to the output.

composition-na

    compute nucleotide frequencies in intervals.

composition-cpg

    compute CpG densities and CpG observed / expected in intervals.

classifier-chipseq

   classify chipseq intervals. Requires a :term:`gff`
   file with genomic annotations (see :doc:`gtf2gff`.)

motif

   Search for a specified motif e.g. using --motif-sequence=TTTT.


Usage
-----

For example, the following command will compute the CpG composition
in the intervals supplied::

   cgat bed2table --counter=composition-cpg < in.bed > out.tsv

If the input is this :term:`bed` formatted file::

  chr19   60118   60120   1       100.0000
  chr19   60171   60173   2       100.0000
  chr19   60182   60184   3       100.0000
  chr19   60339   60341   4       100.0000
  chr19   60375   60377   5       100.0000
  chr19   60110   60118   noCpG   100.0000
  chr19   60118   60119   Conly   100.0000
  chr19   60119   60120   Gonly   100.0000

then the output is the following table:

+------+-----+-----+-----+--------+---------+-----------+----------+
|contig|start|end  |name |score   |CpG_count|CpG_density|CpG_ObsExp|
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60118|60120|1    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60171|60173|2    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60182|60184|3    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60339|60341|4    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60375|60377|5    |100.0000|1        |1.0        |2.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60110|60118|noCpG|100.0000|0        |0.0        |0.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60118|60119|Conly|100.0000|0        |0.0        |0.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+
|chr19 |60119|60120|Gonly|100.0000|0        |0.0        |0.0       |
+------+-----+-----+-----+--------+---------+-----------+----------+

Type::

   cgat bed2table.py --help

for command line help.

Command line options
---------------------

'''
import re
import sys
import collections
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGAT.IOTools as IOTools
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.SequenceProperties as SequenceProperties
import CGAT.Intervals as Intervals
import numpy
import pysam

import CGAT.GeneModelAnalysis as GeneModelAnalysis


class Counter(object):

    def __init__(self, fasta=None, *args, **kwargs):
        self.fasta = fasta

    def update(self, bed):
        self.bed = bed
        self.count(bed)

    def getSegments(self):
        return [self.bed]

    def getHeader(self):
        return '\t'.join(self.headers)


class CounterLength(Counter):
    headers = ['length']

    def count(self, bed):
        self.length = bed.end - bed.start

    def __str__(self):
        return str(self.length)


class CounterOverlap(Counter):

    '''count overlap for each interval in tracks.'''

    def __init__(self, filename, *args, **kwargs):

        assert filename is not None,\
            "please supply filename for CounterOverlap"

        Counter.__init__(self, *args, **kwargs)

        self.filename = filename

        E.info("reading intervals from %s" % self.filename)

        self.index = Bed.readAndIndex(
            IOTools.openFile(self.filename, "r"),
            per_track=True)

        E.info("read intervals for %s tracks" % len(self.index))

        self.tracks = list(self.index.keys())
        self.headers = []
        for track in self.tracks:
            self.headers.extend(["%s_nover" % track, "%s_bases" % track])

    def count(self, bed):
        '''update internal counts.'''

        results = []
        for track in self.tracks:
            try:
                overlaps = [(x[0], x[1])
                            for x in
                            self.index[track][bed.contig]
                            .find(bed.start, bed.end)]
            except KeyError:
                overlaps = []

            results.append((len(overlaps),
                            Intervals.calculateOverlap(
                                [(bed.start, bed.end), ],
                                Intervals.combine(overlaps))))

        self.data = results

    def __str__(self):
        '''output overlap of interval in *bed*'''

        r = []
        for track, result in zip(self.tracks, self.data):
            r.append("\t".join((str(result[0]), str(result[1]))))

        return "\t".join(r)

CounterPeaksResult = collections.namedtuple(
    "CounterPeaksResult", ("length nreads avgval peakval npeaks peakcenter"))


class CounterPeaks(Counter):

    '''compute number of extent of peaks in an interval.'''

    headers = None

    def __init__(self,
                 bamfiles,
                 offsets,
                 control_bamfiles,
                 control_offsets, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs)
        if not bamfiles:
            raise ValueError("supply --bam-file options for readcoverage")

        assert len(offsets) == 0 or len(bamfiles) == len(
            offsets), "number of bamfiles not the same as number of offsets"
        assert len(control_offsets) == 0 or \
            len(control_bamfiles) == len(control_offsets), \
            "number of control bamfiles not the same as number of offsets"

        self.bamfiles = bamfiles
        self.offsets = offsets
        self.control_bamfiles = control_bamfiles
        self.control_offsets = control_offsets

        self.headers = list(CounterPeaksResult._fields)
        if self.control_bamfiles:
            self.headers.extend(
                ["control_%s" % x for x in CounterPeaksResult._fields])

    def _count(self, bed, bamfiles, offsets):
        '''count reads in bed interval.'''

        contig, start, end = bed.contig, bed.start, bed.end

        length = end - start
        try:
            counts = numpy.zeros(length)
        except ValueError as msg:
            raise ValueError("Error negative length obtained: "
                             " message=%s contig=%s, start=%s, end=%s" %
                             (msg, contig, start, end))
        nreads = 0

        if offsets:
            # if offsets are given, shift tags.
            for samfile, offset in zip(bamfiles, offsets):

                shift = offset // 2
                # for peak counting I follow the MACS protocoll,
                # see the function def __tags_call_peak in PeakDetect.py
                # In words
                # Only take the start of reads (taking into account the strand)
                # add d/2=offset to each side of peak and start accumulate
                # counts.
                # for counting, extend reads by offset
                # on + strand shift tags upstream
                # i.e. look at the downstream window
                xstart, xend = max(0, start - shift), max(0, end + shift)

                for read in samfile.fetch(contig, xstart, xend):
                    nreads += 1
                    # some reads are assigned to a contig and position, but
                    # are flagged as unmapped - these might not have an alen
                    # attribute.
                    if read.is_unmapped:
                        continue

                    if read.is_reverse:
                        # rstart = read.pos + read.alen - offset
                        # offset = 2 * shift
                        try:
                            rstart = read.pos + read.alen - offset
                        except TypeError as msg:
                            raise TypeError("Error message =", msg,
                                            "read.pos =", read.pos,
                                            "read.alen =", read.alen,
                                            "offset =", offset,
                                            "query name =", read.qname,
                                            "length of read =", read.rlen)
                    else:
                        rstart = read.pos + shift

                    rend = rstart + shift
                    rstart = max(0, rstart - start)
                    rend = min(length, rend - start)
                    counts[rstart:rend] += 1

        else:
            for samfile in bamfiles:
                for read in samfile.fetch(contig, start, end):
                    nreads += 1
                    rstart = max(0, read.pos - start)
                    rend = min(length, read.pos - start + read.rlen)
                    counts[rstart:rend] += 1

        length = end - start
        avgval = numpy.mean(counts)
        peakval = max(counts)

        # set other peak parameters
        peaks = numpy.array(list(range(0, length)))[counts >= peakval]
        npeaks = len(peaks)
        # peakcenter is median coordinate between peaks
        # such that it is a valid peak in the middle
        peakcenter = start + peaks[npeaks // 2]

        return CounterPeaksResult(length, nreads, avgval,
                                  peakval, npeaks, peakcenter)

    def count(self, bed):
        '''count reads per position.

        If offsets are given, shift tags by offset / 2 and extend
        by offset / 2.
        '''

        self.result = self._count(bed, self.bamfiles, self.offsets)
        if self.control_bamfiles:
            self.control = self._count(
                bed, self.control_bamfiles, self.control_offsets)

    def __str__(self):
        if self.control_bamfiles:
            return "\t".join(map(str, self.result + self.control))
        else:
            return "\t".join(map(str, self.result))


class CounterCompositionNucleotides(Counter):

    headers = SequenceProperties.SequencePropertiesNA().getHeaders()

    def __init__(self, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs)
        self.result_class = SequenceProperties.SequencePropertiesNA
        assert self.fasta, "Counter requires a genomic sequence"

    def count(self, bed):
        s = self.fasta.getSequence(bed.contig, "+", bed.start, bed.end)
        self.result = self.result_class()
        self.result.loadSequence(s)

    def __str__(self):
        return str(self.result)


class CounterMotif(Counter):

    headers = ['motif_counts']

    def __init__(self, motif, *args, **kwargs):
        Counter.__init__(self, *args, **kwargs)
        self.motif = motif

    def count(self, bed):
        s = self.fasta.getSequence(bed.contig, "+", bed.start, bed.end)
        self.result = len([x for x in
                           re.finditer(r'(?=(%s))' % self.motif, s)])

    def __str__(self):
        return str(self.result)


class CounterCompositionCpG(CounterCompositionNucleotides):

    '''compute CpG frequencies as well as nucleotide frequencies.

    Note that CpG density is calculated across the merged exons of a
    transcript. Thus, there might be difference between the CpG on a
    genomic level and on the transrcipt level depending on how many
    genomic CpG are lost across an intron-exon boundary or how many
    transcript CpG are created by exon fusion.

    '''

    headers = SequenceProperties.SequencePropertiesCpg().getHeaders()

    def __init__(self, *args, **kwargs):
        CounterCompositionNucleotides.__init__(self, *args, **kwargs)
        self.result_class = SequenceProperties.SequencePropertiesCpg

    def count(self, bed):

        s = self.fasta.getSequence(bed.contig, "+", bed.start, bed.end)
        self.result = self.result_class()
        self.result.loadSequence(s)


class ClassifierChIPSeq(GeneModelAnalysis.Classifier):

    """classify ChIPSeq intervals based on a reference annotation.

    This assumes the input is a genome annotation derived from an
    ENSEMBL gtf file created with gff2gtf.py.

    In contrast to transcripts, the intervals are fuzzy. Hence the
    classification is based on a mixture of full/partial overlap.

    An interval is classified as:

    cds
       mostly part of a CDS. These would be intervals fully within a CDS exon.
    utr
       mostly part of UTR. These are intervals fully within the UTR of a gene.
    intergenic
       mostly intergenic. These are intervals fully within the
       intergenic region and more than 1kb from the closest exon.
    upstream
       not any of the above and partly upstream of a gene. These are
       intervals that might overlap part of the UTR or the 1kb segment
       before to the 5'-terminal exon of a gene.
    downstream
       not any of the abore and partly downstream of a gene. These are
       intervals that might overlap part of the UTR or the 1kb segment
       after to the 3'-terminal exon of a gene.
    intronic
       not any of the above and partly intronic. Note that these could
       also include promotors of short alternative transcripts that
       skip one or more of the first exons.
    ambiguous
       none of the above

    """

    header = ["is_cds", "is_utr", "is_upstream", "is_downstream",
              "is_intronic", "is_intergenic", "is_flank", "is_ambiguous"]

    # sources to use for classification
    sources = ("", )

    # minimum coverage of a transcript to assign it to a class
    mThresholdMinCoverage = 95

    # full coverage of a transcript to assign it to a class
    mThresholdFullCoverage = 99

    # some coverage of a transcript to assign it to a class
    mThresholdSomeCoverage = 10

    def update(self, bed):

        # convert to a gtf entry
        gtf = GTF.Entry()

        gtf.fromBed(bed)
        gtf.feature = 'exon'
        GeneModelAnalysis.Classifier.update(self, [gtf])

    def count(self):

        for key in self.mKeys:
            self.mCounters[key].update(self.mGFFs)

        def s_min(*args):
            return sum([abs(self.mCounters[x].mPOverlap1)
                        for x in args]) >= self.mThresholdMinCoverage

        def s_excl(*args):
            return sum([abs(self.mCounters[x].mPOverlap1)
                        for x in args]) < (100 - self.mThresholdMinCoverage)

        def s_full(*args):
            return sum([abs(self.mCounters[x].mPOverlap1)
                        for x in args]) >= self.mThresholdFullCoverage

        def s_some(*args):
            return sum([abs(self.mCounters[x].mPOverlap1)
                        for x in args]) >= self.mThresholdSomeCoverage

        self.mIsCDS, self.mIsUTR, self.mIsIntergenic = False, False, False
        self.mIsUpStream = False
        self.mIsDownStream = False
        self.mIsIntronic = False
        self.mIsFlank, self.mIsAmbiguous = False, False

        self.mIsCDS = s_full(":CDS")
        self.mIsUTR = s_full(":UTR", ":UTR3", ":UTR5")
        self.mIsIntergenic = s_full(":intergenic", ":telomeric")

        if not(self.mIsCDS or self.mIsUTR or self.mIsIntergenic):
            self.mIsUpStream = s_some(":5flank", ":UTR5")
            if not self.mIsUpStream:
                self.mIsDownStream = s_some(":3flank", ":UTR3")
                if not self.mIsDownStream:
                    self.mIsIntronic = s_some(":intronic")
                    if not self.mIsIntronic:
                        self.mIsFlank = s_some(":flank")

        self.mIsAmbiguous = not(self.mIsUTR or
                                self.mIsIntergenic or
                                self.mIsIntronic or
                                self.mIsCDS or
                                self.mIsUpStream or
                                self.mIsDownStream or
                                self.mIsFlank)

    def __str__(self):

        def to(v):
            if v:
                return "1"
            else:
                return "0"

        h = [to(x) for x in (self.mIsCDS,
                             self.mIsUTR,
                             self.mIsUpStream,
                             self.mIsDownStream,
                             self.mIsIntronic,
                             self.mIsIntergenic,
                             self.mIsFlank,
                             self.mIsAmbiguous,
                             )]

        for key in self.mKeys:
            h.append(str(self.mCounters[key]))
        return "\t".join(h)


def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome [default=%default].")

    parser.add_option(
        "-b", "--bam-file", dest="bam_files", type="string",
        help="filename with read mapping information. Multiple files can be "
        "submitted in a comma-separated list [default=%default].")

    parser.add_option(
        "--control-bam-file", dest="control_bam_files", type="string",
        help="filename with read mapping information for input/control. "
        "Multiple files can be submitted in a comma-separated list "
        "[default=%default].")

    parser.add_option(
        "--filename-format", dest="filename_format", type="choice",
        choices=("bed", "gff", "gtf"),
        help="format of secondary stream [default=%default].")

    parser.add_option(
        "-c", "--counter", dest="counters", type="choice", action="append",
        choices=("length",
                 "overlap",
                 "peaks",
                 "composition-na",
                 "composition-cpg",
                 "classifier-chipseq",
                 "motif"),
        help="select counters to apply [default=%default].")

    parser.add_option(
        "--motif-sequence", dest="motif_sequence", type="string",
        help="specify a sequence to search for"
        "[default=%default].")

    parser.add_option(
        "-o", "--offset", dest="offsets", type="int", action="append",
        help="tag offsets for tag counting - supply as many as there "
        "are bam-files [default=%default].")

    parser.add_option(
        "--control-offset", dest="control_offsets", type="int",
        action="append",
        help="control tag offsets for tag counting - supply as many as "
        "there are bam-files [default=%default].")

    parser.add_option(
        "-a", "--output-all-fields", dest="all_fields", action="store_true",
        help="output all fields in original bed file, by default only "
        "the first 4 are output [default=%default].")

    parser.add_option(
        "--output-bed-headers", dest="bed_headers", type="string",
        help="supply ',' separated list of headers for bed component "
        "[default=%default].")

    parser.add_option(
        "-f", "--gff-file", dest="filename_gff", type="string",
        action="append", metavar='bed',
        help="filename with extra gff files. The order is important "
        "[default=%default].")

    parser.add_option(
        "--has-header", dest="has_header", action="store_true",
        help="bed file with headers. Headers and first columns are "
        "preserved [default=%default]")

    parser.set_defaults(
        genome_file=None,
        counters=[],
        bam_files=None,
        offsets=[],
        control_bam_files=None,
        control_offsets=[],
        all_fields=False,
        filename_format=None,
        bed_headers=None,
        filename_gff=[],
        has_header=False,
        motif_sequence=None
    )

    (options, args) = E.Start(parser)

    if options.bed_headers is not None:
        bed_headers = [x.strip() for x in options.bed_headers.split(",")]
        if len(bed_headers) < 3:
            raise ValueError("a bed file needs at least three columns")
    else:
        bed_headers = None

    if options.has_header:
        while 1:
            line = options.stdin.readline()
            if not line:
                E.warn("empty bed file with no header")
                E.Stop()
                return
            if not line.startswith("#"):
                break
        bed_headers = line[:-1].split("\t")

    if "motif" in options.counters and not options.motif_sequence:
        raise ValueError("if using motif must specify a motif-sequence")

    # get files
    if options.genome_file:
        fasta = IndexedFasta.IndexedFasta(options.genome_file)
    else:
        fasta = None

    if options.bam_files:
        bam_files = []
        for bamfile in options.bam_files.split(","):
            bam_files.append(pysam.AlignmentFile(bamfile, "rb"))
    else:
        bam_files = None

    if options.control_bam_files:
        control_bam_files = []
        for bamfile in options.control_bam_files.split(","):
            control_bam_files.append(pysam.AlignmentFile(bamfile, "rb"))
    else:
        control_bam_files = None

    counters = []

    for c in options.counters:
        if c == "length":
            counters.append(CounterLength(fasta=fasta,
                                          options=options))

        elif c == "overlap":
            counters.append(CounterOverlap(filename=options.filename_gff[0],
                                           fasta=fasta,
                                           options=options))
            del options.filename_gff[0]
        elif c == "peaks":
            counters.append(CounterPeaks(bam_files,
                                         options.offsets,
                                         control_bam_files,
                                         options.control_offsets,
                                         options=options))
        elif c == "composition-na":
            counters.append(CounterCompositionNucleotides(fasta=fasta,
                                                          options=options))
        elif c == "composition-cpg":
            counters.append(CounterCompositionCpG(fasta=fasta,
                                                  options=options))
        elif c == "classifier-chipseq":
            counters.append(ClassifierChIPSeq(
                filename_gff=options.filename_gff,
                fasta=fasta,
                options=options,
                prefix=None))
            del options.filename_gff[0]

        elif c == "motif":
            counters.append(CounterMotif(fasta=fasta,
                                         motif=options.motif_sequence))

    extra_fields = None

    for bed in Bed.iterator(options.stdin):

        if extra_fields is None:

            # output explicitely given headers
            if bed_headers:
                if len(bed_headers) > bed.columns:
                    raise ValueError(
                        "insufficient columns (%i, expected %i) in %s" %
                        (bed.columns, len(bed_headers), str(bed)))

            else:
                bed_headers = Bed.Headers[:bed.columns]

            options.stdout.write("\t".join(bed_headers))
            options.stdout.write("\t" + "\t".join(
                [x.getHeader() for x in counters]) + "\n")

            extra_fields = list(range(len(bed_headers) - 3))

        for counter in counters:
            counter.update(bed)

        if options.all_fields:
            options.stdout.write(str(bed))
        else:
            options.stdout.write(
                "\t".join([bed.contig,
                           str(bed.start),
                           str(bed.end)] + [bed.fields[x]
                                            for x in extra_fields]))
        for counter in counters:
            options.stdout.write("\t%s" % str(counter))

        options.stdout.write("\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

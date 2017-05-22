'''
Tophat.py - working with tophat/cufflinks output files
======================================================

:Tags: Python

Code
----

'''

import re
import collections
import itertools

from CGAT import IOTools as IOTools


class CuffCompareValues:

    def __init__(self, vals):

        assert len(vals) == 4
        try:
            self.sn = float(vals[0])
        except ValueError:
            self.sn = None

        try:
            self.sp = float(vals[1])
        except ValueError:
            self.sp = None

        try:
            self.fsn = float(vals[2])
        except ValueError:
            self.fsn = None

        try:
            self.fsp = float(vals[3])
        except ValueError:
            self.fsp = None

    def __str__(self):
        return "\t".join(map(IOTools.val2str,
                             (self.sn,
                              self.sp,
                              self.fsn,
                              self.fsp)))

    @classmethod
    def getHeaders(cls):
        return ("sn", "sp", "fsn", "fsp")


class CuffCompareResult:

    def __init__(self):
        (self.baselevel,
         self.exonlevel,
         self.intronlevel,
         self.intronchainlevel,
         self.transcriptlevel,
         self.locuslevel,
         self.missedexons_counts,
         self.missedexons_total,
         self.novelexons_counts,
         self.novelexons_total,
         self.missedintrons_counts,
         self.missedintrons_total,
         self.novelintrons_counts,
         self.novelintrons_total,
         self.missedloci_counts,
         self.missedloci_total,
         self.novelloci_counts,
         self.novelloci_total,
         self.matchingintronchains_count,
         self.matchingloci_count,
         self.query,
         self.query_loci,
         self.query_multi_exon,
         self.reference,
         self.reference_loci,
         self.reference_multi_exon,
         self.loci_multi_exon,
         self.loci_transcripts,
         ) = [None] * 28

    def fromLines(self, lines):
        '''parse from cuffcompare output.'''
        self.is_empty = True

        for line in lines:
            if line.startswith("#"):
                d = line[1:-1].strip()
                if d.startswith("Query"):
                    self.query, self.query_loci, self.query_multi_exon = list(map(int, re.match(
                        "Query mRNAs :\s+(\d+)\s+in\s+(\d+)\s+loci\s+\((\d+)", d).groups()))
                elif d.startswith("Reference"):
                    self.reference, self.reference_loci, self.reference_multi_exon = list(map(int, re.match(
                        "Reference mRNAs :\s+(\d+)\s+in\s+(\d+)\s+loci\s+\((\d+)", d).groups()))
                elif d.startswith("("):
                    self.loci_multi_exon, self.loci_transcripts = re.match(
                        "\((\d+) multi-transcript loci, ~(\S+)", d).groups()
                    self.loci_multi_exon = int(self.loci_multi_exon)
                    self.loci_transcripts = float(self.loci_transcripts)
                continue

            line = line[:-1].strip()
            if not line:
                continue

            try:
                tag, data = line.split(":")
            except ValueError:
                raise ValueError("parsing error in line %s" % line)

            tag = re.sub("\s", "", tag).lower()

            if tag.startswith("novel") or tag.startswith("missed"):

                counts, total = list(map(
                    int, re.match("\s+(\d+)/(\d+)", data).groups()))
                setattr(self, "%s_counts" % tag, counts)
                setattr(self, "%s_total" % tag, total)

            # NICK FIX: This output changed from the previous version of cufflinks. not sure if it means the same
            # but nevertheless the new output will be put into sphinxreport
            elif tag.startswith("matching"):
                counts = int(data)
                setattr(self, "%s_count" % tag, counts)

            elif tag.startswith("total"):
                # stop at total union across all super-loci
                break
            else:
                values = data.strip().split()
                setattr(self, tag, CuffCompareValues(values))
                self.is_empty = False

    def __str__(self):
        return "\t".join(map(str, (
            self.baselevel,
            self.exonlevel,
            self.intronlevel,
            self.intronchainlevel,
            self.transcriptlevel,
            self.locuslevel,
            self.missedexons_counts,
            self.missedexons_total,
            self.novelexons_counts,
            self.novelexons_total,
            self.missedintrons_counts,
            self.missedintrons_total,
            self.novelintrons_counts,
            self.novelintrons_total,
            self.missedloci_counts,
            self.missedloci_total,
            self.novelloci_counts,
            self.novelloci_total,
            self.matchingintronchains_count,
            self.matchingloci_count,
            self.query,
            self.query_loci,
            self.query_multi_exon,
            self.reference,
            self.reference_loci,
            self.reference_multi_exon,
            self.loci_multi_exon,
            self.loci_transcripts
        )))

    @classmethod
    def getHeaders(cls):
        m = ("baselevel", "exonlevel", "intronlevel", "intronchainlevel",
             "transcriptlevel", "locuslevel")
        a = ("missed", "novel")
        b = ("exons", "introns", "loci")
        c = ("counts", "total")
        return [ "%s_%s" % (x, y) for x, y in itertools.product( m, CuffCompareValues.getHeaders() ) ] +\
            [ "%s%s_%s" % (x, y, z) for x, y, z in itertools.product( a, b, c) ] +\
            ["matchingintronchains_count", "matchingloci_count"] +\
            ["query", "query_loci", "query_multi_exon",
             "reference", "reference_loci", "reference_multi_exon",
             "loci", "loci_multi_exon"]


def parseTranscriptComparison(infile):
    '''read cufflinks 1.0.3 output in infile stream.

    returns a two-level dictionary mapping with levels track and contig.
    '''

    firstline = infile.readline()
    result = collections.defaultdict(dict)
    tracks = []

    def __blocker(infile):
        # 'Summary for dataset'
        # comes after 'Genomic sequence'.

        blocks = {}
        contig, block = None, []
        # whether to yield
        y = False
        for line in infile:

            if line.startswith("#> Genomic sequence"):
                if block:
                    blocks[contig] = block
                if y:
                    yield dataset, blocks
                    blocks = {}
                    y = False

                contig = re.match(
                    "#> Genomic sequence: (\S+)", line).groups()[0]
                block = []
                continue
            elif line.startswith("#= Summary for dataset:"):
                if block:
                    blocks[contig] = block
                if y:
                    yield dataset, blocks
                    blocks = {}
                    y = False

                dataset = re.match(
                    "#= Summary for dataset: (\S+)", line).groups()[0]
                contig = "all"
                block = []

                y = True
                continue

            block.append(line)

        blocks[contig] = block
        yield dataset, blocks

    for track, blocks in __blocker(infile):
        print(track)
        tracks.append(track)
        for contig, block in blocks.items():
            r = CuffCompareResult()

            r.fromLines(block)
            result[track][contig] = r

    return tracks, result

Locus = collections.namedtuple(
    "Locus", "locus_id contig strand start end transcript_ids transcripts")
Tracking = collections.namedtuple(
    "Tracking", "transfrag_id locus_id ref_gene_id ref_transcript_id code transcripts")
TranscriptInfo = collections.namedtuple(
    "Transfrag", "gene_id transcript_id fmi fpkm conf_lo conf_hi cov len")


def iterate_tracking(infile):
    '''parse .tracking output file from cuffcompare

    returns iterator with list of loci.
    '''
    for line in infile:
        data = line[:-1].split("\t")
        try:
            transfrag_id, locus_id, transcript_id, code = [
                x.strip() for x in data[:4]]
        except ValueError:
            raise ValueError(
                "parsing error: expected transfrag_id, locus_id, transcript_id, code, got %s" % str(data))

        if transcript_id == "-":
            ref_gene_id, ref_transcript_id = "", ""
        else:
            ref_gene_id, ref_transcript_id = transcript_id.split("|")

        transcripts = []
        for c in data[4:]:
            if c == "-":
                transcripts.append(None)
            else:
                # there can be multiple transcripts
                for cc in c.split(","):
                    try:
                        (gene_id, transcript_id, fmi, fpkm, conf_lo,
                         conf_hi, cov, length) = cc.split("|")
                    except ValueError:
                        raise ValueError("parsing error for field '%s'" % cc)
                    (fpkm, conf_lo, conf_hi, cov) = list(map(
                        float, (fpkm, conf_lo, conf_hi, cov)))
                    if length == "-":
                        length = 0
                    (fmi, length) = list(map(int, (fmi, length)))

                    transcripts.append(TranscriptInfo._make(
                        (gene_id, transcript_id, fmi, fpkm, conf_lo, conf_hi, cov, length)))

        yield Tracking._make((transfrag_id, locus_id, ref_gene_id, ref_transcript_id, code, transcripts))


def iterate_locus(infile):
    '''parse .loci output file from cuffcompare

    returns iterator with list of loci.
    '''

    for line in infile:
        data = line[:-1].split("\t")
        locus_id, pos = data[:2]

        contig, strand, start, end = re.match(
            "(\S+)\[(\S*)\](\d+)-(\d+)", pos).groups()
        start, end = list(map(int, (start, end)))

        # some [] contains the 0 byte, convert to empty field
        if strand not in "+-":
            strand = ""
        transcripts = []
        for c in data[2:]:
            if c == "-":
                transcripts.append(())
            else:
                transcripts.append(c.split(","))

        yield Locus._make((locus_id, contig, strand, start, end, transcripts[1], transcripts[1:]))

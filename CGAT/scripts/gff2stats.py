'''gff2stats.py - count features, etc. in gff file
===============================================

:Tags: Genomics Intervals GFF GTF Summary

Purpose
-------

This script generates summary statistics over features,
source, gene_id and transcript_id in one or more :term:`gff`
or :term:`gtf` formatted files.

Usage
-----

Input is either a gff or gtf file; gtf input must be specified
with the --is-gtf option.

Example::

   python gff2stats.py --is-gtf example.gtf > example_sum.tsv

   cat example.gtf

   19  processed_transcript  exon  6634666509  .  -  .  gene_id "ENSG00000225373"; transcript_id "ENST00000592209" ...
   19  processed_transcript  exon  6052160747  .  -  .  gene_id "ENSG00000225373"; transcript_id "ENST00000592209" ...
   19  processed_transcript  exon  6010560162  .  -  .  gene_id "ENSG00000225373"; transcript_id "ENST00000592209" ...
   19  processed_transcript  exon  6634666416  .  -  .  gene_id "ENSG00000225373"; transcript_id "ENST00000589741" ...

   cat example_sum.tsv

   track  contigs  strands  features  sources  genes  transcripts ...
   stdin  1        2        4         23       2924   12752       ...


The counter used is dependent on the file type.  For a gff file, the implemented counters are:

1. number of intervals per contig, strand, feature and source

For a gtf file, the additional implemented counters are:

1. number of genes, transcripts, single exon transcripts
2. summary statistics for exon numbers, exon sizes, intron sizes and
   transcript sizes

The output is a tab-separated table.

Options
-------

The default action of ``gff2stats`` is to count over contigs, strand,
feature and source.  This assumes the input file is a gff file

There is a single option for this script::

``--is-gtf``
   The input file is gtf format.  The output will therefore
   contain summaries over exon numbers, exon sizes, intron sizes and
   transcript sizes in addition to the the number of genes,
   transcripts and single exon transcripts.

Type::

   python gff2stats.py --help

for command line help.

Command line options
--------------------

'''
import sys
import collections
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Stats as Stats
import CGAT.IOTools as IOTools
import CGAT.Intervals as Intervals


class counter_gff:

    fields = ("contigs", "strands", "features", "sources")

    def __init__(self, iter):
        self.iter = iter

        self.counts_contigs = collections.defaultdict(int)
        self.counts_strands = collections.defaultdict(int)
        self.counts_features = collections.defaultdict(int)
        self.counts_sources = collections.defaultdict(int)

    def __next__(self):

        entry = next(self.iter)

        self.counts_contigs[entry.contig] += 1
        self.counts_features[entry.feature] += 1
        self.counts_sources[entry.source] += 1
        self.counts_strands[entry.strand] += 1

        return entry

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self

    def __str__(self):
        return "\t".join(map(str, (len(self.counts_contigs),
                                   len(self.counts_strands),
                                   len(self.counts_features),
                                   len(self.counts_sources))))


class counter_exons:

    fields = ("genes", "transcripts", "single_exon_transcripts",) +\
        tuple(["exon_count_%s" % x for x in Stats.Summary.fields]) +\
        tuple(["exon_size_%s" % x for x in Stats.Summary.fields]) +\
        tuple(["intron_size_%s" % x for x in Stats.Summary.fields]) +\
        tuple(["transcript_size_%s" % x for x in Stats.Summary.fields])

    def __init__(self, iter):

        self.iter = iter

        self.counts_gene_ids = collections.defaultdict(int)
        self.counts_transcript_ids = collections.defaultdict(int)
        self.counts_exons_per_transcript = collections.defaultdict(list)

    def __next__(self):

        while 1:
            entry = next(self.iter)
            if entry.feature == "exon":
                break

        self.counts_gene_ids[entry.gene_id] += 1
        self.counts_transcript_ids[entry.transcript_id] += 1
        self.counts_exons_per_transcript[
            entry.transcript_id].append((entry.start, entry.end))

        return entry

    def next(self):
        return self.__next__()

    def __iter__(self):
        return self

    def __str__(self):

        single_exon_transcripts = 0
        exons_per_transcript = []
        intron_sizes = []
        transcript_lengths = []
        exon_sizes = []

        for x in list(self.counts_exons_per_transcript.values()):

            x.sort()
            x = Intervals.combine(x)
            transcript_lengths.append(x[-1][1] - x[0][0])

            exons_per_transcript.append(len(x))

            for start, end in x:
                exon_sizes.append(end - start)

            if len(x) == 1:
                single_exon_transcripts += 1
                continue

            last_end = x[0][1]
            for start, end in x[1:]:
                intron_sizes.append(start - last_end)
                last_end = end

        return "\t".join(map(str, (len(self.counts_gene_ids),
                                   len(self.counts_transcript_ids),
                                   single_exon_transcripts,
                                   Stats.Summary(exons_per_transcript),
                                   Stats.Summary(exon_sizes),
                                   Stats.Summary(intron_sizes),
                                   Stats.Summary(transcript_lengths),
                                   )))


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id",
                            usage=globals()["__doc__"])

    parser.add_option("--is-gtf", dest="is_gtf", action="store_true",
                      help="input is gtf.")

    parser.set_defaults(
        is_gtf=False,
    )

    (options, args) = E.Start(parser, add_output_options=True)

    if len(args) == 0:
        files = [options.stdin]
    else:
        files = args

    options.stdout.write("track\t%s" % ("\t".join(counter_gff.fields)))

    if options.is_gtf:
        options.stdout.write("\t%s" % ("\t".join(counter_exons.fields)))
    options.stdout.write("\n")

    for f in files:
        if f == options.stdin:
            infile = f
            options.stdout.write("stdin")
        else:
            infile = IOTools.openFile(f)
            options.stdout.write(f)

        counters = []
        if options.is_gtf:
            iterator = GTF.iterator(infile)
            counters.append(counter_gff(iterator))
            counters.append(counter_exons(counters[0]))
        else:
            iterator = GTF.iterator(infile)
            counters.append(counter_gff(iterator))

        c = counters[-1]
        for x in c:
            pass

        for c in counters:
            options.stdout.write("\t%s" % str(c))
        options.stdout.write("\n")

        if infile != sys.stdin:
            infile.close()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

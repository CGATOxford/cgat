'''
maq2assembly.py - convert output from maq cns2view 
===================================================

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

   python maq2assembly.py --help

Type::

   python maq2assembly.py --help

for command line help.

Command line options
--------------------

'''
import sys
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Stats as Stats

USAGE = """python %s [OPTIONS] input1 input2

convert residue statistics (output from maq cns2view) to
a gff file.

Version: $Id: maq2assembly.py 2781 2009-09-10 11:33:14Z andreas $
""" % sys.argv[0]


def reader(infile):

    for line in infile:
        if line[0] == "#":
            continue
        if line.startswith("contig"):
            break

    last_contig, last_pos, start = None, 0, 0
    reads, qualities = [], []
    for line in infile:
        if line[0] == "#":
            continue
        contig, pos, qual, nreads, cov, bestqual = line[:-1].split("\t")
        pos, qual = int(pos) - 1, int(qual)

        if contig != last_contig or pos - 1 != last_pos:
            if last_contig is not None:
                yield last_contig, start, last_pos + 1, reads, qualities
            start = pos
            last_contig = contig
            reads = []
            qualities = []

        last_pos = pos
        reads.append(nreads)
        qualities.append(qual)

    if last_contig is not None:
        yield last_contig, start, last_pos + 1, reads, qualities


class Builder:

    def __init__(self, genome_fasta, genome_queries, options):

        self.mGenomeFasta = genome_fasta
        self.mQueriesFasta = queries_fasta
        self.options = options

        self.mIdFormat = options.output_format

        if self.options.output_filename_pattern:
            self.mOutFile = open(
                self.options.output_filename_pattern % self.mName, "w")
        else:
            self.mOutFile = self.options.stdout

        self.mId = 0

    def __call__(self, id, *args):

        self.mId = id

        self.mOutputId = self.mIdFormat % self.mId
        self.process(*args)
        self.mOutFile.flush()

    def printHeader(self):
        header = self.getHeader()
        if header:
            self.mOutFile.write(header + "\n")

    def getHeader(self):
        return None


class BuilderGFF(Builder):

    """output gff model of read
    """

    mName = "gff"

    def __init__(self, *args, **kwargs):
        Builder.__init__(self, *args, **kwargs)

    def getHeader(self):
        return None

    def process(self, contig, start, end, reads, qualities):

        entry = GTF.Entry()
        entry.start, entry.end = start, end
        entry.gene_id = self.mIdFormat % id
        entry.transcript_id = entry.gene_id
        entry.contig = contig
        entry.feature = "exon"
        entry.source = "maq"

        read_stats = Stats.Summary(reads)

        entry.score = "%5.2f" % read_stats['mean']

        self.mOutFile.write(str(entry) + "\n")


class BuilderCoverage(Builder):

    """compute residue read coverage over all genomic bases that have at least one aligned match
    """

    mName = "coverage"

    def __init__(self, *args, **kwargs):
        Builder.__init__(self, *args, **kwargs)

    def getHeader(self):
        return "id\tcontig\tstart\tend\tsize\tnmatches\tncovered\t%s" % ("\t".join(Stats.DistributionalParameters().getHeaders()))

    def process(self, contig, start, end, reads, qualities):

        aligned = filter(lambda x: x > 0, reads)
        self.mOutFile.write("%s\t%s\t%i\t%i\t%i\t%i\t%i\t%s\n" % (self.mOutputId,
                                                                  contig, start, end, end -
                                                                  start,
                                                                  len(reads),
                                                                  len(aligned),
                                                                  str(Stats.DistributionalParameters(aligned))))


class BuilderQuality(Builder):

    """compute residue quality over all genomic bases that have at least one aligned match
    """

    mName = "quality"

    def __init__(self, *args, **kwargs):
        Builder.__init__(self, *args, **kwargs)

    def getHeader(self):
        return "id\tcontig\tstart\tend\tsize\tnmatches\tncovered\t%s" % ("\t".join(Stats.DistributionalParameters().getHeaders()))

    def process(self, contig, start, end, reads, qualities):

        self.mOutFile.write("%s\t%s\t%i\t%i\t%i\t%i\t%i\t%s\n" % (self.mOutputId,
                                                                  contig, start, end, end -
                                                                  start,
                                                                  len(reads),
                                                                  len(qualities),
                                                                  str(Stats.DistributionalParameters(qualities))))


class BuilderRegion(Builder):

    mName = "regions"

    def __init__(self, *args, **kwargs):
        Builder.__init__(self, *args, **kwargs)

    def getHeader(self):
        return "id\tcontig\tstart\tend\tsize\tnmatches"

    def process(self, contig, start, end, reads, qualities):
        self.mOutFile.write("%s\t%s\t%i\t%i\t%i\t%i\n" % (
            self.mOutputId, contig, start, end, end - start, len(reads)))
        self.mOutFile.flush()


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: maq2assembly.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-f", "--is-forward-coordinates", dest="forward_coordinates",
                      help="translate to forward coordinates.", action="store_true")

    parser.add_option("-p", "--output-filename-pattern", dest="output_filename_pattern", type="string",
                      help="OUTPUT filename pattern for additional data [%default].")

    parser.add_option("--method", dest="methods", type="choice", action="append",
                      choices=("gff", "coverage", "region", "quality"),
                      help="methods to apply [%default].")

    parser.set_defaults(
        output_format="%08i",
        output_filename_pattern="%s",
        methods=[],
    )

    (options, args) = E.Start(parser)

    ################################################
    ################################################
    ################################################
    # pick a processor
    ################################################
    methods = []

    if len(options.methods) == 0:
        raise "please supply at least one method to apply."

    genome_fasta, queries_fasta = None, None

    for method in options.methods:
        if method == "gff":
            methods.append(BuilderGFF(genome_fasta, queries_fasta, options))
        elif method == "coverage":
            methods.append(
                BuilderCoverage(genome_fasta, queries_fasta, options))
        elif method == "quality":
            methods.append(
                BuilderQuality(genome_fasta, queries_fasta, options))
        elif method == "region":
            methods.append(BuilderRegion(genome_fasta, queries_fasta, options))

    for method in methods:
        method.printHeader()

    ninput, noutput = 0, 0
    id = 0
    for contig, start, end, reads, qualities in reader(options.stdin):

        ninput += 1
        id += 1
        for m in methods:
            m(id, contig, start, end, reads, qualities)

        noutput += 1

    options.stdlog.write("# ninput=%i, noutput=%i\n" % (ninput, noutput))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

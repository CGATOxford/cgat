'''
fasta2unique_kmers.py
=============================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Genomics Sequences FASTA Summary

Purpose
-------

This script takes an input :term:`fasta` file and computes
a the number of unique and non-unique k-mers per fasta entry. The
output is a tab delimited file:


Usage
-----
The fasta must be parsed twice, hence the fasta file name must be
provided rather than passing the fasta on stdin

Example


Options
-------

--input-fasta

--kmer-size
   size of kmer - the larger the more memory required

--subset
   only compute the unique kmers over the first n entries.

Type::

   python fasta2unique_kmerst.py --help

for command line help.

TS - Note, the script is expected to be used with a transcriptome
fasta. If used with a genome fasta this will hold in memory the set of
kmers for each contig in turn which could be considerable. The
countUniqueKmers function would need to be re-written to handle a
genome fasta infile. The intention of this function is to only
consider each unique kmer per fasta entry once so an
iterator/generator by itself will not work without a check to see
whether the kmer has already been seen, which brings us back to a set
operation somewhere.

'''

# TS NOTE: this is a very memory hungry script. On the human
# protein-coding transcriptome (~140,000 transcripts) with 31mer, it
# requires ~ 20G RAM!

import sys
import numpy as np
import resource

import CGAT.FastaIterator as FastaIterator
import CGAT.IOTools as IOTools
import CGAT.Experiment as E


# TS. temp function to monitor usage
def using(point=""):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s mb
           ''' % (point, usage[0], usage[1],
                  (usage[2]*resource.getpagesize())/1000000.0)


class KmerCounter(object):

    def __init__(self):
        ''' initiate dictionary and counters '''

        self.kmer2entry = {}

    def shred(self, seq, k_length):
        ''' shred sequence into kmers and update kmer dictionary '''

        kmers = set([seq[i:i+k_length] for i in range(0, len(seq)-k_length)])

        for kmer in kmers:
            if kmer in self.kmer2entry:
                self.kmer2entry[kmer] = 0
            else:
                self.kmer2entry[kmer] = 1

    def countUniqueKmers(self, seq, k_length):
        ''' shred sequence and count unique and non-unique k-mers'''

        unique_kmers = 0
        non_unique_kmers = 0

        kmers = set([seq[i:i+k_length] for i in range(0, len(seq)-k_length)])

        for kmer in kmers:

            if self.kmer2entry[kmer] == 1:
                unique_kmers += 1
            else:
                non_unique_kmers += 1

        return (unique_kmers, non_unique_kmers)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("--input-fasta", dest="fasta", type="str",
                      help="name of fasta infile")

    parser.add_option("-k", "--kmer-size", dest="kmer", type="int",
                      help="supply kmer length")

    parser.add_option("--subset", dest="subset", type="int",
                      help="only analyse the first x entries")

    parser.set_defaults(
        fasta=None,
        kmer=10,
        subset=None)

    (options, args) = E.Start(parser)

    E.info("%s\n" % using("start"))

    assert options.fasta, "must provide a fasta filename (--input-fasta=)"

    k = KmerCounter()

    Iterator = FastaIterator.iterate(IOTools.openFile(options.fasta))

    # total entries also acts as the index for the entry_id
    total_entries = 0

    options.stdout.write("%s\n" % "\t".join((
        "id", "unique_kmers", "non_unique_kmers", "fraction_unique")))

    # iterate transcripts, shred and identify unique kmers
    E.info("shredding fasta to identify unique kmers")
    for entry in Iterator:
        if total_entries % 1000 == 0:
            E.info("1st shred complete for %i entries" % total_entries)

        if options.subset and total_entries >= options.subset:
            break

        k.shred(entry.sequence.upper(), options.kmer)
        total_entries += 1

    total_entries = 0
    Iterator = FastaIterator.iterate(IOTools.openFile(options.fasta))

    # iterate transcripts, shread and count unique kmers
    E.info("re-shredding fasta to count unique kmers")

    for entry in Iterator:
        if total_entries % 1000 == 0:
            E.info("2nd shred complete for %i entries" % total_entries)

        transcript_id = entry.title.split()[0]

        if options.subset and total_entries >= options.subset:
            break

        total_entries += 1

        unique, non_unique = k.countUniqueKmers(
            entry.sequence.upper(), options.kmer)

        fraction = np.divide(float(unique), (unique + non_unique))

        options.stdout.write("%s\n" % "\t".join(
            map(str, (transcript_id, unique, non_unique, fraction))))

    E.info("found %i kmers" % len(k.kmer2entry))
    E.info("written kmer counts for %i contigs" % total_entries)
    # write footer and output benchmark information.
    E.info("%s\n" % using("end"))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

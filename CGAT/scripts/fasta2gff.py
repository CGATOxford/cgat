'''
fasta2gff.py - create random segments from fasta file
=====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics Sequences Intervals FASTA GFF Transformation

Purpose
-------

This script creates a random sample of "exons" from 
a fasta file.

Usage
-----

Example::

   python fasta2gff.py --help

Type::

   python fasta2gff.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.GTF as GTF

USAGE = """python %s [OPTIONS] 



Version: $Id: fasta2gff.py 2861 2010-02-23 17:36:32Z andreas $
""" % sys.argv[0]


def writeHeader(outfile):
    outfile.write("\t".join(("contig",
                             "nresidues",
                             "ngaps",
                             "nseqregions",
                             "ngapregions",
                             "nA", "nC", "nG", "nT",
                             "nN", "nX", "nO")) + "\n")


def main(argv=None):

    parser = E.OptionParser(
        version="%prog version: $Id: fasta2gff.py 2861 2010-02-23 17:36:32Z andreas $")

    parser.add_option("-g", "--genome-file", dest="genome_file", type="string",
                      help="filename with genome.")

    parser.add_option("-a", "--as-gtf", dest="as_gtf", action="store_true",
                      help="output as gtf.")

    parser.add_option("-f", "--fragment-size", dest="fragment_size", type="int",
                      help="fixed size of fragments [default=%default].")

    parser.add_option("-s", "--sample-size", dest="sample_size", type="int",
                      help="fixed size of fragments.")

    parser.set_defaults(
        as_gtf=False,
        genome_file=None,
        fragment_size=1000,
        sample_size=10000,
        pattern_id="%08i",
    )

    (options, args) = E.Start(parser)

    fasta = IndexedFasta.IndexedFasta(options.genome_file)
    contigs = fasta.getContigSizes()

    if options.as_gtf:
        entry = GTF.Entry()
    else:
        entry = GTF.Entry()

    n = 0
    entry.feature = "exon"
    entry.source = "random"

    for x in range(options.sample_size):

        entry.contig, entry.strand, entry.start, entry.end = fasta.getRandomCoordinates(
            options.fragment_size)

        if entry.strand == "-":
            l = contigs[entry.contig]
            entry.start, entry.end = l - entry.end, l - entry.start

        if options.as_gtf:
            entry.gene_id = options.pattern_id % n
            entry.transcript_id = entry.gene_id

        options.stdout.write(str(entry) + "\n")
        n += 1

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

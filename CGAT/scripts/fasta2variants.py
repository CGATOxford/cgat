'''fasta2variants.py - create sequence variants from a set of sequences
====================================================================

:Tags: Genomics Sequences Variants Protein FASTA Transformation

Purpose
-------

This script reads a collection of sequences in :term:`fasta` format
and outputs a table of possible variants. It outputs for each position
in a protein sequence the number of variants.

If the input sequences are nucleotide coding (CDS) sequences, for each
variant a weight is output indicating the number of times that variant
can occur from single nucleotide changes.

Usage
-----

Example::

    python fasta2variants.py -I CCDS_nucleotide.current.fna.gz -L CDS.log -S CDS.output -c

This will take a CDS file as input, save the log and output files, and
count variants based on single nucleotide changes using the -c option.

Type::

    python fasta2variants.py --help

for command line help.

Compressed (.gz) and various fasta format files (.fasta, .fna) are
accepted. If the -c option is specified and the file is not a CDS
sequence the script will throw an error ('length of sequence
'<input_file>' is not a multiple of 3').

Command line options
--------------------

'''
import sys
import collections

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.FastaIterator as FastaIterator


def main(argv=None):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-c", "--is-cds", dest="is_cds", action="store_true",
                      help="input are cds (nucleotide) sequences [%default]")

    parser.set_defaults(
        is_cds=False,
    )

    (options, args) = E.Start(parser, argv=argv)

    options.stdout.write(
        "snpid\tidentifier\tpos\treference\tvariant\tcounts\tweight\n")

    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    snpid = 0

    for entry in FastaIterator.iterate(options.stdin):
        identifier = entry.title

        if options.is_cds:
            cds_sequence = entry.sequence.upper()
            assert len(cds_sequence) % 3 == 0, \
                "length of sequence '%s' is not a multiple of 3" % entry.title

            sequence = Genomics.translate(cds_sequence)
            weights = []
            for pos, cds_pos in enumerate(range(0, len(cds_sequence), 3)):
                codon = cds_sequence[cds_pos:cds_pos + 3]
                counts = collections.defaultdict(int)
                for x in range(0, 3):
                    rna = codon[x]
                    for na in "ACGT":
                        if na == rna:
                            continue
                        taa = Genomics.translate(
                            codon[:x] + na + codon[x + 1:])
                        counts[taa] += 1
                weights.append(counts)

        else:
            sequence = entry.sequence.upper()
            counts = {}
            for x in alphabet:
                counts[x] = 1
            weights = [counts] * len(sequence)

        for pos, ref in enumerate(sequence):

            if ref not in alphabet:
                continue
            w = weights[pos]
            t = float(sum(w.values()))
            for variant in alphabet:
                if variant == ref:
                    continue
                snpid += 1
                options.stdout.write(
                    "%s\n" % "\t".join(
                        ("%010i" % snpid,
                         identifier,
                         str(pos + 1),
                         ref,
                         variant,
                         "%i" % w[variant],
                         "%6.4f" % (w[variant] / t),
                         )))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

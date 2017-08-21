'''fasta2table.py - analyze sequence composition, codon usage, bias and degeneracy 
====================================================================================

:Tags: Genomics Sequences

Purpose
-------

This script reads a collection of sequences in :term:`fasta` format, computes
various sequence properties and outputs them to a tab-separated file.

Documentation
--------------

The following sequence properties can be reported:

length
   sequence length and number of codons (0 for amino-acid sequence)

sequence
   Full nucleotide / amino-acid sequence

hid
   a hash identifier for the sequence

na
   nucleic acid composition, including GC/AT content

dn
   dinucleotide counts

cpg
   CpG density and CpG observed/expected

gaps
    number of gaps and gapped/ungapped regions in the sequences

aa
   amino acid composition (nucleotide sequence must have length divisible by 3)

degeneracy
    count the number of degenerate sites (nucleotide sequence only,
    sequence must have length divisible by 3)

codons
   codon composition (nucleotide sequence only, sequence must have
   length divisible by 3)

codon-usage
    output codon frequencies for each sequence (nucleotide sequence
    only, sequence must have length divisible by 3)

codon-bias
    output codon bias for each sequence (nucleotide sequence only,
    sequence must have length divisible by 3)

codon-translator
    translate codons for each sequence to their frequency (nucleotide
    sequence only, sequence must have length divisible by 3)

Multiple counters can be calculated at the same by specifying 
--section multiple times.

The script can also process fasta description lines (starting >)
either by splitting each line at the first space and taking only the
first part (--split-fasta-identifier), or by any user-supplied python
regular expression (--regex-identifier).

Usage
-----

Example::

   # View fasta file
   head tests/fasta2table.py/na_test.fasta

   # Count CpG dinucleotides
   cgat fasta2table --sections=length,na --split-fasta-identifier < tests/fasta2table.py/test.fasta > na.tsv

In this example we input a fasta file and compute the sequence composition, i.e.
%C, %G, %A, %T as well for each sequence in the set.

+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|id              |length       |nC     |nG     |nA     |nT     |nN   |nUnk |nGC    |nAT    |nCpG  |pC        |pG        |pA        |pT        |pN        |pUnk      |pGC       |pAT       |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|contig-0        |110529       |26348  |22351  |34036  |27794  |0    |0    |48699  |61830  |4990  |0.238381  |0.202218  |0.307937  |0.251463  |0.000000  |0.000000  |0.440599  |0.559401  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|contig-1000000  |121046       |28432  |22851  |36540  |33223  |0    |0    |51283  |69763  |4814  |0.234886  |0.188779  |0.301869  |0.274466  |0.000000  |0.000000  |0.423665  |0.576335  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|contig-2000000  |61046        |11368  |14668  |16666  |18344  |0    |0    |26036  |35010  |2526  |0.186220  |0.240278  |0.273007  |0.300495  |0.000000  |0.000000  |0.426498  |0.573502  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|contig-3000000  |76120        |17653  |14496  |23591  |20380  |0    |0    |32149  |43971  |3099  |0.231910  |0.190436  |0.309919  |0.267735  |0.000000  |0.000000  |0.422346  |0.577654  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|contig-4000000  |66656        |12970  |14994  |19076  |19616  |0    |0    |27964  |38692  |2911  |0.194581  |0.224946  |0.286186  |0.294287  |0.000000  |0.000000  |0.419527  |0.580473  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|contig-5000000  |81994        |15755  |19061  |23070  |24108  |0    |0    |34816  |47178  |3571  |0.192148  |0.232468  |0.281362  |0.294022  |0.000000  |0.000000  |0.424616  |0.575384  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+
|contig-6000000  |106529       |16997  |21653  |21112  |22258  |0    |0    |38650  |43370  |5135  |0.207230  |0.263997  |0.257401  |0.271373  |0.000000  |0.000000  |0.471227  |0.528773  |
+----------------+-------------+-------+-------+-------+-------+-----+-----+-------+-------+------+----------+----------+----------+----------+----------+----------+----------+----------+

Type::

   cgat fasta2table --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import math

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
import CGAT.SequenceProperties as SequenceProperties
import CGAT.FastaIterator as FastaIterator


def main(argv=None):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-w", "--weights-tsv-file", dest="filename_weights",
        type="string",
        help="filename with codon frequencies. Multiple filenames "
        "can be separated by comma.")

    parser.add_option(
        "-s", "--section", dest="sections", type="choice", action="append",
        choices=("length", "sequence", "hid", "na", "aa", "cpg", "dn",
                 "degeneracy", "gaps",
                 "codons", "codon-usage", "codon-translator", "codon-bias"),
        help="which sections to output [%default]")

    parser.add_option(
        "-t", "--sequence-type", dest="seqtype", type="choice",
        choices=("na", "aa"),
        help="type of sequence: na=nucleotides, aa=amino acids [%default].")

    parser.add_option(
        "-e", "--regex-identifier", dest="regex_identifier", type="string",
        help="regular expression to extract identifier from fasta "
        "description line.")

    parser.add_option(
        "--split-fasta-identifier", dest="split_id",
        action="store_true",
        help="split fasta description line (starting >) and use "
        "only text before first space")

    parser.add_option(
        "--add-total", dest="add_total", action="store_true",
        help="add a row with column totals at the end of the table"
        "[%default]")

    parser.set_defaults(
        filename_weights=None,
        pseudocounts=1,
        sections=[],
        regex_identifier="(.+)",
        seqtype="na",
        gap_chars='xXnN',
        split_id=False,
        add_total=False,
    )

    (options, args) = E.Start(parser, argv=argv)

    rx = re.compile(options.regex_identifier)

    reference_codons = []
    if options.filename_weights:
        options.filename_weights = options.filename_weights.split(",")
        for filename in options.filename_weights:
            if filename == "uniform":
                reference_codons.append(Genomics.GetUniformCodonUsage())
            else:
                reference_codons.append(
                    IOTools.ReadMap(IOTools.openFile(filename, "r"),
                                    has_header=True,
                                    map_functions=(str, float)))

        # print codon table differences
        options.stdlog.write(
            "# Difference between supplied codon usage preferences.\n")
        for x in range(0, len(reference_codons)):
            for y in range(0, len(reference_codons)):
                if x == y:
                    continue
                # calculate KL distance
                a = reference_codons[x]
                b = reference_codons[y]
                d = 0
                for codon, p in list(a.items()):
                    if Genomics.IsStopCodon(codon):
                        continue
                    d += b[codon] * math.log(b[codon] / p)

                options.stdlog.write("# tablediff\t%s\t%s\t%f\n" %
                                     (options.filename_weights[x],
                                      options.filename_weights[y],
                                      d))

    iterator = FastaIterator.FastaIterator(options.stdin)

    def getCounter(section):

        if options.seqtype == "na":
            if section == "length":
                s = SequenceProperties.SequencePropertiesLength()
            elif section == "sequence":
                s = SequenceProperties.SequencePropertiesSequence()
            elif section == "hid":
                s = SequenceProperties.SequencePropertiesHid()
            elif section == "na":
                s = SequenceProperties.SequencePropertiesNA()
            elif section == "gaps":
                s = SequenceProperties.SequencePropertiesGaps(
                    options.gap_chars)
            elif section == "cpg":
                s = SequenceProperties.SequencePropertiesCpg()
            elif section == "dn":
                s = SequenceProperties.SequencePropertiesDN()
            # these sections requires sequence length to be a multiple of 3
            elif section == "aa":
                s = SequenceProperties.SequencePropertiesAA()
            elif section == "degeneracy":
                s = SequenceProperties.SequencePropertiesDegeneracy()
            elif section == "codon-bias":
                s = SequenceProperties.SequencePropertiesBias(reference_codons)
            elif section == "codons":
                s = SequenceProperties.SequencePropertiesCodons()
            elif section == "codon-usage":
                s = SequenceProperties.SequencePropertiesCodonUsage()
            elif section == "codon-translator":
                s = SequenceProperties.SequencePropertiesCodonTranslator()
            else:
                raise ValueError("unknown section %s" % section)
        elif options.seqtype == "aa":
            if section == "length":
                s = SequenceProperties.SequencePropertiesLength()
            elif section == "sequence":
                s = SequenceProperties.SequencePropertiesSequence()
            elif section == "hid":
                s = SequenceProperties.SequencePropertiesHid()
            elif section == "aa":
                s = SequenceProperties.SequencePropertiesAminoAcids()
            else:
                raise ValueError("unknown section %s" % section)
        return s

    # setup totals
    totals = {}
    for section in options.sections:
        totals[section] = getCounter(section)

    options.stdout.write("id")
    for section in options.sections:
        options.stdout.write("\t" + "\t".join(totals[section].getHeaders()))

    options.stdout.write("\n")
    options.stdout.flush()

    s = getCounter("hid")
    s.loadSequence("AAAAAAAAA", "na")

    for cur_record in iterator:

        sequence = re.sub(" ", "", cur_record.sequence).upper()

        if len(sequence) == 0:
            raise ValueError("empty sequence %s" % cur_record.title)

        id = rx.search(cur_record.title).groups()[0]

        if options.split_id is True:
            options.stdout.write("%s" % id.split()[0])
        else:
            options.stdout.write("%s" % id)
        options.stdout.flush()

        for section in options.sections:
            s = getCounter(section)
            s.loadSequence(sequence, options.seqtype)
            totals[section].addProperties(s)

            options.stdout.write("\t" + "\t".join(s.getFields()))

        options.stdout.write("\n")

    if options.add_total:
        options.stdout.write("total")
        for section in options.sections:
            options.stdout.write("\t" + "\t".join(totals[section].getFields()))
        options.stdout.write("\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

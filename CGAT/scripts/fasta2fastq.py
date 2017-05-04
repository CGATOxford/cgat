'''fasta2fastq.py - simulate reads from fasta
=====================================

:Tags: Sequences

Purpose
-------

Simulate illumina sequence reads from a fasta file. The number of
reads per entry is randomly selected from the range given. The primary
use case is expected to be the generation of simulation RNA-Seq reads

For RNA-Seq simulations, the premrna-fraction option allows the user to specify
what fraction of the transcripts originate from pre-mRNA. The user must also
supply a second fasta in the same order for the pre-mRNA
(--infile-premrna-fasta). The simulation assumes all pre-mRNA are full length
which is not likely to be the case for real RNA-Seq.
Note: This may lead to many more reads which align to the mRNA than the
apparent ground truth count. It is therefore best to keep the pre-mRNA fraction
low (recommend 0.01).

--counts-method determine whether the --min-counts and --max--counts methods
relate to the number of reads or copies per entry. If copies, then the number
of reads/read pairs per transcript will be determined by the copy number, the
length of the transcript and the fragment length. E.g for copy number = 2,
transcript length = 1000 bp, read length = 100 bp, paired end, 300bp insert
size, the number of reads will be 2 * (1000 / (( 2 * 100) + 300)) = 4

Options
-------

--output-paired-end
   generate paired-end reads (defaults to single end)

--counts-method
  simulate a ground truth number of reads or copies per transcript
  if reads, then min and max specify the number of reads per entry
  if copies, then min and specify the sequencing depth per entry

--counts-min
   the minimum number of reads to simulate for each fasta entry

--counts-max
   the maximum number of reads to simulate for each fasta entry

--sequence-error-phred
   the sequencing error rate (phred scale)

--output-read-length
   the length of the outputted reads

--output-counts
   filename for counts per fasta entry

--output-quality-format
   the format of the sequence qualities (+33 = Sanger)

--insert-length-mean
   the mean insert length

--insert-length-sd
   the standard deviation for the insert length

--premrna-fraction
   the fraction of reads to simulate from pre-mRNA. Default is 0.
   If set, must provide a pre-mRNA fasta file with:
      --infile-premrna-fasta

If generating paired end reads, the second outfile must be specified with:

--output-fastq2


Usage
-----

Recommend sending logging to separate outfile to keep fastq outfile
clean of comments (see example below)

Example::

   cat transcripts.fa | python fasta2fastq.py
   --output-counts=simulation_counts.tsv -L simulation.log
   > simulation_reads.fastq

Type::

   python fasta2fastq.py --help

for command line help.


Important note for generating reads for simulations
---------------------------------------------------
Currently, the output is non-random, e.g it's in the order of the
fasta input. If you want the fastq to be random pipe the output to
sort -R like so:

   cat transcripts.fa | python fasta2fastq.py
   --output-counts=simulation_counts.tsv -L simulation.log |
   paste - - - - |sort -R | sed 's/\t/\n/g' > simulation_reads_random.fastq

If you're outputting paired end fastqs, you can use the following
command to randomise the order by keep the fastq entris paired:

    paste <(zcat %(fastq1_ordered)s) <(zcat %(fastq2_ordered)s) |
    paste - - - - | sort -R | awk -F'\t' '{OFS="\n"; print $1,$3,$5,$7 >
    "%(fastq1_random)s"; print $2,$4,$6,$8 > "%(fastq2_random)s"}'

'''
import sys
import random
import numpy as np
import collections

import CGAT.Experiment as E
import CGAT.IOTools as IOTools

import CGAT.FastaIterator as FastaIterator


def addSeqErrors(read=None, error_rate=10):
    ''' add sequencing errors to a read.
    Error rates are Phred scaled, so 30 = 1/1000'''

    error_rate = 10**(error_rate/-10.0)

    errors_dict = {"G": ["C", "T", "A"],
                   "C": ["G", "T", "A"],
                   "T": ["C", "G", "A"],
                   "A": ["C", "T", "G"],
                   "N": ["C", "T", "G", "A"]}

    probs = np.random.rand(len(read))
    return "".join([base if prob > error_rate and base != "N"
                    else random.choice(errors_dict[base])
                    for prob, base in zip(probs, read)])


def reverseComp(seq):
    ''' return the reverse complement sequence '''

    comp = {"G": "C",
            "C": "G",
            "A": "T",
            "T": "A",
            "N": "N"}

    return "".join([comp[base] for base in seq[::-1]])


def generateRead(entry, read_length=50, error_rate=40, paired=False,
                 insert_mean=0, insert_sd=1):
    ''' generate a read (or read pair) at random from a fasta entry for
    the given read length with sequencing errors according to error
    rate'''

    if paired:

        position = "not_OK"

        while position != "OK":

            r1_start = random.randint(0, len(entry)-read_length)
            r2_start = (r1_start + read_length +
                        int(np.random.normal(insert_mean, insert_sd)))

            if (r2_start <= (len(entry) - read_length) and r2_start >= r1_start):

                position = "OK"

                read1 = entry[r1_start: r1_start+read_length]
                read2 = reverseComp(
                    entry[r2_start: r2_start+read_length])

                final_read1 = addSeqErrors(read1, error_rate)
                final_read2 = addSeqErrors(read2, error_rate)

                return final_read1, final_read2

    else:
        start = random.randint(0, len(entry)-read_length)
        read = entry[start:start+read_length]

        final_read = addSeqErrors(read, error_rate)

        return final_read



def getTitle(entry):
    ''' return the title for an entry'''
    return entry.title.split()[0]



def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version",
                            usage=globals()["__doc__"])

    parser.add_option(
        "--output-quality-format", dest="q_format", type="int",
        help="sequence quality format, e.g 33 = +33/Sanger"
        "[default=%default].")

    parser.add_option(
        "--output-paired-end", dest="paired", action="store_true",
        help="generate paired end reads [default = %default].")

    parser.add_option(
        "--insert-length-mean", dest="insert_mean", type="float",
        help="mean insert length [default = %default].")

    parser.add_option(
        "--insert-length-sd", dest="insert_sd", type="float",
        help="insert length standard deviation [default = %default].")

    parser.add_option(
        "--counts-method", dest="counts_method", type="choice",
        choices=("reads", "copies"),
        help="simulate a ground truth number of reads per entry or"
        "copies per entry [default = %default].")

    parser.add_option(
        "--counts-min", dest="counts_min", type="float",
        help="minimum number of reads/read pairs per fasta entry"
        "or copies per entry [default = %default].")

    parser.add_option(
        "--counts-max", dest="counts_max", type="float",
        help="maximum number of reads/read pairs per fasta entry "
        "or copies per entry [default = %default].")

    parser.add_option(
        "--output-read-length", dest="read_length", type="int",
        help="read length [default = %default].")

    parser.add_option(
        "--sequence-error-phred", dest="phred", type="int",
        help="phred quality score [default = %default].")

    parser.add_option(
        "--output-counts", dest="output_counts", type="string",
        help="name for counts outfile [default=%default].")

    parser.add_option(
        "--output-fastq2", dest="fastq2_out", type="string",
        help="filename for second fastq outfile [default=%default].")

    parser.add_option(
        "--premrna-fraction", dest="premrna_fraction", type="float",
        help="the fraction of reads to simulate from pre-mRNA"
        "[default= % default].")

    parser.add_option(
        "--infile-premrna-fasta", dest="premrna_fasta", type="string",
        help="filename for pre-mRNA fasta[default=%default].")

    parser.set_defaults(
        q_format=33,
        paired=False,
        insert_mean=0,
        insert_sd=1,
        counts_method="reads",
        counts_min=1,
        counts_max=1,
        read_length=50,
        fastq2_out=None,
        output_counts=None,
        phred=30,
        premrna_fraction=0,
        premrna_fasta=None
    )

    (options, args) = E.Start(parser)

    if options.paired:
        assert options.fastq2_out, ("must specify a second fastq outfile for "
                                    "paired end (--output-fastq2)")
        outf2 = IOTools.openFile(options.fastq2_out, "w")

    if options.premrna_fraction:
        assert options.premrna_fasta, ("must specfify the location of the"
                                       "fasta file for the pre-mRNA")

    # the sequence quality string will always be the same so define here
    sequence_quality = chr(options.q_format + options.phred)
    qual = "".join([sequence_quality] * options.read_length)

    if options.premrna_fraction:
        iterator = FastaIterator.iterate_together(
            options.stdin, IOTools.openFile(options.premrna_fasta))
    else:
        iterator = FastaIterator.FastaIterator(options.stdin)

    # set a cut off of twice the read/pair length for short entries
    if options.paired:
        minimum_entry_length = (
            2 * ((options.read_length * 2) + options.insert_mean))
    else:
        minimum_entry_length = 2 * options.read_length

    c = collections.Counter()
    counts = collections.Counter()
    copies = collections.Counter()

    for f_entry in iterator:

        if options.premrna_fraction:

            assert getTitle(f_entry[0]) == getTitle(f_entry[1]), (
                "entry ids do not match: %s != %s" % (
                    f_entry[0].title, f_entry[1].title))
            entry = f_entry[0]
            pre_entry = f_entry[1]

        else:
            entry = f_entry

        # reject short fasta entries
        if len(entry.sequence) < minimum_entry_length:
            E.info("skipping short transcript: %s length=%i"
                   % (entry.title, len(entry.sequence)))
            c['skipped'] += 1
            continue

        else:
            c['not_skipped'] += 1

        if options.paired:
            fragment_length = (
                (2 * options.read_length) + options.insert_mean)
        else:
            fragment_length = options.read_length

        reads_per_entry = float(len(entry.sequence)) / fragment_length

        if options.counts_method == "reads":
            n_reads = random.randint(options.counts_min,
                                     options.counts_max + 1)

            n_copies = float(n_reads) / reads_per_entry

            if options.premrna_fraction:
                n_reads_pre = int(round(n_reads * options.premrna_fraction))

        elif options.counts_method == "copies":

            # random float [0-1]
            rand = np.random.random_sample()
            n_copies = (options.counts_min +
                        (rand * (options.counts_max - options.counts_min)))

            n_reads = int(round(n_copies * reads_per_entry, 0))

            # as n_reads must be rounded to int, need to redefine n_copies
            n_copies = float(n_reads) / reads_per_entry

            if options.premrna_fraction:
                reads_per_pre_entry = (float(len(pre_entry.sequence)) /
                                       fragment_length)
                n_copies_pre = n_copies * options.premrna_fraction
                n_reads_pre = int(round(n_copies_pre * reads_per_pre_entry, 0))
                # as n_reads_pre must be rounded to int, need to
                # redefine n_copies_pre
                n_copies_pre = float(n_reads_pre) / reads_per_pre_entry

        entry_id = getTitle(entry)

        counts[entry_id] = n_reads
        copies[entry_id] = n_copies

        if "N" in entry.sequence.upper():
            E.warn("fasta entry %s contains unknown bases ('N')" % entry_id)

        for i in range(0, n_reads):

            read = generateRead(entry=entry.sequence.upper(),
                                read_length=options.read_length,
                                error_rate=options.phred,
                                paired=options.paired,
                                insert_mean=options.insert_mean,
                                insert_sd=options.insert_sd)

            if options.paired:
                r1, r2 = read
                h1 = "@%s_%i/1" % (entry_id, i)
                h2 = "@%s_%i/2" % (entry_id, i)

                options.stdout.write("\n".join((h1, r1, "+", qual)) + "\n")
                outf2.write("\n".join((h2, r2, "+", qual)) + "\n")

            else:
                h = "@%s_%i/1" % (entry_id, i)

                options.stdout.write("\n".join((h, read, "+", qual)) + "\n")

        if options.premrna_fraction:
            c['pre_counts'] += n_reads_pre
            c['pre_copies'] += n_copies_pre

            for i in range(0, n_reads_pre):

                read = generateRead(entry=pre_entry.sequence.upper(),
                                    read_length=options.read_length,
                                    error_rate=options.phred,
                                    paired=options.paired,
                                    insert_mean=options.insert_mean,
                                    insert_sd=options.insert_sd)

                if options.paired:
                    r1, r2 = read
                    h1 = "@%s_pre-mRNA_%i/1" % (entry_id, i)
                    h2 = "@%s_pre-mRNA_%i/2" % (entry_id, i)

                    options.stdout.write("\n".join((h1, r1, "+", qual)) + "\n")
                    outf2.write("\n".join((h2, r2, "+", qual)) + "\n")

                else:
                    h = "@%s_pre-mRNA_%i/1" % (entry_id, i)

                    options.stdout.write("\n".join((h, read, "+", qual)) + "\n")

    if options.paired:
        outf2.close()

    with IOTools.openFile(options.output_counts, "w") as counts_out:

        counts_out.write("%s\n" % "\t".join(("id", "read_count", "tpm")))

        sum_copies = sum(copies.values())
        sum_counts = sum(counts.values())

        for entry_id, count in counts.items():
            tpm = 1000000 * (float(copies[entry_id]) / sum_copies)
            counts_out.write(
                "%s\n" % "\t".join(map(str, (entry_id, count, tpm))))

    E.info("Reads simulated for %i fasta entries, %i entries skipped"
           % (c['not_skipped'], c['skipped']))

    E.info("Simulated: %i reads (%i mRNA, %i pre-mRNA), "
           "%f transcripts (%f mRNA, %f pre-mRNA)" % (
               sum_counts + c['pre_counts'], sum_counts, c['pre_counts'],
               sum_copies + c['pre_copies'], sum_copies, c['pre_copies']))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

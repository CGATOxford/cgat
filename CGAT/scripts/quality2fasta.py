'''
quality2fasta.py - convert quality scores to fasta/fastq output
===============================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Convert quality scores to a fasta or fastq formatted file. 
Quality scores are positive and mapped to the range [a-zA-Z0-9] 
starting at 0.

Usage
-----

Example::

   python quality2fasta.py --help

Type::

   python quality2fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import CGAT.Experiment as E
import CGAT.IndexedFasta as IndexedFasta
import CGAT.FastaIterator as FastaIterator


class FromFastaIterator:

    def __init__(self, infile, alphabet="fastq", encoding="phred",
                 default=None):
        """default: set all quality scores to this default value."""

        self.mInputIterator = FastaIterator.FastaIterator(infile)
        self.mOutputIterator = self._iterate()
        self.mNOverFlow = 0
        self.mNUnderFlow = 0
        self.mDefault = default

        # how to convert a phred like score to a character:
        if alphabet == "fastq":
            self.mMapScore2Char = [chr(33 + x) for x in range(0, 93)]
        elif alphabet == "solexa":
            self.mMapScore2Char = [chr(64 + x) for x in range(0, 128)]
        elif alphabet == "printable":
            self.mMapScore2Char = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        else:
            raise "unknown alphabet %s" % alphabet

        if encoding == "phred":
            # no change for phred scores
            self.mMapScore2Score = list(range(0, 93))
        elif encoding == "solexa":
            # solexa encoding
            self.mMapScore2Score = [
                int(10.0 * log(1.0 + 10 ** (x / 10.0)) / log(10) + .499) for x in range(-64, 65)]

        self.mOverFlow = "^"
        self.mUnderFlow = "."
        self.mNInput = 0
        self.mNOutput = 0

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.mOutputIterator)

    def _iterate(self):
        """iterate over muliple files."""

        while 1:
            cur_entry = next(self.mInputIterator)

            self.mNInput += 1
            if self.mDefault:
                values = [self.mDefault] * \
                    len(re.sub("\s", "", cur_entry.sequence))
            else:
                values = list(map(int, re.split(" +", cur_entry.sequence.strip())))

            s = []
            for v in values:

                v = self.mMapScore2Score[v]

                if v < 0:
                    c = self.mUnderflow
                    self.mNUnderFlow += 1
                try:
                    c = self.mMapScore2Char[v]
                except IndexError:
                    c = self.mOverFlow
                    self.mNOverFlows += 1

                s.append(c)

            self.mNOutput += 1
            yield cur_entry.title.strip(), "".join(s)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: quality2fasta.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-f", "--format", dest="format", type="choice",
                      choices=("fasta", ),
                      help="input format [%default].")

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("fasta", "fastq"),
                      help="output format - if fastq is chosen, also supply a sequence file [%default].")

    parser.add_option("-a", "--alphabet", dest="alphabet", type="choice",
                      choices=("fastq", "solexa", "printable"),
                      help="characters to use for quality scores [%default].")

    parser.add_option("-e", "--encoding", dest="encoding", type="choice",
                      choices=("phred", "solexa"),
                      help="encoding of quality scores [%default].")

    parser.add_option("-i", "--build-index", dest="build_index", type="string",
                      help="build an index. Supply the database name [%default].")

    parser.add_option("-s", "--filename-sequences", dest="filename_sequences", type="string",
                      help="input filename with file of sequences in fasta format - sorted in the same way as the quality file [%default].")

    parser.add_option("-d", "--set-to-default", dest="default_value", type="int",
                      help="set all quality codes to the default value. Supply the fasta sequence instead of the quality codes [%default].")

    parser.set_defaults(
        format="fasta",
        output_format="fasta",
        build_index=None,
        filename_sequences=None,
        alphabet="fastq",
        encoding="phred",
        default_value=None,
    )

    (options, args) = E.Start(parser)

    ninput, noutput = 0, 0

    if options.format == "fasta":
        iterator = FromFastaIterator(
            sys.stdin, alphabet=options.alphabet, default=options.default_value)

    if options.output_format == "fasta":

        if options.build_index:
            IndexedFasta.createDatabase(options.build_index,
                                        iterator)
        else:
            while 1:
                try:
                    r = next(iterator)
                except StopIteration:
                    break
                t, s = r
                options.stdout.write(">%s\n%s\n" % (t, s))

    elif options.output_format == "fastq":

        if not options.filename_sequences:
            raise ValueError("please supply a filename with sequences")

        iterator_sequence = FastaIterator.FastaIterator(
           IOTools.openFile(options.filename_sequences, "r"))

        while 1:
            qual, seq = None, None
            try:
                qual = next(iterator)
                seq = next(iterator_sequence)
            except StopIteration:
                if qual and not seq:
                    options.stdlog.write("# sequence file incomplete\n")
                elif seq and not qual:
                    options.stdlog.write("# quality file incomplete\n")

            qt, qs = qual
            st, ss = seq.title, seq.sequence
            assert qt == st, "sequence and quality identifiers incongruent: %s != %s" % (
                qt, st)
            options.stdout.write("@%s\n%s\n+\n%s\n" % (qt, ss, qs))

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, noverflow=%i, nunderflow=%i\n" %
                             (iterator.mNInput,
                              iterator.mNOutput,
                              iterator.mNOverFlow,
                              iterator.mNUnderFlow))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

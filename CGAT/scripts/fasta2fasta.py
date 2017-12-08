'''fasta2fasta.py - operate on sequences
=====================================

:Tags: Sequences

Purpose
-------

perform operations (masking, renaming) on a stream of fasta formatted sequences.

Available edit operations are:

translate
   translate sequences using the standard genetic code.

translate-to-stop
   translate until first stop codon

truncate-at-stop
   truncate sequence at first stop codon

back-translate
   convert nucleotide sequence to peptide sequence
   Requires parameter of second fasta file with peptide sequences.

mark-codons
   adds a space after each codon

apply-map
   rename sequence identifiers from a given map Requires parameter
   with filename of a map. The map is a tab-separated file mapping old
   to new names.

build-map
   rename sequence identifiers numerically and save output in a
   tab-separated file.  Requires parameter with filename of a map. The
   map is a tab-separated file mapping new to old names and will be
   newly created. Any exiting file of the same name will be
   overwritten.

pseudo-codons
   translate, but keep register with codons

interleaved-codons
   mix amino acids and codons

filter
   remove sequence according to certain criteria. For example,
   --method=filter --filter-method=min-length=5  --filter-method=max-length=10

map-codons:

remove-gaps
   remove all gaps in the sequence

mask-stops
   mask all stop codons

mask-seg
   mask sequence by running seg

mask-bias
   mask sequence by running bias

mask-codons
   mask codon sequence given a masked amino acid sequence.
   Requires parameter with masked amino acids in fasta format.

mask-incomplete-codons
   mask codons that are partially masked or gapped

mask-soft
   combine hard-masked (NNN) sequences with unmasked sequences to generate
   soft masked sequence (masked regions in lower case)

remove-stops
   remove stop codons

upper
   convert sequence to upper case

lower
   convert sequence to lower case

reverse-complement
   build the reverse complement

shuffle
   shuffle each sequence

sample
   select a certain proportion of sequences

Parameters are given to the option ``parameters`` in a comma-separated
list in the order that the edit operations are called upon.

Exclusion/inclusion is tested before applying any id mapping.

Usage
-----

Example::

   python fasta2fasta.py --method=translate < in.fasta > out.fasta

Type::

   python fasta2fasta.py --help

for command line help.

Command line options
---------------------

'''
import sys
import string
import re
import random
from future.moves.itertools import zip_longest

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Genomics as Genomics
import CGAT.FastaIterator as FastaIterator
import CGAT.Masker as Masker


def getCodons(sequence, gap_chars="-."):
    """get codons in sequence."""
    codons, codon = [], []
    full_codons, full_codon = [], []
    for x in range(len(sequence)):
        c = sequence[x]
        full_codon.append(c)
        if c not in gap_chars:
            codon.append(c)
        if len(codon) % 3 == 0:
            codons.append(codon)
            full_codons.append(full_codon)
            codon = []
            full_codon = []

    if full_codon:
        full_codons.append(full_codon)
        codons.append(codon)

    return full_codons, codons


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-m", "--method", dest="methods", type="choice", action="append",
        choices=("translate",
                 "translate-to-stop",
                 "truncate-at-stop",
                 "back-translate",
                 "mark-codons",
                 "apply-map",
                 "build-map",
                 "pseudo-codons",
                 "filter",
                 "interleaved-codons",
                 "map-codons",
                 "remove-gaps",
                 "mask-seg",
                 "mask-bias",
                 "mask-codons",
                 "mask-incomplete-codons",
                 "mask-stops",
                 "mask-soft",
                 "remove-stops",
                 "upper",
                 "lower",
                 "reverse-complement",
                 "sample",
                 "shuffle"),
        help="method to apply to sequences.")

    parser.add_option(
        "-p", "--parameters", dest="parameters", type="string",
        help="parameter stack for methods that require one "
        "[default=%default].")

    parser.add_option(
        "-x", "--ignore-errors", dest="ignore_errors", action="store_true",
        help="ignore errors [default = %default].")

    parser.add_option("--sample-proportion", dest="sample_proportion",
                      type="float",
                      help="sample proportion [default = %default].")

    parser.add_option(
        "--exclude-pattern", dest="exclude_pattern", type="string",
        help="exclude all sequences with ids matching pattern "
        "[default = %default].")

    parser.add_option(
        "--include-pattern", dest="include_pattern", type="string",
        help="include only sequences with ids matching pattern "
        "[default = %default].")

    parser.add_option(
        "--filter-method", dest="filter_methods", type="string",
        action="append",
        help="filtering methods to apply "
        "[default = %default].")

    parser.add_option(
        "-t", "--sequence-type", dest="type", type="choice",
        choices=("aa", "na"),
        help="sequence type (aa or na) [%default]. This option determines "
        "which characters to use for masking [default = %default].")

    parser.add_option(
        "-l", "--template-identifier", dest="template_identifier",
        type="string",
        help="template for numerical identifier [default = %default] "
        "for the operation --build-map. A %i is replaced by the position "
        "of the sequence in the file.")

    parser.set_defaults(
        methods=[],
        parameters="",
        type="na",
        aa_mask_chars="xX",
        aa_mask_char="x",
        na_mask_chars="nN",
        na_mask_char="n",
        gap_chars="-.",
        gap_char="-",
        template_identifier="ID%06i",
        ignore_errors=False,
        exclude_pattern=None,
        include_pattern=None,
        sample_proportion=None,
        filter_methods=[],
    )

    (options, args) = E.Start(parser)
    options.parameters = options.parameters.split(",")

    rx_include, rx_exclude = None, None
    if options.include_pattern:
        rx_include = re.compile(options.include_pattern)
    if options.exclude_pattern:
        rx_exclude = re.compile(options.exclude_pattern)

    iterator = FastaIterator.FastaIterator(options.stdin)

    nseq = 0

    map_seq2nid = {}

    if "apply-map" in options.methods:
        map_seq2nid = IOTools.ReadMap(open(options.parameters[0], "r"))
        del options.parameters[0]

    if options.type == "na":
        mask_chars = options.na_mask_chars
        mask_char = options.na_mask_char
    else:
        mask_chars = options.aa_mask_chars
        mask_char = options.aa_mask_char

    if "map-codons" in options.methods:
        map_codon2code = IOTools.ReadMap(open(options.parameters[0], "r"))
        del options.parameters[0]

    if "mask-soft" in options.methods:
        f = options.parameters[0]
        del options.parameters[0]
        hard_masked_iterator = FastaIterator.FastaIterator(open(f, "r"))

    if "mask-codons" in options.methods or "back-translate" in options.methods:

        # open a second stream to read sequences from
        f = options.parameters[0]
        del options.parameters[0]

        other_iterator = FastaIterator.FastaIterator(open(f, "r"))

    ninput, noutput, nerrors, nskipped = 0, 0, 0, 0

    if "sample" in options.methods:
        if not options.sample_proportion:
            raise ValueError("specify a sample proportion")
        sample_proportion = options.sample_proportion
    else:
        sample_proportion = None

    filter_min_sequence_length = None
    filter_max_sequence_length = None
    filter_id_list = None
    for f in options.filter_methods:
        if f.startswith("min-length"):
            filter_min_sequence_length = int(f.split("=")[1])
        elif f.startswith("max-length"):
            filter_max_sequence_length = int(f.split("=")[1])
        elif f.startswith("id-file"):
            filter_id_list = [line[:-1] for line in IOTools.openFile(f.split("=")[1])]

    def raiseIfNotCodon(l, title):
        '''raise ValueError if sequence length l is not divisible by
        3'''

        if l % 3 != 0:
            raise ValueError(
                "length of sequence %s not divisible by 3" % (title))

    while 1:
        try:
            cur_record = next(iterator)
        except StopIteration:
            break

        if cur_record is None:
            break
        nseq += 1
        ninput += 1

        sequence = re.sub(" ", "", cur_record.sequence)
        l = len(sequence)

        if rx_include and not rx_include.search(cur_record.title):
            nskipped += 1
            continue

        if rx_exclude and rx_exclude.search(cur_record.title):
            nskipped += 1
            continue

        if sample_proportion:
            if random.random() > sample_proportion:
                continue

        if not (filter_id_list is None or cur_record.title in filter_id_list):
            nskipped += 1
            continue

        for method in options.methods:

            if method == "translate":
                # translate such that gaps are preserved
                seq = []

                ls = len(re.sub('[%s]' % options.gap_chars, sequence, ""))

                if ls % 3 != 0:
                    msg = "length of sequence %s (%i) not divisible by 3" % (
                        cur_record.title, ls)
                    nerrors += 1
                    if options.ignore_errors:
                        E.warn(msg)
                        continue
                    else:
                        raise ValueError(msg)

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:
                    aa = Genomics.MapCodon2AA(codon)
                    seq.append(aa)

                sequence = "".join(seq)

            elif method == "back-translate":
                # translate from an amino acid alignment to codon alignment
                seq = []

                try:
                    other_record = next(other_iterator)
                except StopIteration:
                    raise ValueError("run out of sequences")

                if cur_record.title != other_record.title:
                    raise "sequence titles don't match: %s %s" % (
                        cur_record.title, other_record.title)

                other_sequence = re.sub(
                    "[ %s]" % options.gap_chars, "", other_record.sequence)

                if len(other_sequence) % 3 != 0:
                    raise ValueError(
                        "length of sequence %s not divisible by 3" %
                        (other_record.title))

                r = re.sub("[%s]" % options.gap_chars, "", sequence)
                if len(other_sequence) != len(r) * 3:
                    raise ValueError(
                        "length of sequences do not match: %i vs %i" %
                        (len(other_sequence), len(r)))

                x = 0
                for aa in sequence:
                    if aa in options.gap_chars:
                        c = options.gap_char * 3
                    else:
                        c = other_sequence[x:x + 3]
                        x += 3
                    seq.append(c)

                sequence = "".join(seq)

            elif method == "pseudo-codons":
                raiseIfNotCodon(l, cur_record.title)
                seq = []

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:

                    aa = Genomics.MapCodon2AA(codon)
                    seq.append(aa)

                sequence = "   ".join(seq)

            elif method == "reverse-complement":
                sequence = sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

            elif method in ("mask-stops", "remove-stops"):
                c = []
                codon = []
                new_sequence = []

                if method == "mask-stops":
                    char = options.na_mask_char
                elif method == "remove-stops":
                    char = options.gap_char

                for x in sequence:

                    if x not in options.gap_chars:
                        codon.append(x.upper())

                    c.append(x)

                    if len(codon) == 3:
                        codon = "".join(codon).upper()
                        # mask all non-gaps
                        if Genomics.IsStopCodon(codon):

                            for x in c:
                                if x in options.gap_chars:
                                    new_sequence.append(x)
                                else:
                                    new_sequence.append(char)
                        else:
                            new_sequence += c

                        c = []
                        codon = []

                new_sequence += c

                sequence = "".join(new_sequence)

            elif method == "mask-soft":
                # Get next hard masked record and extract sequence and length
                try:
                    cur_hm_record = next(hard_masked_iterator)
                except StopIteration:
                    break
                hm_sequence = re.sub(" ", "", cur_hm_record.sequence)
                lhm = len(hm_sequence)
                new_sequence = []

                # Check lengths of unmasked and soft masked sequences the same
                if l != lhm:
                    raise ValueError(
                        "length of unmasked and hard masked sequences not "
                        "identical for record %s" %
                        (cur_record.title))

                # Check if hard masked seq contains repeat (N), if so replace N
                # with lowercase sequence from unmasked version
                if sequence == hm_sequence:
                    pass
                else:
                    for x, y in zip_longest(sequence, hm_sequence):
                        if y == "N":
                            new_sequence += x.lower()
                        else:
                            new_sequence += x.upper()
                sequence = "".join(new_sequence)

            elif method == "map-codons":
                raiseIfNotCodon(l, cur_record.title)
                seq = []

                for codon in (sequence[x:x + 3].upper()
                              for x in range(0, l, 3)):

                    if codon not in map_codon2code:
                        aa = "X"
                    else:
                        aa = map_codon2code[codon]
                    seq.append(aa)

                sequence = "".join(seq)

            elif method == "interleaved-codons":
                raiseIfNotCodon(l, cur_record.title)
                seq = []

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:

                    aa = Genomics.MapCodon2AA(codon)
                    seq.append("%s:%s" % (aa, codon))

                sequence = " ".join(seq)

            elif method == "translate-to-stop":
                seq = []

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:

                    if Genomics.IsStopCodon(codon):
                        break

                    aa = Genomics.MapCodon2AA(codon)
                    seq.append(aa)

                sequence = "".join(seq)

            elif method == "truncate-at-stop":
                seq = []

                for codon in [sequence[x:x + 3] for x in range(0, l, 3)]:

                    if Genomics.IsStopCodon(codon):
                        break
                    seq.append(codon)

                sequence = "".join(seq)

            elif method == "remove-gaps":

                seq = []
                for s in sequence:
                    if s in options.gap_chars:
                        continue
                    seq.append(s)

                sequence = "".join(seq)

            elif method == "upper":
                sequence = sequence.upper()

            elif method == "lower":
                sequence = sequence.lower()

            elif method == "mark-codons":
                raiseIfNotCodon(l, cur_record.title)
                seq = []

                sequence = " ".join([sequence[x:x + 3]
                                     for x in range(0, l, 3)])

            elif method == "apply-map":
                id = re.match("^(\S+)", cur_record.title).groups()[0]
                if id in map_seq2nid:
                    rest = cur_record.title[len(id):]
                    cur_record.title = map_seq2nid[id] + rest

            elif method == "build-map":
                # build a map of identifiers
                id = re.match("^(\S+)", cur_record.title).groups()[0]
                new_id = options.template_identifier % nseq
                if id in map_seq2nid:
                    raise "duplicate fasta entries - can't map those: %s" % id
                map_seq2nid[id] = new_id
                cur_record.title = new_id

            elif method == "mask-bias":
                masker = Masker.MaskerBias()
                sequence = masker(sequence)

            elif method == "mask-seg":
                masker = Masker.MaskerSeg()
                sequence = masker(sequence)

            elif method == "shuffle":
                s = list(sequence)
                random.shuffle(s)
                sequence = "".join(s)

            elif method == "mask-incomplete-codons":
                seq = list(sequence)
                for x in range(0, l, 3):
                    nm = len([x for x in seq[x:x + 3] if x in mask_chars])
                    if 0 < nm < 3:
                        seq[x:x + 3] = [mask_char] * 3
                sequence = "".join(seq)

            elif method == "mask-codons":
                # mask codons based on amino acids given as reference
                # sequences.
                other_record = next(other_iterator)

                if other_record is None:
                    raise ValueError("run out of sequences.")

                if cur_record.title != other_record.title:
                    raise ValueError(
                        "sequence titles don't match: %s %s" %
                        (cur_record.title, other_record.title))

                other_sequence = re.sub(" ", "", other_record.sequence)

                if len(other_sequence) * 3 != len(sequence):
                    raise ValueError(
                        "sequences for %s don't have matching lengths %i - %i" %
                        (cur_record.title, len(other_sequence) * 3,
                         len(sequence)))

                seq = list(sequence)
                c = 0
                for x in other_sequence:
                    if x in options.aa_mask_chars:
                        if x.isupper():
                            seq[c:c + 3] = [options.na_mask_char.upper()] * 3
                        else:
                            seq[c:c + 3] = [options.na_mask_char.lower()] * 3
                    c += 3

                sequence = "".join(seq)

        l = len(sequence)
        if filter_min_sequence_length is not None and \
           l < filter_min_sequence_length:
            nskipped += 1

        if filter_max_sequence_length is not None and \
           l > filter_max_sequence_length:
            nskipped += 1
            continue

        options.stdout.write(">%s\n%s\n" % (cur_record.title, sequence))
        noutput += 1

    if "build-map" in options.methods:
        p = options.parameters[0]
        if p:
            outfile = IOTools.openFile(p, "w")
        else:
            outfile = options.stdout

        outfile.write("old\tnew\n")
        for old_id, new_id in list(map_seq2nid.items()):
            outfile.write("%s\t%s\n" % (old_id, new_id))
        if p:
            outfile.close()

    E.info("ninput=%i, noutput=%i, nskipped=%i, nerrors=%i" %
           (ninput, noutput, nskipped, nerrors))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

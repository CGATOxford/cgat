'''diff_fasta.py - compare contents of two fasta files
===================================================

:Tags: Genomics Sequences FASTA Comparison

Purpose
-------

This script takes two sets of fasta sequences and matches the
identifiers. It then compares the sequences with the same identifiers
and, depending on the output options selected, outputs

   * which sequences are missing
   * which sequences are identical
   * which sequences are prefixes/suffixes of each other

An explanatory field is appended to output sequence identifiers.
An explanation of the different field values is provided in the log.

Options
-------
-s, --correct-gap-shift
   This option will correct shifts in alignment gaps between
   two sequences being compared

-1, --pattern1
   regular expression pattern to extract identifier from in sequence 1

-2, --pattern2
   regular expression pattern to extract identifier from in sequence 2

Depending on the option ``--output-section`` the following are output:

  diff
     identifiers of sequences that are different

  seqdiff
     identifiers of sequences that are different plus sequence

  missed
    identifiers of seqences that are missing from one set or the other

This script is of specialized interest and has been used
in the past to check if ENSEMBL gene models had been
correctly mapped into a database schema.

Usage
-----

Example::

   cat a.fasta | head

   >ENSACAP00000004922
   MRSRNQGGESSSSGKFSKSKPIINTGENQNLQEDAKKKNKSSRKEE ...
   >ENSACAP00000005213
   EEEEDESNNSYLYQPLNQDPDQGPAAVEETAPSTEPALDINERLQA ...
   >ENSACAP00000018122
   LIRSSSMFHIMKHGHYISRFGSKPGLKCIGMHENGIIFNNNPALWK ...

   python diff_fasta.py --output-section=missed --output-section=seqdiff a.fasta b.fasta

   cat diff.out

   # Legend:
   # seqs1:          number of sequences in set 1
   # seqs2:          number of sequences in set 2
   # same:           number of identical sequences
   # diff:           number of sequences with differences
   # nmissed1:       sequences in set 1 that are not found in set 2
   # nmissed2:       sequences in set 2 that are not found in set 1
   # Type of sequence differences
   # first:          only the first residue is different
   # last:           only the last residue is different
   # prefix:         one sequence is prefix of the other
   # selenocysteine: difference due to selenocysteines
   # masked:         difference due to masked residues
   # fixed:          fixed differences
   # other:          other differences


Type::

   python diff_fasta.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import CGAT.FastaIterator as FastaIterator
import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def MapIdentifiers(seqs, pattern):

    rx = re.compile(pattern)

    for k, s in list(seqs.items()):
        try:
            nk = rx.search(k).groups()[0]
        except AttributeError:
            raise ValueError(
                "identifier can not be parsed from '%s' "
                "pattern='%s'" % (k, pattern))

        del seqs[k]
        seqs[nk] = s


def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-s", "--correct-gap-shift", dest="correct_shift",
        action="store_true",
        help="correct gap length shifts in alignments. "
        "Requires alignlib_lite.py [%default]")

    parser.add_option(
        "-1", "--pattern1", dest="pattern1", type="string",
        help="pattern to extract identifier from in identifiers1. "
        "[%default]")

    parser.add_option(
        "-2", "--pattern2", dest="pattern2", type="string",
        help="pattern to extract identifier from in identifiers2. "
        "[%default]")

    parser.add_option(
        "-o", "--output-section", dest="output", type="choice",
        action="append",
        choices=("diff", "missed", "seqdiff"),
        help="what to output [%default]")

    parser.set_defaults(correct_shift=False,
                        pattern1="(\S+)",
                        pattern2="(\S+)",
                        output=[])

    (options, args) = E.Start(parser)

    if len(args) != 2:
        raise ValueError("two files needed to compare.")

    if options.correct_shift:
        try:
            import alignlib_lite
        except ImportError:
            raise ImportError(
                "option --correct-shift requires alignlib_lite.py_ "
                "but alignlib not found")

    seqs1 = dict([
        (x.title, x.sequence) for x in FastaIterator.iterate(
            IOTools.openFile(args[0], "r"))])
    seqs2 = dict([
        (x.title, x.sequence) for x in FastaIterator.iterate(
            IOTools.openFile(args[1], "r"))])

    if not seqs1:
        raise ValueError("first file %s is empty." % (args[0]))
    if not seqs2:
        raise ValueError("second file %s is empty." % (args[1]))

    MapIdentifiers(seqs1, options.pattern1)
    MapIdentifiers(seqs2, options.pattern2)

    nsame = 0
    nmissed1 = 0
    nmissed2 = 0
    ndiff = 0
    ndiff_first = 0
    ndiff_last = 0
    ndiff_prefix = 0
    ndiff_selenocysteine = 0
    ndiff_masked = 0
    nfixed = 0
    found2 = {}

    write_missed1 = "missed" in options.output
    write_missed2 = "missed" in options.output
    write_seqdiff = "seqdiff" in options.output
    write_diff = "diff" in options.output or write_seqdiff

    for k in sorted(seqs1):
        if k not in seqs2:
            nmissed1 += 1
            if write_missed1:
                options.stdout.write("---- %s ---- %s\n" % (k, "missed1"))
            continue

        found2[k] = 1

        s1 = seqs1[k].upper()
        s2 = seqs2[k].upper()
        m = min(len(s1), len(s2))

        if s1 == s2:
            nsame += 1
        else:
            status = "other"

            ndiff += 1

            if s1[1:] == s2[1:]:
                ndiff_first += 1
                status = "first"
            elif s1[:m] == s2[:m]:
                ndiff_prefix += 1
                status = "prefix"
            elif s1[:-1] == s2[:-1]:
                ndiff_last += 1
                status = "last"
            else:
                if len(s1) == len(s2):
                    # get all differences: the first and last residues
                    # can be different for peptide sequences when
                    # comparing my translations with ensembl peptides.
                    differences = []
                    for x in range(1, len(s1) - 1):
                        if s1[x] != s2[x]:
                            differences.append((s1[x], s2[x]))

                    l = len(differences)
                    # check for Selenocysteins
                    if len([x for x in differences if x[0] == "U" or x[1] == "U"]) == l:
                        ndiff_selenocysteine += 1
                        status = "selenocysteine"

                    # check for masked residues
                    elif len([x for x in differences if x[0] in "NX" or x[1] in "NX"]) == l:
                        ndiff_masked += 1
                        status = "masked"

            # correct for different gap lengths
            if options.correct_shift:

                map_a2b = alignlib_lite.py_makeAlignmentVector()

                a, b = 0, 0
                keep = False

                x = 0
                while x < m and not (a == len(s1) and b == len(s2)):
                    try:
                        if s1[a] != s2[b]:
                            while s1[a] == "N" and s2[b] != "N":
                                a += 1
                            while s1[a] != "N" and s2[b] == "N":
                                b += 1

                            if s1[a] != s2[b]:
                                break
                    except IndexError:
                        print("# index error for %s: x=%i, a=%i, b=%i, l1=%i, l2=%i" % (k, x, a, b, len(s1), len(s2)))
                        break

                    a += 1
                    b += 1
                    map_a2b.addPairExplicit(a, b, 0.0)
                    # check if we have reached the end:
                else:
                    keep = True
                    nfixed += 1
                    f = alignlib_lite.py_AlignmentFormatEmissions(map_a2b)
                    print("fix\t%s\t%s" % (k, str(f)))

                if not keep:
                    print("# warning: not fixable: %s" % k)

            if write_diff:
                options.stdout.write("---- %s ---- %s\n" % (k, status))

            if write_seqdiff:
                options.stdout.write("< %s\n> %s\n" % (seqs1[k], seqs2[k]))

    for k in sorted(list(seqs2.keys())):
        if k not in found2:
            nmissed2 += 1
            if write_missed2:
                options.stdout.write("---- %s ---- %s\n" % (k, "missed2"))

    options.stdlog.write("""# Legend:
""")

    E.info("seqs1=%i, seqs2=%i, same=%i, ndiff=%i, nmissed1=%i, nmissed2=%i" %
           (len(seqs1), len(seqs2), nsame, ndiff, nmissed1, nmissed2))

    E.info(
        "ndiff=%i: first=%i, last=%i, prefix=%i, selenocysteine=%i, masked=%i, fixed=%i, other=%i" %
        (ndiff, ndiff_first, ndiff_last, ndiff_prefix,
         ndiff_selenocysteine, ndiff_masked, nfixed,
         ndiff - ndiff_first - ndiff_last - ndiff_prefix -
         ndiff_selenocysteine - ndiff_masked - nfixed))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

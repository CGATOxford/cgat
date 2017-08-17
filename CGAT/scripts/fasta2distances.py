'''
fasta2distances.py - analyze pairs of sequences
===============================================

:Tags: Python

Purpose
-------

This script computes various distances between sequences.

Usage
-----

Example::

   python fasta2distances.py --help

Type::

   python fasta2distances.py --help

for command line help.

Command line options
--------------------

'''
import sys
import math
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.FastaIterator as FastaIterator


def FilterAlignedPairForPositions(seq1, seq2, method):
    """given the method, return a set of aligned sequences
    only containing certain positions.

    Available filters:
    all:        do nothing.
    codon1,codon2,codon3: return 1st, 2nd, 3rd codon positions only.
    d4: only changes within fourfold-degenerate sites
    """

    l1 = len(seq1)
    l2 = len(seq2)

    if method == "all":
        return seq1, seq2
    elif method == "codon1":
        return ("".join([seq1[x] for x in range(0, l1, 3)]),
                "".join([seq2[x] for x in range(0, l2, 3)]))
    elif method == "codon2":
        return ("".join([seq1[x] for x in range(1, l1, 3)]),
                "".join([seq2[x] for x in range(1, l2, 3)]))
    elif method == "codon3":
        return ("".join([seq1[x] for x in range(2, l1, 3)]),
                "".join([seq2[x] for x in range(2, l2, 3)]))
    elif method == "d4":
        s1 = []
        s2 = []
        for x in range(0, l1, 3):
            codon1 = seq1[x:x + 3]
            codon2 = seq2[x:x + 3]
            try:
                aa1, deg11, deg12, deg13 = Genomics.GetDegeneracy(codon1)
                aa2, deg11, deg22, deg23 = Genomics.GetDegeneracy(codon2)
            except KeyError:
                continue
            if aa1 == aa2 and deg13 == 4 and deg23 == 4:
                s1.append(codon1[2])
                s2.append(codon2[2])
        return "".join(s1), "".join(s2)


def CalculateDistanceJC69(info, do_gamma=False, alpha=None):
    """return Jukes-Cantor distance.
    """
    try:
        p = float(info.mNDifferent) / info.mNAligned

        if do_gamma:
            # not done yet
            distance = 0.75 * alpha * (pow(1 - 4 * p / 3, -1 / alpha) - 1)
            variance = p * (1 - p) / (pow(1 - 4 * p / 3, -2 / (alpha + 1)) * L)
        else:
            distance = -0.75 * math.log(1.0 - 4.0 * p / 3.0)
            variance = p * (1.0 - p) / \
                (math.pow(1.0 - 4.0 * p / 3, 2.0) * info.mNAligned)
    except:
        raise ValueError

    return distance, variance


def CalculateDistanceT92(info):
    """
    P,Q: transitions, transversions frequencies
    q: G+C content

    d = -2q(1 - q)loge(1 - P/[2q(1 - q)] - Q) -[1 -2q(1 -q)]loge(1 - 2Q)/2,(4.18)
    V(d) = [c12P + c32Q - (c1P + c3Q)2]/n,(4.19)
    where c1 = 1/(1 - P/[2q(1 - q)] - Q), c2 = 1/(1 - 2Q), c3 = 2q(1 - q)(c1 - c2) + c2, and q is the G+C content

    Note: result is undefined if
        the number of transversions is >= 0.5
        the G+C content is 0

    raises ValueErrors for undefined results
    """
    gc = info.getGCContent()

    # if there are no GC or no AT pairs: result is undefined
    if gc == 0 or gc == 1:
        raise ValueError

    wg = 2.0 * gc * (1.0 - gc)

    P = float(info.mNTransitions) / info.mNAligned
    Q = float(info.mNTransversions) / info.mNAligned

    a1 = 1.0 - P / wg - Q
    if a1 <= 0:
        raise ValueError

    a2 = 1.0 - 2.0 * Q
    if a2 <= 0:
        raise ValueError

    # print a1, a2, wg, gc, "p=", P, "q=", Q, str(info)

    distance = -wg * math.log(a1) - 0.5 * (1.0 - wg) * math.log(a2)

    c1 = 1 / a1
    c2 = 1 / a2
    c3 = wg * (c1 - c2) + c2

    variance = (
        c1 * c1 * P + c3 * c3 * Q - math.pow(c1 * P + c3 * Q, 2.0)) / info.mNAligned

    return distance, variance


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: fasta2distances.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("--filters", dest="filters", type="string",
                      help="Filters to use for filtering sequences [all|codon1|codon2|codon3|d4].")
    parser.add_option("--fields", dest="fields", type="string",
                      help="Fields to output [aligned|nunaligned1|nunaligned2|identical|transitions|transversions|jc69|t92].")

    parser.set_defaults(
        filename_map=None,
        filters="all,codon1,codon2,codon3,d4",
        gap_char="-",
        fields="aligned,unaligned1,unaligned2,identical,transitions,transversions,jc69,t92",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    options.filters = options.filters.split(",")
    options.fields = options.fields.split(",")

    iterator = FastaIterator.FastaIterator(options.stdin)

    headers = ["id1", "id2"]
    for f in options.filters:
        headers += list(["%s_%s" % (f, x) for x in options.fields])

    options.stdout.write("\t".join(headers) + "\n")

    while 1:
        try:
            cur_record = next(iterator)
            if cur_record is None:
                break
            first_record = cur_record
            cur_record = next(iterator)
            if cur_record is None:
                break
            second_record = cur_record

        except StopIteration:
            break

        if len(first_record.sequence) != len(second_record.sequence):
            raise "sequences %s and %s of unequal length" % (
                first_record.title, second_record.title)

        if len(first_record.sequence) % 3 != 0:
            raise "sequence %s not multiple of 3" % first_record.title

        # old: Bio.Alphabet.IUPAC.extended_dna.letters
        alphabet = "ACGT" + options.gap_char

        result = []
        for f in options.filters:

            s1, s2 = FilterAlignedPairForPositions(first_record.sequence,
                                                   second_record.sequence,
                                                   f)

            info = Genomics.CalculatePairIndices(s1, s2, options.gap_char)

            for field in options.fields:

                if field == "aligned":
                    c = "%i" % info.mNAligned
                elif field == "unaligned1":
                    c = "%i" % info.mNUnaligned1
                elif field == "unaligned2":
                    c = "%i" % info.mNUnaligned2
                elif field == "transversions":
                    c = "%i" % info.mNTransversions
                elif field == "transitions":
                    c = "%i" % info.mNTransitions
                elif field == "identical":
                    c = "%i" % info.mNIdentical
                elif field == "jc69":
                    try:
                        c = "%6.4f" % CalculateDistanceJC69(info)[0]
                    except ValueError:
                        c = "nan"
                elif field == "t92":
                    try:
                        c = "%6.4f" % CalculateDistanceT92(info)[0]
                    except ValueError:
                        c = "nan"
                else:
                    raise "Unknown field %s" % field

                result.append(c)

        options.stdout.write("%s\t%s\t%s\n" % (first_record.title,
                                               second_record.title,
                                               "\t".join(result)))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

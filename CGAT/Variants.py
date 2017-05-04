'''
Variants.py - 
======================================================

:Tags: Python

Code
----

'''
import collections

from CGAT import Experiment as E
from CGAT import Genomics as Genomics
from CGAT import NCL as ncl

Variant = collections.namedtuple("Variant",
                                 "pos, reference, genotype")
ExtendedVariant = collections.namedtuple(
    "Variant",
    "start, end, reference, action, has_wildtype, variantseqs")


def updateVariants(variants, lcontig, strand, phased=True):
    '''update variants such that they use same coordinate
    system (and strand) as the transcript

    fixes 1-ness of variants
    '''

    new_variants = []
    is_positive = Genomics.IsPositiveStrand(strand)

    for variant in variants:

        pos = variant.pos
        genotype = bytes(variant.genotype)
        reference = bytes(variant.reference)

        # fix 1-ness of variants
        # pos -= 1

        if len(genotype) == 1:
            variantseqs = list(Genomics.decodeGenotype(genotype))
            has_wildtype = reference in variantseqs
            action = "="
            start, end = pos, pos + 1
        else:

            variantseqs = [x[1:] for x in genotype.split("/")]
            lvariant = max([len(x) for x in variantseqs])
            if not phased:
                variantseqs = [x for x in variantseqs if x]
            has_wildtype = "*" in genotype

            if "+" in genotype and "-" in genotype:
                # both insertion and deletion at position
                # the range is given by the deletion
                # see below for explanations
                if genotype.startswith("+"):
                    action = ">"
                    variantseqs[1] += "-" * (lvariant - len(variantseqs[1]))
                else:
                    action = "<"
                    variantseqs[0] += "-" * (lvariant - len(variantseqs[0]))

                start, end = pos + 1, pos + lvariant + 1

            elif "-" in genotype:
                action = "-"
                # samtools: deletions are after the base denoted by snp.position
                #   * <- deletion at 1
                # 0 1 2 3 4 5 6
                #     - -
                # 6 5 4 3 2 1 0
                # deletion of 2+3 = (2,4)
                # on reverse: (7-4, 7-2) = (3,5)
                start, end = pos + 1, pos + lvariant + 1

                # deletions of unequal length are filled up with "-"
                # This is necessary to deal with negative strands:
                # -at/-atg on the positive strand deletes a t [g]
                # -at/-atg on the negative strand deletes [g] t a
                variantseqs = [x + "-" * (lvariant - len(x))
                               for x in variantseqs]

            elif "+" in genotype:
                action = "+"
                # indels are after the base denoted by position
                # as region use both flanking base so that negative strand
                # coordinates work
                # insertion between position 2 and 3
                #     * <- insection at pos 2
                # 0 1 2i3 4
                # 4 3 2i1 0
                # is insertion between 1 and 2 in reverse
                # including both flanking residues makes it work:
                # (2,3) = (5-3,5-2) = (2,3)
                # but:
                # (2,4) = (5-4,5-2) = (1,3)
                start, end = pos, pos + 2

        # revert strand
        if not is_positive:
            reference = Genomics.complement(reference)
            variantseqs = [Genomics.complement(x.upper()) for x in variantseqs]
            start, end = lcontig - end, lcontig - start

        new_variants.append(ExtendedVariant._make((
            start, end, reference.upper(), action, has_wildtype, variantseqs)))

    return new_variants


def mergeVariants(variants):
    '''merge overlapping variants.

    Overlapping variants occur if there are two deletions
    at the same location:

        WT      ACTG  
        Allele1 -CT-   
        Allele2 ----

    This will be encoded by samtools as (0-based coordinates)::

        0 * -A/ACTG
        3 * -G/-G

    This upsets the re-constitution algoritm.

    This method separates these two variants into two non-overlapping
    variants making use of variable length deletions.

        0 * -A/-A
        1 * ---G/-CTG

    Another case:

        WT      ACTG  
        Allele1 ACT-   
        Allele2 ----

    This will be encoded by samtools as (0-based coordinates)::

        0 * */-ACTG
        3 * -G/*

    This method separates these two as::

        0 * */-ACT
        3 * -G/-G

    '''

    if len(variants) == 0:
        return []

    # sorts by start and then end
    variants.sort()
    merged_variants = []

    def _add(offset, dest, src):
        for x, c in enumerate(src):
            dest[x + offset] = c

    def _split(seq0, seq1):
        # split
        was_0, was_1 = seq0[0] == "-", seq1[0] == "-"
        for x, cc in enumerate(zip(seq0, seq1)):
            is_0, is_1 = cc[0] == "-", cc[1] == "-"
            # yield all changes
            if (is_0 ^ was_0) or (is_1 ^ was_1):
                yield x, was_0, was_1
            was_0, was_1 = is_0, is_1

        yield x + 1, was_0, was_1

    last = variants[0]
    for this in variants[1:]:

        if this.start < last.end and \
                this.action == "-" and \
                last.action == "-":

            E.warn("merging overlapping deletions: %s and %s" %
                   (str(last), str(this)))

            mend = max(last.end, this.end)
            mstart = min(this.start, last.start)
            l = mend - mstart

            seq0 = list("-" * l)
            seq1 = list("-" * l)

            _add(last.start - mstart, seq0, last.variantseqs[0])
            _add(last.start - mstart, seq1, last.variantseqs[1])
            _add(this.start - mstart, seq0, this.variantseqs[0])
            _add(this.start - mstart, seq1, this.variantseqs[1])

            last_x = 0
            n = []
            for x, was_0, was_1 in _split(seq0, seq1):
                if last_x == x:
                    continue

                this = ExtendedVariant._make((mstart + last_x,
                                              mstart + x,
                                              "*",
                                              last.action,
                                              was_0 ^ was_1,
                                              ["".join(seq0[last_x:x]),
                                               "".join(seq1[last_x:x])],
                                              ))
                n.append(this)
                last_x = x

            E.warn("overlapping deletions merged in %i blocks as: %s" %
                   (len(n), list(map(str, n))))
            merged_variants.extend(n[:-1])
            this = n[-1]
        else:
            merged_variants.append(last)

        last = this

    merged_variants.append(last)

    return merged_variants


def indexVariants(variants):
    '''build index of variants for ranged retrieval.'''
    index = ncl.NCL()

    for start, end, reference, action, has_wildtype, variantseqs in variants:
        index.add(start, end, (reference, action, has_wildtype, variantseqs))

    return index


def buildAlleles(sequence, variants, reference_start=0, phased=True):
    '''build alleles for ``sequence`` adding ``variants``.

    Variants are assumed to be in 0-based coordinates on the same strand as the sequence. 
    ``reference_start`` is the position of the first base of ``sequence``. Set to 0, if
    the positions in ``variants`` are relative to ``sequence``.
    '''

    def _delete(allele, del_start, del_end, variant, sequence, startoffset, endoffset, feature_start, feature_end):
        '''little helper: update ``allele`` with a deletion ``del_start:del_end``.
        '''

        # truncate variant according to the feature
        variant = variant[startoffset:len(variant) - endoffset]

        n = variant.count("-")
        if n:
            if variant.startswith("-"):
                del_start += n
                variant = variant[n:]
            else:
                del_end -= n
                variant = variant[:-n]

        # due to gaps, the variant is not actually within the feauture
        if del_start >= del_end:
            return

        refseq = sequence[del_start:del_end].upper()

        assert refseq == variant, \
            'reference base mismatch at deletion: expected %s %s %s, got %s[%i:%i] at feature=%i-%i, variant=%i-%i, relative=%i-%i, del=%i-%i, action=%s' % \
            (sequence[del_start - 10:del_start],
             refseq,
             sequence[del_end:del_end + 10],
             variant, startoffset, len(variant) - endoffset,
             feature_start, feature_end,
             var_start, var_end,
             rel_start, rel_end,
             del_start, del_end,
             action)

        l = del_end - del_start

        # assert len("".join(allele[del_start:del_end])) == l, \
        #     "deletion conflicts with other indels: " \
        #     "got %s[%i:%i] (ref=%s, allele=%s) at feature=%i-%i, variant=%i-%i, relative=%i-%i, del=%i-%i, action=%s" % \
        #     (variant, startoffset, len(variant)-endoffset,
        #      refseq, str(allele[del_start:del_end]),
        #      feature_start, feature_end,
        #      var_start, var_end,
        #      rel_start, rel_end,
        #      del_start, del_end,
        #      action)

        allele[del_start:del_end] = [""] * l

    allele1 = list(sequence.lower())
    allele2 = list(sequence.lower())

    if reference_start is None:
        feature_start = 0
    else:
        feature_start = reference_start

    feature_end = feature_start + len(sequence)

    # main loop: insert variants into allele sequences
    for var_start, var_end, reference, action, has_wildtype, variantseqs in variants:

        # skip variants that are out-of-range
        if var_end <= feature_start or var_start >= feature_end:
            continue

        is_homozygous = len(variantseqs) == 1 and not has_wildtype

        rel_start, rel_end = var_start - feature_start, var_end - feature_start
        startoffset = max(0, feature_start - var_start)
        endoffset = max(0, var_end - feature_end)
        pruned_start, pruned_end = max(
            0, rel_start), min(len(sequence), rel_end)

        if action == "=":

            if E.global_options.loglevel >= 10:
                E.debug("adding SNP at postition %i: reference=%s variants=%s" %
                        (var_start, reference, variantseqs))

            if allele1[rel_start] == "" or allele2[rel_start] == "":
                # these can be cases, where a base is deleted in one allele,
                # but recorded as a homozygous substitution in another allele.
                E.warn("substitution conflicts with a deletion - ignored: %s" %
                       str((var_start, var_end, reference, action, has_wildtype, variantseqs)))
                continue

            assert rel_start >= 0
            assert sequence[rel_start].upper() == reference, \
                'reference base mismatch: expected %s %s %s, got %s at feature=%i-%i, variant=%i-%i, relative=%i-%i, pruned=%i-%i, action=%s' % \
                (sequence[rel_start - 10:rel_start],
                 sequence[rel_start].upper(),
                 sequence[rel_start + 1:rel_start + 10],
                 reference,
                 feature_start, feature_end,
                 var_start, var_end,
                 rel_start, rel_end,
                 pruned_start, pruned_end,
                 action)

            if phased:
                allele1[rel_start] = variantseqs[0] + allele1[rel_start][1:]
                allele2[rel_start] = variantseqs[1] + allele2[rel_start][1:]
            elif is_homozygous:
                allele1[rel_start] = variantseqs[0] + allele1[rel_start][1:]
                allele2[rel_start] = variantseqs[0] + allele2[rel_start][1:]
            else:
                if has_wildtype:
                    if reference == variantseqs[0]:
                        allele2[rel_start] = variantseqs[
                            1] + allele2[rel_start][1:]
                    else:
                        allele2[rel_start] = variantseqs[
                            0] + allele2[rel_start][1:]
                else:
                    allele1[rel_start] = variantseqs[
                        0] + allele1[rel_start][1:]
                    allele2[rel_start] = variantseqs[
                        1] + allele2[rel_start][1:]

        elif action == "-":
            if phased:
                _delete(allele1, pruned_start, pruned_end, variantseqs[0],
                        sequence, startoffset, endoffset, feature_start, feature_end)
                _delete(allele2, pruned_start, pruned_end, variantseqs[1],
                        sequence, startoffset, endoffset, feature_start, feature_end)
            elif is_homozygous:
                _delete(allele1, pruned_start, pruned_end, variantseqs[0],
                        sequence, startoffset, endoffset, feature_start, feature_end)
                _delete(allele2, pruned_start, pruned_end, variantseqs[0],
                        sequence, startoffset, endoffset, feature_start, feature_end)
            else:
                if has_wildtype:
                    _delete(allele2, pruned_start, pruned_end, variantseqs[0],
                            sequence, startoffset, endoffset, feature_start, feature_end)
                else:
                    _delete(allele1, pruned_start, pruned_end, variantseqs[0],
                            sequence, startoffset, endoffset, feature_start, feature_end)
                    _delete(allele2, pruned_start, pruned_end, variantseqs[1],
                            sequence, startoffset, endoffset, feature_start, feature_end)

        elif action == "+":
            # ignore insertions at position -1
            if rel_start < 0:
                continue

            if phased:
                allele1[rel_start] += variantseqs[0].upper()
                allele2[rel_start] += variantseqs[1].upper()
            elif is_homozygous:
                allele1[rel_start] += variantseqs[0].upper()
                allele2[rel_start] += variantseqs[0].upper()
            else:
                if has_wildtype:
                    allele2[rel_start] += variantseqs[0].upper()
                else:
                    allele1[rel_start] += variantseqs[0].upper()
                    allele2[rel_start] += variantseqs[1].upper()

        elif action == ">":
            # indel
            if rel_start >= 0:
                allele1[rel_start] += variantseqs[0].upper()
            _delete(allele2, pruned_start, pruned_end, variantseqs[1],
                    sequence, startoffset, endoffset, feature_start, feature_end)

        elif action == "<":
            # delin
            if rel_start >= 0:
                allele2[rel_start] += variantseqs[1].upper()
            _delete(allele1, pruned_start, pruned_end, variantseqs[0],
                    sequence, startoffset, endoffset, feature_start, feature_end)

    assert len(sequence) == len(allele1)
    assert len(sequence) == len(allele2)

    return (allele1, allele2)


def buildOffsets(variants, phased=True, contig=None):
    '''collect coordinate offsets.

    This methods takes a set of variants and computes
    coordinates offsets based on indels.

    Conflicting variants will be removed. 

    Returns a list of variants, a list of removed variants and a list of offsets.
    '''

    variants.sort()
    offset0, offset1 = 0, 0
    offsets0, offsets1 = [(0, 0)], [(0, 0)]

    def addDeletion(offsets, start, end, offset, variantseq):

        for x, c in enumerate(variantseq):
            if c != "-":
                offset -= 1
                offsets.append((start + x + 1, offset))

        return offset

    last_end = 0

    removed_variants, kept_variants = [], []
    for variant in variants:

        start, end, reference, action, has_wildtype, variantseqs = variant

        if action == "=":
            kept_variants.append(variant)
            # ignore substitutions
            continue

        if last_end >= start:
            removed_variants.append(variant)
            continue

        kept_variants.append(variant)
        is_homozygous = len(variantseqs) == 1 and not has_wildtype

        last_end = end
        if action == "-":
            if phased:
                offset0 = addDeletion(
                    offsets0, start, end, offset0, variantseqs[0])
                offset1 = addDeletion(
                    offsets1, start, end, offset1, variantseqs[1])
            elif is_homozygous:
                offset0 = addDeletion(
                    offsets0, start, end, offset0, variantseqs[0])
                offset1 = addDeletion(
                    offsets1, start, end, offset1, variantseqs[0])
            else:
                if has_wildtype:
                    offset1 = addDeletion(
                        offsets1, start, end, offset1, variantseqs[0])
                else:
                    offset0 = addDeletion(
                        offsets0, start, end, offset0, variantseqs[0])
                    offset1 = addDeletion(
                        offsets1, start, end, offset1, variantseqs[1])

        elif action == "+":

            if phased:
                offset0 += len(variantseqs[0]) - variantseqs[0].count("-")
                offset1 += len(variantseqs[1]) - variantseqs[1].count("-")
            elif is_homozygous:
                offset0 += len(variantseqs[0])
                offset1 += len(variantseqs[0])
            else:
                if has_wildtype:
                    offset1 += len(variantseqs[0])
                else:
                    offset0 += len(variantseqs[0]) - variantseqs[0].count("-")
                    offset1 += len(variantseqs[1]) - variantseqs[1].count("-")
            offsets0.append((end - 1, offset0))
            offsets1.append((end - 1, offset1))

        elif action == ">":
            offset0 += len(variantseqs[0]) - variantseqs[0].count("-")
            offsets0.append((end - 1, offset0))
            offset1 = addDeletion(
                offsets1, start, end, offset1, variantseqs[1])

        elif action == "<":
            offset0 = addDeletion(
                offsets0, start, end, offset0, variantseqs[0])
            offset1 += len(variantseqs[1]) - variantseqs[1].count("-")
            offsets1.append((end - 1, offset1))

        if E.global_options.loglevel >= 8:
            print("# variant=", start, end, reference,
                  action, has_wildtype, variantseqs)
            print("# offsets:", offsets0[-1], offsets1[-1])
        # add offsets - applies to the end of the variant

    return (kept_variants, removed_variants, (offsets0, offsets1))

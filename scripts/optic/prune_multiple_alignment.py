##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
optic/prune_multiple_alignment.py -
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python optic/prune_multiple_alignment.py --help

Type::

   python optic/prune_multiple_alignment.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import warnings
import scipy
import scipy.cluster
import Bio
import Bio.Cluster
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.Mali as Mali
import CGAT.Exons as Exons


# ignore warnings from networkx/matplotlib that a display
# can not be found
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import networkx

USAGE = """python %s [OPTIONS] < exonerate_output > filtered

Prune a nucleotide multiple alignment according to a master sequence.

1. Go in codon steps through the multiple alignment according
to the master sequence.

2. Remove all columns in other sequences, that
        1. fall out of frame
        2. are incomplete codons
        3. are unaligned

Lower case characters are interpreted as unaligned characters, unless
the option --ignore-case is given.
""" % sys.argv[0]


def GetPides(mali, id1, ids, first_res, last_res):
    """calculate pid between sequence id1 and all others in ids."""

    result = []
    sequence1 = mali[id1][first_res:last_res]

    for id2 in ids:
        sequence2 = mali[id2][first_res:last_res]

        # calculate percent identity - gaps treated as separate characters.
        n, t = 0, 0
        for a, b in zip(sequence1, sequence2):
            if a == b:
                n += 1
            t += 1

        result.append(float(n) / t)
    return result


def GetInconsistentClusters(mali, nclusters, first_res, last_res,
                            inconsistencies, options):
    """cluster sequences in mali range into nclusters.

    """
    ##################################################
    # assign consistent entries

    ninc = len(inconsistencies)

    features = scipy.array([[0.0] * ninc] * ninc)

    for x in range(ninc - 1):

        id1, s1, t1, g1, q1 = inconsistencies[x]
        sequence1 = mali[id1][first_res:last_res]

        for y in range(x + 1, ninc):

            id2, s2, t2, g2, q2 = inconsistencies[x]

            if s2 == s1 and g1 == g2:
                weight = 1.0
            else:
                sequence2 = mali[filtered_identifiers[y]][first_res:last_res]

                # calculate percent identity - gaps treated as separate
                # characters.
                n, t = 0, 0
                for a, b in zip(sequence1, sequence2):
                    if a == b:
                        n += 1
                    t += 1

                weight = 1.0 - float(n) / t

                if options.loglevel >= 4:
                    options.stdlog.write(
                        "# weighting: %s\t%s\t%5.2f\n" % (t1, t2, weight))
                    options.stdlog.write(
                        "# %s\n# %s\n" % (sequence1, sequence2))

            features[x][y] = weight
            features[y][x] = weight

    result = Bio.Cluster.treecluster(distancematrix=features)
    assignments = result.cut(nclusters)

    clusters = [[] for x in range(nclusters)]

    for x in range(ninc):
        c = assignments[x]
        clusters[c].append(inconsistencies[x][0])

    return clusters

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def Assign2Clusters(mali, identifiers, clusters, first_res, last_res, options):
    """assign identifiers to clusters.

    Distance is max minimum distance.
    """

    transcripts = {}
    genomes = {}

    new_clusters = [[] for x in range(len(clusters))]
    # remove identifiers that are completely empty:
    for id in identifiers:

        best = 0
        best_pide = 0
        for c in range(len(clusters)):
            pide = min(GetPides(mali, id, clusters[c], first_res, last_res))

            if pide > best_pide:
                best, best_pide = c, pide

        new_clusters[best].append(id)

    for c in range(len(clusters)):
        new_clusters[c] += clusters[c]

    return new_clusters

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def SplitExons(mali, exons, options, separator="|", gap_char="-", masters=None):
    """split exons in mali."""

    identifiers = mali.getIdentifiers()

    ##################################################################
    # compute pairs of transcripts from the same genes
    pairs = {}
    if masters:
        masters = set(masters)
    genes = map(lambda x: (x[0], x[2]), map(
        lambda x: x.split(separator), identifiers))
    for t1 in range(len(genes) - 1):
        if masters and identifiers[t1] not in masters:
            continue
        for t2 in range(t1 + 1, len(genes)):
            if genes[t1] == genes[t2]:
                if genes[t1] not in pairs:
                    pairs[genes[t1]] = []
                pairs[genes[t1]].append((t1, t2))

    ##################################################################
    # translate residue numbers to genomic coordinates
    # add None for unaligned positions
    columns = [[] for x in range(mali.getWidth())]
    positions = []

    for id in identifiers:

        sequence = mali[id]

        ee = exons[id]
        e = 0
        genomic_pos = ee[e].mGenomeFrom
        peptide_pos = 0

        for x in range(len(sequence)):

            c = sequence[x]

            if c in options.gap_chars:
                columns[x].append(None)
            else:
                columns[x].append(genomic_pos)

                # advance on peptide and on genomic coordinates
                peptide_pos += 1
                genomic_pos += 1

                if peptide_pos >= ee[e].mPeptideTo:
                    e += 1
                    if e >= len(ee):
                        x += 1
                        break
                    genomic_pos = ee[e].mGenomeFrom

        for x in range(x, len(sequence)):
            columns[x].append(None)

    # collect inconsistent columns
    # for each inconsistent column remember the maxinum number of
    # inconsistent transcripts per gene
    inconsistent_columns = []
    for c in range(len(columns)):
        column = columns[c]
        # set with ids involved in inconsistencies
        ic = set()

        component_sizes = []
        for gene, pairs_per_gene in pairs.items():

            # The number of clusters is given by the maximum number of
            # consistent groups per gene
            consistency_graph = networkx.Graph()

            is_inconsistent = False
            for g1, g2 in pairs_per_gene:

                if not column[g1]:
                    continue
                if not column[g2]:
                    continue

                consistency_graph.add_node(g1)
                consistency_graph.add_node(g2)

                if column[g1] != column[g2]:
                    if options.loglevel >= 6:
                        options.stdlog.write(
                            "# column %i: inconsistency: %s - %i <---> %s - %i\n" %
                            (c, identifiers[g1], column[g1],
                             identifiers[g2], column[g2]))

                    ic.add(
                        (identifiers[g1],) + tuple(
                            identifiers[g1].split(options.separator)))
                    ic.add(
                        (identifiers[g2],) + tuple(
                            identifiers[g2].split(options.separator)))
                    is_inconsistent = True
                else:
                    consistency_graph.add_edge(g1, g2)

            components = networkx.connected_components(consistency_graph)

            if options.loglevel >= 6:
                if is_inconsistent:
                    options.stdlog.write(
                        "# column %i: inconsistency for gene %s - %s\n" %
                        (c, str(gene), str(components)))

            component_sizes.append(len(components))

        # count maximum transripts per gene
        if not ic:
            continue

        max_counts = max(component_sizes)
        inconsistent_columns.append((c, max_counts, ic))

    if options.loglevel >= 1:
        options.stdlog.write(
            "# found %i inconsistent columns.\n" % len(inconsistent_columns))

    if not inconsistent_columns:
        return

    if options.loglevel >= 5:
        for column, mm, ic in inconsistent_columns:
            options.stdlog.write(
                "# inconsistent column: %i\t%i: %s\n" % (column, mm, str(ic)))

    # convert to segments
    # adjacent segments are joined
    # if there is gap between segments, check, if it
    # it is the same sequence set.
    inconsistent_segments = []
    fc, max_counts, last_ic = inconsistent_columns[0]
    lc = fc

    for column, mm, ic in inconsistent_columns[1:]:

        if last_ic != ic:
            inconsistent_segments.append(
                (fc, lc + 1, max_counts, list(last_ic)))
            fc = column
            lc = column
            max_counts = 0
            last_ic = ic

        lc = column
        last_ic = ic
        max_counts = max(max_counts, mm)

    inconsistent_segments.append((fc, lc + 1, max_counts, list(last_ic)))

    if options.loglevel >= 1:
        options.stdlog.write(
            "# found %i inconsistent segments.\n" % len(inconsistent_segments))

    if options.loglevel >= 4:
        for first_res, last_res, nclusters, ic in inconsistent_segments:
            options.stdlog.write("# inconsistent segment: %i:%i\t%i\t%s\n" % (
                first_res, last_res, nclusters, str(ic)))

    # resolve each inconsistent segment separately
    inconsistent_segments.reverse()

    for first_res, last_res, num_clusters, inconsistencies in inconsistent_segments:

        if last_res - first_res < options.min_segment_length:
            continue

        if options.loglevel >= 1:
            options.stdlog.write("# resolving segment %i:%i with %i clusters.\n" % (
                first_res, last_res, num_clusters))

        # obtain clusters of inconsistent entries
        clusters = GetInconsistentClusters(
            mali, num_clusters, first_res, last_res, inconsistencies, options)

        # assign remaining entries to clusters
        incs = set(map(lambda x: x[0], inconsistencies))
        clusters = Assign2Clusters(
            mali, set(identifiers).difference(incs), clusters, first_res, last_res, options)

        if options.loglevel >= 4:
            for c in range(len(clusters)):
                options.stdlog.write("# cluster %i for range %i:%i: %s\n" % (
                    c, first_res, last_res, ";".join(clusters[c])))

        all = set(identifiers)

        # sanity check
        if sum(map(lambda x: len(x), clusters)) != len(identifiers):
            raise "incomplete clusters"

        # no need to insert gaps into all clusters
        c = 0
        d = last_res - first_res
        for cluster in clusters[:-1]:
            mali.insertColumns(first_res, d, keep_fixed=set(cluster))
            first_res += d
            mali.writeToFile(open("cluster_%i" % c, "w"), format="fasta")
            c += 1

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def GetFrameColumns(mali, master, gap_chars="-."):
    """get columns in frame according to master."""

    columns = []
    try:
        sequence = mali[master]
    except KeyError:
        return columns

    x = 0
    t = 0
    while x < len(sequence):
        c = 0
        codon = []
        while x < len(sequence) and c < 3:
            # and sequence[x] in string.uppercase:
            if sequence[x] not in gap_chars:
                codon.append(x)
                c += 1
                t += 1
            x += 1

        if len(codon) == 3:
            columns.append(codon)

    if t % 3 != 0:
        raise "master %s has length %i, which is not divisible by 3" % (
            master, t)

    return columns

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def AddExonInformation(sequence,
                       exons,
                       gap_char="-",
                       mask_lowercase=False,
                       mask_char="x"):
    """add exon information to an aligned string.

    lowercase characters will be masked.
    """
    new_chars = []
    c = 0
    e = 0
    is_upper = True

    for x in sequence:
        if x == gap_char:
            pass
        else:
            if c == exons[e].mPeptideTo:
                if e < len(exons) - 1:
                    e += 1
                    if is_upper:
                        is_upper = False
                    else:
                        is_upper = True
            c += 1

            if x in string.lowercase:
                x = mask_char

            if is_upper:
                x = string.upper(x)
            else:
                x = string.lower(x)

        new_chars.append(x)

    return string.join(new_chars, "")


def checkCodon(codon, options):

    is_aligned = False
    is_ok = True
    is_all_gaps = True
    for x in codon:
        is_all_gaps = is_all_gaps and x in options.gap_chars

        # a codon will be masked, if it either
        # 1. contains a gap character
        # 2. is an unaligned character, i.e.,
        # exons and masked, or no exons and lowerwase
        residue_is_unaligned = (x in options.gap_chars) or \
                               (not ignore_case and x in string.lowercase) or \
                               (ignore_case and x in mask_chars)
        is_aligned = is_aligned or not residue_is_unaligned
        is_ok = is_ok and not residue_is_unaligned

    return is_ok, is_aligned, is_all_gaps

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/prune_multiple_alignment.py 2654 2009-05-06 13:51:22Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-m", "--master-identifier", dest="master", type="string",
                      help="master sequence.")

    parser.add_option("-p", "--master-pattern", dest="master_pattern", type="string",
                      help="master pattern.")

    parser.add_option("--master-species", dest="master_species", type="string",
                      help="species to use as master sequences.")

    parser.add_option("-t", "--translate", dest="filename_translation", type="string",
                      help="filename on where to store translated sequences.")

    parser.add_option("-e", "--exons-file", dest="filename_exons", type="string",
                      help="filename on where to exon information.")

    parser.add_option("-c", "--mark-codons", dest="mark_codons", action="store_true",
                      help="mark codons.")

    parser.add_option("-i", "--ignore-case", dest="ignore_case", action="store_true",
                      help="ignore case (otherwise: lowercase are unaligned chars).")

    parser.add_option("--remove-stops", dest="remove_stops", action="store_true",
                      help="remove stop codons.")

    parser.add_option("--mask-stops", dest="mask_stops", action="store_true",
                      help="mask stop codons.")

    parser.add_option("--mask-char", dest="mask_char", type="string",
                      help="masking character to use.")

    parser.add_option("-f", "--remove-frameshifts", dest="remove_frameshifts", action="store_true",
                      help="remove columns corresponding to frameshifts.")

    parser.add_option("--mask-master", dest="mask_master", action="store_true",
                      help="columns in master to be removed are masked to keep residue numbering.")

    parser.add_option("-s", "--split-exons", dest="split_exons", action="store_true",
                      help="split columns aligned to different exons in the same gene.")

    parser.add_option("-a", "--target", dest="target", type="choice",
                      choices=("paml", ),
                      help="perform cleaning up for certain targets.")

    parser.set_defaults(
        gap_char="-",
        mask_char="n",
        gap_chars="-.",
        separator="|",
        master=None,
        master_species=None,
        filename_translation=None,
        filename_exons=None,
        master_pattern=None,
        remove_stops=False,
        mark_codons=False,
        mask_unaligned=False,
        split_exons=False,
        remove_frameshifts=False,
        min_segment_length=5,
        ignore_case=False,
        mask_stops=False,
        target=None,
        mask_master=False,
    )

    (options, args) = E.Start(parser)

    if options.target == "paml":
        options.mask_stops = True
        options.mask_char = "n"
        options.remove_frameshifts = True

        if options.loglevel >= 1:
            options.stdlog.write(
                "# setting output to paml : removing frameshifts, masking stops with '%s'.\n" % (options.mask_char))

    # 1. read multiple alignment in fasta format
    mali = Mali.Mali()

    mali.readFromFile(sys.stdin)

    if options.loglevel >= 1:
        options.stdlog.write("# read mali with %i entries.\n" % len(mali))

    if len(mali) == 0:
        raise "empty multiple alignment"

    identifiers = mali.getIdentifiers()

    masters = []
    if options.master:
        masters = options.master.split(",")
    elif options.master_pattern:
        for id in identifiers:
            if re.search(options.master_pattern, id):
                masters.append(id)
    elif options.master_species:
        for id in identifiers:
            if options.master_species == id.split(options.separator)[0]:
                masters.append(id)
    else:
        masters.append(identifiers[0])

    if options.loglevel >= 2:
        options.stdlog.write("# master sequences are: %s\n" % str(masters))
        options.stdlog.flush()

    if options.filename_exons:
        exons = Exons.ReadExonBoundaries(open(options.filename_exons, "r"),
                                         filter=set(identifiers),
                                         from_zero=True)

        if options.loglevel >= 2:
            options.stdlog.write("# read exons %i sequences.\n" % len(exons))
    else:
        exons = {}

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # translate characters to upper/lower case according to exon info.
    ##########################################################################
    if exons:
        for id in identifiers:
            if id in exons:
                mali.getSequence(id).mString = AddExonInformation(
                    mali[id], exons[id], mask_char=options.mask_char)

    elif options.ignore_case:
        # convert all to uppercase
        mali.upper()

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # untangle misaligned exons
    ##########################################################################
    if exons and options.split_exons:

        # first split with masters
        if len(masters) > 0:
            SplitExons(mali, exons, masters=masters, options=options)

            if options.loglevel >= 4:
                mali.writeToFile(open("log_mali1", "w"), format="fasta")

        SplitExons(mali, exons, options)

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # remove frameshifts
    ##########################################################################
    if options.remove_frameshifts:
        out_of_frame_columns = []
        if len(masters) == 1:

            frame_columns = GetFrameColumns(mali, masters[0],
                                            gap_chars=options.gap_chars)

        else:

            columns = []

            for id in masters:
                columns += GetFrameColumns(mali, id,
                                           gap_chars=options.gap_chars)

            if len(columns) == 0:
                columns += GetFrameColumns(mali, identifiers[0],
                                           gap_chars=options.gap_chars)

            # sort all columns by tuple. The "shortest" codon will be first: (1,2,3) before (1,2,100),
            # and (1,2,100) before (1,3,4).
            columns.sort(lambda x, y: cmp((x[0], x[2]), (y[0], y[2])))

            # select codons
            frame_columns = []
            last_codon = columns[0]

            for codon in columns[1:]:
                # skip identical codons
                if codon == last_codon:
                    continue

                # take first (shortest) codon in case of identical first
                # residue
                if codon[0] == last_codon[0]:
                    continue

                # if not overlapping, keep
                if codon[0] > last_codon[2]:
                    frame_columns.append(last_codon)
                else:
                    out_of_frame_columns += last_codon

                # if overlapping, but out of register: skip
                last_codon = codon

            frame_columns.append(last_codon)

        # build set of skipped columns
        frame_set = set()
        for column in frame_columns:
            for c in column:
                frame_set.add(c)

        # columns that contain a master sequence that is out of
        # frame
        out_of_frame_set = set(out_of_frame_columns)
        out_of_frame_set = out_of_frame_set.difference(frame_set)

        if options.loglevel >= 1:
            options.stdlog.write(
                "# found %i/%i columns in frame\n" % (len(frame_columns) * 3, mali.getWidth()))

            if options.loglevel >= 5:
                options.stdlog.write(
                    "# frame columns: %i\n" % (len(frame_columns)))
                x = 0
                for column in frame_columns:
                    options.stdlog.write(
                        "# %i\t%s\n" % (x, ",".join(map(str, column))))
                    x += 1

            if options.loglevel >= 5:
                options.stdlog.write(
                    "# Out-of frame columns with residue of masters: %i\n" % (len(out_of_frame_set)))
                options.stdlog.write(
                    "# %s" % ",".join(map(str, out_of_frame_columns)))

        mask_chars = (
            string.upper(options.mask_char), string.lower(options.mask_char))

        to_delete = []

        ignore_case = exons or options.ignore_case

        for id in identifiers:

            ngaps, nmasked = 0, 0

            sequence = mali.getSequence(id).mString

            if options.loglevel >= 7:
                options.stdlog.write(
                    "# processing sequence %s of length %i with gaps\n" % (id, len(sequence)))

            # treat masters differently if they are only to be masked, not
            # pruned.
            # simple mask all characters that are to skipped
            fragments = []
            nstops, ncodons, naligned = 0, 0, 0

            codon = []
            chars = []

            is_master = id in masters

            for x in range(len(sequence)):
                c = sequence[x]

                # delete columns that do not align to
                # a master.
                if x not in frame_set and x not in out_of_frame_set:
                    continue

                chars.append(c)
                if c not in options.gap_chars:
                    codon.append(c)
                if len(codon) % 3 == 0:
                    codon = "".join(codon)
                    codon_is_ok, codon_is_aligned, codon_is_all_gaps = checkCodon(
                        codon, options)

                    if codon_is_aligned:
                        naligned += 1

                    to_mask = False
                    if codon_is_all_gaps:
                        ngaps += len(chars)
                    elif codon_is_ok:
                        ncodons += 1
                        if string.upper(codon) in ("TAG", "TAA", "TGA"):
                            nstops += 1
                            to_mask = True
                    else:
                        to_mask = True
                        nmasked += 1

                    if to_mask:
                        for i in range(len(chars)):
                            if chars[i] not in options.gap_chars:
                                chars[i] = options.mask_char

                    fragments.append("".join(chars))
                    chars = []
                    codon = []

            # mask incomplete codons at the end
            if chars:
                for i in range(len(chars)):
                    if chars[i] not in options.gap_chars:
                        chars[i] = options.mask_char
                fragments.append("".join(chars))

            s = string.join(fragments, "")
            if options.loglevel >= 1:
                options.stdlog.write("# sequence: %s\tpositions: %i\taligned:%i\tcodons: %i\t stops: %i\tgaps: %i\tnmasked: %i\n" % (
                    id, len(fragments), naligned, ncodons, nstops, ngaps, nmasked))
                options.stdlog.flush()

            # postpone deletion in order to not
            # confuse the iteration of ids
            if naligned == 0:
                options.stdlog.write(
                    "# sequence: %s removed because there are no aligned nucleotides.\n" % id)
                to_delete.append(id)
            elif ncodons == 0:
                options.stdlog.write(
                    "# sequence: %s removed because there are no aligned codons.\n" % id)
                to_delete.append(id)
            else:
                mali.setSequence(id, string.join(fragments, ""))

        for id in to_delete:
            del mali[id]

    for id in identifiers:
        if options.mark_codons:
            a = mali[id]
            f = lambda x: a[x:x + 3]
            s = string.join([f(x) for x in range(0, len(a), 3)], " ")
        else:
            s = mali[id]
        options.stdout.write(">%s\n%s\n" % (id, s))

    if options.filename_translation:
        outfile = open(options.filename_translation, "w")
        for id in mali.keys():
            outfile.write(">%s\n%s\n" %
                          (id, Genomics.TranslateDNA2Protein(mali[id])))
        outfile.close()

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

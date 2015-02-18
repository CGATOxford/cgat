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
optic/orthologs2transcripts.py -
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

   python optic/orthologs2transcripts.py --help

Type::

   python optic/orthologs2transcripts.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import getopt
import tempfile
import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.Exons as Exons
import alignlib_lite
import CGAT.WrapperDialign as WrapperDialign
import CGAT.WrapperDBA as WrapperDBA
# clustalw wrapper not up-to-date
# import CGAT.WrapperClustal as WrapperClustal
import CGAT.AlignedPairs as AlignedPairs

USAGE = """python %s [OPTIONS] < orthologs > genes

Version: $Id: optic/orthologs2transcripts.py 2774 2009-09-10 10:00:32Z andreas $

Convert a list of orthologous sequences into a list of orthologous transcripts.

Introns are aligned by the program. They can be anchored with adjacent exon fragment
(option --extend-introns), which are later removed.

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
--mode=                         [genes|transcripts]. If set to genes, ortholog transcripts are predicted.
--map1, --map2                  map of genes to transcripts for first/second column
--peptides1, --peptides2        peptides of transcripts for first/second column
--cds1, --cds2                  files with coding sequences (exon boundaries)
--transcripts1, --transcripts2  files with transcripts (genomic sequences)
--write-exons=                  write exons (unaligned/full/exons/common/caligned)
                                alignment per exon pairs:
                                        unaligned: equivalent exons in unaligned format
                                        caligned:  equivalent exons aligned on the genome (large gaps
                                                        possible due to introns in case of shifted
                                                        boundaries).
                                alignment over all exons:
                                        common:    concatenation of caligned parts.
                                        exons:     alignment of exons on dna removing all gaps.
                                        full:      alignment on genomic dna of exons with all gaps.
--write-introns=                write introns (unaligned/aligned)
--write-pairs=                  write paired transcript information (all)
--extend-introns=               extend intronic sequences with x residues from adjacent exons
--only-best                     take only best ortholog transcripts
--mask                          mask softmasked sequences
--min-intron-length=            minimum length of intron to align
--max-intron-length=            maximum length of intron to align
--missing-max-missing=          maximum number of missing exons
--missing-min-present=          if there is a missing exon, at least # exons have to be present
--boundaries-max-slippage=      maximum of intron boundaries that they can deviate from reference (in codons)
--boundaries-max-missed=        maximum number of boundaries that can be missed
--boundaries-allow-missed=      if there is a missed boundary, at least # others have to be found.
--compress                      write alignment in compressed form
--min-coverage=                 minimum coverage
--mode-genome1=                 [indexed|all]
--mode-genome2=                 [indexed|all]
--min-exon-size=                minimum exon size (important for exon counts) (in codons).
""" % sys.argv[0]


param_long_options = ["verbose=", "help", "map1=", "map2=",
                      "peptides1=", "peptides2=",
                      "cds1=", "cds2=",
                      "transcripts1=", "transcripts2=",
                      "write-exons=", "write-introns=", "write-pairs=",
                      "extend-introns=", "only-best", "mask",
                      "min-intron-length=", "max-intron-length=",
                      "boundaries-max-missed=",
                      "boundaries-max-slippage=",
                      "boundaries-allow-missed=",
                      "missing-max-missing=", "missing-min-present=",
                      "compress", "min-coverage=", "mode=",
                      "min-exon-size=", "min-alignment-exon-overlap="
                      "mode-genome1=", "mode-genome2=",
                      "version"]


param_short_options = "v:hg:"

param_loglevel = 2

param_mode = "genes"

# format for genomic sequences
param_mode_genome1 = "all"
param_mode_genome2 = "all"

param_genome_file = "genome_%s.fasta"

param_filename_peptides1 = None
param_filename_peptides2 = None
param_filename_map1 = None
param_filename_map2 = None
param_filename_transcripts1 = None
param_filename_transcripts2 = None
param_filename_cds1 = None
param_filename_cds2 = None

# minimum length of an exon (in codons): exons smaller than this are ignored
param_min_exon_size = 3

# minimum overlap required between exon and alignment (in codons)
param_min_alignment_exon_overlap = 3

# maximum number of gaps for ortholog pairs
param_max_gaps = 10

# whether or not to check number of introns
param_check_num_introns = 1

# whether or not to check number of introns
param_check_exon_boundaries = 1

# allowable margin for exon boundaries to vary.
param_boundaries_max_slippage = 3

# allowable number of exon boundaries that are wrong
param_boundaries_max_missed = 0

# above this number of exon boundiares, some boundaries may be missed:
param_boundaries_allow_missed = 0

# maximum intron length for alignment: 0 is no restriction
param_max_intron_length = 0

param_compress = 0

# allowable

# flags of output
param_write_introns = ()
param_write_exons = ()
param_write_pairs = ()

param_extend_introns = 0

param_min_coverage = 80.0

param_only_best = 0

param_mask = None

param_min_intron_length = None
param_max_intron_length = None

param_missing_min_present = 0
param_missing_max_missing = 0

# keep smatrix global, so that data does not get deleted when object goes
# out of scope!
global_substitution_matrix = None

param_min_score_sw = 50

# minimum coverage of orthologous exons.
param_min_exon_coverage = 80

# size of dp matrix after which to switch to tuple
# alignment
param_max_matrix_size = 10000000

# ------------------------------------------------------------


def makeAlignatorDNA(type="EMBOSS"):
    """make alignator with DNA substitution matrix.

    EMBOSS style matrix:
    identity = 5
    mismatch = -4
    gop = -16
    gep = -4
    """

    matrix = []

    if type == "EMBOSS":
        # build empty matrix
        for m in range(21):
            matrix.append([-4] * 21)

        # A <-> A
        matrix[0][0] = 5

        # C <-> C
        matrix[1][1] = 5

        # G <-> G
        matrix[5][5] = 5

        # T <-> T
        matrix[16][16] = 5

        gop = -16.0
        gep = -4.0
    else:
        raise "unknown MATRIX."

    handle_tmpfile, filename_tmpfile = tempfile.mkstemp()
    for m in matrix:
        os.write(handle_tmpfile, string.join(map(str, m), "\t") + "\n")
    os.close(handle_tmpfile)

    smatrix = alignlib_lite.readSubstitutionMatrixAA(filename_tmpfile)

    return smatrix, gop, gep

# ------------------------------------------------------------


def makeSubstitutionMatrix(type="EMBOSS"):
    """make alignator with DNA substitution matrix.

    EMBOSS style matrix:
    identity = 5
    mismatch = -4
    gop = -16
    gep = -4

    ClustalW style matrix:
    match = 1 mismatch = 0
    gop = -10 gep = -0.1
    """

    matrix = []

    if type == "emboss":
        match = 5
        mismatch = -4
        gop = -16.0
        gep = -4.0
    elif type == "blastn":
        match = 1
        mismatch = -3
        gop = -5.0
        gep = -2.0
    elif type == "clustal":
        match = 1
        mismatch = 0
        gop = -10.0
        gep = -0.1
    elif type == "iub":
        match = 1.9
        mismatch = 0
        gop = -10.0
        gep = -0.1
    elif type == "trans":
        match = 1
        mismatch = -5
        gop = -10.0
        gep = -0.1
    else:
        raise "unknown MATRIX."

    # build empty matrix
    for m in range(21):
        matrix.append([mismatch] * 21)

    # A <-> A
    matrix[0][0] = match
    # C <-> C
    matrix[1][1] = match
    # G <-> G
    matrix[5][5] = match
    # T <-> T
    matrix[16][16] = match

    if type == "trans":
        # A <-> G
        matrix[0][5] = matrix[5][0] = match
        # C <-> T
        matrix[1][16] = matrix[16][1] = match

    handle_tmpfile, filename_tmpfile = tempfile.mkstemp()
    for m in matrix:
        os.write(handle_tmpfile, string.join(map(str, m), "\t") + "\n")
    os.close(handle_tmpfile)

    smatrix = alignlib_lite.readSubstitutionMatrixAA(filename_tmpfile)
    os.remove(filename_tmpfile)
    return smatrix, gop, gep

# ------------------------------------------------------------


def RemoveRedundantTranscripts(transcripts, peptides):
    """remove redundant entries. Also check for presence.
    """
    sequences = []
    new = []

    new.append(transcripts[0])
    try:
        sequences.append(peptides[transcripts[0][0]])
    except KeyError:
        return new

    for t in transcripts[1:]:
        try:
            n = peptides[t[0]]
        except KeyError:
            continue
        for s in sequences:
            if s == n:
                break
        else:
            new.append(t)
            sequences.append(n)

    return new

# ------------------------------------------------------------


def ReadMapGene2Transcripts(infile, peptides, filter=None):
    """read map of gene to transcripts.
    """
    map_gene2transcripts = {}
    for line in infile:
        if line[0] == "#":
            continue
        try:
            gene_id, sbjct_from, sbjct_to, sbjct_strand, sbjct_token, prediction_id = line[
                :-1].split("\t")
        except ValueError:
            continue

        if filter and gene_id not in filter:
            continue
        if gene_id not in map_gene2transcripts:
            map_gene2transcripts[gene_id] = []
        map_gene2transcripts[gene_id].append(
            (prediction_id, sbjct_token, sbjct_strand, sbjct_from, sbjct_to))

    for key in map_gene2transcripts.keys():
        map_gene2transcripts[key] = RemoveRedundantTranscripts(
            map_gene2transcripts[key], peptides)

    return map_gene2transcripts

# ------------------------------------------------------------


def WriteExons(token1, peptide1, cds1, transcript1,
               token2, peptide2, cds2, transcript2,
               peptide_map_a2b):

    if param_loglevel >= 3:
        for cd in cds1:
            print "#", str(cd)
        for cd in cds2:
            print "#", str(cd)
        print "# peptide_map_a2b", str(alignlib_lite.AlignmentFormatExplicit(peptide_map_a2b))
        sys.stdout.flush()

    dna_map_a2b = Genomics.AlignmentProtein2CDNA(peptide_map_a2b,
                                                 cds1, cds2)

    if len(cds1) != len(cds2):
        if param_loglevel >= 4:
            print ""  # WARNING: different number of exons!"

    seq1 = alignlib_lite.makeSequence(transcript1)
    seq2 = alignlib_lite.makeSequence(transcript2)
    tmp_map_a2b = alignlib_lite.makeAlignmentVector()

    dialign = WrapperDialign.Dialign("-n")
    dialignlgs = WrapperDialign.Dialign("-n -it -thr 2 -lmax 30 -smin 8")
    dba = WrapperDBA.DBA()
    # clustal = WrapperClustal.Clustal()

    matrix, gop, gep = global_substitution_matrix
    alignator_nw = alignlib_lite.makeAlignatorDPFullDP(
        alignlib_lite.ALIGNMENT_GLOBAL, gop, gep, matrix)
    alignator_sw = alignlib_lite.makeAlignatorDPFullDP(
        alignlib_lite.ALIGNMENT_LOCAL, gop, gep, matrix)

    # concatenated alignments for exons:
    # 1: only the common parts
    ali_common1 = ""
    ali_common2 = ""

    e1, e2 = 0, 0
    while cds1[e1].mGenomeTo <= dna_map_a2b.getRowFrom():
        e1 += 1
    while cds2[e2].mGenomeTo <= dna_map_a2b.getColFrom():
        e2 += 1

    nskipped, nerrors = 0, 0

    if param_loglevel >= 5:
        nmapped = 0
        for x in range(dna_map_a2b.getRowFrom(), dna_map_a2b.getRowTo() + 1):
            if dna_map_a2b.mapRowToCol(x) >= 0:
                nmapped += 1
        print "# nmapped=", nmapped
        print str(alignlib_lite.AlignmentFormatEmissions(dna_map_a2b))

    # declare alignments used
    map_intron_a2b = alignlib_lite.makeAlignmentVector()

    result = Exons.CompareGeneStructures(
        cds1, cds2, map_cmp2ref=peptide_map_a2b)

    if param_loglevel >= 2:
        print result.Pretty("#")

    nskipped_exons, nskipped_introns = 0, 0

    last_e1, last_e2 = None, None

    for link in result.mEquivalences:

        if link.mCoverage <= param_min_exon_coverage:
            nskipped_exons += 1
            continue

        e1, e2 = link.mId1, link.mId2

        c1 = cds1[e1]
        c2 = cds2[e2]
        exon_fragment1 = transcript1[c1.mGenomeFrom:c1.mGenomeTo]
        exon_fragment2 = transcript2[c2.mGenomeFrom:c2.mGenomeTo]

        #######################################################################
        # write unaligned exons
        if param_write_exons:
            pair = AlignedPairs.UnalignedPair()

            pair.mCategory = "exon"
            pair.mToken1 = token1
            pair.mId1 = e1 + 1
            pair.mNum1 = len(cds1)
            pair.mLen1 = len(exon_fragment1)
            pair.mSequence1 = exon_fragment1
            pair.mToken2 = token2
            pair.mId2 = e2 + 1
            pair.mNum2 = len(cds2)
            pair.mLen2 = len(exon_fragment2)
            pair.mSequence2 = exon_fragment2
            pair.mFrom1, pair.mTo1 = c1.mGenomeFrom, c1.mGenomeTo,
            pair.mFrom2, pair.mTo2 = c2.mGenomeFrom, c2.mGenomeTo,

            print str(pair)
            sys.stdout.flush()

        #######################################################################
        # build alignment for overlap of both exons
# tmp_map_a2b.clear()
# alignlib_lite.copyAlignment( tmp_map_a2b, dna_map_a2b,
# c1.mGenomeFrom + 1, c1.mGenomeTo )

# if param_loglevel >= 5:
# print "# alignment: %i-%i" % (c1.mGenomeFrom + 1, c1.mGenomeTo)
# for x in alignlib_lite.writeAlignmentTable( tmp_map_a2b ).split("\n"):
# print "#", x
# if tmp_map_a2b.getLength() == 0:
# if param_loglevel >= 1:
# print "# WARNING: empty alignment between exon %i (from %i to %i) and exon %i" % \
#                       (e1,c1.mGenomeFrom + 1, c1.mGenomeTo, e2)
# print "## peptide_map_a2b", peptide_map_a2b.getRowFrom(), peptide_map_a2b.getRowTo(),\
#                       peptide_map_a2b.getColFrom(), peptide_map_a2b.getColTo(), \
# Alignlib.writeAlignmentCompressed(peptide_map_a2b)
# print "## dna_map_a2b", dna_map_a2b.getRowFrom(), dna_map_a2b.getRowTo(),\
#                       dna_map_a2b.getColFrom(), dna_map_a2b.getColTo(), \
# Alignlib.writeAlignmentCompressed(dna_map_a2b)
# for cd in cds1: print "##", str(cd)
# for cd in cds2: print "##", str(cd)
#             nerrors += 1
# continue
#         data = map(lambda x: x.split("\t"), alignlib_lite.writePairAlignment( seq1, seq2, tmp_map_a2b  ).split("\n"))
# if "caligned" in param_write_exons :
# print "exon\tcaligned\t%s\t%i\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s" % ( token1, e1,
#                                                                                token2, e2,
#                                                                                data[0][0], data[0][2],
#                                                                                data[1][0], data[1][2],
# data[0][1], data[1][1] )
#         ali_common1 += data[0][1]
#         ali_common2 += data[1][1]
        #######################################################################
        # write alignment of introns for orthologous introns
        # orthologous introns are between orthologous exons
        if param_write_introns:

            if last_e1 is not None:
                if e1 - last_e1 != 1 or e2 - last_e2 != 1:
                    nskipped_introns += 1
                else:
                    pair = AlignedPairs.UnalignedPair()

                    intron_from1 = cds1[e1 - 1].mGenomeTo
                    intron_to1 = cds1[e1].mGenomeFrom
                    intron_from2 = cds2[e2 - 1].mGenomeTo
                    intron_to2 = cds2[e2].mGenomeFrom

                    intron_fragment1 = transcript1[intron_from1:intron_to1]
                    intron_fragment2 = transcript2[intron_from2:intron_to2]

                    if len(intron_fragment1) == 0 or len(intron_fragment2) == 0:
                        print "## ERROR: empty intron fragments: %i-%i out of %i and %i-%i out of %i." %\
                              (intron_from1, intron_to1, len(transcript1),
                               intron_from2, intron_to2, len(transcript2))
                        continue

                    pair.mCategory = "intron"
                    pair.mToken1 = token1
                    pair.mId1 = e1 + 1
                    pair.mNum1 = len(cds1) - 1
                    pair.mLen1 = len(intron_fragment1)
                    pair.mFrom1 = intron_from1
                    pair.mTo1 = intron_to1
                    pair.mSequence1 = intron_fragment1
                    pair.mToken2 = token2
                    pair.mId2 = e2 + 1
                    pair.mNum1 = len(cds2) - 1
                    pair.mLen2 = len(intron_fragment2)
                    pair.mFrom2 = intron_from2
                    pair.mTo2 = intron_to2
                    pair.mSequence2 = intron_fragment2

                    if (param_min_intron_length and len(intron_fragment1) < param_min_intron_length) or \
                            (param_min_intron_length and len(intron_fragment2) < param_min_intron_length) or \
                            (param_max_intron_length and len(intron_fragment1) > param_max_intron_length) or \
                            (param_max_intron_length and len(intron_fragment2) > param_max_intron_length):
                        if param_loglevel >= 1:
                            print "# skipped: fragment lengths out of bounds for: %s\t%s\t%s\t%s\t%i\t%i" %\
                                  (token1, e1, token2, e2,
                                   len(intron_fragment1),
                                   len(intron_fragment2))
                            sys.stdout.flush()
                            nskipped += 1

                    print str(pair)

# else:
#                         anchored_from1 = intron_from1 - param_extend_introns
#                         anchored_to1 = intron_to1 + param_extend_introns
#                         anchored_from2 = intron_from2 - param_extend_introns
#                         anchored_to2 = intron_to2 + param_extend_introns

#                         anchored_fragment1 = transcript1[anchored_from1:anchored_to1]
#                         anchored_fragment2 = transcript2[anchored_from2:anchored_to2]

# for method in param_write_introns:

# if param_loglevel >= 2:
# print "## aligning with method %s" % method
# sys.stdout.flush

# map_intron_a2b.clear()

# if method == "unaligned":

#                                 from1, to1, ali1, from2, to2, ali2 = 0, 0, intron_fragment1, 0, 0, intron_fragment2

# elif method in ("dialigned", "dbaligned", "clusaligned", "dialignedlgs"):

#                                 tmp_intron_a2b = alignlib_lite.makeAlignmentVector()

# if param_loglevel >= 1:
# print "# aligning with method %s two fragments of length %i and %i" % (method,
# len(anchored_fragment1),
# len(anchored_fragment2))
# sys.stdout.flush()

# if method == "dialigned":
#                                     result = dialign.Align( anchored_fragment1, anchored_fragment2, tmp_intron_a2b )
# elif method == "dialignedlgs":
#                                     result = dialignlgs.Align( anchored_fragment1, anchored_fragment2, tmp_intron_a2b )
# elif method == "dbaligned":
#                                     result = dba.Align( anchored_fragment1, anchored_fragment2, tmp_intron_a2b )
# elif method == "clusaligned":
#                                     result = clustal.Align( anchored_fragment1, anchored_fragment2, tmp_intron_a2b )
# if not result or result.getLength() == 0:
# if param_loglevel >= 1:
# print "# Error: empty intron alignment"
# sys.stdout.flush()
#                                     nerrors += 1
# continue
#                                 tmp_intron_a2b.moveAlignment( anchored_from1, anchored_from2 )
# alignlib_lite.copyAlignment( map_intron_a2b, tmp_intron_a2b,
#                                                        intron_from1 + 1, intron_to1,
# intron_from2 + 1, intron_to2 )
# elif method == "nwaligned":
#                                 seq1.useSegment( cds1[e1-1].mGenomeTo + 1, cds1[e1].mGenomeFrom )
#                                 seq2.useSegment( cds2[e2-1].mGenomeTo + 1, cds2[e2].mGenomeFrom )
#                                 alignator_nw.Align( seq1, seq2, map_intron_a2b )
# seq1.useFullLength()
# seq2.useFullLength()
# elif method == "swaligned":
#                                 seq1.useSegment( cds1[e1-1].mGenomeTo + 1, cds1[e1].mGenomeFrom )
#                                 seq2.useSegment( cds2[e2-1].mGenomeTo + 1, cds2[e2].mGenomeFrom )
#                                 alignlib_lite.performIterativeAlignment( map_intron_a2b, seq1, seq2, alignator_sw, param_min_score_sw )
# seq1.useFullLength()
# seq2.useFullLength()
# else:
#                                 raise "unknown method %s" % method
# if map_intron_a2b.getLength() > 0:
# if param_compress:
#                                     from1, to1 = map_intron_a2b.getRowFrom(), map_intron_a2b.getRowTo()
#                                     from2, to2 = map_intron_a2b.getColFrom(), map_intron_a2b.getColTo()
#                                     ali1, ali2 = Alignlib.writeAlignmentCompressed( map_intron_a2b )
# else:
# data = map(lambda x: x.split("\t"),
# alignlib_lite.writePairAlignment( seq1, seq2, map_intron_a2b  ).split("\n"))
# if len(data) < 2:
#                                         data=[ ( 0, "", 0), (0, "", 0)]
#                                     from1, ali1, to1 = data[0]
#                                     from2, ali2, to2 = data[1]
# print string.join(map(str, ("intron",
# method,
#                                                         token1, e1, len(cds1) - 1, len(intron_fragment1),
#                                                         token2, e2, len(cds2) - 1, len(intron_fragment2),
# map_intron_a2b.getNumGaps(),
# map_intron_a2b.getLength(),
#                                                         map_intron_a2b.getLength() - map_intron_a2b.getNumGaps(),
#                                                         from1, to1, ali1,
#                                                         from2, to2, ali2,
#                                                         intron_from1, intron_to1,
# intron_from2, intron_to2)), "\t")
# sys.stdout.flush()
        last_e1, last_e2 = e1, e2

    ##########################################################################
    # write concatenated exons
# for method in param_write_exons:
# if method == "common":
# print "exon\tcommon\t%s\t%i\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s" % ( token1, 0,
#                                                                            token2, 0,
#                                                                            0, 0,
#                                                                            0, 0,
# ali_common1, ali_common2 )
# elif method == "exons":
# Write full alignment without gaps.
# This will not care about exon boundaries and gaps.
# data = map(lambda x: x.split("\t"),
# alignlib_lite.writePairAlignment( seq1, seq2, dna_map_a2b  ).split("\n"))

# try:
#                 from1, s1, to1, from2, s2, to2 = data[0] + data[1]
# except ValueError:
#                 from1, to1, from2, to2 = 0, 0, 0, 0
#                 s1, s2 = "", ""
#                 nerrors += 1
# except IndexError:
#                 from1, to1, from2, to2 = 0, 0, 0, 0
#                 s1, s2 = "", ""
#                 nerrors += 1

# if from1:
# if len(s1) != len(s2):
# print "# WARNING: alignment of different lengths: %i and %i" % (len(s1), len(s2))
#                     nerrors += 1
#                     from1, to1, from2, to2 = 0, 0, 0, 0
#                     s1, s2 = "", ""
# else:
#                     a1, a2 = [], []
# for x in range( min(len(s1), len(s2)) ):
# if s1[x] != "-" and s2[x] != "-":
#                             a1.append( s1[x] )
#                             a2.append( s2[x] )
#                     s1 = string.join(a1, "")
#                     s2 = string.join(a2, "")

# print "exon\texons\t%s\t%i\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s" % ( (token1, 0,
#                                                                              token2, 0,
#                                                                              from1, to1,
#                                                                              from2, to2,
# s1, s2 ) )
# elif method == "full":
# write full alignment (do not care about exon boundaries)
# data = map(lambda x: x.split("\t"),
# alignlib_lite.writePairAlignment( seq1, seq2, dna_map_a2b  ).split("\n"))
#             if len(data) < 2: data=[ ( 0, "", 0), (0, "", 0)]
# print "exon\tfull\t%s\t%i\t%s\t%i\t%s\t%s\t%s\t%s\t%s\t%s" % ( token1, 0,
#                                                                            token2, 0,
#                                                                            data[0][0], data[0][2],
#                                                                            data[1][0], data[1][2],
# data[0][1], data[1][1] )

    if param_loglevel >= 3:
        print "# skipped_exons=%i, skipped_introns=%i" % (nskipped_exons, nskipped_introns)

    return nerrors, nskipped

# ------------------------------------------------------------


def GetOrthologTranscripts(transcripts1, peptides1, cds1,
                           transcripts2, peptides2, cds2):
    """sort out ortholog relationships between
    transcripts of orthologous genes.

    Orthologs have:
        the same number of exons        
        compatible intron/exon boundaries

    For the remaining transcript pairs, take reciprocal bet hits.

    I see the following:
    0: 0(100%), 1: 0(94%), 2: 0,1(100%)
    0: 0(100%), 1: 0,1,2(100%)

    Selecting 1-0 first, would result in a suboptimal match, because one transcript
    is longer than the other, while matching up 0-0 and 2-1 would be better.

    Objective function: it is the maximal matching/assignment problem. Use greedy
    implementation instead. Assign as much as possible according to descending weights.
    """

    alignator = alignlib_lite.makeAlignatorDPFull(
        alignlib_lite.ALIGNMENT_LOCAL, -10.0, -2.0)

    # for long sequence: use dot alignment with tuple size of three
    dottor = alignlib_lite.makeAlignatorTuples(3)
    alignator_dots = alignlib_lite.makeAlignatorDotsSquared(
        param_gop, param_gep, dottor)

    seqs1 = map(lambda x: alignlib_lite.makeSequence(
        peptides1[x[0]]), transcripts1)
    seqs2 = map(lambda x: alignlib_lite.makeSequence(
        peptides2[x[0]]), transcripts2)

    if param_loglevel >= 4:
        print "# building sequence 1"
    for i in range(len(seqs1)):
        if transcripts1[i][0] not in cds1:
            if param_loglevel >= 4:
                print "# %s not found" % transcripts1[i][0]

    if param_loglevel >= 4:
        print "# building sequence 2"

    for i in range(len(seqs2)):
        if transcripts2[i][0] not in cds2:
            if param_loglevel >= 4:
                print "# %s not found" % transcripts1[i][0]

    if param_loglevel >= 4:
        print "# all-vs-all alignment"

    # do all versus all alignment
    alis1 = []
    alis2 = []
    for i in range(len(seqs1)):
        alis1.append([])
    for i in range(len(seqs2)):
        alis2.append([])

    if param_loglevel >= 3:

        print "#################################"

        for i in range(len(seqs1)):
            for cd in cds1[transcripts1[i][0]]:
                print "#", str(cd)
        print "# versus"
        for i in range(len(seqs2)):
            for cd in cds2[transcripts2[i][0]]:
                print "#", str(cd)
        sys.stdout.flush()

    weights = {}
    for i in range(len(seqs1)):
        prediction_id1, sbjct_token1, sbjct_strand1, sbjct_from1, sbjct_to1 = transcripts1[
            i]

        for j in range(len(seqs2)):
            prediction_id2, sbjct_token2, sbjct_strand2, sbjct_from2, sbjct_to2 = transcripts2[
                j]
            map_a2b = alignlib_lite.makeAlignmentVector()

            m = seqs1[i].getLength() * seqs2[j].getLength()

            if param_loglevel >= 3:
                print "# Starting alignment of pair (%i,%i) of lengths %s:%i and %s:%i" %\
                      (i, j, prediction_id1, seqs1[
                       i].getLength(), prediction_id2, seqs2[j].getLength())
                sys.stdout.flush()

            if m > param_max_matrix_size:
                # switch to tuple alignment if sequences are too large
                if param_loglevel >= 2:
                    print "# WARNING: sequences are of length %i and %i: switching to dot alignment." % (seqs1[i].getLength(), seqs2[j].getLength())
                    sys.stdout.flush()

                alignator_dots.align(map_a2b, seqs1[i], seqs2[j])
            else:
                alignator.align(map_a2b, seqs1[i], seqs2[j])

            coverage_a = 100.0 * \
                (map_a2b.getRowTo() - map_a2b.getRowFrom() + 1) / \
                seqs1[i].getLength()
            coverage_b = 100.0 * \
                (map_a2b.getColTo() - map_a2b.getColFrom() + 1) / \
                seqs2[j].getLength()

            # get copy of cds, but only those overlapping with alignment
            c1 = Exons.GetExonsRange(cds1[prediction_id1],
                                     (map_a2b.getRowFrom() - 1) * 3,
                                     (map_a2b.getRowTo()) * 3 + 1,
                                     full=False,
                                     min_overlap=param_min_alignment_exon_overlap,
                                     min_exon_size=param_min_exon_size)
            c2 = Exons.GetExonsRange(cds2[prediction_id2],
                                     (map_a2b.getColFrom() - 1) * 3,
                                     (map_a2b.getColTo()) * 3 + 1,
                                     full=False,
                                     min_overlap=param_min_alignment_exon_overlap,
                                     min_exon_size=param_min_exon_size)

            # check exon boundaries, look at starts, skip first exon
            def MyMap(a, x):
                while x <= a.getRowTo():
                    c = a.mapRowToCol(x)
                    if c:
                        return c
                    x += 1
                else:
                    return 0

            mapped_boundaries = map(
                lambda x: MyMap(map_a2b, x.mPeptideFrom / 3 + 1), c1[1:])
            mapped_boundaries.sort()
            reference_boundaries = map(
                lambda x: x.mPeptideFrom / 3 + 1, c2[1:])
            reference_boundaries.sort()

            nmissed_cmp2ref = Exons.CountMissedBoundaries(
                mapped_boundaries, reference_boundaries, param_boundaries_max_slippage)
            nmissed_ref2cmp = Exons.CountMissedBoundaries(
                reference_boundaries, mapped_boundaries, param_boundaries_max_slippage)

            min_nmissed = min(nmissed_cmp2ref, nmissed_ref2cmp)

            # set is_ok for the whole thing
            # no intron: is ok
            is_ok = 0
            if (len(c1) == 1 and len(c2) == 1):
                is_ok = 1
            else:
                # allow for missed boundaries, if param_boundaries_allow_missed
                # > 0
                if min_nmissed == 0:
                    is_ok = 1
                else:
                    if param_boundaries_allow_missed and \
                            len(mapped_boundaries) >= param_boundaries_allow_missed and \
                            min_nmissed <= param_boundaries_max_missed:
                        is_ok = 1

            cc = min(coverage_a, coverage_b)
            if cc >= param_min_coverage:
                is_ok_coverage = 1
            else:
                is_ok_coverage = 0

            # check for missing introns
            is_ok_exons = 1
            if abs(len(c1) - len(c2)) != 0:
                if param_missing_max_missing:
                    if ((abs(len(c1) - len(c2)) > param_missing_max_missing) or
                            (min(len(c1), len(c2)) < param_missing_min_present)):
                        is_ok_exons = 0
                else:
                    is_ok_exons = 0

            if param_loglevel >= 3:
                print "# i=", i, "li=", len(c1), "j=", j, "lj=", len(c2), \
                      "boundaries_ok=", is_ok, \
                      "nexons_ok=", is_ok_exons, \
                      "missed_c2r=", nmissed_cmp2ref, \
                      "missed_r2c=", nmissed_ref2cmp, \
                      "min_cov=", cc, \
                      "mapped=", mapped_boundaries, \
                      "reference=", reference_boundaries

                print "#", string.join(map(str, (alignlib_lite.AlignmentFormatEmissions(map_a2b),
                                                 map_a2b.getNumGaps(), coverage_a, coverage_b)), "\t")
                sys.stdout.flush()

            # dump out pairs
            for method in param_write_pairs:
                if method == "all":
                    print string.join(map(str, (
                        "pair", method,
                        prediction_id1,
                        prediction_id2,
                        sbjct_token1, sbjct_strand1, sbjct_from1, sbjct_to1, seqs1[
                            i].getLength(),
                        sbjct_token2, sbjct_strand2, sbjct_from2, sbjct_to2, seqs2[
                            j].getLength(),
                        map_a2b.getRowFrom(), map_a2b.getRowTo(), row_ali,
                        map_a2b.getColFrom(), map_a2b.getColTo(), col_ali,
                        map_a2b.getNumGaps(), coverage_a, coverage_b,
                        nmissed_cmp2ref, mapped_boundaries,
                        nmissed_ref2cmp, reference_boundaries,
                        i, j, len(c1), len(c2), cc, is_ok, is_ok_exons, is_ok_coverage)), "\t")
                elif method == "alignment":
                    print string.join(map(str, (
                        "pair", method,
                        prediction_id1, prediction_id2,
                        map_a2b.getRowFrom(), map_a2b.getRowTo(), row_ali,
                        map_a2b.getColFrom(), map_a2b.getColTo(), col_ali,
                        map_a2b.getNumGaps(), coverage_a, coverage_b)), "\t")
                elif method == "location":
                    print string.join(map(str, (
                        "pair", method,
                        prediction_id1,
                        prediction_id2,
                        sbjct_token1, sbjct_strand1, sbjct_from1, sbjct_to1, seqs1[
                            i].getLength(),
                        sbjct_token2, sbjct_strand2, sbjct_from2, sbjct_to2, seqs2[j].getLength())), "\t")
            if not is_ok_exons:
                if param_loglevel >= 4:
                    print "# rejected %i and %i: too many exons difference." % (i, j)
                continue

            if param_check_exon_boundaries:
                if not is_ok:
                    continue

            if cc < param_min_coverage:
                continue

            if cc in weights:
                weights[cc] = []

            alis1[i].append((coverage_a, j))
            alis2[j].append((coverage_b, i))

            weights[cc].append((i, j, map_a2b))

    # sort out alignments
    ww = weights.keys()
    ww.sort()
    ww.reverse()

    pairs = []
    assigned1 = {}
    assigned2 = {}

    if param_loglevel >= 3:
        print "# alis1=", alis1
        print "# alis2=", alis2
        print "# --------------------------------------"

    for w in ww:
        for i, j, map_a2b in weights[w]:
            if i not in assigned1 and j not in assigned2:
                pairs.append((transcripts1[i], transcripts2[j], w, map_a2b))
                assigned1[i] = 1
                assigned2[j] = 1
        if len(assigned1) == len(transcripts1):
            break
        if len(assigned2) == len(transcripts2):
            break

    return pairs


def PruneHash(hash, *args):

    for k, vv in hash.items():
        new_vv = []
        for v in vv:
            found = 1
            for arg in args:
                found = found and v[0] in arg
            if found:
                new_vv.append(v)
            else:
                if param_loglevel >= 4:
                    print "# eliminated transcripts for %s: %s" % (k, v)
        hash[k] = new_vv

    return hash


def ReadTranscriptsAndCds(transcript_ids1, transcript_ids2):

    if param_loglevel >= 1:
        print "# reading %i left and %i right transcripts" % (len(transcript_ids1), len(transcript_ids2))
        sys.stdout.flush()
    if param_loglevel >= 1:
        print "# reading exon boundaries."
        sys.stdout.flush()

    cds1 = Exons.ReadExonBoundaries(
        open(param_filename_cds1, "r"), filter=transcript_ids1, reset=True)
    cds2 = Exons.ReadExonBoundaries(
        open(param_filename_cds2, "r"), filter=transcript_ids2, reset=True)

    if param_loglevel >= 1:
        print "# read %i left and %i right cds" % (len(cds1), len(cds2))
        sys.stdout.flush()

    if param_loglevel >= 2:
        if len(cds1) != len(transcript_ids1):
            print "# missed in left:  %s" % ":".join(set(transcript_ids1.keys()).difference(cds1.keys()))
        if len(cds2) != len(transcript_ids2):
            print "# missed in right: %s" % ":".join(set(transcript_ids2.keys()).difference(cds2.keys()))

    if param_loglevel >= 1:
        print "# reading genomic sequences."
        sys.stdout.flush()

    transcripts1 = {}
    if param_filename_transcripts1:
        if param_mode_genome1 == "indexed":
            transcripts1 = Genomics.ParseFasta2HashFromIndex(
                param_filename_transcripts1, filter=transcript_ids1)
        else:
            transcripts1 = Genomics.ReadGenomicSequences(open(param_filename_transcripts1, "r"),
                                                         do_reverse=0,
                                                         filter=transcript_ids1,
                                                         mask=param_mask)
    transcripts2 = {}
    if param_filename_transcripts2:
        if param_mode_genome2 == "indexed":
            transcripts2 = Genomics.ParseFasta2HashFromIndex(
                param_filename_transcripts2, filter=transcript_ids2)
        else:
            transcripts2 = Genomics.ReadGenomicSequences(open(param_filename_transcripts2, "r"),
                                                         do_reverse=0,
                                                         filter=transcript_ids2,
                                                         mask=param_mask)
    if param_loglevel >= 1:
        print "# read %i left and %i right transcript sequences" % (len(transcripts1), len(transcripts2))
        sys.stdout.flush()

    return transcripts1, transcripts2, cds1, cds2

# ------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    try:
        optlist, args = getopt.getopt(
            sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o == "--map1":
            param_filename_map1 = a
        elif o == "--map2":
            param_filename_map2 = a
        elif o == "--cds1":
            param_filename_cds1 = a
        elif o == "--cds2":
            param_filename_cds2 = a
        elif o == "--peptides1":
            param_filename_peptides1 = a
        elif o == "--peptides2":
            param_filename_peptides2 = a
        elif o == "--transcripts1":
            param_filename_transcripts1 = a
        elif o == "--transcripts2":
            param_filename_transcripts2 = a
        elif o == "--write-exons":
            param_write_exons = filter(lambda x: x != "", a.split(","))
        elif o == "--write-introns":
            param_write_introns = filter(lambda x: x != "", a.split(","))
        elif o == "--write-pairs":
            param_write_pairs = filter(lambda x: x != "", a.split(","))
        elif o == "--extend-introns":
            param_extend_introns = int(a)
        elif o == "--min-intron-length":
            param_min_intron_length = int(a)
        elif o == "--max-intron-length":
            param_max_intron_length = int(a)
        elif o == "--only-best":
            param_only_best = 1
        elif o == "--mask":
            param_mask = 1
        elif o == "--compress":
            param_compress = 1
        elif o == "--mode":
            param_mode = a
        elif o == "--boundaries-max-slippage":
            param_boundaries_max_slippage = int(a)
        elif o == "--boundaries-max-missed":
            param_boundaries_max_missed = int(a)
        elif o == "--boundaries-allow-missed":
            param_boundaries_allow_missed = int(a)
        elif o == "--missing-max-missing":
            param_missing_max_missing = int(a)
        elif o == "--missing-min-present":
            param_missing_min_present = int(a)
        elif o == "--min-coverage":
            param_min_coverage = int(a)
        elif o == "--min-exon-size":
            param_min_exon_size = int(a)
        elif o == "--min-alignment-exon-overlap":
            param_min_alignment_exon_overlap = int(a)

    if len(args) > 0:
        print USAGE, "no arguments required."
        sys.exit(2)

    print E.GetHeader()
    print E.GetParams()
    sys.stdout.flush()

    global_substitution_matrix = makeSubstitutionMatrix("emboss")

    # get a list of used transcripts first
    pairs = map(lambda x: x[:-1].split("\t")[:2],
                filter(lambda x: x[0] != "#", sys.stdin.readlines()))

    if param_loglevel >= 1:
        print "# read %i pairs" % len(pairs)
        print "# reading peptide sequences."
        sys.stdout.flush()

    peptides1 = Genomics.ReadPeptideSequences(
        open(param_filename_peptides1, "r"))
    peptides2 = Genomics.ReadPeptideSequences(
        open(param_filename_peptides2, "r"))

    if param_loglevel >= 1:
        print "# received %i left and %i right peptide sequences" % (len(peptides1), len(peptides2))
        sys.stdout.flush()

    ninput, noutput, npairs, total_nerrors, total_nskipped = 0, 0, 0, 0, 0

    if param_mode == "genes":

        gene_ids1 = {}
        gene_ids2 = {}
        for gene_id1, gene_id2 in pairs:
            gene_ids1[gene_id1] = 1
            gene_ids2[gene_id2] = 1

        if param_loglevel >= 1:
            print "# received %i left and %i right gene_ids" % (len(gene_ids1), len(gene_ids2))
            sys.stdout.flush()

        if param_loglevel >= 1:
            print "# reading maps of genes to transcripts."
            sys.stdout.flush()

        map_gene2transcripts1 = ReadMapGene2Transcripts(
            open(param_filename_map1, "r"), peptides1, gene_ids1)
        map_gene2transcripts2 = ReadMapGene2Transcripts(
            open(param_filename_map2, "r"), peptides2, gene_ids2)

        transcript_ids1 = {}
        for x in map_gene2transcripts1.values():
            for y in x:
                transcript_ids1[y[0]] = 1
        transcript_ids2 = {}
        for x in map_gene2transcripts2.values():
            for y in x:
                transcript_ids2[y[0]] = 1

        transcripts1, transcripts2, cds1, cds2 = ReadTranscriptsAndCds(
            transcript_ids1, transcript_ids2)

        if param_loglevel >= 1:
            print "# reading has finished."
            sys.stdout.flush()

        # prune gene2transcripts, so that it only contains entries present
        PruneHash(map_gene2transcripts1,
                  transcripts1, peptides1, cds1)
        PruneHash(map_gene2transcripts2,
                  transcripts2, peptides2, cds2)

        for gene_id1, gene_id2 in pairs:

            ninput += 1

            if param_loglevel >= 1:
                print "# processing %s and %s" % (gene_id1, gene_id2)

            t1 = map_gene2transcripts1[gene_id1]
            t2 = map_gene2transcripts2[gene_id2]

            pairs = GetOrthologTranscripts(t1, peptides1, cds1,
                                           t2, peptides2, cds2)

            nerrors, nskipped = 0, 0
            for transcript1, transcript2, weight, map_a2b in pairs:

                q1 = transcript1[0]
                q2 = transcript2[0]

                print string.join(map(str, ("transcript",
                                            string.join(
                                                map(str, transcript1), "\t"),
                                            len(cds1[q1]), cds1[q1][
                                                0].mGenomeFrom, cds1[q1][-1].mGenomeTo,
                                            string.join(
                                                map(str, transcript2), "\t"),
                                            len(cds2[q2]), cds2[q2][
                                                0].mGenomeFrom, cds2[q2][-1].mGenomeTo,
                                            )), "\t")
                npairs += 1

                if map_a2b:
                    if q1 in transcripts1 and q2 in transcripts2:
                        nerrors, nskipped = WriteExons(q1, peptides1[q1], cds1[q1], transcripts1[q1],
                                                       q2, peptides2[q2], cds2[
                                                           q2], transcripts2[q2],
                                                       map_a2b)

                total_nerrors += nerrors
                total_nskipped += nskipped

                if param_only_best:
                    break

            if param_loglevel >= 3:
                print "# in1=%i, in2=%i, selected=%i, missed=%i, errors=%i, skipped=%i" % (len(t1), len(t2),
                                                                                           len(pairs),
                                                                                           min(len(t1), len(
                                                                                               t2)) - len(pairs),
                                                                                           nerrors, nskipped)

            print string.join(map(str, ("gene",
                                        gene_id1, gene_id2,
                                        len(t1), len(t2), len(pairs),
                                        min(len(t1), len(t2)) - len(pairs),
                                        nerrors, nskipped)), "\t")

            noutput += 1

    elif param_mode == "transcripts":

        transcript_ids1 = {}
        transcript_ids2 = {}
        for x, y in pairs:
            transcript_ids1[x] = 1
            transcript_ids2[y] = 1

        transcripts1, transcripts2, cds1, cds2 = ReadTranscriptsAndCds(
            transcript_ids1, transcript_ids2)

        if param_loglevel >= 1:
            print "# reading has finished."
            sys.stdout.flush()

        alignator = alignlib_lite.makeAlignatorDPFull(
            alignlib_lite.ALIGNMENT_LOCAL, -10.0, -2.0)

        for q1, q2 in pairs:

            ninput += 1

            if param_loglevel >= 1:
                print "# processing %s and %s" % (q1, q2)

            if q1 in transcripts1 and q2 in transcripts2:
                map_a2b = alignlib_lite.makeAlignmentVector()

                alignator.align(map_a2b,
                                alignlib_lite.makeSequence(peptides1[q1]),
                                alignlib_lite.makeSequence(peptides2[q2]))

                if map_a2b.getLength() == 0:
                    if param_loglevel >= 1:
                        print "# Alignment failed between %s and %s" % (q1, q2)
                        sys.stdout.flush()

                    ntotal_errors += 1
                    continue

                nerrors, nskipped = WriteExons(q1, peptides1[q1], cds1[q1], transcripts1[q1],
                                               q2, peptides2[q2], cds2[
                                                   q2], transcripts2[q2],
                                               map_a2b)

                total_nerrors += nerrors
                total_nskipped += nskipped
                noutput += 1

    print "# ninput=%i, noutput=%i, npairs=%i, nerrors=%i, nskipped=%i" % (ninput, noutput, npairs, total_nerrors, total_nskipped)
    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

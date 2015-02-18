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
optic/analyze_orthology_pairwise.py - 
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

   python optic/analyze_orthology_pairwise.py --help

Type::

   python optic/analyze_orthology_pairwise.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import time
import copy
import gzip
import CGAT.Experiment as E
import CGAT.Orthologs as Orthologs
import CGAT.Genomics as Genomics
import csv
import scipy
import scipy.stats
import numpy
import alignlib_lite
import CGAT.AlignedPairs as AlignedPairs
import CGAT.TreeTools as TreeTools
import CGAT.BlastAlignments as BlastAlignments

"""program $Id: optic/analyze_orthology_pairwise.py 2781 2009-09-10 11:33:14Z andreas $

analyse orthology sets.

"""


def AnalyseOrphans(orphans, outfile,
                   all_genes_this,
                   aligned_genes_this,
                   assigned_genes_this,
                   aligned_genes_other,
                   assigned_genes_other,
                   map_query2best,
                   map_transcript2location_other,
                   map_transcript2gene_other,
                   options):

    categories = ('nopreds', 'alternative',
                  'missed-class', 'missed-input', 'missed-replaced',
                  'missed-suboptimal',
                  'missed-link', 'missed-assignment',
                  'unknown')

    sums = {}
    for x in categories:
        sums[x] = {'all': 0}

    fields = ['geneid', 'ngenes', 'ntranscripts', 'status', 'nmatches', 'bestpid', 'bestcoverage', 'genes',
              'best_id', 'best_gid', 'best_query', 'best_pide', 'best_coverage', 'best_class']

    outfile.write("\t".join(fields) + "\n")
    writer = csv.DictWriter(outfile,
                            fields,
                            dialect=options.csv_dialect,
                            lineterminator=options.csv_lineterminator,
                            extrasaction='ignore')

    f = 0
    for g, t in orphans.items():

        CountOrphans(outfile, g, t,
                     sums,
                     dbhandle, writer,
                     all_genes_this,
                     aligned_genes_this,
                     assigned_genes_this,
                     aligned_genes_other,
                     assigned_genes_other,
                     map_query2best=map_query2best,
                     map_transcript2location_other=map_transcript2location_other,
                     map_transcript2gene_other=map_transcript2gene_other,
                     options=options)

        f += 1
        if f > 40:
            break

    outfile.write("# summary of orphans in %s versus %s\n" %
                  (options.schema1, options.schema2))
    outfile.write("category\tall\t%s\n" % "\t".join(options.quality_priority))
    for k in categories:
        outfile.write("%s\t%i" % (k, sums[k]['all']))
        for x in options.quality_priority:
            if x in sums[k]:
                outfile.write("\t%i" % sums[k][x])
            else:
                outfile.write("\t0")

        outfile.write("\n")

    outfile.write("# ninput=%i, nskipped=%i, ngenes=%i\n" %
                  (ninput, nskipped, len(orphans)))


def CountOrphans(outfile,
                 geneid, transcripts, sums,
                 dbhandle, writer,
                 all_genes,
                 aligned_genes_this,
                 assigned_genes_this,
                 aligned_genes_other,
                 assigned_genes_other,
                 map_query2best,
                 map_transcript2location_other,
                 map_transcript2gene_other,
                 options):
    """analyse a list of orphans.
    """

    input_classes = options.input_classes

    row = {}
    row['geneid'] = geneid

    genes = Orthologs.GetGenes(transcripts)

    row['ngenes'] = len(genes)
    row['ntranscripts'] = len(transcripts)
    row['genes'] = ",".join(genes.keys())

    # if map_query2best is defined use it.
    if map_query2best:
        # get best matching gene
        best_match = None
        best = 0
        for x in genes.keys():
            if x not in map_query2best:
                continue
            (nmatches, best_query_coverage, best_pidentity,
             prediction_id, gene_id, query_token, quality,
             query_coverage, pidentity,
             sbjct_token, sbjct_strand,
             sbjct_genome_from, sbjct_genome_to) = map_query2best[x]

            v = query_coverage * pidentity

            if v >= best:
                best = v
                best_match = x

        status = None

        if best_match:
            (nmatches, best_query_coverage, best_pidentity,
             prediction_id, gene_id, query_token, quality,
             query_coverage, pidentity,
             sbjct_token, sbjct_strand,
             sbjct_genome_from, sbjct_genome_to) = map_query2best[best_match]

            row['nmatches'] = nmatches
            row['bestpid'] = best_pidentity
            row['bestcoverage'] = best_query_coverage

            # -----------------------------------------------------------
            is_alternative = False
            is_ingraph = False
            has_match = False
            for g in genes.keys():
                if g in assigned_genes:
                    is_alternative = True
                if g in aligned_genes:
                    is_ingraph = True

            if is_alternative:
                status = "alternative"
            else:
                row['best_id'] = prediction_id
                row['best_gid'] = gene_id
                row['best_query'] = best_match
                row['best_class'] = quality
                row['best_coverage'] = query_coverage
                row['best_pide'] = pidentity

                # check if assigned to other gene
                if prediction_id in map_transcript2gene:
                    g = map_transcript2gene[prediction_id]
                    if g in assigned_genes_other:
                        status = "missed-suboptimal"

                if not status:
                    # check if best prediction overlaps with assignments to
                    # other genes
                    overlaps = []

                    for transcript, location in map_transcript2location.items():
                        if location[0] == sbjct_token and location[1] == sbjct_strand and \
                                min(location[3], sbjct_genome_to) - max(location[2], sbjct_genome_from) > 0:
                            overlaps.append(transcript)

                    is_replaced = False
                    if overlaps:
                        alternatives = {}
                        for o in overlaps:
                            if o in map_transcript2gene:
                                g = map_transcript2gene[o]
                                if g in assigned_genes_other:
                                    alternatives[g] = 1
                        if len(alternatives) > 0:
                            is_replaced = True
                            outfile.write("# alternative assignments available for %s: %s\n" % (
                                prediction_id, ";".join(alternatives)))

                    if is_replaced:
                        status = "missed-replaced"
                    elif quality not in input_classes:
                        status = "missed-class"
                    elif not is_ingraph:
                        status = "missed-input"
                    elif str(gene_id) in aligned_genes[best_match]:
                        status = "missed-assignment"
                    else:
                        status = "missed-link"

                if not status:
                    status = "unknown"

                x = sums[status]
                if quality not in x:
                    x[quality] = 0
                x[quality] += 1

        else:
            row['nmatches'] = 0
            status = "nopreds"

    else:
        # prediction information available

        # -----------------------------------------------------------
        is_alternative = False
        is_ingraph = False
        is_replaced = False
        has_match = False
        for g in genes.keys():
            if g in assigned_genes_this:
                is_alternative = True
            if g in aligned_genes_this:
                is_ingraph = True

        if is_alternative:
            status = "alternative"
        else:
            if is_replaced:
                status = "missed-replaced"
            elif quality not in input_classes:
                status = "missed-class"
            elif not is_ingraph:
                status = "missed-input"
            else:
                status = "missed-link"

        if not status:
            status = "unknown"

            x = sums[status]
            if quality not in x:
                x[quality] = 0
            x[quality] += 1

    row['status'] = status
    sums[status]['all'] += 1

    writer.writerow(row)

# ------------------------------------------------------------------------


def PrintCrossAssignments(outfile,
                          assigned_genes, found,
                          schema, clusters,
                          first=True,
                          map_transcript2location={},
                          map_contig2junk={},
                          separator="|"):
    """print genes with more than one transcript and where the transcripts
    are assigned to different orthologous genes.
    """

    vals = []
    out = []
    nalternative_transcripts = 0
    njunk = 0
    ndisjoint = 0
    noverlap = 0
    for g in found.keys():

        l = len(assigned_genes[g])
        if l > 1:
            # check for alternative transcript clusters
            # all other clusters must only contain a single gene
            # and it must always be the same gene
            is_alternative_transcripts = True
            other_genes = {}
            for x in assigned_genes[g]:
                t1, t2, g1, g2, w = clusters[x]
                if first:
                    if len(g2) > 1:
                        is_alternative_transcripts = False
                    for gg in g2:
                        if gg not in other_genes:
                            other_genes[gg] = []
                        other_genes[gg] += t2
                else:
                    if len(g1) > 1:
                        is_alternative_transcripts = False
                    for gg in g1:
                        if gg not in other_genes:
                            other_genes[gg] = []
                        other_genes[gg] += t1

            if len(other_genes) > 1:
                is_alternative_transcripts = False

            if is_alternative_transcripts:
                l = 1
                nalternative_transcripts += 1
            else:
                # check if genes are on assembled chromosomes:
                is_junk = False

                if map_contig2junk and map_transcript2location:
                    non_junk_genes = {}
                    for x, ttt in other_genes.items():
                        for tt in ttt:
                            xs, xt, xg, xq = tt.split(separator)
                            if map_transcript2location[xt][0] not in map_contig2junk:
                                non_junk_genes[xg] = 1

                    # only one gene without junk location, all others with junk
                    # location
                    if len(non_junk_genes) <= 1:
                        is_junk = True

                if is_junk:
                    l = 1
                    njunk += 1

        vals.append(l)
        if l > 1:
            status = "?"

            # check transcripts for genes in other genome and see if they
            # overlap.
            if map_transcript2location:
                if CheckAllOverlap(other_genes, map_transcript2location):
                    status = "overlap"
                    noverlap += 1
                else:
                    status = "disjoint"
                    ndisjoint += 1

            out.append((l, g, status, other_genes))

    outfile.write("# %s: histogram of cross-assigned genes (%i alternative transcript clusters fused)\n" %
                  (schema, nalternative_transcripts))
    outfile.write("bin\tcounts")

    h = scipy.stats.histogram2(vals, range(1, 20))
    for x in range(1, 20):
        outfile.write("%i\t%i\n" % (x, h[x - 1]))

    outfile.write("# %s: cross-assigned genes and transcripts\n" % schema)
    out.sort()
    for x, g, status, o in out:
        outfile.write("%i\t%s\t%s\t%s\n" % (x, g, status, ";".join(o.keys())))
        for a, b in o.items():
            outfile.write("\t\t\t\t%s\t%s\n" % (a, ";".join(b)))

    outfile.write("# %s: nalternatives=%i, junk=%i, ncross=%i, ndisjoint=%i, noverlap=%i\n" % (schema,
                                                                                               nalternative_transcripts, njunk,
                                                                                               len(out),
                                                                                               ndisjoint, noverlap))


def GetAssignments(genes, map_transcript2location, separator="|"):

    assignments = []
    # get maximum range per gene
    for g, vv in genes.items():

        max_to = 0
        min_from = 0
        for v in vv:
            s, t, g, q = v.split(separator)
            (sbjct_token, sbjct_strand, sbjct_genome_from,
             sbjct_genome_to) = map_transcript2location[t]
            if min_from == 0:
                min_from = sbjct_genome_from
            else:
                min_from = min(min_from, sbjct_genome_from)
            max_to = max(max_to, sbjct_genome_to)

        assignments.append((sbjct_token, sbjct_strand, min_from, max_to))

    assignments.sort()

    return assignments

# ------------------------------------------------------------------------


def CheckOverlap(genes, map_transcript2location,
                 separator="|"):
    """check if there is overlap between paralogs.

    returns false if at least one pair of transcripts per gene is overlapping.
    """

    assignments = GetAssignments(genes, map_transcript2location, separator)

    for x in range(1, len(assignments)):
        if assignments[x - 1][0] == assignments[x][0] and \
           assignments[x - 1][1] == assignments[x][1] and \
           assignments[x - 1][3] > assignments[x][2]:
            return False

    return True

# ------------------------------------------------------------------------


def CheckAllOverlap(genes, map_transcript2location,
                    separator="|"):
    """check if there is overlap between paralogs.

    returns true if all pairs overlap.
    """

    assignments = GetAssignments(genes, map_transcript2location, separator)

    for x in range(1, len(assignments)):
        if assignments[x - 1][0] != assignments[x][0] or \
           assignments[x - 1][1] != assignments[x][1] or \
           assignments[x - 1][3] < assignments[x][2]:
            return False

    return True

# ------------------------------------------------------------------------


def FilterTranscriptsByPriority(outfile, chunk):
    """given a chunk of overlapping transcripts, keep only
    those of the highest quality gene."""

    chunk.sort()

    best_gene = chunk[0][1].mGene

    filtered = []
    for quality, transcript in chunk:
        if transcript.mGene == best_gene:
            filtered.append(transcript)
        else:
            outfile.write("# removed overlapping prediction: %s\n" %
                          str(transcript))

    return filtered
# ------------------------------------------------------------------------


def RemoveRedundancy(outfile,
                     transcripts,
                     map_transcript2location,
                     map_quality2priority,
                     separator="|"):
    """remove overlapping predictions from different genes (conservatively)
    """

    if len(transcripts) == 0:
        return transcripts

    assignments = []

    # get all ranges for transcripts
    for tt in transcripts:
        try:
            (sbjct_token, sbjct_strand,
             sbjct_genome_from, sbjct_genome_to) = map_transcript2location[tt.mTranscript]
        except KeyError:
            raise "key %s not in map_transcript2location, examples are: %s" % (
                tt.mTranscript, map_transcript2location.keys()[:10])

        assignments.append(
            (sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to, tt))

    # cluster transcripts by overlap
    assignments.sort()
    chunk = []
    last_token, last_strand, last_from, last_to, last_tt = assignments[0]
    chunk.append((map_quality2priority[last_tt.mQuality], last_tt))

    new_transcripts = []
    for this_token, this_strand, this_from, this_to, this_tt in assignments[1:]:
        if last_token != this_token or \
                last_strand != this_strand or \
                this_from > last_to:
            new_transcripts += FilterTranscriptsByPriority(outfile, chunk)
            chunk = []
            last_to = this_to
            last_token, last_strand = this_token, this_strand

        last_to = max(last_to, this_to)
        chunk.append((map_quality2priority[this_tt.mQuality], this_tt))

    new_transcripts += FilterTranscriptsByPriority(outfile, chunk)

    return new_transcripts

# -----------------------------------------------------------------------------


def CountScores(outfile, t1, t2,
                graph, separator="|",
                cds1={}, cds2={},
                options={}):
    """count scores between t1 and t2 in graph.

    return lists of scores between t1 and within clusters and
    the number of not found vertices and links
    """

    between_scores = []
    within_scores = []
    nmissed_vertices = 0
    nmissed_links_within = 0
    nmissed_links_between = 0

    for tt1 in t1:
        if tt1 in graph:
            # within links
            for tt2 in t1:
                if tt1 == tt2:
                    continue
                for tt3, score in graph[tt1]:
                    if tt3 == tt2:
                        within_scores.append(score)
                        if options.loglevel >= 3:
                            outfile.write("# %s\t%s\t%6.4f\n" %
                                          (tt1, tt2, score))
                        break
                else:
                    score = 0.0
                    # check if same genes
                    # check if cds are identical
                    if tt1 in cds1 and tt2 in cds1 and \
                       (cds1[tt1] == cds1[tt2] or
                            cds1[tt1] in cds1[tt2] or
                            cds1[tt2] in cds1[tt1]):
                        within_scores.append(score)
                        if options.loglevel >= 3:
                            outfile.write(
                                "# %s\t%s\t%6.4f\tsame sequence\n" % (tt1, tt2, score))
                    else:
                        # check if same gene
                        xs, xt, xg, xq = tt1.split(separator)
                        ys, yt, yg, yq = tt2.split(separator)
                        if xg != yg:
                            if options.loglevel >= 3:
                                outfile.write(
                                    "# %s\t%s\t%6.4f\tlink not found\n" % (tt1, tt2, score))
                            nmissed_links_within += 1
                        else:
                            if options.loglevel >= 3:
                                outfile.write(
                                    "# %s\t%s\t%6.4f\tsame gene\n" % (tt1, tt2, score))
            # between links
            for tt2 in t2:
                for tt3, score in graph[tt1]:
                    if tt3 == tt2:
                        between_scores.append(score)
                        if options.loglevel >= 3:
                            outfile.write("# %s\t%s\t%6.4f\n" %
                                          (tt1, tt2, score))
                        break
                else:
                    score = 0.0
                    # check if cds are identical
                    if tt1 in cds1 and tt2 in cds2 and \
                       (cds1[tt1] == cds2[tt2] or
                            cds1[tt1] in cds2[tt2] or
                            cds2[tt2] in cds1[tt1]):
                        between_scores.append(score)
                        if options.loglevel >= 3:
                            outfile.write(
                                "# %s\t%s\t%6.4f\tsame sequence\n" % (tt1, tt2, score))
                    else:
                        if options.loglevel >= 3:
                            outfile.write(
                                "# %s\t%s\t%6.4f\tlink not found\n" % (tt1, tt2, score))
                        nmissed_links_between += 1
        else:
            if len(t1) != 1:
                nmissed_vertices += 1

    return between_scores, within_scores, nmissed_vertices, nmissed_links_within, nmissed_links_between

# -----------------------------------------------------------------------------


def IsLineageSpecificDuplication(outfile, t1, t2, graph, separator="|",
                                 cds1={}, cds2={},
                                 options={}):
    """check wether duplication of t2 is lineage specific.
    """

    # get scores between genes within a species and between species
    between_scores1, within_scores1, nmissed_vertices1, nmissed_links_within1, nmissed_links_between1 = \
        CountScores(outfile, t1, t2, graph, separator, cds1, cds2, options)
    between_scores2, within_scores2, nmissed_vertices2, nmissed_links_within2, nmissed_links_between2 = \
        CountScores(outfile, t2, t1, graph, separator, cds2, cds1, options)

    outfile.write("# missed: vertices1=%i, vertices2=%i, links_within1=%i, links_between1=%i, links_within2=%i,links_between2=%i\n" %
                  (nmissed_vertices1, nmissed_vertices2,
                   nmissed_links_within1, nmissed_links_between1,
                   nmissed_links_within2, nmissed_links_between2))

    between_scores = between_scores1 + between_scores2

    if len(between_scores) > 0:
        min_between = min(between_scores)
    else:
        min_between = 0

    if len(within_scores1) > 0:
        max_within1 = max(within_scores1)
    else:
        max_within1 = 0
    if len(within_scores2) > 0:
        max_within2 = max(within_scores2)
    else:
        max_within2 = 0

    outfile.write("# within1=%6.4f: %s\n" %
                  (max_within1, ";".join(map(str, within_scores1))))
    outfile.write("# within2=%6.4f: %s\n" %
                  (max_within2, ";".join(map(str, within_scores2))))
    outfile.write("# between=%6.4f: %s\n" %
                  (min_between, ";".join(map(str, between_scores))))

    return max_within1 < min_between, max_within2 < min_between

#-------------------------------------------------------------------------


def GetAlignedPairs(genes, cds):
    """get nucleotide differences for a set of genes.

    Each gene can be represented by several transcripts.

    Return min/max distance.
    """

    matrix, gop, gep = Genomics.makeSubstitutionMatrix("emboss")

    alignator = alignlib_lite.makeAlignatorDPFull(
        alignlib_lite.ALIGNMENT_LOCAL, gop, gep, matrix)

    #-------------------------------------------
    def MyAlignFunction(s1, s2, map_a2b):
        alignator.align(
            map_a2b, alignlib_lite.makeSequence(s1), alignlib_lite.makeSequence(s2))

    gg = genes.keys()

    map_a2b = alignlib_lite.makeAlignmentVector()

    pairs = {}

    for x in range(len(gg) - 1):
        for y in range(x + 1, len(gg)):
            for tx in genes[gg[x]]:
                for ty in genes[gg[y]]:

                    if tx not in cds or ty not in cds:
                        continue

                    p = AlignedPairs.AlignedPair(
                        sequence1=cds[tx], sequence2=cds[ty])
                    p.Align(MyAlignFunction, verbose=0)

                    if tx not in pairs:
                        pairs[tx] = {}
                    if ty not in pairs:
                        pairs[ty] = {}
                    pairs[tx][ty] = p
                    pairs[ty][tx] = p

    return pairs

# ------------------------------------------------------------------------


def CheckLocations(genes,
                   map_transcript2location,
                   do_strict=False,
                   map_contig2junk={},
                   map_contig2chromosome={},
                   max_local_duplication=100000,
                   separator="|"):
    """check the location for a set of genes and transcripts.

    local:              local gene cluster.
    non-local:          non-local gene clusters (but on the same contig).
    junk:               all but one gene lying on junk chromosomes.
    trans:              transposition to other chromosome
    muller:             transposition to other mueller element of chromosome
    """

    # count number of assignments per contig
    contigs = {}
    all_contigs = {}
    ngenes = 0
    map_gene2location = {}

    for g, vv in genes.items():

        max_to = 0
        min_from = 0
        for v in vv:

            if v.mTranscript not in map_transcript2location:
                sbjct_token = "dummy"
                sbjct_strand = "+"
                sbjct_genome_from = min_from
                sbjct_genome_to = max_to
                continue

            (sbjct_token, sbjct_strand, sbjct_genome_from,
             sbjct_genome_to) = map_transcript2location[v.mTranscript]

            if min_from == 0:
                min_from = sbjct_genome_from
            else:
                min_from = min(min_from, sbjct_genome_from)
            max_to = max(max_to, sbjct_genome_to)

        all_contigs[sbjct_token] = 1

        map_gene2location[g] = (sbjct_token, sbjct_strand, min_from, max_to)

        # ignore junk chromosomes
        if sbjct_token in map_contig2junk:
            continue

        if sbjct_token not in contigs:
            contigs[sbjct_token] = []

        contigs[sbjct_token].append((sbjct_strand, min_from, max_to))
        ngenes += 1

    # only one gene is left, all other must have been on junk chromosomes
    if len(genes) > 1 and ngenes <= 1:
        return "junk", all_contigs.keys(), map_gene2location

    if do_strict and ngenes < len(genes):
        return "junk", all_contigs.keys(), map_gene2location

    if ngenes == 0:
        return "junk", all_contigs.keys(), map_gene2location

    if ngenes == 1:
        return "single", all_contigs.keys(), map_gene2location

    status = "unknown"

    if len(contigs) > 1:
        #######################################################################
        # duplications to different chromosomes
        mapped = {}
        for x in contigs.keys():

            if x in map_contig2chromosome:
                z = map_contig2chromosome[x]
            else:
                z = x

            if z not in mapped:
                mapped[z] = 0
            mapped[z] += 1

        if len(mapped) == 1:
            status = "muller"
        else:
            status = "trans"
    else:
        #######################################################################
        # duplications on the same chromosome
        sbjct_token, residues = contigs.items()[0]

        residues.sort()
        distances = []
        # overlap on same chromosomes and strands
        for x in range(1, len(residues)):
            if residues[x][0] == residues[x - 1][0]:
                if residues[x][1] < residues[x - 1][2]:
                    status = "overlap"
                    break
                else:
                    distances.append(residues[x][1] - residues[x - 1][2])
        else:
            for d in distances:
                if d > max_local_duplication:
                    status = "nonlocal"
                    break
            else:
                status = "local"

    return (status, all_contigs.keys(), map_gene2location)

# ------------------------------------------------------------------------


def GetBestQuality(transcripts,
                   quality_priority=[]):
    """get best quality for a set of transcripts."""

    found = {}

    for v in transcripts:
        found[v.mQuality] = 1

    for x in quality_priority:
        if x in found:
            return x

    return ""

# ------------------------------------------------------------------------


def GetQualities(genes,
                 quality_priority=[]):
    """get best quality code from a set of genes."""
    qualities = []
    for g, transcripts in genes.items():

        q = GetBestQuality(transcripts, quality_priority)
        qualities.append(q)

    return qualities

# ------------------------------------------------------------------------


def UpdateCountsForCluster(outfile,
                           cluster_id, schema,
                           counts,
                           genes, transcripts,
                           trees,
                           map_transcript2location,
                           ortholog_distance,
                           options):
    """analyse duplications for a given cluster

    ortholog_distance: distance for normalization."""

    if options.loglevel >= 4:
        options.stdlog.write("# updating counts for cluster %i: %i genes and %i transcripts\n" %
                             (cluster_id, len(genes), len(transcripts)))
        options.stdlog.flush()

    if options.max_duplications and len(transcripts) > options.max_duplications:
        return "transposon", [], []

    # count number of contigs that genes are lying on
    status, locations, gene2location = CheckLocations(genes, map_transcript2location,
                                                      do_strict=False,
                                                      map_contig2junk=options.map_contig2junk,
                                                      map_contig2chromosome=options.map_contig2chromosome,
                                                      max_local_duplication=options.max_local_duplication,
                                                      separator=options.separator)

    qualities = GetQualities(genes, options.quality_priority)

    if status not in counts:
        counts[status] = 0

    counts[status] += 1
    counts["all"] += 1

    outfile.write("species\t%i\t%s\t%s\t%s\t%s\t%s\n" % (cluster_id, schema,
                                                         status,
                                                         ",".join(
                                                             genes.keys()),
                                                         ",".join(qualities),
                                                         ",".join(locations)))

    # print results: members of clusters
    for g, vv in genes.items():

        for v in vv:
            if v.mTranscript in map_transcript2location:
                (sbjct_token, sbjct_strand, sbjct_genome_from,
                 sbjct_genome_to) = map_transcript2location[v.mTranscript]
            else:
                sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to = "dummy", "0", 0, 0

            outfile.write("members\t%i\t%s\t%s\n" %
                          (cluster_id,
                           schema,
                           "\t".join([str(v), ] + map(str, (sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to)))))

    distances = []
    duplications = []

    str_transcripts = map(str, transcripts)

    for tree in trees:

        time_t1 = time.time()
        # check if tree is monophyletic for all transcripts in the genes:
        is_monophyletic = TreeTools.IsMonophyleticForTaxa(
            tree, str_transcripts)

        if options.loglevel >= 4:
            options.stdlog.write(
                "# checked monophyly in %5.2f seconds.\n" % (time.time() - time_t1))
            options.stdlog.flush()

        if is_monophyletic:

            time_t1 = time.time()

            branchpoints = TreeTools.CountBranchPoints(tree, str_transcripts)

            if options.loglevel >= 4:
                options.stdlog.write(
                    "# counted branchpoints in %5.2f seconds.\n" % (time.time() - time_t1))
                options.stdlog.flush()

            if branchpoints:

                for str_children, height, branchlength in branchpoints:

                    if len(str_children) <= 1:
                        continue

                    children = map(
                        lambda x: Orthologs.Transcript(x), str_children)

                    distances.append(height)

                    is_pseudogene = False

                    for child in children:
                        if child.mQuality in options.pseudogenes:
                            is_pseudogene = True
                            break

                    if is_pseudogene:
                        duplication_status = "pseudo"
                    else:
                        duplication_status = "functional"

                    this_genes = Orthologs.GetGenes(children)
                    this_status, this_locations, this_gene2location = CheckLocations(this_genes,
                                                                                     map_transcript2location,
                                                                                     do_strict=True,
                                                                                     map_contig2junk=options.map_contig2junk,
                                                                                     map_contig2chromosome=options.map_contig2chromosome,
                                                                                     max_local_duplication=options.max_local_duplication,
                                                                                     separator=options.separator)

                    temp_tree = copy.deepcopy(tree)
                    TreeTools.PruneTree(temp_tree, str_children)

                    duplications.append(
                        (cluster_id, this_status, duplication_status, height, children))

                    if ortholog_distance > 0:
                        rel_height = "%6.4f" % (height / ortholog_distance)
                    else:
                        rel_height = "NaN"

                    locs = []
                    for child in children:
                        if child.mGene in this_gene2location:
                            locs.append(
                                "%s:%s:%s:%i:%i" % ((child.mGene,) + this_gene2location[child.mGene]))

                    outfile.write("branchpoint\t%i\t%s\t%s\t%s\t%f\t%s\t%s\t%s\t%s\t%s\n" % (cluster_id,
                                                                                             schema,
                                                                                             this_status,
                                                                                             duplication_status,
                                                                                             height,
                                                                                             rel_height,
                                                                                             ";".join(
                                                                                                 this_locations),
                                                                                             ";".join(
                                                                                                 str_children),
                                                                                             ";".join(
                                                                                                 locs),
                                                                                             TreeTools.Tree2Newick(
                                                                                                 temp_tree),
                                                                                             ))

                if ortholog_distance > 0:
                    rel_heights = map(
                        lambda x: x / ortholog_distance, distances)
                else:
                    rel_heights = "NaN"

                outfile.write("branchpoints\t%i\t%s\t%s\t%s\t%s\n" % (cluster_id,
                                                                      schema,
                                                                      status,
                                                                      ";".join(
                                                                          map(str, distances)),
                                                                      ";".join(map(str, rel_heights))))

    return status, distances, duplications

# ------------------------------------------------------------------------


def PrintResultsDuplicationsDistances(outfile,
                                      categories,
                                      histogram_data,
                                      options):
    """write histograms of duplication distances."""

    ###################################
    # construct and write histograms
    num_bins = 100
    bins = map(lambda x: float(x) / 20.0, range(0, num_bins))
    histograms1 = {}
    histograms2 = {}
    vals0 = []
    vals1 = []
    for key, vals in histogram_data.items():
        if key not in categories:
            continue
        h = scipy.stats.histogram2(vals[0], bins)
        histograms1[key] = h
        h = scipy.stats.histogram2(vals[1], bins)
        histograms2[key] = h
        vals0 += vals[0]
        vals1 += vals[1]

    h0 = scipy.stats.histogram2(vals0, bins)
    h1 = scipy.stats.histogram2(vals1, bins)

    outfile.write("# duplications - all histograms for %s and %s\n" %
                  (options.schema1, options.schema2))
    outfile.write("bin\t('sum','sum')\t%s\n" %
                  "\t\t".join(map(str, categories)))
    for b in range(0, num_bins):
        outfile.write("%5.2f" % bins[b])
        outfile.write("\t%i\t%i" % (h0[b], h1[b]))
        for x in categories:
            if x in histograms1 and x in histograms2:
                outfile.write("\t%i\t%i" %
                              (histograms1[x][b], histograms2[x][b]))
            else:
                outfile.write("\t0\t0")

        outfile.write("\n")
    outfile.write("total")
    outfile.write(
        "\t%i\t%i" % (reduce(lambda x, y: x + y, h0), reduce(lambda x, y: x + y, h0)))
    for x in categories:
        if x in histograms1 and x in histograms2:
            outfile.write("\t%i\t%i" % (reduce(lambda x, y: x + y, histograms1[x]),
                                        reduce(lambda x, y: x + y, histograms2[x])))
        else:
            outfile.write("\t0\t0")
    outfile.write("\n")


def GetLocationKeys():
    return ("local", "nonlocal", "muller", "trans", "junk", "unknown", "overlap", "single")


def GetLocationHash():
    """return hash with location counts initialized to 0."""
    h = {}
    for x in GetLocationKeys():
        h[x] = 0
    return h


def GetFunctionKeys():
    return ("functional", "pseudo")


def GetFunctionHash():
    """return hash with location counts initialized to 0."""
    h = {}
    for x in GetFunctionKeys():
        h[x] = 0
    return h

# ------------------------------------------------------------------------


def PrintResultsDuplicationsPairs(outfile, cluster_id, schema, members, duplications, ortholog_distance):
    """print pairs of duplicated genes."""

    results_locations = GetLocationHash()
    results_functions = GetFunctionHash()

    for cluster_id, location_status, duplication_status, height, children in duplications:
        results_locations[location_status] += 1
        results_functions[duplication_status] += 1

    for member in members:
        outfile.write("%s\t%i\t%s\t%i\t%s\t%s\t%s\n" %
                      ("pairs",
                       cluster_id, schema, len(members), str(member),
                       "\t".join(
                           map(str, [results_locations[x] for x in GetLocationKeys()])),
                       "\t".join(map(str, [results_functions[x] for x in GetFunctionKeys()]))))

# ------------------------------------------------------------------------


def PrintResultsDuplicationsType(outfile, results):

    kk = results.keys()
    kk.sort()
    outfile.write("\t".join(kk) + "\n")
    outfile.write("\t".join(map(lambda x: str(results[x]), kk)) + "\n")


# ------------------------------------------------------------------------
def AnalyseDuplications(outfile,
                        orthologs,
                        trees,
                        graph_genes, graph_kaks,
                        map_transcript2location1, map_transcript2location2,
                        schema1, schema2,
                        max_local_duplication=100000,
                        separator="|",
                        tablename_predictions="predictions",
                        map_contig2chromosome={},
                        map_contig2junk={},
                        cds1={}, cds2={},
                        quality_priority=[],
                        options={}):

    outfile.write("pairs\tcluster_id\tschema\tnmembers\tid\t%s\t%s\n" %
                  ("\t".join(GetLocationKeys()), "\t".join(GetFunctionKeys())))

    outfile.write("## duplications of %s to %s\n" % (schema1, schema2))

    lines = []

    results1 = {"all": 0}
    results2 = {"all": 0}
    all_duplications = []
    histogram_data = {("all", "all"): [[], []]}
    nskipped = 0
    cluster_id = 0
    ninput = 0
    nstrict = 0
    ndegenerate = 0
    for cluster_id in range(len(orthologs)):

        transcripts1, transcripts2, g1, g2, weights = orthologs[cluster_id]

        t1 = map(str, transcripts1)
        t2 = map(str, transcripts2)

        ninput += 1

        # skip 1:1 clusters
        if len(g1) == 1 and len(g2) == 1:
            nstrict += 1
            continue

        ndegenerate += 1

        outfile.write("##################################################\n")
        outfile.write("summary\t%i\t%i\n" %
                      (cluster_id, len(trees[cluster_id])))

        ######################################################################
        # calculate distances between orthologs
        # and build subtrees for quicker processing.

        time_t1 = time.time()

        orthologs_distances = []
        subtrees = []

        for tree in trees[cluster_id]:

            # tree.display()
            # print t1 + t2

            node = TreeTools.Reroot(tree, t1 + t2)
            subtree = TreeTools.GetSubtree(tree, node)
            subtrees.append(subtree)
            # subtree.display()
            # print "distances1", TreeTools.GetDistancesBetweenTaxa( subtree, t1, t2 )
            # print "distances2", TreeTools.GetDistancesBetweenTaxa( tree, t1,
            # t2 )
            orthologs_distances += TreeTools.GetDistancesBetweenTaxa(
                subtree, t1, t2)

        if options.loglevel >= 4:
            options.stdlog.write(
                "# collected distances between orthologs in %5.2f seconds.\n" % (time.time() - time_t1))
            options.stdlog.flush()

        orthologs_distances = map(lambda x: x[2], orthologs_distances)

        if len(orthologs_distances) == 0:
            outfile.write(
                "# skipped %i: no ortholog distances.\n" % (cluster_id))
            nskipped += 1
            continue

        ortholog_distance = numpy.mean(orthologs_distances)
        outfile.write("orthologs\t%i\t%i\t%6.4f\t%6.4f\t%6.4f\n" %
                      (len(orthologs_distances),
                       cluster_id,
                       ortholog_distance,
                       numpy.std(orthologs_distances),
                       numpy.median(orthologs_distances)))

        ######################################################################
        # get branchpoints and update counts
        time_t1 = time.time()
        status1, distances1, duplications1 = UpdateCountsForCluster(outfile,
                                                                    cluster_id, options.schema1,
                                                                    results1,
                                                                    g1, t1,
                                                                    subtrees,
                                                                    map_transcript2location1,
                                                                    ortholog_distance,
                                                                    options)

        if options.loglevel >= 4:
            options.stdlog.write("# %i: UpdateCountsForCluster1: %5.2f seconds.\n" % (
                cluster_id, time.time() - time_t1))
            options.stdlog.flush()

        time_t1 = time.time()
        status2, distances2, duplications2 = UpdateCountsForCluster(outfile,
                                                                    cluster_id, options.schema2,
                                                                    results2,
                                                                    g2, t2,
                                                                    subtrees,
                                                                    map_transcript2location2,
                                                                    ortholog_distance,
                                                                    options)

        if options.loglevel >= 4:
            options.stdlog.write("# %i: UpdateCountsForCluster1: %5.2f seconds.\n" % (
                cluster_id, time.time() - time_t1))
            options.stdlog.flush()

        if len(distances1) == 0 and len(distances2) == 0:
            outfile.write("# skipped %i : %s=%i %s=%i\n" % (
                cluster_id, status1, len(distances1), status2, len(distances2)))
            nskipped += 1
            continue

        key = (status1, status2)

        if key not in histogram_data:
            histogram_data[key] = [[], []]

        histogram_data[key][0] += distances1
        histogram_data[key][1] += distances2

        histogram_data[("all", "all")][0] += distances1
        histogram_data[("all", "all")][1] += distances2

        time_t1 = time.time()
        PrintResultsDuplicationsPairs(
            outfile, cluster_id, schema1, transcripts1, duplications2, ortholog_distance)
        PrintResultsDuplicationsPairs(
            outfile, cluster_id, schema2, transcripts2, duplications1, ortholog_distance)

        if options.loglevel >= 4:
            options.stdlog.write(
                "# %i: Output: %5.2f seconds.\n" % (cluster_id, time.time() - time_t1))
            options.stdlog.flush()

    outfile.write("#" * 30 + "\n")
    outfile.write("# summary table for duplication type in %s: %s to %s\n" % (
        schema1, schema1, schema2))
    PrintResultsDuplicationsType(outfile, results1)

    outfile.write("#" * 30 + "\n")
    outfile.write("# summary table for duplication type in %s: %s to %s\n" % (
        schema2, schema1, schema2))
    PrintResultsDuplicationsType(outfile, results2)

    ###################################
    # duplication times
    kk = histogram_data.keys()
    kk.sort()
    for key in kk:
        status1, status2 = key
        outfile.write("duplication\t%s\t%s\t%s\t%s\t%i\t%i\t%s\t%s\n" %
                      (schema1, schema2,
                       status1, status2,
                       len(histogram_data[key][0]),
                       len(histogram_data[key][1]),
                       ",".join(
                           map(lambda x: "%6.4f" % x, histogram_data[key][0])),
                       ",".join(
                           map(lambda x: "%6.4f" % x, histogram_data[key][1])),
                       ))

    outfile.write("# histogram of duplications: all categories: %s to %s, skipped=%i\n" % (
        options.schema1, options.schema2, nskipped))
    PrintResultsDuplicationsDistances(outfile, kk, histogram_data, options)

    outfile.write("# histogram of duplications: selected categories: %s to %s, skipped=%i\n" % (
        options.schema1, options.schema2, nskipped))
    kk = (("single", "local"), ("local", "single"), ("local", "local"))
    PrintResultsDuplicationsDistances(outfile, kk, histogram_data, options)

    outfile.write("# input=%i, direct=%i, degenerate=%i, skipped=%i\n" %
                  (ninput, nstrict, ndegenerate, nskipped))

# ------------------------------------------------------------------------


def GetPrediction2Location(dbhandle, schema,
                           tablename_predictions="overview",
                           tablename_genes="genes",
                           use_genes=False):

    map_transcript2location = {}

    if use_genes:
        statement = """
        SELECT g.gene_id, p.sbjct_token, p.sbjct_strand,
        MIN(export_sbjct_genome_from), MAX(export_sbjct_genome_to)
        FROM %s.%s AS p, %s.%s AS g
        WHERE p.prediction_id = g.prediction_id
        GROUP BY g.gene_id, p.sbjct_token, p.sbjct_strand
        """ % \
            (schema, tablename_predictions,
             schema, tablename_genes)
    else:
        statement = """
        SELECT prediction_id, sbjct_token, sbjct_strand,
        export_sbjct_genome_from, export_sbjct_genome_to
        FROM %s.%s 
        """ % \
            (schema, tablename_predictions)

    cc = dbhandle.cursor()
    cc.execute(statement)
    rr = cc.fetchall()
    cc.close()

    for prediction_id, sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to in rr:
        map_transcript2location[str(prediction_id)] = (
            sbjct_token, sbjct_strand, sbjct_genome_from, sbjct_genome_to)

    return map_transcript2location

# ------------------------------------------------------------------------


def GetRepeats(dbhandle, schema, repeats_queries):
    """build a list of predictions which are repeats."""

    ##########################################################################
    # get all best matches
    statement = """
    SELECT DISTINCT prediction_id
    FROM %s.predictions
    WHERE query_token IN ('%s')""" % \
        (schema, "','".join(repeats_queries))

    cc = dbhandle.cursor()
    cc.execute(statement)
    rr = cc.fetchall()
    cc.close()

    repeats_list = {}

    for r in rr:
        repeats_list[r[0]] = 1

    return repeats_list

# ------------------------------------------------------------------------


def GetBestMatches(dbhandle, schema):

    ##########################################################################
    # get all best matches
    statement = """
    SELECT prediction_id,
    gene_id, query_token, class,
    query_coverage, pidentity, 
    sbjct_token, sbjct_strand,
    full_sbjct_genome_from, full_sbjct_genome_to
    FROM %s.overview
    WHERE is_best_prediction = TRUE
    ORDER BY query_token ASC, query_coverage * pidentity DESC""" % \
        (schema)

    cc = dbhandle.cursor()
    cc.execute(statement)
    rr = cc.fetchall()
    cc.close()

    best = 0
    best_match = None
    # get best match for all possible genes

    map_query2best = {}
    best_pide, best_coverage, best_identity, nmatches = 0, 0, 0, 0
    last_query_token = None
    for r in rr:

        (prediction_id, gene_id,
         query_token, quality,
         query_coverage, pidentity,
         sbjct_token, sbjct_strand,
         sbjct_genome_from, sbjct_genome_to,
         ) = r

        if last_query_token != query_token:
            if best_match:
                map_query2best[last_query_token] = (
                    nmatches, best_pide, best_coverage) + tuple(best_match)

            last_query_token = query_token
            best, best_pide, best_coverage, best_identity, nmatches = 0, 0, 0, 0, 0
            best_match = None

        v = query_coverage * pidentity

        if v > best:
            best = v
            best_match = r[:]
            best_match[0] = str(best_match[0])
            best_match[1] = str(best_match[1])

        if query_coverage > best_coverage:
            best_coverage = query_coverage
        if pidentity > best_identity:
            best_identity = pidentity
        nmatches += 1

    if best_match:
        map_query2best[last_query_token] = (
            nmatches, best_pide, best_coverage) + tuple(best_match)

    return map_query2best

# ------------------------------------------------------------------------


def WriteOrthologs(outfile, schema, genes, transcripts):
    """write list of ortholog status for individual genes/transcripts."""

    for g, s in genes:
        outfile.write("%s\tgene\t%s\t%s\n" % (schema, str(g), s))

    for t, s in transcripts:
        outfile.write("%s\ttranscript\t%s\t%s\n" % (schema, str(t), s))

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


def ReadOrphans(infile, options):
    """read orphans from a file.
    """

    ninput, nskipped = 0, 0

    orphans1, orphans2 = {}, {}
    for line in infile:

        if line[0] == "#":
            continue
        elif line[0] == ">":
            cluster_id = re.match(">cluster# (\d+)", line[:-1]).groups()[0]
            continue
        elif line[0] == "\n":
            continue

        ninput += 1

        transcripts = map(
            lambda x: Orthologs.Transcript(x), line[:-1].strip().split("\t"))

        if len(transcripts) == 0:
            nskipped += 1
            continue

        transcripts1 = []
        transcripts2 = []

        # sort transcripts by schema
        for t in transcripts:
            if t.mSchema == options.schema1:
                transcripts1.append(t)
            elif t.mSchema == options.schema2:
                transcripts2.append(t)

        def addToOrphans(orphans, transcripts):
            genes = Orthologs.GetGenes(transcripts)

            for g in genes:
                if g not in orphans:
                    orphans[g] = []
                orphans[g] += genes[g]

        addToOrphans(orphans1, transcripts1)
        addToOrphans(orphans2, transcripts2)

    if options.loglevel >= 1:
        options.stdlog.write("# orphans: read %i/%i orphaned genes: ninput=%i, nskipped=%i\n" % (len(orphans1),
                                                                                                 len(orphans2),
                                                                                                 ninput, nskipped))

    return orphans1, orphans2

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/analyze_orthology_pairwise.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-e", "--schema", dest="schema", type="string",
                      help="schema.")

    parser.add_option("--test", dest="test", type="int",
                      help="test run - stop after reading # links.")

    parser.add_option("--filename-cds", dest="filename_cds", type="string",
                      help="input filename with sequence information - determines the input set.")

    parser.add_option("--filename-kaks", dest="filename_kaks", type="string",
                      help="input filename with kaks results.")

    parser.add_option("--filename-blast", dest="filename_blast", type="string",
                      help="input filename with blast results.")

    parser.add_option("--filename-bitscore-distances", dest="filename_bitscore_distances", type="string",
                      help="filename with bitscore distances.")

    parser.add_option("--filename-orphans", dest="filename_orphans", type="string", action="append",
                      help="input filename with orphans from orthology assignment procedure.")

    parser.add_option("--filename-orthologs", dest="filename_orthologs", type="string", action="append",
                      help="input filename with orthologs from orthology assignment procedure.")

    parser.set_defaults(
        schema="dmel_vs_dmel2",
        filename_cds=None,
        filename_blast=None,
        filename_bitscore_distances=None,
        filename_orphans=[],
        filename_orthologs=[],
        threshold_max_bitdist=0.80,
        threshold_max_ks=5.0,
        threshold_filter_classes="CG,SG,PG,RG,CP,SP,PP",
        test=None)

    (options, args) = E.Start(parser)

    if options.threshold_filter_classes:
        options.threshold_filter_classes = options.threshold_filter_classes.split(
            ",")

    if len(options.filename_orthologs) != len(options.filename_orphans):
        raise "please specify the same number of files for orhtologs and orphans."

    ##########################################################################
    # get identifiers from fasta file
    all_identifiers = Genomics.ReadPeptideSequences(
        open(options.filename_cds, "r"))
    all_transcripts = (all_identifiers.keys())

    all_genes = {}
    for id in all_transcripts:
        t = Orthologs.Transcript(id)
        if t.mGene not in all_genes:
            all_genes[t.mGene] = []
        all_genes[t.mGene].append(t)

    if options.loglevel >= 1:
        options.stdlog.write("# read %i transcripts and %i genes.\n" % (
            len(all_transcripts), len(all_genes)))

    ##########################################################################
    # check how many have matches in BLAST
    if options.loglevel >= 1:
        options.stdlog.write("# reading blast information ...")
        options.stdlog.flush()

    # record best match in other genome
    best_matches_blast = {}
    ninput = 0
    for line in gzip.open(options.filename_blast, "r"):
        if line[0] == "#":
            continue
        ninput += 1
        link = BlastAlignments.Link()
        link.Read(line)
        t1 = Orthologs.Transcript(link.mQueryToken)
        t2 = Orthologs.Transcript(link.mSbjctToken)

        if t1.mSchema == options.schema and t2.mSchema != options.schema:
            pass
        elif t2.mSchema == options.schema and t1.mSchema != options.schema:
            t1, t2 = t2, t1
            link.Reverse()
        else:
            continue

        if t1.mGene not in best_matches_blast:
            best_matches_blast[t1.mGene] = {'coverage': 0, 'pid': 0}

        b = best_matches_blast[t1.mGene]
        coverage = 100.0 * \
            (link.mQueryTo - link.mQueryFrom + 1) / link.mQueryLength
        if b["coverage"] < coverage:
            b["coverage"] = coverage
        if b['pid'] < link.mPercentIdentity:
            b['pid'] = link.mPercentIdentity

        if options.test and ninput > options.test:
            break

    if options.loglevel >= 1:
        options.stdlog.write(" finished.\n")
        options.stdlog.flush()

    if options.loglevel >= 1:
        options.stdlog.write(
            "# read best matches for %i genes.\n" % (len(best_matches_blast)))

    ##########################################################################
    # read bitscore alignments
    if options.loglevel >= 1:
        options.stdlog.write("# reading bitscore information ...")
        options.stdlog.flush()

    # record best match in other genome
    best_matches_bitdist = {}

    ninput = 0
    for line in gzip.open(options.filename_bitscore_distances, "r"):
        if line[0] == "#":
            continue
        ninput += 1

        query_token, sbjct_token, distance = line[:-1].split("\t")[:3]

        distance = float(distance)

        t1 = Orthologs.Transcript(query_token)
        t2 = Orthologs.Transcript(sbjct_token)

        if t1.mSchema == options.schema and t2.mSchema != options.schema:
            pass
        elif t2.mSchema == options.schema and t1.mSchema != options.schema:
            t1, t2 = t2, t1
            link.Reverse()
        else:
            continue

        if t1.mGene not in best_matches_bitdist:
            best_matches_bitdist[t1.mGene] = 100.0

        if best_matches_bitdist[t1.mGene] > distance:
            best_matches_bitdist[t1.mGene] = distance

        if options.test and ninput > options.test:
            break

    if options.loglevel >= 1:
        options.stdlog.write(" finished.\n")
        options.stdlog.flush()

    if options.loglevel >= 1:
        options.stdlog.write(
            "# read bitscore distances for %i genes.\n" % (len(best_matches_bitdist)))

    ##########################################################################
    # read kaks information
    if options.loglevel >= 1:
        options.stdlog.write("# reading kaks information ...")
        options.stdlog.flush()

    # record best match in other genome
    best_matches_ks = {}
    best_matches_ka = {}

    ninput = 0
    for line in open(options.filename_kaks, "r"):
        if line[0] == "#":
            continue
        if line.startswith("seq1\tseq2"):
            continue
        ninput += 1
        query_token, sbjct_token, ka, ks = line[:-1].split("\t")[:4]
        ka, ks = float(ka), float(ks)

        t1 = Orthologs.Transcript(query_token)
        t2 = Orthologs.Transcript(sbjct_token)

        if t1.mSchema == options.schema and t2.mSchema != options.schema:
            pass
        elif t2.mSchema == options.schema and t1.mSchema != options.schema:
            t1, t2 = t2, t1
            link.Reverse()
        else:
            continue

        if t1.mGene not in best_matches_ks:
            best_matches_ks[t1.mGene] = 999.0

        if best_matches_ks[t1.mGene] > ks:
            best_matches_ks[t1.mGene] = ks

        if t1.mGene not in best_matches_ka:
            best_matches_ka[t1.mGene] = 999.0

        if best_matches_ka[t1.mGene] > ka:
            best_matches_ka[t1.mGene] = ka

        if options.test and ninput > options.test:
            break

    if options.loglevel >= 1:
        options.stdlog.write(" finished.\n")
        options.stdlog.flush()

    if options.loglevel >= 1:
        options.stdlog.write(
            "# read ka and ks distances for %i genes.\n" % (len(best_matches_ks)))

    genes = all_genes.keys()
    genes.sort()

    ##########################################################################
    # read orphan information
    orphans = []
    for filename_orphans in options.filename_orphans:
        oo = {}
        for line in open(filename_orphans, "r"):
            if line[0] == "#":
                continue
            if line[0] == ">":
                continue
            tt = map(lambda x: Orthologs.Transcript(x), line[:-1].split("\t"))
            for t in tt:
                if t.mSchema != options.schema:
                    continue
                if t.mGene not in oo:
                    oo[t.mGene] = []
                oo[t.mGene].append(t)
        orphans.append(oo)

        if options.loglevel >= 1:
            options.stdlog.write(
                "# read orphans for %i genes from %s.\n" % (len(oo), filename_orphans))

    ##########################################################################
    # read ortholog information
    orthologs = []
    for filename_orthologs in options.filename_orthologs:
        oo = {}
        for line in open(filename_orthologs, "r"):
            if line[0] == "#":
                continue
            if line[0] == ">":
                continue
            tt = map(
                lambda x: Orthologs.Transcript(x), line[:-1].split("\t")[1:])
            for t in tt:
                if t.mSchema != options.schema:
                    continue
                if t.mGene not in oo:
                    oo[t.mGene] = []
                oo[t.mGene].append(t)
        orthologs.append(oo)

        if options.loglevel >= 1:
            options.stdlog.write(
                "# read orthologs for %i genes from %s.\n" % (len(oo), filename_orthologs))

    options.stdout.write("""# Index for stage information.
# 1: passed BLAST       - gene has BLAST information.
# 2: passed ALI         - gene has passed alignment quality threshold.
# 3: passed CLASS       - gene has passed class filter.
# 4: passed ORTHOLOG_KA - gene is assigned to orthologs on ka level.
# 5: passed KAKS        - gene has KAKS information.
# 6: passed ORTHOLOG_KS - gene is assigned to orthologs on ks level.
""")

    stages = {}
    for g in genes:
        stage = 0
        if g in best_matches_blast:
            # passed blast
            stage += 1
            if g in best_matches_bitdist and \
               best_matches_bitdist[g] <= options.threshold_max_bitdist:
                # passed bitdist
                stage += 1
                found = False
                for t in all_genes[g]:
                    if t.mQuality in options.threshold_filter_classes:
                        found = True
                        break
                if found:
                    # passed bitdist
                    stage += 1
                    if g in orthologs[0]:
                        # passed orthologs_ka
                        stage += 1
                        if g in best_matches_ks and best_matches_ks[g] < options.threshold_max_ks:
                            # passed kaks
                            stage += 1
                            if g in orthologs[1]:
                                # passed orthologs_ks
                                stage += 1
        stages[g] = stage

    #################################################################
    #################################################################
    #################################################################
    # Compute stage information
    options.stdout.write(
        "gene\tntranscripts\tstage\tblast_pid\tblast_cov\tbitdist\tka\tks")

    for x in range(len(options.filename_orphans)):
        options.stdout.write("\torphans%i\torthologs%i" % (x, x))

    options.stdout.write("\n")

    for g in genes:

        blast_pid, blast_cov, bitdist, ks, ka = "na", "na", "na", "na", "na"
        ntranscripts = "%i" % len(all_genes[g])
        stage = "%i" % stages[g]
        if g in best_matches_blast:
            blast_pid, blast_cov = "%5.2f" % best_matches_blast[g][
                "pid"], "%5.2f" % best_matches_blast[g]["coverage"]

        if g in best_matches_bitdist:
            bitdist = "%5.2f" % best_matches_bitdist[g]
        if g in best_matches_ka:
            ka = "%6.4f" % best_matches_ka[g]
        if g in best_matches_ks:
            ks = "%6.4f" % best_matches_ks[g]

        options.stdout.write("\t".join((g, ntranscripts, stage, blast_pid, blast_cov, bitdist,
                                        ka, ks)))

        for x in range(len(options.filename_orphans)):
            norphans, northologs = "na", "na"
            if g in orphans[x]:
                norphans = "%i" % len(orphans[x][g])
            if g in orthologs[x]:
                northologs = "%i" % len(orthologs[x][g])

            options.stdout.write("\t%s\t%s" % (norphans, northologs))

        options.stdout.write("\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

################################################################################
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
#################################################################################
'''
mali_evaluate.py - 
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

   python mali_evaluate.py --help

Type::

   python mali_evaluate.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import getopt
import math

USAGE="""python %s [OPTIONS] < exonerate_output > filtered

Evaluate a multiple alignment. Remove gaps if so wished.

Version = $Id: mali_evaluate.py 2782 2009-09-10 11:40:29Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-o, --output-pattern            pattern for output multiple alignments
-p, --master-pattern=           pattern identifying sequences of the query species from identifier
-s, --species-pattern=          pattern identifying species from identifier
-e, --exons=                    exon information for evaluation
-c, --cluster                   cluster by overlap
-f, --remove-fragments          remove fragments
-p, --prefix=                   prefix for cluster identifier
--min-overlap-percent           minimum percent overlap between sequence paris
--min-overlap-residues          minimum overlap in residues between sequence pairs
--min-coverage-percent          minimum percent coverage between sequence paris
--min-coverage-residues         minimum coverage in residues between sequence pairs
--components                    filename with components to be analyses separately in the multiple alignment
""" % sys.argv[0]

param_long_options=["verbose=", "help", "file-output=",
                    "min-overlap-residues=", "min-overlap-percent=",
                    "min-coverage-residues=", "min-coverage-percent=",
                    "components=",
                    "species-pattern=", "master-pattern=", "output-pattern=",
                    "exons=", "cluster", "remove-fragments", "prefix=",
                    "version"
                    ]

param_short_options="v:ho:s:p:e:"

param_loglevel = 1

param_gap_char = "-"
param_mask_char = "x"

param_output_pattern = None

param_filename_exons = None
param_filename_components = None

param_min_overlap_residues = 10
param_min_overlap_percent  = 0.5

param_min_coverage_residues = 10
param_min_coverage_percent  = 0.5


param_max_exons_difference = 1

param_species_pattern = "^([^|:@]+)[^|:@]"
param_master_pattern = None

param_threshold_splipping_exon_boundary = 9
param_evaluate_single_exons_min_coverage = 80
param_evaluate_min_percent_exon_identity = 70

## clustering parameters
param_clustering_max_distance = 20
param_max_compatibility = 50

## permissiveness for merging clusters, set to < 0 for disallowing merging
param_cluster_merging_max_distance = 2

param_do_cluster = False

param_remove_fragments = False

param_pattern_prefix = ""

import CGAT.Experiment as E
import CGAT.MaliIO as MaliIO
import scipy
import CGAT.Exons as Exons
import alignlib
import CGAT.Genomics as Genomics
import numpy

##--------------------------------------------------------------------------------------- 
def GetOverlapMatrix( mali, identifiers, gap_char = "-" ):
    """get overlap matrix for all pairwise combinations in mali.
    """
    wmali = len(identifiers)
    lmali = len(mali[identifiers[0]])

    coverages = numpy.zeros( wmali )

    matrix_overlap = numpy.zeros( (wmali, wmali) )
    matrix_total = numpy.ones( (wmali, wmali) )    
    
    for x in range(wmali):
        a = mali[identifiers[x]]

        ncoverage = 0
        for i in range(lmali):
            if a[i] != param_gap_char:
                  ncoverage += 1
                
        coverages[x] = ncoverage 

        for y in range(x+1, wmali):
            
            noverlap = 0
            
            b = mali[identifiers[y]]
            ntotal = 0
            for i in range(lmali):
                if a[i] != gap_char and b[i] != gap_char and \
                       a[i] in string.uppercase and b[i] in string.uppercase:
                    noverlap += 1
                if a[i] != gap_char or b[i] != gap_char:
                    ntotal += 1
                    
            matrix_total[x][y] = ntotal
            matrix_total[y][x] = ntotal
            matrix_overlap[x][y] = noverlap
            matrix_overlap[y][x] = noverlap

    return matrix_overlap, matrix_total, coverages

##--------------------------------------------------------------------------------------- 
def WriteOverlap( mali, identifiers, prefix = "" ):

    """check multiple alignment for overlap.

    This subroutine dumps out various summary parameters.
    """

    if len(identifiers) <= 1: return []
    
    wmali = len(identifiers)
    lmali = len(mali[identifiers[0]])

    ## check for non-overlapping pairs
    nnon_overlapping = 0
    nfragments = 0
    pairs = []
    ppairs = []
    coverages = []
    fragments = []

    for x in range(wmali):
        a = mali[identifiers[x]]

        ncoverage = 0
        for i in range(lmali):
            if a[i] != param_gap_char:
                  ncoverage += 1
                
        coverages.append( ncoverage )

        if ncoverage < param_min_coverage_residues or (float(ncoverage) / lmali) < param_min_overlap_percent:
            if param_loglevel == 3:
                print "# fragment: %s\t%i\t%i\t%.2f" % (identifiers[x], ncoverage, lmali, float(ncoverage) / lmali )
            
            nfragments += 1
            fragments.append( identifiers[x] )
                
        for y in range(x+1, wmali):
            
            if x == y: continue
            
            noverlap = 0
            
            b = mali[identifiers[y]]
            ntotal = 0
            for i in range(lmali):
                if a[i] != param_gap_char and b[i] != param_gap_char:
                    noverlap += 1
                if a[i] != param_gap_char or b[i] != param_gap_char:
                    ntotal += 1
                    
            pairs.append( noverlap )
            ppairs.append( float(noverlap) / ntotal )
            if param_loglevel >= 5:
                print "# overlap: %s\t%s\t%i\t%i\t%.2f" % (identifiers[x], identifiers[y], noverlap, ntotal, float(noverlap) / ntotal )
            
            if noverlap < param_min_overlap_residues or (float(noverlap) / ntotal) < param_min_overlap_percent:
                nnon_overlapping += 1
            
                if param_loglevel == 3:
                    print "# failed overlap: %s\t%s\t%i\t%i\t%.2f" % (identifiers[x], identifiers[y], noverlap, ntotal, float(noverlap) / ntotal )

    pcoverages = map( lambda x: float(x) / lmali, coverages)    

    print "%s\tfailed\t" % prefix +\
          string.join(map( str, (nfragments, wmali,
                                 "%.2f" % (float(nfragments) / wmali),                                               
                                 nnon_overlapping, wmali * (wmali - 1) / 2,
                                 "%.2f" % (float(nnon_overlapping) / (wmali * (wmali - 1) / 2))
                                 )),"\t")

    print "%s\tpairs\t" % prefix +\
          string.join(map(str,
                          (min(pairs), max(pairs), "%.2f" % scipy.mean(pairs), scipy.median(pairs), "%.2f" % scipy.std(pairs))), "\t")
                                 
    print "%s\tppairs\t" % prefix +\
          string.join(map(lambda x: "%.2f" % x, (min(ppairs),
                                                 max(ppairs),
                                                 scipy.mean(ppairs),
                                                 scipy.median(ppairs),
                                                 scipy.std(ppairs))),"\t")
    print "%s\tcov\t" % prefix +\
          string.join( map(str, 
                           (min(coverages),
                            max(coverages),
                            "%.2f" % scipy.mean(coverages),
                            scipy.median(coverages),
                            "%.2f" % scipy.std(coverages))), "\t")
                                 
    print "%s\tpcov\t" % prefix +\
          string.join(map(lambda x: "%.2f" % x, (min(pcoverages),
                                                 max(pcoverages),
                                                 scipy.mean(pcoverages),
                                                 scipy.median(pcoverages),
                                                 scipy.std(pcoverages))), "\t")

    return fragments

##---------------------------------------------------------------------------------------
def GetIntronPositions( seq, gap_char = "-" ):
    """get intron positions."""

    last_upper = None

    transitions = []
    last_char = 0
    
    for x in range(len(seq)):
        
        if seq[x] == gap_char: continue

        last_char = x
        
        if last_upper == None:
            transitions.append(x)
        elif seq[x] in string.uppercase and not last_upper:
            transitions.append(x)
        elif seq[x] in string.lowercase and last_upper:
            transitions.append(x)
            
        last_upper = seq[x] in string.uppercase

    transitions.append( last_char )

    return transitions

##---------------------------------------------------------------------------------------
def MyMap( map, r, distance = 3, offset = 0):
    """returns
    -1, if no mapping possible with distance of 3 residues.
    0:  if mapping is out of bounds
    """
    if r < map.getRowFrom() or r > map.getRowTo():
        return 0

    b = min(1, r - distance)
    a = r
    while a >= b:
        c = map.mapRowToCol(a)
        if c != 0:
            return c + offset
        a -= 1

    b = max(map.getRowTo(), r + distance )
    a = r

    while a <= b:
        c = map.mapRowToCol(a)
        if c != 0: return c + offset
        a += 1

    return -1

def WriteGeneStructureCorrespondence( mali, identifiers, exons, param_master_pattern, gap_char = "-" , prefix = "" ):
    """split multiple alignment into clusters of orthologous transcripts.

    Orthologous transcripts are defined by similarity of gene structure to
    query sequences.

    Also: return matrix of gene structure compatibility

    0   : perfect compatibility (exact match)

    ratio of missed exon boundaries to total exon boundaries.
    
    100 : no compatibility
    """

    wmali = len(identifiers)
    lmali = len(mali[identifiers[0]])

    matrix_compatibility = numpy.zeros( (wmali, wmali) )

    if len(identifiers) == 0: return
    wmali = len(identifiers)
    lmali = len(mali[identifiers[0]])

    nok = 0
    nperfect   = 0

    ntotal_exons = 0
    nidentical_exons = 0
    nskipped_exons = 0

    ref_nok = 0
    ref_nperfect = 0

    ref_ntotal_exons = 0
    ref_nidentical_exons = 0
    ref_nskipped_exons = 0
    ref_ntotal = 0
    
    rx = re.compile( param_master_pattern )

    ## list of number of exons
    anexons = []

    ## exons in reference
    ref_nexons = 0
    x = 0
    for key1 in identifiers:

        seq = mali[key1]
        
        matches = []
        unassigned = []

        is_perfect = False

        anexons.append( len (exons[key1]) )
        if rx.search( key1 ):
            ref_nexons = len(exons[key1] )

        y = 0
        for key2 in identifiers:

            if key2 == key1: continue
            
            if param_loglevel >= 3:
                print "#############################################"
                print "# comparing %s to %s" % (key1, key2)
            
            mref = 0
            mcmp = 0
            
            seq_master = mali[key2]
            ref_exons = exons[key2]

            map_cmp2ref = MaliIO.getMapFromMali( seq, seq_master, gap_char )

            ## map exon boundaries to reference sequence
            cmp_exons = []

            if param_loglevel >= 5:
                print alignlib.writeAlignataTable(map_cmp2ref)
            
            for e in exons[key1]:
                ne = e.GetCopy()
                ne.mPeptideFrom = MyMap( map_cmp2ref, e.mPeptideFrom + 1, 3, -1 )
                ne.mPeptideTo   = MyMap( map_cmp2ref, e.mPeptideTo, 3, 0 )
                cmp_exons.append(ne)

            ## massage boundaries for terminal exons:
            if cmp_exons[0].mPeptideFrom <= 0: cmp_exons[0].mPeptideFrom = ref_exons[0].mPeptideFrom
            if cmp_exons[-1].mPeptideTo  <= 0: cmp_exons[-1].mPeptideTo  = ref_exons[-1].mPeptideTo
                
            if param_loglevel >= 4:
                for e in exons[key1]:
                    print "# exon", str(e)
                    
            if param_loglevel >= 3:                    
                for e in cmp_exons:
                    print "# exon", str(e)
                for e in ref_exons:
                    print "# exon", str(e)

            ## do exon comparison
            comparison = Exons.CompareGeneStructures( cmp_exons, 
                                                      ref_exons,
                                                      threshold_min_pide = 0,
                                                      threshold_slipping_exon_boundary = param_threshold_splipping_exon_boundary )


            if param_loglevel >= 3:
                print comparison.Pretty( prefix = "# EVAL: ")

            ## analyse results
            min_nexons = min(len(cmp_exons), len(ref_exons))                
            max_nexons = max(len(cmp_exons), len(ref_exons))

            similarity = (max_nexons - comparison.mNumIdenticalExons) * (abs( comparison.mNumDifferenceExons))
            
            is_perfect = False
            is_ok = False
            status = []

            # non-equivalent exon pairs
            ne = len(cmp_exons) - comparison.mNumIdenticalExons - comparison.mNumSkippedExons

            is_perfect = False
            is_ok = False
            if comparison.mNumIdenticalExons == 0:
                # F: complete and utter failure, no excuses
                status.append( "F" )
            else:
                if ne == 0:
                    # P: perfect conservation
                    status.append( "=" )
                    is_ok = True
                    is_perfect = True
                elif ne == min_nexons - comparison.mNumSkippedExons:
                    # D: completely different predictions
                    status.append( "D" )
                elif ne in (1,2):
                    # A: almost conserved
                    status.append( "A" )
                    is_ok = True
                elif ne > 2:
                    # M : mostly conserved (in case of long proteins that is good enough).
                    if (100 * comparison.mNumIdenticalExons) / max_nexons > param_evaluate_min_percent_exon_identity:
                        status.append( "M" )
                    else:
                    # S : spuriously conserved
                        status.append( "S" )                    
                else:
                    # U: unconserved
                    status.append( "U" )

            if len(cmp_exons) > len(ref_exons):
                status.append( ">" )
            elif len(ref_exons) < len(cmp_exons):
                status.append( "<" )
            else:
                status.append( "=" )

            if min_nexons == max_nexons and min_nexons == 1:
                status.append( "S" )
            elif min_nexons == 1 and max_nexons == 2:
                status.append( "s")
            elif min_nexons == 2 and max_nexons == 2:
                status.append( "D" )
            elif min_nexons == 2 and max_nexons > 2:
                status.append( "d" )
            elif min_nexons == max_nexons:
                status.append( "M" )
            elif min_nexons > 2 and max_nexons > 2:
                status.append( "m" )                
            else:
                status.append( "U")
                
            status = string.join( status, "")

            structure_compatibility = 100
            
            if is_ok:
                nok += 1
                structure_compatibility = 100 - 100 * (comparison.mNumIdenticalExons + comparison.mNumSkippedExons) / len(cmp_exons)
            if is_perfect:
                nperfect += 1
                structure_compatibility = 0
                
            if abs(comparison.mNumDifferenceExons) > param_max_exons_difference:
                compatibility_value = 100
            else:
                compatibility_value = structure_compatibility

            t = comparison.mNumRefBoundaries + comparison.mNumCmpBoundaries

            if t == 0:
                compatibility_value = 0
            else:
                compatibility_value = 100 * (comparison.mNumMissedRefBoundaries + comparison.mNumMissedCmpBoundaries) / t
                
            matrix_compatibility[x][y] = compatibility_value
            
            nidentical_exons += comparison.mNumIdenticalExons
            nskipped_exons   += comparison.mNumSkippedExons
            ntotal_exons     += len(cmp_exons)

            if param_loglevel >= 2:
                print "%s\tgenepair\t%s\t%s\t%s\t%i\t%i\t%i\t%s" % (prefix, key1, key2, status, compatibility_value,
                                                                    len(cmp_exons), len(ref_exons), str(comparison))

            ## comparison to reference: count separately:
            if rx.search( key2 ):
                ref_nidentical_exons += comparison.mNumIdenticalExons
                ref_nskipped_exons   += comparison.mNumSkippedExons
                ref_ntotal_exons     += len(cmp_exons)
                if is_ok: ref_nok += 1
                if is_perfect: ref_nperfect += 1
                ref_ntotal += 1

            y += 1
            
        x += 1
                
    ntotal = wmali * ( wmali - 1)

    print "%s\tallstructure\t%i\t%i\t%i\t%6.4f\t%6.4f\t%i\t%i\t%i\t%6.4f\t%6.4f" % (prefix,
                                                                                    ntotal, nperfect, nok,
                                                                                    float(nperfect) / ntotal, float(nok) / ntotal,
                                                                                    ntotal_exons, nidentical_exons, nskipped_exons, 
                                                                                    float(nidentical_exons) / ntotal_exons,
                                                                                    float(nidentical_exons + nskipped_exons) / ntotal_exons)

    if ref_ntotal > 0:
        if ref_ntotal_exons == 0:
            raise "no exons in reference : ref_ntotal_exons = 0, ref_ntotal = %i" % (ref_ntotal)
        
        print "%s\trefstructure\t%i\t%i\t%i\t%6.4f\t%6.4f\t%i\t%i\t%i\t%6.4f\t%6.4f" % (prefix,
                                                                                        ref_ntotal, ref_nperfect, ref_nok,
                                                                                        float(ref_nperfect) / ref_ntotal, float(ref_nok) / ref_ntotal,
                                                                                        ref_ntotal_exons, ref_nidentical_exons, ref_nskipped_exons, 
                                                                                        float(ref_nidentical_exons) / ref_ntotal_exons,
                                                                                        float(ref_nidentical_exons + ref_nskipped_exons) / ref_ntotal_exons)
        
    print "%s\tnexons\t%i\t%i\t" % (prefix,
                                    len(anexons), ref_nexons) +\
                                    string.join(map(lambda x: "%.2f" % x, (min(anexons),
                                                                           max(anexons),
                                                                           scipy.mean(anexons),
                                                                           scipy.median(anexons),
                                                                           scipy.std(anexons))), "\t")

    return matrix_compatibility

##---------------------------------------------------------------------------------------                                                         
def oldClusterMatrixClosestDistance( matrix_weights,
                                  identifiers,
                                  master_pattern,
                                  species_pattern,
                                  matrix_compatibility = None):
    """cluster identifiers by their coverage.

    """

    rxm = re.compile( master_pattern )
    rxs = re.compile( species_pattern )

    # positions of masters in identifier list
    master_ids = []
    # members of each clusters
    clusters= {}
    # set of species per cluster
    species = {}
    # unassigned identifiers (no masters)
    unassigned = {}

    # initialize
    for x in range(len(identifiers)):
        i = identifiers[x]
        g = rxs.search( i )
        if g:
            s = g.groups()[0]
            species[i] = {}
            
        if rxm.search( i ):
            master_ids.append( x )
            clusters[i] = []
        else:
            unassigned[i] = 1

    if len(master_ids) == 0:
        clusters = { None : identifiers }
        return clusters, []

    ## merge clusters
    skip = {}
    new = []
    for x in range(0, len(master_ids)-1):
        xm = master_ids[x]
        xi = identifiers[xm]
        if xi in clusters:
            new.append( xm )            
            for y in range(x+1, len(master_ids)):
                ym = master_ids[y]
                yi = identifiers[ym]
                if yi not in clusters: continue
                if matrix_compatibility[xm][ym] <= param_cluster_merging_max_distance:
                    if param_loglevel >= 3:
                        print "# merging %s with %s" % (xi, yi )
                    clusters[xi] += clusters[yi] + [yi]
                    del clusters[yi]

    master_ids = new

    # build edges
    edges = []
    best_matches = {}
    # if masters are identical, the first one gets all the
    # the members. 
    for x in range(len(identifiers)):
        
        i = identifiers[x]
        g = rxs.search( i )
        if g:
            s = g.groups()[0]
        else:
            s = None
            
        if i in clusters: continue
        
        for m in master_ids:
            
            c = matrix_weights[x][m]
            
            if matrix_compatibility:
                co = max(matrix_compatibility[x][m], matrix_compatibility[m][x])
            else:
                co = 0
                
            if c <= param_clustering_max_distance and co <= param_max_compatibility:
                edges.append( (c, x, m, s, co) )
                    

        # best: best entry per identifier with compatible gene structure
        best = None
        best_m = None

        if param_loglevel >= 5:
            print "# pair:", i, identifiers[best_m], matrix_weights[x][best_m]
            
        for m in master_ids:
            if param_loglevel >= 5:
                print "# pair:", i, identifiers[m], matrix_weights[x][m]

            if matrix_compatibility:
                co = max(matrix_compatibility[x][m], matrix_compatibility[m][x])
            else:
                co = 0

            if matrix_compatibility and co <= param_max_compatibility:
                if best == None or matrix_weights[x][m] < best:
                    best = matrix_weights[x][m]
                    best_m = m

        if best != None and best <= param_clustering_max_distance:
            best_matches[i] = (best_m, best)
            
    # go through edge list by increasing distance and add identifiers
    # to clusters where they don't have another member of the same species.
    edges.sort()
    for coverage, x, m, s, compatibility in edges:

        cm = identifiers[m]
        i = identifiers[x]

        if i not in unassigned: continue
        
        if s and s not in species[cm]:
            species[cm][s] = 1
            clusters[cm].append( i )
            del unassigned[i]
            if param_loglevel >= 4:
                print "# assigning %s to %s with weight %f and compatibility %f" % ( i, cm, coverage, compatibility )

    ## assign remaining to best cluster
    for i in unassigned.keys():
        if i in best_matches:
            bm, coverage = best_matches[i]
            cm = identifiers[bm]
            if param_loglevel >= 4:            
                print "# assigning %s to %s as best match with weight %f" % ( i, cm, coverage )            
            clusters[cm].append( i )
            del unassigned[i]
            
    
    # assign masters to clusters. If they have no members,
    # assign to best match
    reps = []
    mems = []
    for m in master_ids:
        i = identifiers[m]
        if len(clusters[i]) > 0:
            clusters[i].append( i )
            reps.append( m )
        else:
            mems.append( m )

    for x in mems:
        i = identifiers[x]
        best = matrix_weights[x][reps[0]]
        best_m = reps[0]
        for m in reps[1:]:
            if matrix_weights[x][m] < best:
                best = matrix_weights[x][m]
                best_m = m

        if best <= param_clustering_max_distance:
            if param_loglevel >= 4:
                print "# assigning %s to %s with weight %f" % ( i, identifiers[best_m], best)
            clusters[identifiers[best_m]].append(i)
            del clusters[i]
        else:
            unassigned[i] = 1
            clusters[i].append( i )

#     # if two masters are identical to another, merge
#     skip = {}
#     for x in range(0, len(reps)-1):
#         xi = reps[x]
#         for y in range( x+1, len(reps)):
#             yi = reps[y]
#             if y in skip: continue
#             if matrix_weights[xi][yi] <= param_cluster_merging_max_distance:
#                 if param_loglevel >= 3:
#                     print "# merging %s with %s" % (identifiers[xi], identifiers[yi] )
#                 clusters[identifiers[xi]] += clusters[identifiers[yi]]
#                 del clusters[identifiers[yi]]
#                 skip[y] = 1
            
    return clusters, unassigned.keys()

##---------------------------------------------------------------------------------------                                                         
def ClusterMatrixClosestDistance( matrix_weights,
                                  identifiers,
                                  master_pattern,
                                  species_pattern,
                                  matrix_compatibility = None):
    """cluster identifiers by their coverage.

    """

    rxm = re.compile( master_pattern )
    rxs = re.compile( species_pattern )

    # positions of masters in identifier list
    master_ids = []
    # members of each clusters
    clusters= {}
    # set of species per cluster
    species = {}
    # unassigned identifiers (no masters)
    unassigned = {}

    # initialize, each cluster contains one master
    for x in range(len(identifiers)):
        i = identifiers[x]
        g = rxs.search( i )
        if g:
            s = g.groups()[0]
            species[i] = {}
            
        if rxm.search( i ):
            master_ids.append( x )
            clusters[i] = [ (i, x, "", 0.0, True) ] 
        else:
            unassigned[i] = 1

    if len(master_ids) == 0:
        clusters = { None : identifiers }
        return clusters, []

    if param_loglevel >= 2:
        print "# CLUSTERING: at start: %i/%i clusters for %i unassigned" % (len(clusters), len(master_ids), len(unassigned))
        
    ## merge clusters if masters are similar
    skip = {}
    new = []
    for x in range(0, len(master_ids)):
        xm = master_ids[x]
        xi = identifiers[xm]
        if xi in clusters:
            new.append( xm )            
            for y in range(x+1, len(master_ids)):
                if y == x: continue
                ym = master_ids[y]
                yi = identifiers[ym]
                if yi not in clusters: continue
                if matrix_weights[xm][ym] <= param_cluster_merging_max_distance:
                    if param_loglevel >= 3:
                        print "# merging %s with %s" % (xi, yi )
                    clusters[xi] += clusters[yi]
                    del clusters[yi]

    if param_loglevel >= 2:
        print "# CLUSTERING: after merging: %i clusters for %i unassigned" % (len(clusters), len(unassigned))

    master_ids = new

    # assign matches to best cluster
    for x in range(len(identifiers)):
        
        id = identifiers[x]
        g = rxs.search( id )
        if g:
            species = g.groups()[0]
        else:
            species = None
            
        if id not in unassigned: continue
        
        # best: best entry per identifier with compatible gene structure
        best = None
        best_m = None

        for m in master_ids:
            if param_loglevel >= 5:
                print "# pair:", id, identifiers[m], matrix_weights[x][m]

            if best == None or matrix_weights[x][m] < best:
                best = matrix_weights[x][m]
                best_m = m
        
        if best != None and best <= param_clustering_max_distance:
            if param_loglevel >= 3:
                print "# assigning %s to %s with weight %f" % (id, identifiers[best_m], best)
            clusters[identifiers[best_m]].append( (id, x, species, best, False) )
            del unassigned[id]

    if param_loglevel >= 2:
        print "# CLUSTERING: after assignment: %i clusters for %i unassigned" % (len(clusters), len(unassigned))

    # for each cluster sort according to compatibility and keep best
    for m in master_ids:

        new = []
        to_sort = []
        master_id = identifiers[m]
        for id, index, species, weight, is_master in clusters[master_id]:
            if is_master:
                new.append( id )
            else:
                to_sort.append( (species, matrix_compatibility[index][m], weight, id) )

        # this sorts first by species and then by compatibility in increasing order
        to_sort.sort()
            
        last_species = None
        
        for species, compatibility, weight, id in to_sort:

            if last_species == species:
                if param_loglevel >= 3:
                    print "# cluster: %s: removing %s at weight %f and compatibility %f" % (master_id, id, weight, compatibility)
                unassigned[id] = 1                
                continue
            else:
                if param_loglevel >= 3:
                    print "# cluster: %s: keeping %s at weight %f and compatibility %f" % (master_id, id, weight, compatibility)
                
            last_species = species
            new.append( id )
            
        clusters[master_id] = new

    if param_loglevel >= 2:
        print "# CLUSTERING: after compatiblitity: %i clusters for %i unassigned" % (len(clusters), len(unassigned))
        
    return clusters, unassigned.keys()


##---------------------------------------------------------------------------------------
def WriteSpeciesCoverage( identifiers, species_pattern, prefix = ""):
    """write number of species present in identifiers."""

    species = {}
    nunknown = 0
    rx = re.compile( species_pattern )
    
    for i in identifiers:
        x = rx.search(i)
        if x:
            s = x.groups()[0]
        else:
            nunknown += 1
            s = "unknown"
            
        if s not in species: species[s] = 0
        species[s] += 1

    if len(species) == 0:
        max_species = 0
    else:
        max_species = max(species.values())
    print "%s\tspecies\t%i\t%i\t%i\t%i" % (prefix,
                                           len(species),
                                           max_species,
                                           len(filter( lambda x: x >= 2, species.values())),
                                           nunknown )


##---------------------------------------------------------------------------------------
def WriteCodonSummary( mali, identifiers, frame_columns, prefix = "", gap_char = "-"):
    """write codon summary."""
    
    new_mali = {}
    aligned = []
    codons = []
    stops = []
    nclean = 0
    total_no_stops = 0
    for key, seq in core_mali.items():
        new_mali[key], naligned, ncodons, nstops = MaliIO.getCodonSequence( seq, frame_columns, param_gap_char, remove_stops = True)
        aligned.append( naligned )
        codons.append( ncodons )
        stops.append( nstops )
        if nstops == 0:
            total_no_stops += 1
        if naligned == ncodons and nstops == 0:
            nclean += 1
            
    print "%s\tcodons\t%i\t%i\t" % (prefix, nclean, total_no_stops) +\
          string.join(map(lambda x: "%.2f" % x, (min(aligned),
                                                 max(aligned),
                                                 scipy.mean(aligned),
                                                 scipy.median(aligned),
                                                 scipy.std(aligned))),"\t") + "\t" +\
          string.join(map(lambda x: "%.2f" % x, (min(codons),
                                                 max(codons),
                                                 scipy.mean(codons),
                                                 scipy.median(codons),
                                                 scipy.std(codons))),"\t") + "\t" +\
          string.join(map(lambda x: "%.2f" % x, (min(stops),
                                                 max(stops),
                                                 scipy.mean(stops),
                                                 scipy.median(stops),
                                                 scipy.std(stops))),"\t")
                                                 
        
    return new_mali

##---------------------------------------------------------------------------------------
def WriteRadius( mali, identifiers, prefix = "", gap_char = "-"):
    """write percent identities in pairwise comparisons both for nucleotide acids and amino acids."""

    pides_na = []
    seq_aa = []
    for x in range(0,len(identifiers)):

        seq_aa.append( Genomics.TranslateDNA2Protein(mali[identifiers[x]]))
                       
        for y in range(x+1, len(identifiers)):
            if x == y: continue
            pides_na.append( MaliIO.getPercentIdentity( mali[identifiers[x]], mali[identifiers[y]], gap_char) )

    pides_aa = []
    for x in range(0,len(identifiers)-1):
        for y in range(x+1, len(identifiers)):        
            pides_aa.append( MaliIO.getPercentIdentity( seq_aa[x], seq_aa[y], gap_char) )            
            
    print "%s\tpide\t%i\t" % (prefix, len(pides_na)) +\
          string.join(map(lambda x: "%.2f" % x, (min(pides_na),
                                                 max(pides_na),
                                                 scipy.mean(pides_na),
                                                 scipy.median(pides_na),
                                                 scipy.std(pides_na))),"\t") + "\t" +\
          string.join(map(lambda x: "%.2f" % x, (min(pides_aa),
                                                 max(pides_aa),
                                                 scipy.mean(pides_aa),
                                                 scipy.median(pides_aa),
                                                 scipy.std(pides_aa))),"\t") 
    
            
##------------------------------------------------------------
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)
        elif o in ( "--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-o", "--file-output"):
            param_filename_output = a
        elif o == "--min-overlap-residues":
            param_min_overlap_residues = int(a)
        elif o == "--min-overlap-percent":
            param_min_overlap_percent = int(a)
        elif o == "--min-coverage-residues":
            param_min_coverage_residues = int(a)
        elif o == "--min-coverage-percent":
            param_min_coverage_percent = int(a)
        elif o in ("-m", "--master-pattern"):
            param_master_pattern = a
        elif o in ("-s", "--species-pattern"):
            param_species_pattern = a
        elif o in ("-e", "--exons"):
            param_filename_exons = a
        elif o in ("-c", "--cluster"):
            param_do_cluster = True
        elif o in ("-f", "--remove-fragments"):
            param_remove_fragments = True
        elif o in ("-p", "--prefix"):
            param_prefix = a
        elif o == "--components":
            param_filename_components = a

    print E.GetHeader()
    print E.GetParams()
            
    ## 1. read multiple alignment in fasta format
    all_mali, all_identifiers = MaliIO.readFasta( sys.stdin )

    if len(all_identifiers) == 0:
        raise "alignment is empty."

    if param_loglevel >= 1:
        print "# read mali with %i entries." % len(all_identifiers)

    if param_filename_components:

        infile = open(param_filename_components, "r")
        components = {}
        for line in infile:
            if line[0] == "#": continue
            if line[0] == ">": continue
            a, b = line[:-1].split("\t")[:2]
            if b not in components:
                components[b] = []
            components[b].append(a)
        
        if param_loglevel >= 1:
            print "# read %i components." % len(components)
    else:
        components = { 'all': all_identifiers }

    if param_filename_exons:
        exons = Exons.ReadExonBoundaries( open(param_filename_exons, "r"), filter=all_mali )
        if param_loglevel >= 2:
            print "# read %i exons." % len(exons)
    else:
        exons = {}

    print "# PREFIX\tsummary\tNSEQUENCES\tNASSIGNED\tNCLUSTERS\tNASSIGNED\tUNASSIGNED"
    print "# PREFIX\tcluster\tNMEMBERS\tMEMBERS"
    print "# PREFIX\tfragments\tNFRAGMENTS\tFRAGMENTS"

    print "# PREFIX\tpide\tNPAIRS\tNAMIN\tNAMAX\tNAMEAN\tNAMEDIAN\tNASTDDEV\tAAMIN\tAAMAX\tAAMEAN\tAAMEDIAN\tAASTDDEV"
    
    print string.join( ("# PREFIX", "codons", 
                        "NCLEAN","NNOSTOPS",
                        "ALIGNED_MIN","ALIGNED_MAX","ALIGNED_MEAN","ALIGNED_MEDIAN","ALIGNED_STDDEV",
                        "CODONS_MIN","CODONS_MAX","CODONS_MEAN","CODONS_MEDIAN","CODONS_STDDEV",
                        "STOPS_MIN","STOPS_MAX","STOPS_MEAN","STOPS_MEDIAN","STOPS_STDDEV"), "\t")

    print string.join( ("# PREFIX", "allstructure", 
                        "NTOTAL","NPERFECT", "NOK",
                        "PPERFECT", "POK",
                        "NTOTAL_EXONS", "NIDENTICAL_EXONS", "NSKIPPED_EXONS",
                        "PIDENTICAL_EXONS", "PIDENTICAL_EXONS"), "\t")

    print string.join( ("# PREFIX", "refstructure", 
                        "NTOTAL","NPERFECT", "NOK",
                        "PPERFECT", "POK",
                        "NTOTAL_EXONS", "NIDENTICAL_EXONS", "NSKIPPED_EXONS",
                        "PIDENTICAL_EXONS", "PIDENTICAL_EXONS"), "\t")

    print string.join( ("# PREFIX", "exons", 
                        "NDATA","NREF_EXONS",
                        "NEXONS_MIN","NEXONS_MAX","NEXONS_MEAN","NEXONS_MEDIAN","NEXONS_STDDEV"), "\t")

    print string.join( ("# PREFIX", "species", 
                        "NSPECIES","SPECIES_MAX",
                        "MAX_PER_SPECIES","UNKNOWN"), "\t")

    print string.join( ("# PREFIX", "failed",
                        "NFAILED_SEQS","NTOTAL_SEQS", "PFAILED_SEQS",
                        "NFAILED_PAIRS","NTOTAL_PAIRS", "PFAILED_PAIRS"), "\t")

    print string.join( ("# PREFIX", "npairs",
                        "NPAIRS_MIN","NPAIRS_MAX","NPAIRS_MEAN","NPAIRS_MEDIAN","NPAIRS_STDDEV"), "\t")
                        
    print string.join( ("# PREFIX", "ppairs",
                        "PPAIRS_MIN","PPAIRS_MAX","PPAIRS_MEAN","PPAIRS_MEDIAN","PPAIRS_STDDEV"), "\t")

    print string.join( ("# PREFIX", "cov",
                        "COV_MIN","COV_MAX","COV_MEAN","COV_MEDIAN","COV_STDDEV"), "\t")

    print string.join( ("# PREFIX", "pcov",
                        "PCOV_MIN","PCOV_MAX","PCOV_MEAN","PCOV_MEDIAN","PCOV_STDDEV"), "\t")

    print string.join( ("# PREFIX", "genepair",
                        "STATUS", "COMPATIBILITY", "CMP_NEXONS", "REF_NEXONS", Exons.ComparisonResult().GetHeader()), "\t")
                        
    for key, identifiers in components.items():

        prefix = "%s%s_all" % (param_prefix, key)
        
        ## 1. remove gaps in multiple alignment
        mali = MaliIO.removeGappedColumns( MaliIO.getSubset( all_mali, identifiers), param_gap_char )

        ## 1. remove gaps in multiple alignment
        if param_do_cluster:
            matrix_overlap, matrix_total, coverages = GetOverlapMatrix( mali, identifiers)

            matrix_coverage = 100 - ((100 * matrix_overlap) / matrix_total)

            if param_loglevel >= 3:
                print "# overlap\n", matrix_overlap
                print "# total\n", matrix_total
                print "# coverage\n", matrix_coverage
                print "# coverages\n", coverages

            matrix_compatibility = WriteGeneStructureCorrespondence( mali, identifiers, exons, param_master_pattern, prefix = prefix )

            if param_loglevel >= 3:
                print "# compatibility\n", matrix_compatibility

            clusters, unassigned = ClusterMatrixClosestDistance( matrix_coverage, identifiers,
                                                                 param_master_pattern, param_species_pattern,
                                                                 matrix_compatibility)
        else:
            clusters = { 'all': identifiers }

        cluster_id = 0

        print "%s%s\tsummary\t%i\t%i\t%i\t%i\t%s" % (param_prefix, key,
                                                     len(all_identifiers),
                                                     len(all_identifiers) - len(unassigned),
                                                     len(clusters),
                                                     len(unassigned),                                                
                                                     string.join(unassigned, "\t"))

        for k, v in clusters.items():

            cluster_id += 1
            prefix = "%s%s_%i" % (param_prefix, key, cluster_id)
            
            if param_loglevel >= 1:
                print "###########################################################"

            if len(v) == 0:
                print "# cluster %s contains no sequence" % (prefix)
                continue
                
            if len(v) == 1:
                print "# cluster %s contains only one sequence %s" % (prefix, v[0])
                continue
                
            # 1. remove gaps in multiple alignment
            ungapped_mali = MaliIO.removeGappedColumns( MaliIO.getSubset( mali, v ), param_gap_char )

            # write overlap parameters
            fragments = WriteOverlap( ungapped_mali, v, prefix = prefix )

            not_fragments = filter( lambda x: x not in fragments, v )

            print "%s\tcluster\t%s\t%i\t%s" % (prefix, k, len(not_fragments), string.join(not_fragments, "\t"))

            if len(fragments) > 0:
                print "%s\tfragments\t%i\t%s" % (prefix, len(fragments), string.join(fragments, "\t"))

            if len(not_fragments) == 1:
                continue

            if param_remove_fragments:
                core_mali = MaliIO.removeGappedColumns( MaliIO.getSubset( ungapped_mali, not_fragments ), param_gap_char )
            else:
                core_mali = ungapped_mali

            # check how many exons are recovered.
            WriteGeneStructureCorrespondence( core_mali, not_fragments, exons, param_master_pattern, prefix = prefix )

            # write coverage of genomes
            WriteSpeciesCoverage( not_fragments, species_pattern = param_species_pattern, prefix = prefix )

            # get frame columns
            frame_columns = MaliIO.getFrameColumnsForMasterPattern( core_mali, not_fragments, param_master_pattern, param_gap_char )

            # translate multiple alignment
            codon_mali = WriteCodonSummary( core_mali, not_fragments, frame_columns, prefix, param_gap_char )
            
            # write coverage of genomes
            WriteRadius( codon_mali, not_fragments, prefix = prefix )

            if param_output_pattern:
                filename = param_output_pattern % cluster_id
                outfile = open( filename, "w")
                for id in core_mali.keys():
                    outfile.write( ">%s\n%s\n" % (id, core_mali[id]))
                outfile.close()

    print E.GetFooter()
    
    

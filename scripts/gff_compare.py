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
gff_compare.py - compare two gene sets
======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script benchmarks two gene predictions. Predictions are given as gff files.

This script is transcript aware but not gene aware, i.e. if there are several
transcripts for a gene, overcounting happens.

Usage
-----

Example::

   python gff_compare.py --help

Type::

   python gff_compare.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import re
import optparse

USAGE="""python %s [OPTIONS] target reference

""" % sys.argv[0]

import CGAT.Experiment as E
import CGAT.PredictionParser as PredictionParser
import numpy

##------------------------------------------------------------------------
def Old():
    
    ref_nucleotides    = numpy.zeros( max_end - min_start, numpy.int0 )
    target_nucleotides = numpy.zeros( max_end - min_start, numpy.int0 )

    # set range overlapping by reference
    for x in targets:
        target_nucleotides[x.start - min_start:x.end - min_start] = 1

    for x in references:
        ref_nucleotides[x.start - min_start:x.end - min_start] = 1

    not_ref = numpy.logical_not(ref_nucleotides)
    not_target = numpy.logical_not(target_nucleotides)
    
    v = numpy.logical_and( ref_nucleotides, target_nucleotides) 
    tp = numpy.dot( v, v )

    v = numpy.logical_and( ref_nucleotides, not_target )
    fn = numpy.dot( v, v )    
    
    v = numpy.logical_and( not_ref, target_nucleotides) 
    fp = numpy.dot( v, v )

    v = numpy.logical_and( not_ref, not_target)
    tn = numpy.dot( v, v )

##------------------------------------------------------------------------
def GetCounts( nucleotides, val ):
    """count number of occurances of val in array.
    
    can't find a way without running into memory troubles.
    """
    t = 0
    for x in nucleotides:
        if x == val: t += 1
    return t

##------------------------------------------------------------------------
def AnalyseOverlaps( references, targets ):
    """calculate nucleotide overlap.

    All references and targets should be on the same chromsome.
    Name and strand are not checked.
    """

    if len(references) == 0 or len(targets) == 0:
        return( (0,0,0,0) )
    
    targets_start = min( map(lambda x: x.start, targets) )
    targets_end   = max( map(lambda x: x.end, targets) )

    refs_start = min( map(lambda x: x.start, references) )
    refs_end   = max( map(lambda x: x.end, references) )
    
    min_start = min( targets_start, refs_start )
    max_end   = max( targets_end, refs_end)

    nucleotides    = numpy.zeros( max_end - min_start, numpy.int0 )
    
    # set ranges: bit 0: target, bit 1: reference
    for x in targets[:100]:
        nucleotides[x.start - min_start:x.end - min_start] = 1

    for x in references[:100]:
        nucleotides[x.start - min_start:x.end - min_start] = \
                             numpy.array( [2] * (x.end - x.start), numpy.int0 )

    tp = GetCounts( nucleotides, 3 )
    fp = GetCounts( nucleotides, 1 )
    tn = GetCounts( nucleotides, 0 )
    fn = GetCounts( nucleotides, 2 )
    
    return (tp, fp, tn, fn )

##------------------------------------------------------------------------
def CalculateSpecificitySensitivity( tp, fp, tn, fn ):

    if fp + tp == 0:
        spec = 1.0
    else:
        spec = float(tp) / (fp + tp)

    if tp + fn == 0:
        sens = 0.0
    else:
        sens = float(tp) / (tp + fn)
        
    return spec, sens

##------------------------------------------------------------------------
def CalculateCorrelationCoefficient( tp, fp, tn, fn ):

    x = ( ( tp + fn ) * ( tn + fp ) * ( tp + fp ) * ( tn + fn ) )
    if x == 0:
        return 0.0
    else:
        return (tp * tn - fn * fp ) / x
               

##------------------------------------------------------------------------
def GetFirstOverlaps( gffs, index ):
    """get next overlapping gffs from array starting at index.
    """
    
    this_name = gffs[index].mName
    this_strand = gffs[index].strand
    l = len(gffs)
    min_start = gffs[index].start
    max_end = gffs[index].end

    overlaps = [gffs[index]]
    index += 1
    
    while index < l:
        if this_name != gffs[index].mName or \
           this_strand != gffs[index].strand or\
           max_end < gffs[index].start:
            return overlaps, index, min_start, max_end

        max_end = max( max_end, gffs[index].end )
        overlaps.append( gffs[index] )
        index += 1
        
    return overlaps, index, min_start, max_end        

##------------------------------------------------------------------------
def CountMatchesPerGene( gffs,
                         rx_gene,
                         rx_other,
                         write = (),
                         outfile = sys.stdout ):
    """use status information in gffs to check for completely matches/extra genes.
    """

    # map of gene to matching exons
    # match contains tuple of:
    # nexons, nextra, complete matches to other_gene, partial matches to other gene
    map_gene2matches = {}
    for x in gffs:

        try:
            gene = rx_gene.search( x.mInfo ).groups()[0]
        except AttributeError:
            print "# ERROR: could not find gene identifier in %s" % x.mInfo
            sys.stdout.flush()
            sys.exit(1)
         
        # total, extra, matches(per other gene), partial(per other gene)
        if gene not in map_gene2matches:
            map_gene2matches[gene] = [0, 0, {}, {}] 

        info = map_gene2matches[gene]
        info[0] += 1

        try:
            if x.mStatus == "extra":
                info[1] += 1
            elif x.mStatus in ("match", "partial"):

                for y in x.mMatches:
                    other_gene = rx_other.search( y.mInfo ).groups()[0] 
                    if other_gene not in info[2]: info[2][other_gene] = 0
                    if other_gene not in info[3]: info[3][other_gene] = 0

                    if GTF.Identity( x, y, max_slippage=options.max_exon_slippage ):
                        info[2][other_gene] += 1
                        info[3][other_gene] += 1                    
                    elif GTF.HalfIdentity( x, y, max_slippage=options.max_exon_slippage ):
                        info[3][other_gene] += 1
        except AttributeError:
            print "Programming ERROR: no mStatus for ", str(x)
            
    total_extra = 0
    total_match = 0
    total_partial_match = 0
    total = len(map_gene2matches)
    for gene, info in map_gene2matches.items():

        l, nextra, matches, partials = info
        
        if nextra == l:
            total_extra += 1
            if "missed" in write:
                outfile.write( "missed\t%s\n" % (gene,))                                        
        else:
            for x, y in matches.items():
                if y == l:
                    total_match += 1
                    if "match" in write:
                        outfile.write( "match\t%s\t%s\n" % (gene, x))                        
                    break
                        
            else:
                for x,y in partials.items():            
                    if y == l:
                        total_partial_match += 1
                        if "partial" in write:
                            outfile.write( "partial\t%s\t%s\n" % (gene, x))
                        break

    return total, total_match, total_partial_match, total_extra

##------------------------------------------------------------------------
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: gff_compare.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"] )

    parser.add_option("-f", "--write-full", dest="write_full",
                      help="write full gff entries.", action="store_true"  )
    parser.add_option("-e", "--write-matched-exons", dest="write_matched_exons", 
                      help="write matched exons.", action="store_true"  )
    parser.add_option("-o", "--write-missed-exons", dest="write_missed_exons", action="store_true",
                      help="write missed exons."  )
    parser.add_option("-g", "--write-missed-genes", dest="write_missed_genes", action="store_true",
                      help="write missed genes."  )
    parser.add_option("-r", "--regex-reference", dest="regex_reference", type="string",
                      help="regular expression mapping exon to transcript in reference."  )
    parser.add_option("-t", "--regex-target", dest="regex_target", type="string",
                      help="regular expression mapping exon to transcript in target."  )
    parser.add_option("--no-nucleotides", dest="do_nucleotides", action="store_false",
                      help="skip nucleotide benchmark."  )
    parser.add_option("--no-exons", dest="do_exons", action="store_false",
                      help="skip exon benchmark."  )
    parser.add_option("--no-genes", dest="do_genes", action="store_false",
                      help="skip gene benchmark."  )
    parser.add_option("--out-filename-pattern", dest="outfile_pattern", type="string",
                      help="output filename pattern for extra info (%s will be substituted with reference,target)."  )
    
    parser.set_defaults(
        remove_redundancy = False,
        max_exon_slippage = 9,
        write_missed_exons = False,
        write_matched_exons = False,
        write_missed_genes = False,
        write_wrong_exons = False,
        write_wrong_genes = False,
        do_nucleotides = True,
        do_exons = True,
        do_genes = True,
        regex_reference = None,
        regex_target = None,
        outfile_pattern = "%s.info",
        )

    (options, args) = E.Start( parser )

    if len(args) != 2:
        print USAGE
        print "two arguments required"
        sys.exit(1)

    input_filename_target, input_filename_reference = args

    if options.loglevel >= 1:
        print "# target entries from %s" % input_filename_target
        print "# reading target entries ...",
        sys.stdout.flush()
        
    gff_targets = GTF.readFromFile( open( input_filename_target, "r" ) )
    
    if options.loglevel >= 1:
        print "finished: %i" % (len(gff_targets))
        sys.stdout.flush()

    if options.loglevel >= 1:
        print "# reference entries from %s" % input_filename_reference
        print "# reading reference entries ...",
        sys.stdout.flush()    

    gff_references = GTF.readFromFile( open( input_filename_reference, "r" ) )

    if options.loglevel >= 1:
        print "finished: %i" % (len(gff_references))
        sys.stdout.flush()

    if options.remove_redundancy:
        gff_targets = GTF.CombineOverlaps( gff_targets )
        gff_references = GTF.CombineOverlaps( gff_references )    

        if options.loglevel >= 1:    
            print "# after filtering: targets=%i, references=%i" % (len(gff_targets), len(gff_references))

    #################################################################################
    ## sort exons
    if options.loglevel >= 1:
        print "# sorting exons ...",
        sys.stdout.flush()
        
    gff_targets.sort( lambda x,y: cmp( (x.mName, x.strand, x.start, x.end),
                                       (y.mName, y.strand, y.start, y.end) ) )
    

    gff_references.sort( lambda x,y: cmp( (x.mName, x.strand, x.start, x.end),
                                          (y.mName, y.strand, y.start, y.end) ) )
    

    ntargets = len(gff_targets)
    nreferences = len(gff_references)

    if options.loglevel >= 1:
        print "finished"
        sys.stdout.flush()
        
    #################################################################################
    # get nucleotide level accuracy
    # process each fragment separately
    if options.do_nucleotides:
        print """##############################################################################
# nucleotide accuracy
#
# TP: true postives  : nucleotides to be predicted coding in both target and reference.
# FN: true negatives : nucleotides to be predicted non-coding in both target and reference.
# FP: false positives: nucleotides to be predicted non-coding in target but not in reference.
# FN: false negatives: nucleotides to be predicted non-coding in reference but not in target.
#
# Sensitivity: TP / TP + FN = correct exons / actual exons
# Specificity: TP / TP + FP
# CC: correlation coefficient:  (tp*tn-fn*fp)/((tp+fn)*(tn+fp)*(tp+fp)*(tn+fn))"""
    
        headers = ("contig", "strand", "tp", "fp", "tn", "fn", "sp", "sn", "cc" )
    
        print "\t".join(headers)
        
        first_r, first_t = 0, 0
        r, t = 0, 0

        ttp, tfp, ttn, tfn = 0, 0, 0, 0

        ## this only works, if all contigs in reference are present in target.
        while r < nreferences and t < ntargets:

            this_name = gff_references[r].mName
            this_strand = gff_references[r].strand

            ## get all in references
            while r < nreferences and \
                      gff_references[r].mName == this_name and \
                      gff_references[r].strand == this_strand:
                r+= 1

            ## skip over extra contigs in target
            while t < ntargets and \
                      (gff_targets[t].mName != this_name or \
                       gff_targets[t].strand != this_strand):
                t+= 1
            first_t = t

            ## get all in targets
            while t < ntargets and \
                      gff_targets[t].mName == this_name and \
                      gff_targets[t].strand == this_strand:                  
                t+= 1
                
            tp, fp, tn, fn =  AnalyseOverlaps( gff_references[first_r:r],
                                               gff_targets[first_t:t] )

            spec, sens = CalculateSpecificitySensitivity( tp, fp, tn, fn )
            cc = CalculateCorrelationCoefficient( tp, fp, tn, fn )
            print "%s\t%s\t%i\t%i\t%i\t%i\t%5.2f\t%5.2f\t%5.2f" % (this_name, this_strand, tp, fp, tn, fn, spec, sens, cc )

            ttp += tp
            tfp += fp
            ttn += tn
            tfn += fn
            first_r, first_t = r, t

        spec, sens = CalculateSpecificitySensitivity( ttp, tfp, ttn, tfn )
        cc = CalculateCorrelationCoefficient( ttp, tfp, ttn, tfn )
        print "%s\t%s\t%i\t%i\t%i\t%i\t%5.2f\t%5.2f\t%5.2f" % ("all", "all", ttp, tfp, ttn, tfn, spec, sens, cc )

        sys.stdout.flush()
        
    #################################################################################
    if options.do_exons or options.do_genes:

        print """##############################################################################
# exon accuracy (see Burset & Guigo (1996))
#
# TP: true postives  : exons matching both exon boundaries exactly (a deviation of
#                         maximal %i nucleotides is allowed).
# FN: true negatives : not counted.
# FP: false positives: exons in target not overlapping with in reference.
# FN: false negatives: exons overlapping, but with deviant boundaries and exons
#                         in reference not overlapping with any target exon.
#
# ME: missed exons: exons in reference not overlapping with exons in target
#                 (add to false negative count).
# WE: wrong exons: exons in target not overlapping with exons in reference
#                 (add to false positive count).
#
# Sensitivity: TP / (TP + FN) = correct exons / actual exons
# Specificity: TP / (TP + FP) = correct exons / predicted exons
# CC: (sensitivity + specificity) / 2
#
# Values are calculated where both exon boundaries of a target exon have to be
# correct (category "full") and where only one exon boundary has to be correct (catgory "half").""" %\
        (options.max_exon_slippage)

        headers = ("category", "contig", "strand", "tp", "fp", "tn", "fn", "sp", "sn", "cc", "me", "we", "me", "we" )

        print "\t".join(headers)

        r, t = 0, 0
        next_r, next_t = r, t

        # strict false positves/negatives
        tp, fp, tn, fn = 0, 0, 0, 0
        ttp, tfp, ttn, tfn = 0, 0, 0, 0
        # partial false positives/negatives
        ptp, pfp, ptn, pfn = 0, 0, 0, 0
        tptp, tpfp, tptn, tpfn = 0, 0, 0, 0

        # missed and wrong exons
        missed_exons, wrong_exons = 0, 0
        tmissed_exons, twrong_exons = 0, 0    

        ## Flag set, if partial overlap in previous pair
        last_partial_overlap = False
        ## Flag set, if partial overlap and reference was last increased
        last_increased_ref = False

        while r < nreferences and t < ntargets:

            this_name = gff_references[r].mName
            this_strand = gff_references[r].strand

            ## get overlap segments
            if next_r == r:
                ref_overlaps, next_r, ref_start, ref_end = GetFirstOverlaps( gff_references, r )
            if next_t == t:
                target_overlaps, next_t, target_start, target_end = GetFirstOverlaps( gff_targets, t )        

            if options.loglevel >= 3:
                print "########################################################"
                for x in ref_overlaps:
                    print "#", str(x)
                for x in target_overlaps:            
                    print "#", str(x)

            do_summary = False
            ## check strand switch in reference
            if next_r < nreferences and \
                   (this_name != gff_references[next_r].mName or \
                    this_strand != gff_references[next_r].strand):
                if options.loglevel >= 3:
                    print "# target advance"
                do_summary = True

                last_increased_ref = False
                last_partial_overlap = False
                    
                ## advance in target until next name is found
                next_name = gff_references[next_r].mName
                next_strand = gff_references[next_r].strand            
                while next_t < ntargets and \
                          next_name != gff_targets[next_t].mName or \
                          next_strand != gff_targets[next_t].strand:
                    fp += 1
                    pfp += 1
                    target_overlaps, next_t, target_start, target_end = GetFirstOverlaps( gff_targets, next_t )

                for x in gff_targets[t:next_t]: x.mStatus = "extra"
                for x in gff_references[r:next_r]: x.mStatus = "extra"                
                
                r, t = next_r, next_t
            ## check strand switch in target
            elif next_t < ntargets and \
                     (this_name != gff_targets[next_t].mName or \
                      this_strand != gff_targets[next_t].strand):
                    ## advance in reference until next name is found
                if options.loglevel >= 3:
                    print "# reference advance"
                do_summary = True

                last_increased_ref = False
                last_partial_overlap = False
                
                next_name = gff_targets[next_t].mName
                next_strand = gff_targets[next_t].strand            
                while next_r < nreferences and \
                          next_name != gff_references[next_r].mName or \
                          next_strand != gff_references[next_r].strand:
                    fn += 1
                    pfn += 1
                    reference_overlaps, next_r, references_start, references_end = GetFirstOverlaps( gff_references, next_r )

                for x in gff_targets[t:next_t]: x.mStatus = "extra"                
                for x in gff_references[r:next_r]: x.mStatus = "extra"

                r, t = next_r, next_t
            ## otherwise                
            else:

                ref_status, target_status = None, None
                
                if options.loglevel >= 3:
                    print "# same chromosome"
                        
                ## overlap between segments            
                if min(ref_end, target_end) - max(ref_start, target_start) > 0:

                    # clear flags
                    last_increased_ref = False
                    last_partial_overlap = False
                    found = False
                    
                    for rr in ref_overlaps:
                        xfound = False
                        for tt in target_overlaps:
                            if GTF.Identity( rr, tt, max_slippage=options.max_exon_slippage ):
                                xfound = True
                                break
                        if xfound:
                            found = True
                            break
                        
                    if found:
                        ref_status = "match"
                        target_status = "match"
                        tp += 1
                        ptp += 1
                        if options.write_matched_exons:
                            print "############# matching exons ###########################" 
                            for x in ref_overlaps:
                                print "#", str(x)
                            for x in target_overlaps:            
                                print "#", str(x)
                    else:
                        fn += 1

                        ## check for one-sided matches
                        for rr in ref_overlaps:
                            xfound = False
                            for tt in target_overlaps:
                                if GTF.HalfIdentity( rr, tt, max_slippage=options.max_exon_slippage ):
                                    xfound = True
                                    break
                            if xfound:
                                found = True
                                break

                        if found:
                            ptp += 1
                            code = "partial"
                            ref_status = "partial"
                            target_status = "partial"
                        else:
                            pfn += 1
                            code = "complete"
                            ref_status = "mismatch"
                            target_status = "mismatch"
                            
                        if options.write_missed_exons:
                            print "############# %s non-overlapping exons ###########################" % code
                            for x in ref_overlaps:
                                print "#", str(x)
                            for x in target_overlaps:            
                                print "#", str(x)
                                
                    #################################################################
                    ## r, t = next_r, next_t                                                        
                    if ref_end == target_end:
                        r, t = next_r, next_t
                    elif ref_end < target_end:
                        r = next_r
                        last_increased_ref = True
                        last_partial_overlap = True
                    else:
                        t = next_t
                        last_increased_ref = False
                        last_partial_overlap = True
                        
                ## non-overlap between segments                
                else:

                    if ref_end < target_start:
                        
                        # for non-overlap, check whether there was partial overlap before
                        # and reference was not increased.
                        # if there was, just increment reference, but do not count.

                        if not (last_partial_overlap and not last_increased_ref):
                        
                            if options.write_missed_exons:
                                print "############# missed exon ###########################" 
                                for x in ref_overlaps:
                                    print "#", str(x)
                            missed_exons += 1
                            fn += 1
                            pfn += 1
                            ref_status = "extra"
                            
                        r = next_r
                        
                    else:
                        
                        # for non-overlap, check whether there was partial overlap before
                        # and target was not increased.
                        # if there was, just increment target, but do not count.

                        if not (last_partial_overlap and last_increased_ref):
                            if options.write_wrong_exons:
                                print "############# wrong exon ###########################" 
                                for x in target_overlaps:
                                    print "#", str(x)
                            
                            wrong_exons += 1
                            fp += 1
                            pfp += 1
                            target_status = "extra"
                            
                        t = next_t
                        
                    last_partial_overlap = False

                if options.loglevel >= 3:
                    print "# ref_status=%s, target_status=%s" % (ref_status, target_status)
                    
                if ref_status:
                    for rr in ref_overlaps:
                        rr.mStatus = ref_status
                        
                    if ref_status in ("match", "partial") and options.do_genes:
                        for rr in ref_overlaps:                            
                            rr.mMatches = target_overlaps

                if target_status:
                    for tt in target_overlaps:
                        tt.mStatus = target_status
                    
                    if target_status in ("match", "partial") and options.do_genes:
                        for tt in target_overlaps:
                            tt.mMatches = ref_overlaps

            if do_summary or r >= nreferences or t >= ntargets:
                ttp += tp
                tfp += fp
                ttn += tn
                tfn += fn

                tptp += ptp
                tpfp += pfp
                tptn += ptn
                tpfn += pfn

                tmissed_exons += missed_exons
                twrong_exons += wrong_exons

                if tp + fn != 0:
                    pmissed_exons = "%5.2f" % (float(missed_exons) / (tp + fn))
                else:
                    pmissed_exons = "0"

                if tp + fp != 0:
                    pwrong_exons = "%5.2f" % (float(wrong_exons) / (tp + fp))
                else:
                    pwrong_exons = "na"
                    
                spec, sens = CalculateSpecificitySensitivity( tp, fp, tn, fn )
                cc = (spec + sens) / 2.0
                print "full\t%s\t%s\t%i\t%i\t%i\t%i\t%5.2f\t%5.2f\t%5.2f\t%i\t%i\t%s\t%s" % \
                      (this_name, this_strand,
                       tp, fp, tn, fn,
                       spec, sens, cc, 
                       missed_exons, wrong_exons,
                       pmissed_exons, pwrong_exons )

                spec, sens = CalculateSpecificitySensitivity( ptp, pfp, ptn, pfn )
                cc = (spec + sens) / 2.0
                print "half\t%s\t%s\t%i\t%i\t%i\t%i\t%5.2f\t%5.2f\t%5.2f\t%i\t%i\t%s\t%s" %\
                      (this_name, this_strand,
                       ptp, pfp, ptn, pfn,
                       spec, sens, cc,
                       missed_exons, wrong_exons,
                       pmissed_exons, pwrong_exons )

                tp, fp, tn, fn = 0, 0, 0, 0
                ptp, pfp, ptn, pfn = 0, 0, 0, 0
                missed_exons, wrong_exons = 0, 0

        if t < ntargets:
            for x in gff_targets[t:ntargets]: x.mStatus = "extra"
        if r < nreferences:            
            for x in gff_references[r:nreferences]: x.mStatus = "extra"

        spec, sens = CalculateSpecificitySensitivity( ttp, tfp, ttn, tfn )
        cc = (spec + sens) / 2.0
        print "full\t%s\t%s\t%i\t%i\t%i\t%i\t%5.2f\t%5.2f\t%5.2f\t%i\t%i\t%5.2f\t%5.2f" % \
              ("all", "all", ttp, tfp, ttn, tfn,
               spec, sens, cc,
               tmissed_exons, twrong_exons,
               float(tmissed_exons) / (ttp + tfn),
               float(twrong_exons) / (ttp + tfp ) )

        spec, sens = CalculateSpecificitySensitivity( tptp, tpfp, tptn, tpfn )
        cc = (spec + sens) / 2.0
        print "half\t%s\t%s\t%i\t%i\t%i\t%i\t%5.2f\t%5.2f\t%5.2f\t%i\t%i\t%5.2f\t%5.2f" %\
              ("all", "all", tptp, tpfp, tptn, tpfn,
               spec, sens, cc,
               tmissed_exons, twrong_exons,
               float(tmissed_exons) / (ttp + tfn),
               float(twrong_exons) / (ttp + tfp ) )             

    if options.do_genes and \
           options.regex_reference and \
           options.regex_target:

        print """##############################################################################
# gene accuracy (see Burset & Guigo (1996))
#
# TP: true postives  : genes with exons matching both exon boundaries exactly (a deviation of
#                         maximal %i nucleotides is allowed).
# FN: true negatives : not counted.
# FP: false positives: genes in target not overlapping with in reference.
# FN: false negatives: genes overlapping, but with deviant exon boundaries and genes
#                         in reference not overlapping with any target exon.
#
# MG: missed genes: exons in reference not overlapping with exons in target
#                 (add to false negative count).
# WG: wrong genes: exons in target not overlapping with exons in reference
#                 (add to false positive count).
#
# Sensitivity: TP / (TP + FN) = correct genes / actual genes
# Specificity: TP / (TP + FP) = correct genes / predicted genes
# CC: (sensitivity + specificity) / 2
#
# Values are calculated where both exon boundaries of a target exon have to be
# correct (category "full") and where only one exon boundary has to be correct (catgory "half").""" %\
        (options.max_exon_slippage)

        out_options = []
        if options.write_missed_genes:
            out_options.append("missed")

        if options.loglevel >= 2:
            print "# counting matches for reference."
            sys.stdout.flush()
            
        (ref_total, ref_match, ref_partial, ref_extra) =\
                    CountMatchesPerGene( gff_references,
                                         re.compile(options.regex_reference),
                                         re.compile(options.regex_target),
                                         write = out_options,
                                         outfile = open(options.outfile_pattern % "reference", "w") )                                         

        if options.loglevel >= 2:
            print "# counting matches for target."
            sys.stdout.flush()
        
        (target_total, target_match, target_partial, target_extra) =\
                       CountMatchesPerGene( gff_targets,
                                            re.compile(options.regex_target),
                                            re.compile(options.regex_reference),
                                            write = out_options,
                                            outfile = open(options.outfile_pattern % "target", "w") )        

        if options.loglevel >= 1:
            print "# reference: genes=%6i, matches=%6i, partial=%6i, extra=%6i" % \
                  (ref_total, ref_match, ref_partial, ref_extra)
            print "# target   : genes=%6i, matches=%6i, partial=%6i, extra=%6i" % \
                  (target_total, target_match, target_partial, target_extra)            
        
        headers = ("category", "tp", "fp", "tn", "fn", "sp", "sn", "cc", "mg", "wg", "mg", "wg" )        
        print "\t".join(headers)

        tp = ref_match
        fp = target_extra
        tn = 0
        fn = ref_total - ref_match
        wrong_genes = target_extra
        missed_genes = ref_extra 

        spec, sens = CalculateSpecificitySensitivity( tp, fp, tn, fn )
        cc = (spec + sens) / 2.0

        if tp + fp == 0: fp = nreferences
        
        print "full\t%i\t%i\t%i\t%i\t%5.2f\t%5.2f\t%5.2f\t%i\t%i\t%5.2f\t%5.2f" % \
              (tp, fp, tn, fn,
               spec, sens, cc,
               missed_genes, wrong_genes,
               float(missed_genes) / (tp + fn),
               float(wrong_genes) / (tp + fp) )

        tp = ref_match + ref_partial
        fp = target_extra
        tn = 0
        fn = ref_total - ref_match - ref_partial
        wrong_genes = target_extra
        missed_genes = ref_extra 
        
        spec, sens = CalculateSpecificitySensitivity( tp, fp, tn, fn )
        cc = (spec + sens) / 2.0
        print "half\t%i\t%i\t%i\t%i\t%5.2f\t%5.2f\t%5.2f\t%i\t%i\t%5.2f\t%5.2f" % \
              (tp, fp, tn, fn,
               spec, sens, cc,
               missed_genes, wrong_genes,
               float(missed_genes) / (tp + fn),
               float(wrong_genes) / (tp + fp) )

    E.Stop()    



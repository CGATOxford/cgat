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
codemls2tsv.py - analyze results from site-specific codeml run
=======================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

The input is either:

   * filenames with a set of files of related codeml runs

summary-numbers:

jalview: annotation file for Jalview

Usage
-----

Example::

   python codemls2tsv.py --help

Type::

   python codemls2tsv.py --help

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
import optparse
import time

import CGAT.Genomics as Genomics
import CGAT.Experiment as E
import CGAT.WrapperCodeML as WrapperCodeML
import CGAT.Stats as Stats
import CGAT.Mali as Mali
import alignlib

def selectPositiveSites( results, selection_mode, options, mali = None ):
    """returns sites, which are consistently estimated to be positively selected.

    Depending on the option selection_mode, various sites are selected:

    'all': all positive sites are returned
    'consistent': only positive sites that are positive in all models and runs
    'emes': only sites that are > 0.9 in one model and at least > 0.5 in all other models

    If mali is given, positions that are not fully aligned are removed.

    """

    ## filter and extract functions
    if selection_mode == "emes":
        filter_f = lambda x: x.mProbability >= 0.5 and x.mOmega >= options.filter_omega
    else:
        filter_f = lambda x: x.mProbability >= options.filter_probability and x.mOmega >= options.filter_omega
        
    extract_f = lambda x: x.mResidue

    ## maximum significance per site (for emes)
    max_per_site = {}

    total_sites = set()

    first = True
    
    for result in results:

        for model in options.models:
        
            sites = result.mSites[model]
            
            s1, s2 = set(), set()
            if "neb" in options.analysis:
                s1 = set(map( extract_f, filter( filter_f, sites.mNEB.mPositiveSites)))
                for x in filter( filter_f, sites.mNEB.mPositiveSites):
                    if x.mResidue not in max_per_site: max_per_site[x.mResidue] = 0
                    max_per_site[x.mResidue] = max( x.mProbability, max_per_site[x.mResidue] )
                
            if "beb" in options.analysis:
                s2 = set(map( extract_f, filter( filter_f, sites.mBEB.mPositiveSites)))
                for x in filter( filter_f, sites.mBEB.mPositiveSites):
                    if x.mResidue not in max_per_site: max_per_site[x.mResidue] = 0
                    max_per_site[x.mResidue] = max( x.mProbability, max_per_site[x.mResidue] )


            s = s1.union( s2 )
            
            if first:
                total_sites = s
                first = False
            else:
                if selection_mode == "all":
                    total_sites = total_sites.union( s )
                elif selection_mode == "consistent":
                    total_sites = total_sites.intersection( s )
                elif selection_mode == "emes":
                    total_sites = total_sites.intersection( s )


    if selection_mode == "emes":
        if options.loglevel >= 2:
            options.stdlog.write("# before EMES filtering %i positive sites: mode %s, P>%5.2f\n" % (len(total_sites), selection_mode, 0.5) )
        
        # filter according to emes: maximum significance larger than 0.9
        total_sites = set( filter( lambda x: max_per_site[x] > 0.9, total_sites) )

        if options.loglevel >= 2:
            options.stdlog.write("# after EMES filtering %i positive sites: mode %s, P>%5.2f\n" % (len(total_sites), selection_mode, 0.9) )

    else:
        if options.loglevel >= 2:
            options.stdlog.write("# extracted %i positive sites: mode %s, P>%5.2f\n" % (len(total_sites), selection_mode, options.filter_probabiltiy) )

    if mali and options.filter_mali:
        if options.filter_mali == "gaps":
            nfiltered = 0
            mali_length = mali.getLength()

            column_data = map(lambda x: Mali.MaliData( x, gap_chars = "Nn", mask_chars = "-."), mali.getColumns())
            new_sites = set()

            for x in total_sites:

                ## PAML uses one-based coordinates
                column = column_data[x-1]

                if column.mNChars != mali_length:
                    nfiltered += 1
                    if options.loglevel >= 3:
                        options.stdlog.write("# rejected position %i due to mali\n" % x )
                    continue

                new_sites.add( x )

            total_sites = new_sites

        if options.loglevel >= 2:
            options.stdlog.write("# after MALI filtering %i positive sites\n" % (len(total_sites) ))

    return total_sites, max_per_site

def mapSites2Codons( total_sites, max_per_site ):
    """map sites to codons and return codons
    """
    
    max_per_codon = {}
    total_codons = set()
        
    for site in total_sites:
        codon, pos = divmod( site - 1, 3 )
        codon += 1
        total_codons.add( codon )
        if codon not in max_per_codon:
            max_per_codon[codon] = 0
        max_per_codon[codon] = max(max_per_codon[codon], max_per_site[site] )
            
    return total_codons, max_per_codon

def convertMali2Mali( mali ):
    """convert a mali to a profile."""

    new_mali = alignlib.makeMultipleAlignment()
    for id in mali.getIdentifiers():
        s = alignlib.makeAlignatumFromString( mali[id] )
        s.thisown = 0
        new_mali.addAlignatum( s )

    return new_mali

if __name__ == "__main__":
    
    parser = E.OptionParser( version = "%prog version: $Id: codemls2tsv.py 2781 2009-09-10 11:33:14Z andreas $" )

    parser.add_option("--methods", dest="methods", type="choice", action="append",
                      choices=("summary-numbers", "jalview",
                               "positive-site-table", "positive-site-list",
                               "count-positive-sites" ),
                      help="methods for analysis."  )

    parser.add_option("--selection-mode", dest="selection_mode", type="choice", 
                      choices=("all", "consistent", "emes" ),
                      help="how to select positive sites."  )

    parser.add_option("--prefix", dest="prefix", type="string",
                      help="prefix for rows.")

    parser.add_option("--pattern-input-filenames", dest="pattern_input_filenames", type="string",
                      help="input pattern."  )

    parser.add_option("--filter-probability", dest="filter_probability", type="float",
                      help="threshold for probability above which to include positive sites [default=%default]."  )

    parser.add_option("--filter-omega", dest="filter_omega", type="float",
                      help="threshold for omega above which to include positive sites [default=%default]."  )

    parser.add_option("--models", dest="models", type="string",
                      help="restrict output to set of site specific models."  )

    parser.add_option("--analysis", dest="analysis", type="string",
                      help="restrict output to set of analysis [beb|neb]."  )

    parser.add_option("--significance-threshold", dest="significance_threshold", type="float",
                      help="significance threshold for log-likelihood test."  )

    parser.add_option("--filter-mali", dest="filter_mali", type="choice",
                      choices=("none", "gaps"),
                      help="filter by mali to remove gapped positions.")

    parser.add_option("--filename-mali", dest="filename_mali", type="string",
                      help="filename with multiple alignment used for calculating sites - used for filtering" )

    parser.add_option("--filename-map-mali", dest="filename_map_mali", type="string",
                      help="filename with multiple alignment to map sites onto." )

    parser.add_option("--jalview-titles", dest="jalview_titles", type="string",
                      help="comma separated list of jalview annotation titles." )

    parser.add_option("--jalview-symbol", dest="jalview_symbol", type="string",
                      help="symbol to use in jalview." )

    parser.set_defaults(
        methods = [],
        prefix = None,
        filter_probability = 0,
        filter_omega = 0,
        models = "",
        analysis = "",
        significance_threshold = 0.05,
        selection_mode = "consistent",
        filename_mali = None,
        filename_map_mali = None,
        jalview_symbol = "*",
        jalview_titles = "",
        filter_mali = None,
        )

    (options, args) = E.Start( parser )

    if options.jalview_titles:
        options.jalview_titles = options.jalview_titles.split(",")
    else:
        options.jalview_titles = args    
        
    options.models = options.models.split(",")
    options.analysis = options.analysis.split(",")
    
    for a in options.analysis:
        if a not in ("beb", "neb"):
            raise "unknown analysis section: '%s', possible values are 'beb' and/or 'neb'" % a

    for a in options.models:
        if a not in ("8", "2", "3"):
            raise "unknown model: '%s', possible values are 2, 3, 8" % a
    
    codeml = WrapperCodeML.CodeMLSites()

    ## filter and extract functions
    filter_f = lambda x: x.mProbability >= options.filter_probability and x.mOmega >= options.filter_omega
    extract_f = lambda x: x.mResidue

    ## read multiple results
    results = []
    ninput, noutput, nskipped = 0, 0, 0

    headers = []
    for f in args:
        ninput += 1
        try:
            results.append(codeml.parseOutput( open(f, "r").readlines() ))
        except WrapperCodeML.UsageError:
            if options.loglevel >= 1:
                options.stdlog.write("# no input from %s\n" % f)
            nskipped += 1
            continue
        noutput += 1
        headers.append( f )
        
    ## map of nested model (key) to more general model
    map_nested_models = { '8' : '7',
                          '2' : '1',
                          '3' : '0' }

    if options.filename_mali:
        mali = Mali.Mali()
        mali.readFromFile( open(options.filename_mali, "r") )
    else:
        mali = None


    ###############################################################
    ###############################################################
    ###############################################################
    ## use multiple alignment to map residues to a reference mali
    ## or a sequence.
    ###############################################################        
    if options.filename_map_mali:
        
        if not mali:
            raise "please supply the input multiple alignment, if residues are to be mapped."

        ## translate the alignments
        def translate( s ):
            sequence = s.mString
            seq = []
            for codon in [ sequence[x:x+3] for x in range(0, len(sequence), 3) ]:
                aa = Genomics.MapCodon2AA( codon )
                seq.append( aa )

            s.mString = "".join(seq)


        tmali = Mali.Mali()
        tmali.readFromFile( open(options.filename_mali, "r" ) )
        tmali.apply( translate )
        
        tmap_mali = Mali.Mali()
        tmap_mali.readFromFile( open(options.filename_map_mali, "r") )

        if tmap_mali.getAlphabet() == "na":
            tmap_mali.apply( translate )
        
        map_old2new = alignlib.makeAlignmentVector()

        mali1 = alignlib.makeProfileFromMali( convertMali2Mali( tmali ) )

        if tmap_mali.getLength() == 1:
            
            s = tmap_mali.values()[0].mString
            mali2 = alignlib.makeSequence( s )
            ## see if you can find an identical subsequence and then align to thisD
            for x in tmali.values():
                if s in re.sub( "[- .]+", "", x.mString):
                    mali1 = alignlib.makeSequence( x.mString )
                    break
        else:
            mali2 = alignlib.makeProfileFromMali( convertMali2Mali( tmap_mali ) )        

        alignator = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_LOCAL, -10.0, -2.0 )
        alignator.align( map_old2new, mali1, mali2 )

        consensus = tmap_mali.getConsensus()
        
        if options.loglevel >= 4:
            options.stdlog.write( "# alphabet: %s\n" % tmap_mali.getAlphabet() )
            options.stdlog.write( "# orig  : %s\n" % tmali.getConsensus() )
            options.stdlog.write( "# mapped: %s\n" % consensus )
            options.stdlog.write( "# alignment: %s\n" % map_old2new.Write())
    else:
        map_old2new = None

    for method in options.methods:
        
        if method == "summary-numbers":

            options.stdlog.write( \
"""# Numbers of positive sites.
#
# The consistent row/column contains positive sites that are significant
# (above thresholds for probability and omega) for all models/analysis
# that have been selected (label: cons).
#
# The log-likelihood ratio test is performed for model pairs, depending
# on the output chosen.
# Significance threshold: %6.4f
# The pairs are 8 versus 7 and 2 versus 1 and 3 versus 0.
#
""" % options.significance_threshold )
            
            ## write header
            if options.prefix: options.stdout.write("prefix\t")
            
            options.stdout.write("method\tnseq\t" )
            h = []
            for model in options.models:
                for analysis in options.analysis:
                    h.append("%s%s" % (analysis, model) )
                h.append("p%s" % (model))
                h.append("df%s" % (model))
                h.append("chi%s" % (model))
                h.append("lrt%s" % (model))
                
            options.stdout.write("\t".join(h) )
            options.stdout.write("\tcons\tpassed\tfilename\n")

            nmethod = 0

            consistent_cols = [ None for x in range( len(options.analysis) ) ]
            passed_tests = {}
            for m in options.models: passed_tests[m] = 0
            
            for result in results:

                row_consistent = None

                if options.prefix:
                    options.stdout.write("%s" % (options.prefix ))
                
                options.stdout.write("%i" % nmethod)
                options.stdout.write("\t%i" % (result.mNumSequences ))
                                         
                npassed = 0
                
                for model in options.models:

                    sites = result.mSites[model]

                    ## do significance test
                    full_model, null_model = model, map_nested_models[model]
                    
                    lrt = Stats.doLogLikelihoodTest(
                        result.mSites[full_model].mLogLikelihood, 
                        result.mSites[full_model].mNumParameters, 
                        result.mSites[null_model].mLogLikelihood, 
                        result.mSites[null_model].mNumParameters, 
                        options.significance_threshold )

                    x = 0
                    for analysis in options.analysis:
                        
                        if analysis == "neb":
                            s = set(map( extract_f, filter( filter_f, sites.mNEB.mPositiveSites)))
                            
                        elif analysis == "beb":
                            s = set(map( extract_f, filter( filter_f, sites.mBEB.mPositiveSites)))                            
                            
                        options.stdout.write("\t%i" % ( len(s) ) )

                        if not lrt.mPassed:
                            s = set()
                            
                        if row_consistent == None:
                            row_consistent = s
                        else:
                            row_consistent = row_consistent.intersection( s )

                        if consistent_cols[x] == None:
                            consistent_cols[x] = s
                        else:
                            consistent_cols[x] = consistent_cols[x].intersection(s)

                        x += 1

                    if lrt.mPassed:
                        c = "passed"
                        passed_tests[model] += 1
                        npassed += 1
                    else:
                        c = "failed"
                        
                    options.stdout.write("\t%5.2e\t%i\t%5.2f\t%s" %\
                                         (lrt.mProbability, 
                                          lrt.mDegreesFreedom, 
                                          lrt.mChiSquaredValue, 
                                          c))
                    
                options.stdout.write("\t%i\t%i\t%s\n" % (len(row_consistent), npassed, headers[nmethod]))
                
                nmethod += 1
                
            if options.prefix:
                options.stdout.write("%s\t" % options.prefix )

            options.stdout.write("cons" )
            
            row_consistent = None
            total_passed = 0
            for model in options.models:

                x = 0
                
                for analysis in options.analysis:
                    
                    s = consistent_cols[x]
                    if s == None:
                        s = set()
                        
                    options.stdout.write("\t%i" % (len(s)))

                    if row_consistent == None:
                        row_consistent = s 
                    else:
                        row_consistent = row_consistent.intersection( s )

                    x += 1
                    
                options.stdout.write("\tna\t%i" % passed_tests[model] )
                total_passed += passed_tests[model]
                
            options.stdout.write("\t%i\t%i\n" % (len(row_consistent), total_passed) )

        elif method == "jalview":

            options.stdout.write("JALVIEW_ANNOTATION\n" )
            options.stdout.write("# Created: %s\n\n" % (time.asctime(time.localtime(time.time()))))

            l = 1
            x = 0
            for result in results:
                
                sites, significance = selectPositiveSites( [result], options.selection_mode, options, mali )

                codes = [""] * result.mLength

                if len(sites) == 0: continue

                for site in sites:
                    codes[site-1] = options.jalview_symbol

                options.stdout.write("NO_GRAPH\t%s\t%s\n" % (options.jalview_titles[x], "|".join( codes ) ))
                x += 1
                
        elif method == "count-positive-sites":

            sites, significance = selectPositiveSites( results, options.selection_mode, options, mali )

            options.stdout.write( "%i\n" % (len(sites) ))
            
        elif method in ( "positive-site-table", ):

            sites, significance = selectPositiveSites( results, options.selection_mode, options, mali )

            headers = ["site", "P"]
            if map_old2new:
                headers.append( "mapped" )
                headers.append( "Pm" )                
                
            options.stdout.write( "\t".join( headers ) + "\n" )

            sites = list(sites)
            sites.sort()
            nmapped, nunmapped = 0, 0
            for site in sites:
                values = [site, "%6.4f" % significance[site]]

                if map_old2new:
                    r = map_old2new.mapRowToCol(site)
                    if r == 0:
                        values.append("na")
                        values.append("")
                        nunmapped += 1
                        if options.loglevel >= 2:
                            options.stdlog.write("# unmapped residue: %i\n" % site )
                    else:
                        values.append(r)
                        values.append(consensus[r-1])
                        nmapped += 1
                        
                options.stdout.write( "\t".join( map(str, (values ))) + "\n" )

            if options.loglevel >= 1:
                options.stdlog.write("# sites: ninput=%i, noutput=%i, nskipped=%i\n" % (len(sites), nmapped, nunmapped))
                
    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped))
    
    E.Stop()

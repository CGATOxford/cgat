################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
genelist_analysis.py -analyse gene lists (GO, etc.)
===================================================

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

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E
import CGAT.GO as GO
import CGAT.IOTools as IOTools
import CGAT.Stats as Stats

def computeFDR( all_results, 
                qvalue_method = "storey" ):
    '''compute FDR.
    
    update GOResult structure with field .fdr
    '''
    
    # flatten all_results
    results = []
    for key, data in all_results.iteritems():
        results.extend( data.mResults.values() )
    
    observed_min_pvalues = [ min(x.mProbabilityOverRepresentation,
                                 x.mProbabilityUnderRepresentation) for x in results ]
    
    if qvalue_method == "storey":

        # compute fdr via Storey's method
        fdr_data = Stats.doFDR( observed_min_pvalues, vlambda = 0.1)

        E.info( "estimated proportion of true null hypotheses = %6.4f" % fdr_data.mPi0 ) 

        if fdr_data.mPi0 < 0.1:
            E.warn( "estimated proportion of true null hypotheses is less than 10%%  (%6.4f)" % fdr_data.mPi0) 
            
        for result, qvalue in zip( results, fdr_data.mQValues ):
            result.fdr = qvalue

    elif options.qvalue_method == "empirical":
        assert options.sample > 0, "requiring a sample size of > 0"
        raise NotImplementedError("empirical needs some work" )

def outputOntologyResults( all_results, go2infos, options ):
    '''output ontology results

    The following tables will be created:
    tests - overview of individual gene list results
    all - full output
    '''
    
    headers = ["ontology",
               "genelist",
               "code" ] +\
               GO.GOResult().getHeaders() +\
               GO.GOInfo().getHeaders()

    if options.fdr: 
        has_fdr = True
        headers += ["fdr"]

    outfile = options.stdout

    outfile.write("\t".join(headers) + "\n" )

    for testpair, data in all_results.iteritems():
        ontology, genelist = testpair
        for key, result in data.mResults.iteritems():
            code = GO.GetCode( result )

            n = go2infos[ontology].get( key, GO.GOInfo() )
            outfile.write( "\t".join( (ontology, genelist, code, str(result), str(n))) )
            if has_fdr: outfile.write("\t%e" % result.fdr)
            outfile.write("\n")

    
def doOntologyAnalysis( gene_lists, options ):
    '''do ontology analysis - requires options.filename_assignments to be set.'''
    
    
    E.info( "reading association of categories and genes from %s" % (options.filename_assignments) )
    gene2gos, go2infos = GO.ReadGene2GOFromFile( IOTools.openFile(options.filename_assignments ) )
    E.info( "read %i ontologies" % (len(gene2gos)) )

    #############################################################
    ## sort out which ontologies to test
    ontologies = options.ontology

    # test all if none specified
    if not ontologies: ontologies = gene2gos.keys()

    all_results = {}

    for ontology in ontologies:
        gene2go, go2info = gene2gos[ontology], go2infos[ontology]

        if len(go2info) == 0:
            E.warn( "could not find information for terms - could be mismatch between ontologies")

        ngenes, ncategories, nmaps = GO.CountGO( gene2go )        
        E.info( "%s: ontology assignments: %i genes mapped to %i categories (%i maps)" % 
                (ontology, ngenes, ncategories, nmaps) )

        for geneset, x in gene_lists.iteritems():
            foreground, background = x
            
            E.debug( "working on %s - %s" % (ontology, geneset))

            E.info( "%s - %s: (unfiltered) foreground=%i, background=%i" % (ontology, 
                                                                            geneset, 
                                                                            len(foreground), 
                                                                            len(background)))
            
            results = GO.AnalyseGO( gene2go, foreground, background )

            if len(results.mSampleGenes) == 0:
                E.warn( "%s - %s: no genes with GO categories - analysis aborted" % (ontology, geneset ) )
                continue

            # add sampling at this point for empirical FDR
            
            all_results[ (ontology,geneset) ] = results
    
    if options.fdr:
        E.info("computing the FDR with method %s" % options.qvalue_method )
        computeFDR( all_results, 
                    qvalue_method = options.qvalue_method )
    
    outputOntologyResults( all_results,
                           go2infos,
                           options )

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    parser = E.OptionParser( version = "%prog version: $Id: GO.py 2883 2010-04-07 08:46:22Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-s", "--species", dest="species", type="string",
                      help="species to use [default=%default]." )

    parser.add_option("-i", "--slims", dest="filename_slims", type="string",
                      help="filename with GO SLIM categories [default=%default].")

    parser.add_option( "-g", "--genes", dest="filename_genes", type="string",
                       help="filename with genes to analyse [default=%default]." )

    parser.add_option( "-f", "--format", dest="filename_format", type="choice",
                       choices = ("list", "matrix" ),
                       help="filename format [default=%default]." )

    parser.add_option( "-b", "--background", dest="filename_background", type="string",
                       help="filename with background genes to analyse [default=%default]." )

    parser.add_option( "-o", "--sort-order", dest="sort_order", type="choice",
                       choices=("fdr", "pover", "ratio" ),
                       help="output sort order [default=%default]." )

    parser.add_option( "--ontology", dest="ontology", type="choice", action="append",
                       choices=("biol_process","cell_location","mol_function", "mgi" ),
                       help="go ontologies to analyze. Ontologies are tested separately."
                       " [default=%default]." )

    parser.add_option( "-t", "--threshold", dest="threshold", type="float",
                       help="significance threshold [>1.0 = all ]. If --fdr is set, this refers to the fdr, otherwise it is a cutoff for p-values." )

    parser.add_option ("--filename-dump", dest="filename_dump", type="string",
                       help="dump GO category assignments into a flatfile [default=%default]." )

    parser.add_option ("--filename-ontology", dest="filename_ontology", type="string",
                       help="filename with ontology in OBO format [default=%default]." )

    parser.add_option ( "--filename-assignments", dest="filename_assignments", type="string",
                        help="read ontology assignments from a flatfile [default=%default]." )

    parser.add_option ( "--sample", dest="sample", type="int",
                        help="do sampling (with # samples) [default=%default]." )
    
    parser.add_option ( "--filename-output-pattern", dest = "output_filename_pattern", type="string",
                        help="pattern with output filename pattern (should contain: %(go)s and %(section)s ) [default=%default]")

    parser.add_option ( "--output-filename-pattern", dest = "output_filename_pattern", type="string",
                        help="pattern with output filename pattern (should contain: %(go)s and %(section)s ) [default=%default]")
    
    parser.add_option ( "--fdr", dest="fdr", action="store_true",
                       help="calculate and filter by FDR [default=%default]." )

    parser.add_option ( "--go2goslim", dest="go2goslim", action="store_true",
                       help="convert go assignments in STDIN to goslim assignments and write to STDOUT [default=%default]." )

    parser.add_option ( "--gene-pattern", dest = "gene_pattern", type="string",
                        help="pattern to transform identifiers to GO gene names [default=%default].")

    parser.add_option( "--filename-map-slims", dest="filename_map_slims", type="string",
                       help="write mapping between GO categories and GOSlims [default=%default].")

    parser.add_option( "--get-genes", dest="get_genes", type="string",
                       help="list all genes in the with a certain GOID [default=%default]." )

    parser.add_option( "--strict", dest="strict", action="store_true",
                       help="require all genes in foreground to be part of background. "
                       "If not set, genes in foreground will be added to the background [default=%default]." )

    parser.add_option("-q", "--qvalue-method", dest="qvalue_method", type="choice",
                      choices = ( "empirical", "storey" ),
                      help="method to perform multiple testing correction by controlling the fdr [default=%default]."  )

    parser.add_option( "--qvalue-lambda", dest="qvalue_lambda", type="float",
                       help="fdr computation: set lambda to fixed value [default=%default]."  )

    # parser.add_option( "--qvalue-pi0-method", dest="qvalue_pi0_method", type="choice",
    #                    choices = ("smoother", "bootstrap" ),
    #                    help="fdr computation: method for estimating pi0 [default=%default]."  )
    
    parser.set_defaults( species = None,
                         filename_genes = "-",
                         filename_format = "list",
                         filename_ontology = None,
                         filename_assignments = None,
                         filename_background = None,
                         filename_slims = None,
                         ontology = [],
                         filename_dump = None,
                         sample = 0,
                         fdr = False,
                         output_filename_pattern = None,
                         threshold = 0.05,
                         filename_map_slims = None,
                         gene_pattern = None,
                         sort_order = "ratio",
                         get_genes = None,
                         strict = False,
                         qvalue_method = "empirical",
                         qvalue_lambda = None,
                         )

    (options, args) = E.Start( parser, argv = argv )

    # collect all gene lists
    # gene lists are tuples of fg/bg.
    # If bg is None, the full list of genes is 
    # taken as the background.
    gene_lists = {}
    bg = None

    if options.filename_format == "list":
        gene_lists[ "default"] = (GO.ReadGeneList( options.filename_genes, 
                                                   gene_pattern = options.gene_pattern ),
                                  None )

    elif options.filename_format == "matrix":
        bg, genes = GO.ReadGeneLists( options.filename_genes,
                                  gene_pattern = options.gene_pattern ) 

        # use default background        
        gene_lists = dict( [ (x, (y,None)) for x,y in genes.iteritems() ] )

    if bg: E.info("gene list: %20s = %i genes" % ("default background", len(bg)))

    E.info( "read %i gene list from %s" % (len(gene_lists), options.filename_genes) )
    if len(gene_lists) == 0:
        raise ValueError( "no gene lists input")

    # set the default background for those tests without explicit background
    for key in gene_lists.keys():
        if gene_lists[key][1] == None:
            gene_lists[key] = (gene_lists[key][0], bg )

    # check for inconsistencies
    for key, l in gene_lists.iteritems():
        foreground, background = l

        missing = set(foreground).difference( set(background))

        if options.strict:
            assert len(missing) == 0, \
                "%s: %i genes in foreground but not in background: %s" % (key, len(missing), str(missing))
        else:
            if len(missing) != 0:
                E.warn( "%s: %i genes in foreground that are not in background - added to background of %i" %\
                            (key, len(missing), len(background)) )
            background.extend( missing )


    for key, l in gene_lists.iteritems():
        lfg = len(l[0])
        lbg =len(l[1])

        E.debug( "gene list: %20s = %i genes / %i background = %5.2f %%" % 
                 (key, lfg, lbg, 100.0 * lfg / lbg ) )


    ## read gene list assignments
    if options.filename_assignments:
        doOntologyAnalysis( gene_lists, options )


    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

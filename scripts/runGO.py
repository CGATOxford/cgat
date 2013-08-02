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
GO.py - compute GO enrichment from gene lists
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Usage
-----

The script ``GO.py`` will test for enrichment or depletion 
of GO categories within a gene list.

The script uses a hypergeometric test to check if a particular
GO category is enriched in a foreground set with respect
to a background set. Multiple testing is controlled by 
computing a empirical false discovery rate using a sampling 
procedure.

A GO analysis proceeds in three steps:

   1. building gene to GO assignments
   2. create one or more gene lists with foreground and background
   3. run one or more GO analyses for each of the foreground gene lists

This script analyses multiple gene lists in parallel when a matrix of
gene lists is provided. If multiple gene lists are provided, the FDR is 
controlled per gene list and not overall. However, my intuition is that
if the number of tests is large the results should be comparable as if 
the FDR was controlled globally, though I have no proof for this.

Building gene to GO assignments
+++++++++++++++++++++++++++++++

The easiest way to obtain a map from gene identifiers to GO assignments
is to down download GO assignments from the ENSEMBL database. The command
below will download go assignments for the human gene set
and save it in the file :file:`gene2go.data`::

   python GO.py 
      --filename-dump=gene2go.data
      --host=ensembldb.ensembl.org
      --user=anonymous 
      --database=homo_sapiens_core_54_36p
      --port=5306 
   > gene2go.log

In order to use GOslim categories, an additional mapping step needs to be performed.
The sequence of commands is::

    wget http://www.geneontology.org/GO_slims/goslim_goa.obo
    wget http://www.geneontology.org/ontology/gene_ontology.obo
    map2slim -outmap go2goslim.map goslim_goa.obo gene_ontology.obo
    python GO.py 
                --go2goslim 
                --filename-ontology=gene_ontology.obo 
                --slims=go2goslim.map 
                --log=goslim.log 
        < gene2go.data > gene2goslim.data

The first two commands obtain GOslim information. 
`map2slim <http://search.cpan.org/~cmungall/go-perl/scripts/map2slim>`_
is part of Chris Mungall's `go-perl <http://search.cpan.org/~cmungall/go-perl/>`_ module 
and the last command converts the gene-to-GO assignment into gene-to-GOSlim assignments.

The gene-to-GO mapping can be constructed any other way. It is simply 
a table of tab-separated values::

   go_type gene_id go_id   description     evidence
   biol_process    ENSG00000151729 GO:0000002      mitochondrial genome maintenance        NA
   biol_process    ENSG00000025708 GO:0000002      mitochondrial genome maintenance        NA
   biol_process    ENSG00000115204 GO:0000002      mitochondrial genome maintenance        NA
   ...

Building gene lists
+++++++++++++++++++

GO requires a list of genes to test for enrichment. This list is simply
a table with one column of gene identifiers. For example::

   gene_id
   ENSG00000116586
   ENSG00000065809
   ENSG00000164048
   ENSG00000115137
   ENSG00000121210

Alternatively, the gene list can be a multi-column table such as::

   gene_id             dataset1    dataset2
   ENSG00000116586     1           0
   ENSG00000065809     0           0
   ENSG00000164048     1           0
   ENSG00000115137     1           1
   ENSG00000121210     0           1

In this case, enrichment is computed for multiple datasets at once. Make sure
to add the ``%(set)s`` place holder to ``--filename-output-pattern``.

If no background is given, all genes that have GO assignments will constitute
the background. 

Statistics
++++++++++

Enrichment is computed using the hypergeometric test. 

.. todo::
    * apply filtering
    * more stats
    * more FDR

Running the GO analysis
+++++++++++++++++++++++

The command below runs a GO analysis, computing an FDR using 10.000 samples::

    python GO.py 
        --filename-input=gene2go.data
        --genes=foreground
        --background=background 
        --sample=10000
        --fdr 
        --filename-ontology=gene_ontology.obo
        --output-filename-pattern='result/%(set)s.%(go)s.%(section)s' 
   > go.log

The output will be stored in the directory :file:`result` and output files will be
created according to the pattern ``<set>.<go>.<section>``. ``<set>`` is the gene set
that is analysed, ``<go>`` is one of ``biol_process``, ``mol_function`` and ``cell_location``.
``<section>`` denotes the file contents. Files output are:

+------------+----------------------------------------------+
|``section`` | contents                                     |
+------------+----------------------------------------------+
|samples     |sampling statistics                           |
+------------+----------------------------------------------+
|overall     |table with full results                       |
+------------+----------------------------------------------+
|results     |table with only the significant results       |
+------------+----------------------------------------------+
|parameters  |input and sampling parameters                 |
+------------+----------------------------------------------+
|fg          |assigments for genes in the foreground set    |
+------------+----------------------------------------------+

Other options
+++++++++++++

The script can accept other ontologies than just GO ontologies. 

Code
----

'''
import os, sys, string, re, getopt, time, optparse, math, tempfile, subprocess, random
import collections

import scipy
import scipy.stats
import scipy.special
import numpy
import CGAT.Stats as Stats
import CGAT.Database as Database
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.CSV as CSV
from CGAT.GO import *

##---------------------------------------------------------------------------    
def main():

    parser = E.OptionParser( version = "%prog version: $Id: GO.py 2883 2010-04-07 08:46:22Z andreas $", usage = globals()["__doc__"])

    dbhandle = Database.Database()
    
    parser.add_option("-s", "--species", dest="species", type="string",
                      help="species to use [default=%default]." )

    parser.add_option("-i", "--slims", dest="filename_slims", type="string",
                      help="filename with GO SLIM categories [default=%default].")

    parser.add_option( "-g", "--genes", dest="filename_genes", type="string",
                       help="filename with genes to analyse [default=%default]." )

    parser.add_option( "-b", "--background", dest="filename_background", type="string",
                       help="filename with background genes to analyse [default=%default]." )

    parser.add_option( "-m", "--minimum-counts", dest="minimum_counts", type="int",
                       help="minimum count - ignore all categories that have fewer than # number of genes"
                            " [default=%default]." )

    parser.add_option( "-o", "--sort-order", dest="sort_order", type="choice",
                       choices=("fdr", "pvalue", "ratio" ),
                       help="output sort order [default=%default]." )

    parser.add_option( "--ontology", dest="ontology", type="string", action="append",
                       help="go ontologies to analyze. Ontologies are tested separately."
                       " [default=%default]." )

    parser.add_option( "-t", "--threshold", dest="threshold", type="float",
                       help="significance threshold [>1.0 = all ]. If --fdr is set, this refers to the fdr, otherwise it is a cutoff for p-values." )

    parser.add_option ("--filename-dump", dest="filename_dump", type="string",
                       help="dump GO category assignments into a flatfile [default=%default]." )

    parser.add_option ("--filename-gene2name", dest="filename_gene2name", type="string",
                       help="optional filename mapping gene identifiers to gene names [default=%default]." )

    parser.add_option ("--filename-ontology", dest="filename_ontology", type="string",
                       help="filename with ontology in OBO format [default=%default]." )

    parser.add_option ( "--filename-input", dest="filename_input", type="string",
                       help="read GO category assignments from a flatfile [default=%default]." )

    parser.add_option ( "--sample", dest="sample", type="int",
                        help="do sampling (with # samples) [default=%default]." )

    parser.add_option ( "--filename-output-pattern", "--output-filename-pattern", 
                        dest = "output_filename_pattern", type="string",
                        help="pattern with output filename pattern (should contain: %(go)s and %(section)s ) [default=%default]")

    parser.add_option ( "--fdr", dest="fdr", action="store_true",
                        help="calculate and filter by FDR [ReadGene2GOFromFiledefault=%default]." )
    
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
                      choices = ( "empirical", "storey", "BH" ),
                      help="method to perform multiple testing correction by controlling the fdr [default=%default]."  )


    parser.add_option( "--pairwise", dest="compute_pairwise", action="store_true",
                       help="compute pairwise enrichment for multiple gene lists. "
                       "[default=%default]." )

    # parser.add_option( "--qvalue-lambda", dest="qvalue_lambda", type="float",
    #                   help="fdr computation: lambda [default=%default]."  )

    # parser.add_option( "--qvalue-pi0-method", dest="qvalue_pi0_method", type="choice",
    #                    choices = ("smoother", "bootstrap" ),
    #                    help="fdr computation: method for estimating pi0 [default=%default]."  )
    
    parser.set_defaults( species = None,
                         filename_genes = "-",
                         filename_background = None,
                         filename_slims = None,
                         minimum_counts = 0,
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
                         pairs_min_observed_counts = 3,
                         compute_pairwise = False,
                         filename_gene2name = None
                         )

    (options, args) = E.Start( parser, add_mysql_options = True )

    if options.go2goslim:
        convertGo2Goslim( options )
        E.Stop()
        sys.exit(0)

    if options.fdr and options.sample == 0:
        E.warn( "fdr will be computed without sampling" )
        
    #############################################################
    ## dump GO
    if options.filename_dump:
        # set default orthologies to GO
        if not options.ontology:
            options.ontology = ["biol_process", "mol_function", "cell_location"] 

        E.info( "dumping GO categories to %s" % (options.filename_dump) )

        dbhandle.Connect( options )
            
        outfile = IOTools.openFile( options.filename_dump, "w", create_dir = True )
        DumpGOFromDatabase( outfile,
                            dbhandle,
                            options )
        outfile.close()
        E.Stop()
        sys.exit(0)

    #############################################################
    ## read GO categories from file
    if options.filename_input:
        E.info( "reading association of categories and genes from %s" % (options.filename_input) )
        infile = IOTools.openFile(options.filename_input)
        gene2gos, go2infos = ReadGene2GOFromFile( infile )
        infile.close()

    if options.filename_gene2name:
        E.info("reading gene identifier to gene name mapping from %s" % options.filename_gene2name)
        infile = IOTools.openFile( options.filename_gene2name)
        gene2name = IOTools.readMap( infile, has_header = True )
        infile.close()
        E.info("read %i gene names for %i gene identifiers" % (len(set(gene2name.values())),
                                                               len(gene2name)))
    else:
        gene2name = None

    #############################################################
    ## read GO ontology from file
    if options.filename_ontology:
        E.info( "reading ontology from %s" % (options.filename_ontology) )
        
        infile = IOTools.openFile(options.filename_ontology)
        ontology = readOntology( infile )
        infile.close()
        
        def _g():
            return collections.defaultdict( GOInfo )
        go2infos = collections.defaultdict( _g );

        ## substitute go2infos
        for go in ontology.values():
            go2infos[go.mNameSpace][go.mId] = GOInfo( go.mId,
                                                      go_type = go.mNameSpace,
                                                      description = go.mName )

    #############################################################
    ## get foreground gene list
    input_foreground, genelists = ReadGeneLists( options.filename_genes, 
                                                 gene_pattern = options.gene_pattern )

    E.info( "read %i genes for forground in %i gene lists" % (len(input_foreground), len(genelists)) )

    #############################################################
    ## get background
    if options.filename_background:
        
        # nick - bug fix: background is the first tuple element from ReadGeneLists
        input_background = ReadGeneLists( options.filename_background, 
                                          gene_pattern = options.gene_pattern )[0]
        E.info( "read %i genes for background" % len(input_background) )
    else:
        input_background = None

    #############################################################
    ## sort out which ontologies to test
    if not options.ontology: 
        if options.filename_input:
            options.ontology = gene2gos.keys()

    E.info( "found %i ontologies: %s" % (len(options.ontology), options.ontology))

    outfile_summary = options.stdout
    outfile_summary.write( "\t".join( (
                "genelist",
                "ontology",
                "significant",
                "threshold",
                "ngenes",
                "ncategories",
                "nmaps",
                "nforegound",
                "nforeground_mapped",
                "nbackground",
                "nbackground_mapped",
                "nsample_counts",
                "nbackground_counts",
                "psample_assignments",
                "pbackground_assignments") ) + "\n" )

    #############################################################
    ## get go categories for genes
    for test_ontology in options.ontology:

        # store results for aggregate output of multiple gene lists
        all_results = []
        all_significant_results = []
        all_genelists_with_results = []

        E.info( "working on ontology %s" % test_ontology )
        #############################################################
        ## get/read association of GO categories to genes
        if options.filename_input:
            gene2go, go2info = gene2gos[test_ontology], go2infos[test_ontology]
        else:
            E.info( "reading data from database ..." )

            dbhandle.Connect( options )
            gene2go, go2info = ReadGene2GOFromDatabase( dbhandle,
                                                        test_ontology,
                                                        options.database, options.species )

            E.info( "finished" )

        if len(go2info) == 0:
            E.warn( "could not find information for terms - could be mismatch between ontologies")

        ngenes, ncategories, nmaps, counts_per_category = CountGO( gene2go )        
        E.info( "assignments found: %i genes mapped to %i categories (%i maps)" % (ngenes, ncategories, nmaps) )

        if options.minimum_counts > 0:
            to_remove = set([ x for x,y in counts_per_category.iteritems() if y < options.minimum_counts ])
            E.info("removing %i categories with less than %i genes" % (len(to_remove), options.minimum_counts ) )
            removeCategories( gene2go, to_remove )

            ngenes, ncategories, nmaps, counts_per_category = CountGO( gene2go )        
            E.info( "assignments after filtering: %i genes mapped to %i categories (%i maps)" % (ngenes, ncategories, nmaps) )

        for genelist_name, foreground in genelists.iteritems():

            msgs = []
            E.info("processing %s with %i genes" % (genelist_name, len(foreground)))
            ##################################################################
            ##################################################################
            ##################################################################
            ## build background - reconcile with foreground
            ##################################################################
            if input_background == None:
                background = list(gene2go.keys())
            else:
                background = list(input_background)

            # nick - bug-fix backgorund included the foreground in a tuple. 
            # background is the first tuple element 
            missing = foreground.difference(set(background))

            if options.strict:
                assert len(missing) == 0, \
                    "%i genes in foreground but not in background: %s" % (len(missing), str(missing))
            else:
                if len(missing) != 0:
                    E.warn( "%i genes in foreground that are not in background - added to background of %i" %\
                                (len(missing), len(background)) )

                background.extend(missing)

            E.info( "(unfiltered) foreground=%i, background=%i" % (len(foreground), len(background)))

            #############################################################
            ## sanity checks:            
            ## are all of the foreground genes in the dataset
            ## missing = set(genes).difference( set(gene2go.keys()) )
            ## assert len(missing) == 0, "%i genes in foreground set without GO annotation: %s" % (len(missing), str(missing))

            #############################################################            
            ## read GO slims and map GO categories to GO slim categories
            if options.filename_slims:
                go_slims = GetGOSlims( IOTools.openFile(options.filename_slims, "r") )

                if options.loglevel >=1:
                    v = set()
                    for x in go_slims.values():
                        for xx in x: v.add(xx)
                    options.stdlog.write( "# read go slims from %s: go=%i, slim=%i\n" %\
                                              ( options.filename_slims,
                                                len(go_slims), 
                                                len( v ) ))



                if options.filename_map_slims:
                    if options.filename_map_slims == "-":
                        outfile = options.stdout
                    else:
                        outfile=IOTools.openFile(options.filename_map_slims, "w" )

                    outfile.write( "GO\tGOSlim\n" )
                    for go, go_slim in go_slims.items():
                        outfile.write("%s\t%s\n" % (go, go_slim))

                    if outfile != options.stdout:
                        outfile.close()

                gene2go = MapGO2Slims( gene2go, go_slims, ontology = ontology )

                if options.loglevel >=1:
                    ngenes, ncategories, nmaps, counts_per_category = CountGO( gene2go )
                    options.stdlog.write( "# after go slim filtering: %i genes mapped to %i categories (%i maps)\n" % (ngenes, ncategories, nmaps) )

            #############################################################
            ## Just dump out the gene list
            if options.get_genes:
                fg, bg, ng = [], [], []

                for gene, vv in gene2go.items():
                    for v in vv:
                        if v.mGOId == options.get_genes:
                            if gene in genes:
                                fg.append( gene )
                            elif gene in background:
                                bg.append( gene )
                            else:
                                ng.append( gene )

                ## skip to next GO class
                if not (bg or ng): continue

                options.stdout.write( "# genes in GO category %s\n" % options.get_genes )
                options.stdout.write( "gene\tset\n" )
                for x in fg: options.stdout.write("%s\t%s\n" % ("fg", x))
                for x in bg: options.stdout.write("%s\t%s\n" % ("bg", x))           
                for x in ng: options.stdout.write("%s\t%s\n" % ("ng", x))                       

                E.info( "nfg=%i, nbg=%i, nng=%i" % (len(fg), len(bg), len(ng) ))

                E.Stop()
                sys.exit(0)

            #############################################################
            outfile = getFileName( options, 
                                   go = test_ontology,
                                   section = 'foreground',
                                   set = genelist_name )

            outfile.write ("gene_id\n%s\n" % ("\n".join( sorted( foreground) ) ) )
            if options.output_filename_pattern:
                outfile.close()

            outfile = getFileName( options, 
                                   go = test_ontology,
                                   section = 'background',
                                   set = genelist_name )

            outfile.write ("gene_id\n%s\n" % ("\n".join( sorted( background[0]) ) ) )
            if options.output_filename_pattern:
                outfile.close()

            #############################################################
            ## do the analysis
            go_results = AnalyseGO( gene2go, foreground, background )

            if len(go_results.mSampleGenes) == 0:
                E.warn( "%s: no genes with GO categories - analysis aborted" % genelist_name)
                continue

            pairs = go_results.mResults.items()

            #############################################################
            ## calculate fdr for each hypothesis
            if options.fdr:
                fdrs, samples, method  = computeFDRs( go_results, options, test_ontology )
                for x,v in enumerate(pairs):
                    v[1].mQValue = fdrs[v[0]][0]
            else:
                fdrs, samples, method = {}, None
                
            msgs.append( "fdr=%s" % method)

            if options.sort_order == "fdr":
                pairs.sort( lambda x, y: cmp(x[1].mQValue, y[1].mQValue))           
            elif options.sort_order == "ratio":
                pairs.sort( lambda x, y: cmp(x[1].mRatio, y[1].mRatio))
            elif options.sort_order == "pvalue":
                pairs.sort( lambda x, y: cmp(x[1].mPValue, y[1].mPValue))

            #############################################################
            #############################################################
            #############################################################
            ## output the full result
            outfile = getFileName( options, 
                                   go = test_ontology,
                                   section = 'overall',
                                   set = genelist_name )

            outputResults( outfile, pairs, go2info, options, fdrs = fdrs, samples = samples )

            if options.output_filename_pattern:
                outfile.close()

            #############################################################
            #############################################################
            #############################################################
            # filter significant results and output
            filtered_pairs = selectSignificantResults( pairs, fdrs, options )

            nselected = len(filtered_pairs)
            nselected_up = len( [x for x in filtered_pairs if x[1].mRatio > 1 ] )
            nselected_down = len( [x for x in filtered_pairs if x[1].mRatio < 1 ] )
            
            assert nselected_up + nselected_down == nselected

            outfile = getFileName( options, 
                                   go = test_ontology,
                                   section = 'results',
                                   set = genelist_name )

            outputResults( outfile, 
                           filtered_pairs, 
                           go2info, 
                           options,
                           fdrs = fdrs, 
                           samples = samples )
            
            if options.output_filename_pattern:
                outfile.close()

            #############################################################
            #############################################################
            #############################################################
            # save results for multi-gene-list analysis
            all_results.append( pairs )
            all_significant_results.append( filtered_pairs )
            all_genelists_with_results.append( genelist_name )

            #############################################################
            #############################################################
            #############################################################
            ## output parameters
            ngenes, ncategories, nmaps, counts_per_category = CountGO( gene2go )

            outfile = getFileName( options, 
                                   go = test_ontology,
                                   section = 'parameters',
                                   set = genelist_name )

            nbackground = len(background)
            if nbackground == 0:
                nbackground = len(go_results.mBackgroundGenes)

            outfile.write( "# input go mappings for gene list '%s' and category '%s'\n" % (genelist_name, test_ontology ))
            outfile.write( "parameter\tvalue\tdescription\n" )
            outfile.write( "mapped_genes\t%i\tmapped genes\n" % ngenes )
            outfile.write( "mapped_categories\t%i\tmapped categories\n" % ncategories )
            outfile.write( "mappings\t%i\tmappings\n" % nmaps )
            outfile.write( "genes_in_fg\t%i\tgenes in foreground\n" % len(foreground) )
            outfile.write( "genes_in_fg_with_assignment\t%i\tgenes in foreground with GO assignments\n" % (len(go_results.mSampleGenes)) )
            outfile.write( "genes_in_bg\t%i\tinput background\n" % nbackground )
            outfile.write( "genes_in_bg_with_assignment\t%i\tgenes in background with GO assignments\n" % (len(go_results.mBackgroundGenes)) )
            outfile.write( "associations_in_fg\t%i\tassociations in sample\n"     % go_results.mSampleCountsTotal )
            outfile.write( "associations_in_bg\t%i\tassociations in background\n" % go_results.mBackgroundCountsTotal )
            outfile.write( "percent_genes_in_fg_with_association\t%s\tpercent genes in sample with GO assignments\n" % (IOTools.prettyPercent( len(go_results.mSampleGenes) , len(foreground), "%5.2f" )))
            outfile.write( "percent_genes_in_bg_with_associations\t%s\tpercent genes background with GO assignments\n" % (IOTools.prettyPercent( len(go_results.mBackgroundGenes), nbackground, "%5.2f" )))
            outfile.write( "significant\t%i\tsignificant results reported\n" % nselected )
            outfile.write( "significant_up\t%i\tsignificant up-regulated results reported\n" % nselected_up )
            outfile.write( "significant_down\t%i\tsignificant up-regulated results reported\n" % nselected_down )
            outfile.write( "threshold\t%6.4f\tsignificance threshold\n" % options.threshold )        

            if options.output_filename_pattern:
                outfile.close()

            if outfile_summary:
                outfile_summary.write( "\t".join( map(str, ( \
                                genelist_name,
                                test_ontology,
                                nselected,
                                options.threshold,
                                ngenes,
                                ncategories,
                                nmaps,
                                len(foreground),
                                len(go_results.mSampleGenes),
                                nbackground,
                                len(go_results.mBackgroundGenes),
                                go_results.mSampleCountsTotal,
                                go_results.mBackgroundCountsTotal,
                                IOTools.prettyPercent( len(go_results.mSampleGenes) , len(foreground), "%5.2f" ),
                                IOTools.prettyPercent( len(go_results.mBackgroundGenes), nbackground, "%5.2f" ),
                                ",".join( msgs) ) ) ) + "\n" )

            #############################################################
            #############################################################
            #############################################################
            ## output the fg patterns
            outfile = getFileName( options, 
                                   go = test_ontology,
                                   section = 'withgenes',
                                   set = genelist_name )

            outputResults( outfile, pairs, go2info, options, 
                           fdrs = fdrs, 
                           gene2go = gene2go,
                           foreground = foreground,
                           gene2name = gene2name )

            if options.output_filename_pattern:
                outfile.close()

        if len(genelists) > 1:

            ######################################################################
            ######################################################################
            ######################################################################
            ## output various summary files
            ## significant results
            outputMultipleGeneListResults( all_significant_results, 
                                           all_genelists_with_results, 
                                           test_ontology, 
                                           go2info,
                                           options,
                                           section = 'significant')

            ## all results
            outputMultipleGeneListResults( all_results, 
                                           all_genelists_with_results, 
                                           test_ontology, 
                                           go2info,
                                           options,
                                           section = 'all')

            
            if options.compute_pairwise:
                pairwiseGOEnrichment( all_results,
                                      all_genelists_with_results,
                                      test_ontology,
                                      go2info,
                                      options )

    E.Stop()


if __name__ == "__main__":
    sys.exit( main() )

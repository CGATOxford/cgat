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
xrate_blocks.py - 
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

   python xrate_blocks.py --help

Type::

   python xrate_blocks.py --help

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
import math
import time
import random
import types
import tempfile

USAGE="""python %s [OPTIONS]

apply blocks models to sequence data
""" % sys.argv[0]

import CGAT.Experiment as Experiment
import CGAT.Mali as Mali
import CGAT.Genomics as Genomics
import CGAT.RateEstimation as RateEstimation
import CGAT.TreeTools as TreeTools
import CGAT.IOTools as IOTools

from XGram.Generator.Prebuilt import DNA
from XGram.Model import Annotation
import XGram.Run
import Bio.Data.CodonTable

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def writeModel( grammar, section, options):
    """write a model to output file."""
    if section in options.write or "all" in options.write:
        outfile = open( options.output_pattern % section, "w" )
        outfile.write( "%s\n" % grammar.getGrammar())
        outfile.close()

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def prepareGrammar( options ):
    """prepare grammar for custom grammars."""
    
    num_blocks = options.num_blocks

    labels = string.letters.upper()
    annotate_terminals = {}
    for x in range(num_blocks):
        annotations = []
        key = "B%i" % x
        annotations.append( Annotation( row = "STATE",
                                        column = key,
                                        label = labels[x % len(labels)] ))

        annotate_terminals[ key ] = annotations

    input_model = DNA.buildModel( substitution_model = "gtr",
                                  num_blocks = num_blocks,
                                  grammar_type = options.grammar_type,
                                  shared_frequencies = False,
                                  shared_rates = False,
                                  annotate_terminals = annotate_terminals,
                                  )

    rate = 0.2
    for x in range(num_blocks):
        for param in ("alpha", "beta", "gamma", "delta", "theta", "epsilon"):
            p = "B%i_%s" % (x,param)
            input_model.mGrammar.removeParameter( p )
            input_model.mGrammar.addParameter( (p, rate), is_explicit = True )
        rate += 0.1

    grammar = input_model.mGrammar.mRules
    pseudononterminals = dict( [ ( ("NT_B%i*" % x,), ("NT_B%i" % x,) ) for x in range(num_blocks) ] )
    
    prob_same = options.probability_block
    prob_diff = (1.0 - prob_same) / len(pseudononterminals)

    for source in pseudononterminals.keys():
        mapped_source = pseudononterminals[source]
        for target,rule in grammar[source].items():
            if target == mapped_source:
                rule.mRate = (rule.mRate[0], prob_same )
            else:
                rule.mRate = (rule.mRate[0], prob_diff )

    writeModel( input_model, "input", options )

    return input_model

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def prepareMali( mali, tree, map_old2new, options ):
    
    ids = mali.getIdentifiers()

    fh, filename = tempfile.mkstemp()

    os.close(fh)
    outfile = open(filename, "w" )
    
    if tree:
        tree.rescaleBranchLengths( 1.0 )
        tree_options = "#=GF NH %s" % tree.to_string( branchlengths_only=True, format="nh")
    elif mali.getNumSequences() == 2:
        tree_options = "#=GF NH (%s:1.0)%s;" % tuple(map_old2new.values())
    else:
        raise "Please supply a tree."

    mali.writeToFile( outfile, 
                      format="stockholm",
                      write_ranges = False,
                      options = ( tree_options, ) )
    outfile.close()
    
    return filename

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def run( input_model, filename, train, options ):
    
    t1 = time.time()

    xgram = XGram.XGram()

    if options.xrate_min_increment:
        xgram.setMinIncrement( options.xrate_min_increment )

    if train:
        result = xgram.train( input_model, filename )
    else:
        result = xgram.annotate( input_model, filename )

    if options.dump:
        options.stdlog.write( "### input multiple alignment\n" )
        options.stdlog.write( "".join(open(filename).readlines()))
        options.stdlog.write( "### output xrate data\n" )
        options.stdlog.write( "".join(result.mData) )
        options.stdlog.write( "### output xrate log \n" )
        options.stdlog.write( "".join(result.mLog) )

    t2 = time.time()

    trained_model = result.getModel()

    writeModel( trained_model, "trained", options )

    return result

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def trainOnMali( input_model, filename_mali, options ):
    """train a grammar using an annotated mali."""

    result = run( input_model, filename_mali, True, options )

    return result

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def trainMali( mali, options ):
    """train a grammar on a multiple alignment."""

    ## remove empty columns and masked columns
    if options.clean_mali:
        mali.mGapChars = mali.mGapChars + ("n", "N")
        mali.removeGaps( minimum_gaps = 1, frame=1 )
    
    length = mali.getNumColumns()

    input_model = prepareGrammar( options )

    for id in mali.getIdentifiers():
        if options.separator in id:
            species = id.split(options.separator)[0]
            mali.rename( id, species )

    map_new2old = mali.mapIdentifiers()
    map_old2new = IOTools.getInvertedDictionary( map_new2old, make_unique = True )
    
    ids = mali.getIdentifiers()

    if options.input_filename_tree:
        nexus = TreeTools.Newick2Nexus( open(options.input_filename_tree,"r") )
        tree = nexus.trees[0]
        try:
            tree.relabel( map_old2new, warn = True )
        except KeyError, msg:
            raise KeyError( "names in mali and tree are not congruent: %s" % msg )
    else:
        tree = None

    filename_mali = prepareMali( mali, tree, map_old2new, options )
    
    result = trainOnMali( input_model, filename_mali, options )

    options.stdout.write( "%f\t%i" % (result.getLogLikelihood(), result.getNumIterations() ) )
    outputRates( result, options )
    outputAnnotations( result, options )
    
    options.stdout.write("\n")

    return result

##-----------------------------------------------------------------------------------
def outputAnnotations( result, options ):
    """output the annotations in the model."""

    mali = Mali.Mali()

    mali.readFromFile( result.getData(), format="stockholm" )
    annotation = mali.getAnnotation( "STATE" )
    
    l,c,f = 0, None, []
    for x in annotation:
        if x != c:
            if c: f.append( "%s:%i" % (c,l) )
            c = x
            l = 0
        l += 1
        
    f.append( "%s:%i" % (c,l) )

    options.stdout.write( "\t%s" % ",".join(f))

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def outputRates( result, options ):
    """output rates in a grammar."""

    trained_model = result.getModel()

    pis = trained_model.evaluateTerminalFrequencies()
    matrices = trained_model.evaluateRateMatrix()
    terminals = pis.keys()

    for terminal in terminals:
        Q, distance = RateEstimation.getDistanceGTR( pis[terminal], matrices[terminal] )
        options.stdout.write("\t%s" % (options.value_format % distance ) )
        
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def annotateMali( input_model, filename_mali, options ):
    """annotate a multiple alignment using a trained grammar."""

    result = run( input_model, filename_mali, True, options )

    return result

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------

if __name__ == "__main__":
    
    parser = optparse.OptionParser( version = "%prog version: $Id: xrate_blocks.py 2781 2009-09-10 11:33:14Z andreas $" )

    parser.add_option("--input-filename-tree", dest="input_filename_tree", type="string",
                       help="""filename with pyhlogenetic tree.""" )
     
    parser.add_option("--dump", dest="dump", action="store_true",
                      help="dump raw output [%default].")

    parser.add_option("-m", "--method", dest="method", type="choice",
                      choices=("train", "annotate" ),
                      help="method to apply [%default]."  )

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal", "stockholm", "phylip" ),
                      help="input format of multiple alignment"  )

    parser.add_option( "--no-clean", dest="clean_mali", action="store_false",
                      help="do not clean multiple alignment before submitting to codeml. It might take too long for very large sequences.")

    parser.add_option("-w", "--write", dest="write", type="choice", action="append",
                      choices=("input_fixed", "trained_fixed", "input_variable", "trained_variable", "all" ),
                      help="output sections to write."  )

    parser.add_option("-o", "--output-pattern", dest="output_pattern", type="string",
                      help="output pattern for output files." )

    parser.add_option( "--xrate-min-increment", dest="xrate_min_increment", type = "float",
                       help="minimum increment to stop iteration in xrate." )    

    parser.add_option( "--num-blocks", dest="num_blocks", type = "int",
                       help="number of blocks to create in grammar." )    

    parser.add_option("--replicates", dest="replicates", type="int",
                      help="in benchmarking mode expect ## replicates." )

    parser.set_defaults(
        input_format="stockholm",
        input_filename_tree=None, 
        dump = False,
        clean_data = False,
        gap_chars = "-.",
        mask_chars = "nN",
        clean_mali = True,
        method = "train",
        insert_frequencies = False,
        fix_frequencies = False,
        write = [],
        output_pattern = "%s.eg",
        value_format = "%6.4f",
        xrate_min_increment = 0.000001,
        separator = "|",
        single_omega = False,
        shared_frequencies = False,
        shared_rates = False,
        block_size = None,
        num_blocks = 2,
        replicates = None,
        grammar_type = "multiple-blocks",
        probability_block = 0.9,
        )

    (options, args) = Experiment.Start( parser )

    if options.replicates != None:
        # read a sequence collection with possible duplicate names
        # used for benchmarking
        mali = Mali.SequenceCollection()
    else:
        mali = Mali.Mali()
    
    mali.readFromFile( sys.stdin, format = options.input_format )

    if options.method == "train":
        processMali = trainMali
        options.stdout.write("logL\tniterations\t%s\tannotation\n" %\
                             ("\t".join( ["rate%s" % x for x in range( options.num_blocks) ] ) ) )

    elif options.method == "annotate":
        raise "unimplemented"
        processMali = annotateMali

    if options.replicates != None:
        ids = mali.getIdentifiers()
        assert( len(ids) % options.replicates == 0 ) 
        s = len(ids) / options.replicates
        for x in range( 0, len(ids), s ):
            m = Mali.Mali()
            for id in ids[x:x+s]:
                m.addEntry( mali.getEntry( id ) )

            processMali( m, options )
    else:
        processMali( mali, options )

    Experiment.Stop()

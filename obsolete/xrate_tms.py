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
xrate_tms.py - 
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

   python xrate_tms.py --help

Type::

   python xrate_tms.py --help

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

analyze a protein alignment using a model for TM proteins.

""" % sys.argv[0]

import CGAT.Experiment as Experiment
import CGAT.Mali as Mali
import CGAT.Genomics as Genomics
import CGAT.RateEstimation as RateEstimation
import CGAT.TreeTools as TreeTools
import CGAT.IOTools as IOTools

from XGram.Generator.Prebuilt import Codons
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
def prepareGrammar( xgram, mali, tree, map_old2new, blocks, options ):
    """prepare grammar for custom grammars."""
    
    labels = map( lambda x: x[1], blocks )
    nblocks = len(blocks)
    
    annotate_terminals = {}
    for x in range(len(labels)):
        annotations = []
        key = []

        for c in range( 0,3 ):
            t = "B%i_COD%i" % (x, c)
            key.append(t)
            annotations.append( Annotation( row = "STATE",
                                            column = t,
                                            label = labels[x] ))
            
        annotate_terminals[ tuple(key) ] = annotations

    input_model = Codons.buildCodonML( codon_model = "f3x4-fourproducts",
                                       num_blocks = nblocks,
                                       grammar_type = "linear-blocks",
                                       annotate_terminals=annotate_terminals,
                                       shared_frequencies = options.shared_frequencies,
                                       shared_rates = False,
                                       )

    ## manually share rates between blocks
    if options.shared_rates == "kappa":
        for c in range( 0, nblocks):
            input_model.renameParameter( "B%i_Ri" % c, "Ri" )
            input_model.renameParameter( "B%i_Rv" % c, "Rv" )
    elif options.shared_rates == "kappa-ds":
        for c in range( 0, nblocks):
            input_model.renameParameter( "B%i_Ri" % c, "Ri" )
            input_model.renameParameter( "B%i_Rv" % c, "Rv" )
            input_model.renameParameter( "B%i_Rs" % c, "Rs" )
    elif options.shared_rates == "omega":
        for c in range( 0, nblocks):
            input_model.renameParameter( "B%i_Rs" % c, "Rs" )
            input_model.renameParameter( "B%i_Rn" % c, "Rn" )
    elif options.shared_rates == "omega-ds":
        for c in range( 0, nblocks):
            input_model.renameParameter( "B%i_Rv" % c, "Rv" )
            input_model.renameParameter( "B%i_Rs" % c, "Rs" )
            input_model.renameParameter( "B%i_Rn" % c, "Rn" )
    elif options.shared_rates == "ds":
        for c in range( 0, nblocks):
            input_model.renameParameter( "B%i_Rs" % c, "Rs" )
    elif options.shared_rates == "all":
        for c in range( 0, nblocks):
            input_model.renameParameter( "B%i_Rv" % c, "Rv" )
            input_model.renameParameter( "B%i_Rs" % c, "Rs" )
            input_model.renameParameter( "B%i_Rn" % c, "Rn" )
            input_model.renameParameter( "B%i_Ri" % c, "Ri" )

    writeModel( input_model, "input", options )
    
    ids = mali.getIdentifiers()

    fh, filename = tempfile.mkstemp()

    os.close(fh)
    outfile = open(filename, "w" )
    
    ## clip mali by supplied blocks
    mali.clipByAnnotation( "STATE", "".join(labels))

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
    
    ## prefix, code
    if options.shared_frequencies:
        frequency_codes = ( ("", ""), )
    else:
        frequency_codes = blocks
        
    if options.insert_frequencies:
        for prefix, code in frequency_codes:
            temp_mali = mali.getClone()
            temp_mali.clipByAnnotation( "STATE", code )
            RateEstimation.setFrequencies( input_model, temp_mali, prefix )
            
    if options.fix_frequencies:
        for prefix, code in frequency_codes:
            for char in ('a', 'c', 'g', 't'):
                for x in (0, 1, 2):
                    param = "%sp%s%i" % (prefix, char, x)
                    input_model.mGrammar.moveVariableToConst( param )

    writeModel( input_model, "input", options )
    
    t1 = time.time()

    result = xgram.train( input_model, filename )

    if options.dump:
        options.stdlog.write( "".join(result.mData) )
        options.stdlog.write( "".join(result.mLog) )
        mali.writeToFile( options.stdlog, 
                          format="stockholm",
                          write_ranges = False,
                          options = (tree_options,))

    t2 = time.time()
    
    trained_model = result.getModel()

    writeModel( trained_model, "trained", options )

    return result, mali, ids

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def processMali( mali, options ):

    ncols = mali.getNumColumns()

    if ncols == 0:
        raise "refusing to process empty alignment."

    ## add annotation of states
    if options.block_size != None:
        if options.block_size < 1:
            size = int( float( ncols ) / 3.0 * options.block_size) * 3
        else:
            size = int( options.block_size ) * 3
        
        size = min( size, ncols )
        mali.addAnnotation( "STATE", "N" * size + "C" * (ncols - size))
            
    ## remove gene ids
    for id in mali.getIdentifiers():
        if options.separator in id:
            species = id.split(options.separator)[0]
            mali.rename( id, species )

    map_new2old = mali.mapIdentifiers()
    map_old2new = IOTools.getInvertedDictionary( map_new2old, make_unique = True )
    
    ids = mali.getIdentifiers()
    xgram = XGram.XGram()

    if options.xrate_min_increment:
        xgram.setMinIncrement( options.xrate_min_increment )

    ninput, noutput, nskipped = 0, 0, 0

    # remove empty columns and masked columns
    if options.clean_mali:
        mali.mGapChars = mali.mGapChars + ("n", "N")
        mali.removeGaps( minimum_gaps = 1, frame=3 )

    if options.input_filename_tree:
        nexus = TreeTools.Newick2Nexus( open(options.input_filename_tree,"r") )
        tree = nexus.trees[0]
        tree.relabel( map_old2new )
    else:
        tree = None

    annotation = mali.getAnnotation( "STATE" )
    chars = set(list(annotation))
    for c in chars:
        assert c in ("N", "C"), "unknown annotation %s: only 'N' and 'C' are recognized"
    if len(chars) == 1:
        if options.loglevel >= 1:
            options.stdlog.write("# WARNING: only a single block" )
        blocks = ( ("B0_", chars[0]), )
    else:
        blocks = ( ("B0_", "N"), 
                   ("B1_", "C") )
    
    result, mali, ids = prepareGrammar( xgram, mali, tree, map_old2new, blocks, options )

    trained_model = result.getModel()

    pis, matrices = RateEstimation.getRateMatrix( trained_model )

    annotation = mali.getAnnotation( "STATE" )

    for block, code in blocks :

        terminals = ( "%sCOD0" % block,
                      "%sCOD1" % block,
                      "%sCOD2" % block )
        
        pi = pis[terminals]

        if options.shared_rates == "all":
            rate_prefix_rs = ""
            rate_prefix_rn = ""
            rate_prefix_ri = ""
            rate_prefix_rv = ""
        elif options.shared_rates == "kappa":
            rate_prefix_rs = block
            rate_prefix_rn = block
            rate_prefix_ri = ""
            rate_prefix_rv = ""
        elif options.shared_rates == "kappa-ds":
            rate_prefix_rs = ""
            rate_prefix_rn = block
            rate_prefix_ri = ""
            rate_prefix_rv = ""
        elif options.shared_rates == "omega":
            rate_prefix_rs = ""
            rate_prefix_rn = ""
            rate_prefix_ri = block
            rate_prefix_rv = block
        elif options.shared_rates == "omega-ds":
            rate_prefix_rs = ""
            rate_prefix_rn = ""
            rate_prefix_ri = block
            rate_prefix_rv = ""
        elif options.shared_rates == "ds":
            rate_prefix_rs = ""
            rate_prefix_rn = block
            rate_prefix_ri = block
            rate_prefix_rv = block
        else:
            rate_prefix_rs = block
            rate_prefix_rn = block
            rate_prefix_ri = block
            rate_prefix_rv = block
        
        if options.shared_frequencies:
            frequency_prefix = ""
        else:
            frequency_prefix = block

        rs = trained_model.mGrammar.getParameter( '%sRs' % rate_prefix_rs )
        rn = trained_model.mGrammar.getParameter( '%sRn' % rate_prefix_rn )
        ri = trained_model.mGrammar.getParameter( '%sRi' % rate_prefix_ri )
        rv = trained_model.mGrammar.getParameter( '%sRv' % rate_prefix_rv )    

        nchars = annotation.count( code )

        msg = "iter=%i Rs=%6.4f Rn=%6.4f Ri=%6.4f Rv=%6.4f" % ( result.getNumIterations(), rs, rn, ri, rv )
        
        try:
            Q, t = RateEstimation.getQMatrix( pi,
                                              Rsi=rs * ri,
                                              Rsv=rs * rv,
                                              Rni=rn * ri,
                                              Rnv=rn * rv )
            avg_omega = (rs + rn) / 2.0
            Q0, t0 = RateEstimation.getQMatrix( pi,
                                                Rsi = ri * avg_omega,
                                                Rsv = rv * avg_omega,
                                                Rni = ri * avg_omega,
                                                Rnv = rv * avg_omega )

            avg_kappa = (ri + rv) / 2.0
            Q1, t1 = RateEstimation.getQMatrix( pi,
                                                Rsi = rs * avg_kappa,
                                                Rsv = rs * avg_kappa,
                                                Rni = rn * avg_kappa,
                                                Rnv = rn * avg_kappa )

            rI, rV, rS, rN = RateEstimation.countSubstitutions( pi, Q )
            rI0, rV0, rS0, rN0 = RateEstimation.countSubstitutions( pi, Q0 )    
            rI1, rV1, rS1, rN1 = RateEstimation.countSubstitutions( pi, Q1 )    

            dS = rS / (3 * rS0) * t
            dN = rN / (3 * rN0) * t

            o_kappa = options.value_format % ( rI / rI0 * rV0 / rV )
            o_omega = options.value_format % (dN / dS)

            o_dn = options.value_format % dN
            o_ds = options.value_format % dS
            o_rn = options.value_format % rN
            o_rs = options.value_format % rS
            o_rn0 = options.value_format % rN0
            o_rs0 = options.value_format % rS0
            o_t = options.value_format % t
            o_t0 = options.value_format % t0

        except ZeroDivisionError:

            o_kappa = "na"
            o_omega = "na"
            o_dn = "na"
            o_ds = "na"
            o_rn = "na"
            o_rs = "na"
            o_rn0 = "na"
            o_rs0 = "na"
            o_t = "na"
            o_t0 = "na"
            Q = None
            msg = "insufficient data to estimate rate matrix."
        
        options.stdout.write( "\t".join( map(str, (
                        code, block,
                        o_dn, o_ds, o_omega,
                        "na", "na", "na", "na",
                        o_kappa, 
                        result.getLogLikelihood(),
                        "na",
                        nchars ))))

        if options.with_rho:
            options.stdout.write( "\t" + "\t".join( map(str, (o_rn, o_rs, o_t,
                                                              o_rn0, o_rs0, o_t0 ))))
            
        options.stdout.write( "\t%s\n" %  msg )


##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------

if __name__ == "__main__":
    
    parser = optparse.OptionParser( version = "%prog version: $Id: xrate_tms.py 2781 2009-09-10 11:33:14Z andreas $" )

    parser.add_option("--input-filename-tree", dest="input_filename_tree", type="string",
                       help="""filename with pyhlogenetic tree.""" )
     
    parser.add_option("--dump", dest="dump", action="store_true",
                      help="dump raw output.")
    
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

    parser.add_option("--xrate-shared-frequencies", dest="shared_frequencies", action="store_true",
                      help="frequencies are shared between each chain." )

    parser.add_option("--xrate-insert-frequencies", dest="insert_frequencies", action="store_true",
                      help="estimate codon frequencies from input." )

    parser.add_option("--xrate-fix-frequencies", dest="fix_frequencies", action="store_true",
                      help="set initial frequencies to const." )

    parser.add_option("--xrate-uniform-frequencies", dest="insert_frequencies", action="store_false",
                      help="use uniform codon frequencies." )

    parser.add_option("--xrate-estimate-frequencies", dest="fix_frequencies", action="store_false",
                      help="estimate nucleotide frequencies." )

    parser.add_option("--xrate-shared-rates", dest="shared_rates", type="choice",
                      choices=("none", "all", "kappa", "omega", "ds", "kappa-ds", "omega-ds" ),
                      help="""rates are shared between each chain. 
none:      no rates are shared.
all:       all rates are shared. 
kappa:     only kappa rates ri and rv are shared.
kappa-ds:  fix kappa and ds.
""" )

    parser.add_option("--block-size", dest="block_size", type="float",
                      help="set size of first block (in codons)." )

    parser.add_option( "--xrate-min-increment", dest="xrate_min_increment", type = "float",
                       help="minimum increment to stop iteration in xrate." )    

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
        method = "paml",
        insert_frequencies = False,
        fix_frequencies = False,
        write = [],
        output_pattern = "%s.eg",
        value_format = "%6.4f",
        xrate_min_increment = 0.000001,
        with_rho = True,
        separator = "|",
        single_omega = False,
        shared_frequencies = False,
        shared_rates = False,
        block_size = None,
        replicates = None,
        )

    (options, args) = Experiment.Start( parser )

    if options.replicates != None:
        # read a sequence collection with possible duplicate names
        # used for benchmarking
        mali = Mali.SequenceCollection()
    else:
        mali = Mali.Mali()
    
    mali.readFromFile( sys.stdin, format = options.input_format )

    options.stdout.write("seq1\tseq2\tdN\tdS\tdNdS\tN\tS\tdN_err\tdS_err\tkappa\tlnL\ttau\tlen" )

    if options.with_rho:
        options.stdout.write("\trN\trS\tt\trN0\trS0\tt0" )

    options.stdout.write( "\terror_str\n")

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

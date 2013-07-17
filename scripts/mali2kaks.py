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
mali2kaks.py - rate analysis with PAML
======================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script will do pairwise rate analysis on a set of multiply aligned
sequences in :term:`fasta` format.

PAML bails if sequences in a multiple alignment are not overlapping. Thus this
code checks for each pair the overlap. If a pair is not overlapping, it goes
into pairwise mode.

Usage
-----

Example::

   python mali2kaks.py --help

Type::

   python mali2kaks.py --help

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
import tempfile
import subprocess
import optparse
import time
import math

from types import *

import CGAT.Experiment as E
import CGAT.WrapperCodeML as WrapperCodeML
import CGAT.Mali as Mali
import CGAT.Genomics as Genomics
import CGAT.RateEstimation as RateEstimation
import CGAT.IOTools as IOTools
import CGAT.TreeTools as TreeTools

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
def printPairs( pairs, mali, map_new2old, options ):
    """print pairs form codeml."""
    noutput = 0
    for pair in pairs:
        options.stdout.write( "\t".join( map(str, (mali.getEntry(pair.mName2).mId,
                                                   mali.getEntry(pair.mName1).mId,
                                                   pair.mKa, pair.mKs, pair.mKaks,
                                                   pair.mN, pair.mS, "na", "na",
                                                   pair.mKappa, pair.mLogLikelihood,
                                                   pair.mTau ))))
        
        if options.with_rho:
            options.stdout.write( "\t" + "\t".join( map(str, (pair.mRn, pair.mRs, pair.mBranchLength,
                                                              pair.mRn0, pair.mRs0, "na") ) ) )

        if options.with_counts:
            info = Genomics.CalculatePairIndices( mali[pair.mName1], mali[pair.mName2] )
            options.stdout.write( "\t%s" % (str(info)) )

        options.stdout.write( "\t" + pair.mError + "\n" )
        options.stdout.flush()
        noutput += 1
        
    return noutput
##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
def runCodeML( mali, tree, has_non_overlaps, pairs, map_new2old, options ):
    """setup codeml wrapper.

    Sets options and returns a wrapper.
    """

    ids = mali.getIdentifiers()

    ## setup codeml
    codeml_options = {}

    if options.seqtype == "codon":
        codeml_options["seqtype"] = "1"
    elif options.seqtype == "aa":        
        codeml_options["seqtype"] = "2"        
    elif options.seqtype == "trans":
        codeml_options["seqtype"] = "3"

    if options.clean_data:
        codeml_options[ "cleandata" ] = options.clean_data

    if options.omega != None:
        codeml_options[ "omega" ] = str(options.omega)

    if options.kappa != None:
        codeml_options[ "kappa" ] = str(options.kappa)

    if options.fix_kappa:
        codeml_options[ "fix_kappa" ] = "1"

    if options.fix_omega:
        codeml_options[ "fix_omega" ] = "1"

    if options.codon_frequencies != None:
        c = options.codon_frequencies.upper()  
        if c == "UNIFORM":
            a = "0"
        elif c == "F1X4":
            a = "1"
        elif c == "F3X4":
            a = "2"
        elif c == "F61":
            a = "3"
        else:
            a = options.codon_frequencies
        codeml_options[ "CodonFreq" ] = a

    if options.paml_method != None:
        codeml_options[ "paml_method" ] = str(options.method)

    if options.optimization_threshold != None:
        codeml_options[ "Small_Diff" ] = str(options.optimization_threshold)

    ninput, noutput, nskipped = 0, 0, 0
    tstart = time.time()
    
    if pairs and (options.pairwise or has_non_overlaps):
        wrapper = WrapperCodeML.CodeMLPairwise()
        
        ## do pairwise run
        result = WrapperCodeML.CodeMLResultPairs()

        ntotal = (len(ids) * (len(ids) -1 ) ) / 2

        for x,y in pairs:
            m1 = mali.getSequence( ids[x] )                        
            ninput += 1

            temp_mali = Mali.Mali()
            m2 = mali.getSequence( ids[y] )

            temp_mali.addSequence( ids[x], m1.mFrom, m1.mTo, m1.mString )
            temp_mali.addSequence( ids[y], m2.mFrom, m2.mTo, m2.mString )

            ## remove empty columns and masked columns
            if options.clean_mali:
                temp_mali.mGapChars = temp_mali.mGapChars + ("n", "N")
                temp_mali.removeGaps( minimum_gaps = 1, frame=3 )

            if temp_mali.getWidth() < options.min_overlap:
                if options.loglevel >= 1:
                    options.stdlog.write("# pair %s-%s: not computed because only %i residues overlap\n" % (mali.getEntry(ids[x]).mId,
                                                                                                            mali.getEntry(ids[y]).mId,
                                                                                                            temp_mali.getWidth()) )
                nskipped += 1
                continue

            sub_result = wrapper.Run( temp_mali, options = codeml_options, dump=options.dump )
            result.mPairs += sub_result.mPairs

            if options.loglevel >= 1 and ninput % options.report_step == 0:
                options.stdlog.write( "# pairwise computation: %i/%i -> %i%% in %i seconds.\n" % (ninput, ntotal, 100.0 * ninput / ntotal, time.time() - tstart ) )
                options.stdlog.flush()

            noutput += printPairs( sub_result.mPairs,
                                   mali,
                                   map_new2old,
                                   options )

            options.stdout.flush()

        if options.loglevel >= 1:
            options.stdlog.write("# pairwise computation: ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped))
            options.stdlog.flush()
                
    else:
        wrapper = WrapperCodeML.CodeML()

        result = wrapper.Run( mali, tree = tree, options = codeml_options, dump=options.dump )

        result_pairs = WrapperCodeML.CodeMLResultPairs()
        result_pairs.fromResult( result )
        noutput += printPairs( result_pairs.mPairs,
                               mali,
                               map_new2old,
                               options )

        l = mali.getLength()
        if options.loglevel >= 1:        
            options.stdlog.write("# input=%i, npairs=%i, noutput=%i\n" % (l, l * (l - 1 ) / 2, len(result_pairs.mPairs) ) )


##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def checkMatrix( pi, matrix ):
    """check if matrix does fulfill the requirements:
    diagonal sums to -1
    row sum is equal to diagonal entry.
    piQij = pjQji
    """
    codon_table = Bio.Data.CodonTable.standard_dna_table
    codons = codon_table.forward_table.keys()

    ## Matrix needs to be normalized so that branch lenghts in a tree
    ## are effectively measured as expected numbers of nucleotide substituions
    ## per codon
    t = 0.0
    for codon_i in codons:
        t += pi[codon_i] * matrix[codon_i][codon_i]
    # print "check Q diagonal (should be -1): %f" % t        
    assert( "%5.2f" % t == "%5.2f" % -1.0 )
    
    for codon_i in codons:
        t = 0
        for codon_j in codons:
            if codon_i == codon_j: continue
            t += matrix[codon_i][codon_j]
        assert( "%5.2f" % t == "%5.2f" % -matrix[codon_i][codon_i] )
        assert( "%5.2f" % (pi[codon_i] * matrix[codon_i][codon_j]) == "%5.2f" % (pi[codon_j] * matrix[codon_j][codon_i] ) )
        # print "check row %s (should be %f): %f" % (codon_i, -matrix[codon_i][codon_i], t)

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def getSequencesFromStk( filename):
    """quick parse through stk file.
    
    Only gets sequences."""
    sequences = {}
    for line in open(filename, "r"):
        if line[0] == '#': continue
        if line[:2] == '//': break
        id, s = re.split('\s+', line[:-1])
        if id not in sequences:
            sequences[id] = []
        sequences[id].append(s)
    
    for id, s in sequences.items():
        sequences[id] = "".join(s)
        
    return sequences

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def setFrequencies( model, data):
    """set frequencies in a model according to those observed in data (F3X4).
    """
    sequences = getSequencesFromStk( data )

    frequencies = Codons.getFrequenciesPerCodonPosition( sequences.values() )

    ## build a dummy grammar to insert frequencies
    dummy_grammar = XGram.Model.Grammar()
    for x in range(0,3):
        params = []
        for a in ('a', 'c', 'g', 't'):
            params.append( ("p%s%i" % (a.lower(), x), frequencies[x][a.upper()]) )
        dummy_grammar.addVariable( params )
        
    model.mGrammar.copyParameters( dummy_grammar, 
                                   ignore_missing = True)

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def setFixedRates( model, rates):
    """set rates in a model to specified values and fix them."""

    for x in range(len(rates)):
        k = "K%i" % x
        model.mGrammar.setParameter( k, rates[x] )
        model.mGrammar.moveVariableToConst( k )

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def writeModel( grammar, section, options):
    """write a model to output file."""
    if section in options.write or "all" in options.write:
        outfile = open( options.output_pattern % section, "w" )
        outfile.write( "%s\n" % grammar.getGrammar())
        outfile.close()

def prepareGrammar( xgram, mali, options ):
    """prepare grammar for custom grammars."""

    ids = mali.getIdentifiers()

    fh, filename = tempfile.mkstemp()
    os.close(fh)
    outfile = open(filename, "w" )
    mali.writeToFile( outfile, format="stockholm",
                      write_ranges = False,
                      options = ("#=GF NH (%s:1.0)%s;" % tuple(ids),) )
    outfile.close()
    
    if options.xrate_model == "sn":
        infile = open(XGram.PATH_DATA + "/sn.eg", "r")
        input_model = XGram.Parser.parseGrammar( infile.readlines() )                

    elif options.xrate_model == "akaksgc":
        infile = open(XGram.PATH_DATA + "/akaksgc.eg", "r")
        input_model = XGram.Parser.parseGrammar( infile.readlines() )

    elif options.xrate_model in ("f3x4-two", "f3x4-four", "f3x4-fourproducts" ):
        input_model = Codons.buildCodonML(codon_model = options.xrate_model,
                                          fix_kappa = options.fix_kappa,
                                          fix_omega = options.fix_omega )

    if options.xrate_model in ("ef3x4-four",):

        sequences = getSequencesFromStk( filename )
        frequencies = Codons.getFrequenciesPerCodonPosition( sequences.values() )

        codon_frequencies = {}
        if options.xrate_insert_frequencies:
            for c1 in ('A', 'C', 'G', 'T'):
                for c2 in ('A', 'C', 'G', 'T'):
                    for c3 in ('A', 'C', 'G', 'T'):
                        codon = "".join( (c1,c2,c3) )
                        if not Genomics.IsStopCodon(codon):
                            codon_frequencies[codon] = frequencies[0][c1] * frequencies[1][c2] * frequencies[2][c3]

            total = sum( codon_frequencies.values())
            for k,v in codon_frequencies.items():
                codon_frequencies[k] /= total
        else:
            for c1 in ('A', 'C', 'G', 'T'):
                for c2 in ('A', 'C', 'G', 'T'):
                    for c3 in ('A', 'C', 'G', 'T'):
                        codon = "".join( (c1,c2,c3) )
                        codon_frequencies[codon] = 1/61.0

        input_model = Codons.buildCodonML(codon_model = "codons-four",
                                          codon_frequencies = codon_frequencies,
                                          fix_kappa = options.fix_kappa,
                                          fix_omega = options.fix_omega )

    else:

        if options.xrate_insert_frequencies:
            setFrequencies( input_model, filename )

        if options.xrate_fix_frequencies:
            for char in ('a', 'c', 'g', 't'):
                for x in (0, 1, 2):
                    param = "p%s%i" % (char, x)
                    input_model.mGrammar.moveVariableToConst( param )

    if options.dump:
        options.stdlog.write( "## input model:\n%s\n" % input_model.getGrammar() )
                    
    writeModel( input_model, "input", options )
    
    t1 = time.time()

    result = xgram.train( input_model, filename )
    
    t2 = time.time()
    
    trained_model = result.getModel()

    if options.dump:
        options.stdlog.write( "## trained model:\n%s\n" % trained_model.getGrammar() )

    writeModel( trained_model, "trained", options )

    return result, mali, ids

def evaluateGrammar( trained_model ):
    
    ## retrieve the terminal frequencies for all codons
    xpi = trained_model.evaluateTerminalFrequencies()[('COD0', 'COD1', 'COD2')]

    pi = {}
    for codon, f in xpi.items():
        pi["".join(codon).upper()] = f

    ## retrieve the rate matrix for all codons
    xmatrix = trained_model.evaluateRateMatrix()[('COD0', 'COD1', 'COD2')]
    matrix = {}
    for codon1, v in matrix.items():
        x = {}
        for codon2, f in matrix[codon1].items():
            x["".join(codon2).upper()] = v
        matrix["".join(codon1).upper()] = x

    return pi, matrix

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
def runXrateSN( xgram, mali, options ):
    """run xrate using Ians sn.eg grammar."""

    result, mali, ids = prepareGrammar( xgram, mali, options )
    trained_model = result.getModel()
    
    pi, matrix = evaluateGrammar( trained_model )
    
    def getQMatrix( pi, k, s, n):
        """build a q matrix.

        Diagonal elements are set to the negative of the row sums.
        The matrix is normalized such that trace of the matrix is -1.
        """

        codons = Bio.Data.CodonTable.standard_dna_table.forward_table.keys()

        Q = initializeQMatrix( codons )

        trace = 0.0
        for codon_i in codons:
            row_sum = 0.0
            for codon_j in codons:
                if codon_i == codon_j: continue

                is_single, is_synonymous, is_transition = RateEstimation.evaluateCodonPair( codon_i, codon_j )
                
                if not is_single: continue

                if is_synonymous:
                    if is_transition:
                        v = s
                    else:
                        v = s * k
                else:
                    if is_transition:
                        v = n
                    else:
                        v = n * k

                v *= pi[codon_j]
                Q[codon_i][codon_j] = v
                row_sum += v
                
            Q[codon_i][codon_i] = -row_sum
            trace += pi[codon_i] * row_sum

        for codon_i in codons:
            for codon_j in codons:
                Q[codon_i][codon_j] /= trace
                
        return Q, trace

    s = trained_model.mGrammar.getParameter( 's' )
    n = trained_model.mGrammar.getParameter( 'n' )    
    k = trained_model.mGrammar.getParameter( 'k' )    
    not_k = trained_model.mGrammar.getParameter( 'not_k' )    

    Q, t = getQMatrix( pi, k, s, n )
    Q0, t0 = getQMatrix( pi, k, 1, 1 )

    ri, rv, rS, rN = countSubstitutions( pi, Q )
    ri0, rv0, rS0, rN0 = countSubstitutions( pi, Q0 )    

    kappa = ri / rv
    dS = rS / (3 * rS0) * t
    dN = rN / (3 * rN0) * t
    
    if s == None or n == None:
        o_dn, o_ds, o_omega = "na", "na", "na"
        o_rn, o_rn0, o_rs, o_rs0 = "na", "na", "na", "na"
        o_t, o_t0 = "na", "na"
        o_kappa = "na",
        msg = "estimated rate parameters are zero"
    else:
        o_omega = options.value_format % (n / s)
        o_dn = options.value_format % dN
        o_ds = options.value_format % dS
        o_rn = options.value_format % rN
        o_rs = options.value_format % rS
        o_rn0 = options.value_format % rN0
        o_rs0 = options.value_format % rS0
        o_t = options.value_format % t
        o_t0 = options.value_format % t0
        o_kappa = options.value_format % kappa
        msg = "iter=%i s=%6.4f n=%6.4f k=%6.4f ~k=%6.4f" % (result.getNumIterations(), s, n, k, not_k)

    options.stdout.write( "\t".join( map(str, (mali.getEntry(ids[0]).mId,
                                               mali.getEntry(ids[1]).mId,
                                               o_dn, o_ds, o_omega,
                                               "na", "na", "na", "na",
                                               o_kappa, result.getLogLikelihood(),
                                               "na",
                                               o_rn, o_rs, o_t,
                                               o_rn0, o_rs0, o_t0 ))))

    if options.with_counts:
        info = Genomics.CalculatePairIndices( mali[ids[0]], mali[ids[1]] )
        options.stdout.write( "\t%s" % (str(info)) )

    options.stdout.write( "\t%s\n" % msg )

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
def runXrateAKaKsGc( xgram, mali, options ):
    """run xrate using Ians sn.eg grammar."""

    result, mali, ids = prepareGrammar( xgram, mali, options )
    trained_model = result.getModel()

    rsi = trained_model.mGrammar.getParameter( 'Rsi' )
    rsv = trained_model.mGrammar.getParameter( 'Rsv' )
    rni = trained_model.mGrammar.getParameter( 'Rni' )
    rnv = trained_model.mGrammar.getParameter( 'Rnv' )    
    msg = "iter=%i Rsi=%6.4f Rsv=%6.4f Rni=%6.4f Rnv=%6.4f" % (result.getNumIterations(), rsi, rsv, rni, rnv )
        
    outputXRateResult( mali, result, rsi, rsv, rni, rnv, msg )
        
##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
def runXrateF3X4( xgram, mali, options ):
    """run F3X4 type of models using my grammars.
    """

    result, mali, ids = prepareGrammar( xgram, mali, options )

    trained_model = result.getModel()

    if options.xrate_model == "f3x4-fourproducts":
            i_rs = trained_model.mGrammar.getParameter( 'Rs' )
            i_rn = trained_model.mGrammar.getParameter( 'Rn' )
            i_ri = trained_model.mGrammar.getParameter( 'Ri' )
            i_rv = trained_model.mGrammar.getParameter( 'Rv' )    
            msg = "iter=%i Rs=%6.4f Rn=%6.4f Ri=%6.4f Rv=%6.4f" % (result.getNumIterations(), i_rs, i_rn, i_ri, i_rv )
            rsi = i_rs * i_ri
            rsv = i_rs * i_rv
            rni = i_rn * i_ri
            rnv = i_rn * i_rv
    elif options.xrate_model == "f3x4-two":
            i_kappa = trained_model.mGrammar.getParameter( 'kappa' )
            i_omega = trained_model.mGrammar.getParameter( 'omega' )
            msg = "iter=%i kappa=%6.4f omega=%6.4f" % (result.getNumIterations(), i_kappa, i_omega)
            rsi = i_kappa 
            rsv = 1.0
            rni = i_kappa * i_omega
            rnv = i_omega
    else:
        if options.fix_omega:
            i_kappa = trained_model.mGrammar.getParameter( 'kappa' )
            i_not_kappa = trained_model.mGrammar.getParameter( 'not_kappa' )
            msg = "iter=%i kappa=%6.4f not_kappa=%6.4f" % (result.getNumIterations(), i_kappa, i_not_kappa)
            rsi = i_kappa
            rsv = i_not_kappa
            rni = i_kappa
            rnv = i_not_kappa
        elif options.fix_kappa:
            i_omega = trained_model.mGrammar.getParameter( 'omega' )
            i_not_omega = trained_model.mGrammar.getParameter( 'not_omega' )
            rsi = i_omega
            rsv = i_omega
            rni = i_not_omega        
            rnv = i_not_omega
            msg = "iter=%i omega=%6.4f not_omega=%6.4f" % (result.getNumIterations(), i_omega, i_not_omega)        
        else:
            rsi = trained_model.mGrammar.getParameter( 'Rsi' )
            rsv = trained_model.mGrammar.getParameter( 'Rsv' )
            rni = trained_model.mGrammar.getParameter( 'Rni' )
            rnv = trained_model.mGrammar.getParameter( 'Rnv' )    
            msg = "iter=%i Rsi=%6.4f Rsv=%6.4f Rni=%6.4f Rnv=%6.4f" % (result.getNumIterations(), rsi, rsv, rni, rnv )
        
    outputXRateResult( mali, result, rsi, rsv, rni, rnv, msg )

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
def outputXRateResult( mali, result, rsi, rsv, rni, rnv, msg ):
    """output the results of running the Xrate four parameter grammar.
    """
    ids = mali.getIdentifiers()
    
    pi, matrix = RateEstimation.getRateMatrix( result.getModel(), 
                                               terminals=('COD0', 'COD1', 'COD2'))

    if rsi == None:
        o_dn, o_ds, o_omega = "na", "na", "na"
        o_rn, o_rn0, o_rs, o_rs0 = "na", "na", "na", "na"
        o_t, o_t0 = "na", "na"
        o_N, o_S = "na", "na"
        o_kappa = "na",
        msg = "estimated rate parameters are zero"
    else:
        Q, t = RateEstimation.getQMatrix( pi,
                                          Rsi=rsi,
                                          Rsv=rsv,
                                          Rni=rni,
                                          Rnv=rnv )

        ## get rate matrix as if omega was set to 1
        Q0, t0 = RateEstimation.getQMatrix( pi,
                                            Rsi = (rsi + rni) / 2.0,
                                            Rsv = (rsv + rnv) / 2.0,
                                            Rni = (rsi + rni) / 2.0,
                                            Rnv = (rsv + rnv) / 2.0 )

        ## get rate matrix as if kappa was set to 1
        Q1, t1 = RateEstimation.getQMatrix( pi,
                                            Rsi = (rsi + rsv) / 2.0,
                                            Rsv = (rsi + rsv) / 2.0,
                                            Rni = (rni + rnv) / 2.0,
                                            Rnv = (rni + rnv) / 2.0 )

        rI, rV, rS, rN = RateEstimation.countSubstitutions( pi, Q )
        rI0, rV0, rS0, rN0 = RateEstimation.countSubstitutions( pi, Q0 )    
        rI1, rV1, rS1, rN1 = RateEstimation.countSubstitutions( pi, Q1 )    

        # 64.0/61.0 results from the fact that xrate does not normalize
        # the terminals 
        dS = rS / (3 * rS0) * t 
        dN = rN / (3 * rN0) * t 
        
        o_omega = options.value_format % (dN / dS)
        o_dn = options.value_format % dN
        o_ds = options.value_format % dS
        o_rn = options.value_format % rN
        o_rs = options.value_format % rS
        o_rn0 = options.value_format % rN0
        o_rs0 = options.value_format % rS0
        o_t = options.value_format % t
        o_t0 = options.value_format % t0
        o_S = options.value_format % (mali.getNumColumns() * rS0)
        o_N = options.value_format % (mali.getNumColumns() * rN0)

        ## kappa is given normalized by sites like omega
        o_kappa = options.value_format % (rI / rI1 * rV1 / rV)

        ## kappa1 is given by the ratio of the rates NOT normalized by the sites.
        msg += " rI/rV=%f rI0/rV0=%f kappa1=%s" % (rI/rV, 
                                                   rI0/rV0,
                                                   options.value_format % ( (rsi + rni) / (rsv + rnv) ) )
    
    options.stdout.write( "\t".join( map(str, (mali.getEntry(ids[0]).mId,
                                               mali.getEntry(ids[1]).mId,
                                               o_dn, o_ds, o_omega,
                                               o_N, o_S, "na", "na",
                                               o_kappa, result.getLogLikelihood(),
                                               "na" ))))

    if options.with_rho:
        options.stdout.write( "\t" + "\t".join( map(str, (o_rn, o_rs, o_t,
                                                          o_rn0, o_rs0, o_t0 ))))

    if options.with_counts:
        info = Genomics.CalculatePairIndices( mali[ids[0]], mali[ids[1]] )
        options.stdout.write( "\t%s" % (str(info)) )

    options.stdout.write( "\t%s\n" % msg )
    options.stdout.flush()

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
def runXrate( mali, has_non_overlaps, pairs, map_old2new, options ):
    """run xrate on a multiple alignment."""

    ids = mali.getIdentifiers()

    xgram = XGram.XGram()
    if options.xrate_min_increment:
        xgram.setMinIncrement( options.xrate_min_increment )

    ninput, noutput, nskipped = 0, 0, 0

    ## do pairwise run
    for x,y in pairs:
        m1 = mali.getSequence( ids[x] )            
        ninput += 1
        temp_mali = Mali.Mali()
        m2 = mali.getSequence( ids[y] )

        temp_mali.addSequence( m1.mId, m1.mFrom, m1.mTo, m1.mString )
        temp_mali.addSequence( m2.mId, m2.mFrom, m2.mTo, m2.mString )

        ## remove empty columns and masked columns
        if options.clean_mali:
            temp_mali.mGapChars = temp_mali.mGapChars + ("n", "N")
            temp_mali.removeGaps( minimum_gaps = 1, frame=3 )

        if temp_mali.getWidth() < options.min_overlap:
            if options.loglevel >= 1:
                options.stdlog.write("# pair %s-%s: not computed because only %i residues overlap\n" % (mali.getEntry(ids[x]).mId,
                                                                                                        mali.getEntry(ids[y]).mId,
                                                                                                        temp_mali.getWidth()) )

            nskipped += 1
            continue

        if options.xrate_model in ( "sn", ):
            runXrateSN( xgram, temp_mali, options )
        elif options.xrate_model in ( "akaksgc" ):
            runXrateAKaKsGc( xgram, temp_mali, options )            
        else:
            runXrateF3X4( xgram, temp_mali, options )

        if options.loglevel >= 1 and ninput % options.report_step == 0:
            options.stdlog.write( "# pairwise computation: %i/%i -> %i%% in %i seconds.\n" % (ninput, ntotal, 100.0 * ninput / ntotal, time.time() - tstart ) )
            options.stdlog.flush()
            
        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# pairwise computation: ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped))
        options.stdlog.flush()
    
##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
def processMali( mali, options ):

    map_new2old = mali.mapIdentifiers()
    ids = mali.getIdentifiers()
        
    invalid_chars = options.gap_chars + options.mask_chars
    
    has_non_overlaps = False

    pairs = []

    if options.iteration == "all-vs-all":
        for x in range(len(ids)):
            for y in range(0,x):
                pairs.append( (x,y) )
    elif options.iteration == "first-vs-all":
        for y in range(1, len(ids)):
            pairs.append( (0,y) )
    elif options.iteration == "pairwise":
        if len(ids) % 2 != 0:
            raise "uneven number of sequences (%i) not compatible with --iteration=pairwise" % len(ids)
        for x in range(0, len(ids), 2):
            pairs.append( (x, x+1) )
    elif options.iteration == "tree":
        pairs = []
    else:
        raise "unknown iteration mode: %s" % (options.iteration)

    if options.remove_stops:
        for id, entry in mali.items():
            s = entry.mString.upper()
            fragments = []
            for x in range(0,len(s), 3):
                codon = s[x:x+3]
                if Genomics.IsStopCodon( codon ):
                    codon = "NNN"
                    
                fragments.append(codon)
                    
            entry.mString = "".join(fragments)

    for x, y in pairs:
        noverlap = 0
        for a, b in zip( mali[ids[x]],mali[ids[y]] ):
            if a not in invalid_chars and b not in invalid_chars:
                noverlap += 1
                if noverlap >= options.min_overlap:
                    break
        else:
            has_non_overlaps = True
            break
        
    if options.tree:
        tree = TreeTools.Newick2Nexus( options.tree ).trees[0]
        map_old2new = IOTools.getInvertedDictionary( map_new2old, make_unique = True )        
        tree.relabel( map_old2new )
    else:
        tree = None

    if options.method == "paml":
        runCodeML( mali, tree, has_non_overlaps, pairs, map_new2old, options )

    elif options.method == "xrate":
        runXrate( mali, has_non_overlaps, pairs, map_new2old, options )


##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------

if __name__ == "__main__":
    
    parser = E.OptionParser( version = "%prog version: $Id: mali2kaks.py 2781 2009-09-10 11:33:14Z andreas $" )

    parser.add_option( "--set-omega", dest="omega", type="float",
                      help="initial omega value.")

    parser.add_option( "--set-kappa", dest="kappa", type="float",
                       help="initial kappa value.")

    parser.add_option( "--fix-kappa", dest="fix_kappa", action = "store_true",
                       help="do not estimate kappa.")

    parser.add_option( "--fix-omega", dest="fix_omega", action = "store_true",
                       help="do not estimate omega.")

    parser.add_option( "--set-codon-frequencies", dest="codon_frequencies", type="choice",
                       choices=("uniform", "fequal", "f3x4", "f1x4", "f61" ),
                       help="set codon frequencies.")
    
    parser.add_option( "--set-method", dest="paml_method", type="int",
                      help="set paml optimization method [0|1].")
    
    parser.add_option( "--set-sequence-type", dest="seqtype", type="choice",
                       choices = ("codon", "aa", "trans"),
                       help="sequence type.")

    parser.add_option( "--set-clean-data", dest="clean_data", type="choice",
                       choices = ("0", "1"),
                       help="PAML should cleanup data:  0=only gaps within pair are removed, 1=columns in the mali with gaps are removed.")

    parser.add_option("--dump", dest="dump", action="store_true",
                      help="dump raw output [%default].")
    
    parser.add_option( "--set-optimization-threshold", dest="optimization_threshold", type="string",
                      help="set paml optimization threshold [%default].")
    
    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal", "stockholm", "phylip" ),
                      help="input format of multiple alignment [%default]."  )

    parser.add_option( "--pairwise", dest="pairwise", action="store_true",
                      help="force pairwise comparison [%default].")

    parser.add_option( "--iteration", dest="iteration", type="choice",
                       choices=("all-vs-all", "first-vs-all", "pairwise", "tree" ),
                       help="iteration mode [%default]." )

    parser.add_option( "--no-clean", dest="clean_mali", action="store_false",
                      help="do not clean multiple alignment before submitting to codeml. It might take too long for very large sequences.")

    parser.add_option( "--method", dest="method", type="choice",
                       choices=("paml", "xrate"),
                       help = "choose method for rate computation [%default]" )

    parser.add_option( "--xrate-model", dest="xrate_model", type="choice", 
                      choices=("f3x4-two", "f3x4-four", "sn", "akaksgc", "ef3x4-four", "f3x4-fourproducts" ),
                      help="models to use [%default]."  )

    parser.add_option("-w", "--write", dest="write", type="choice", action="append",
                      choices=("input_fixed", "trained_fixed", "input_variable", "trained_variable", "all" ),
                      help="output sections to write [%default]."  )

    parser.add_option("-o", "--output-pattern", dest="output_pattern", type="string",
                      help="output pattern for output files [%default]." )

    parser.add_option("--xrate-insert-frequencies", dest="xrate_insert_frequencies", action="store_true",
                      help="estimate codon frequencies from input [%default]." )

    parser.add_option("--xrate-uniform-frequencies", dest="xrate_insert_frequencies", action="store_false",
                      help="use uniform codon frequencies [%default]." )

    parser.add_option("--xrate-fix-frequencies", dest="xrate_fix_frequencies", action="store_true",
                      help="set initial frequencies to const [%default]." )

    parser.add_option("--xrate-estimate-frequencies", dest="xrate_fix_frequencies", action="store_false",
                      help="estimate nucleotide frequencies [%default]." )

    parser.add_option("--xrate-fix-rates", dest="fix_rates", type="string",
                      help="""fix rates to specified values. Note that the number of rates has to match the ones
in the model. Provide values in a comma-separated list [%default].""")

    parser.add_option( "--xrate-min-increment", dest="xrate_min_increment", type = float,
                       help="minimum increment to stop iteration in xrate [%default]." )    

    parser.add_option( "--min-overlap", dest="min_overlap", type="int",
                       help="minimum overlap between a sequence pair in residues [%default].")
    
    parser.add_option( "--with-rho", dest="with_rho", action="store_true",
                      help="output rho values (substitution rates per codon). This requires a patched version of PAML [%default]."  )

    parser.add_option( "--with-counts", dest="with_counts", action="store_true",
                      help="output counts of aligned positions, transitions and transversions [%default]."  )

    parser.add_option( "--remove-stops", dest="remove_stops", action="store_true",
                      help="remove stop codons [%default]."  )

    parser.add_option("--replicates", dest="replicates", type="int",
                      help="in benchmarking mode expect ## replicates [%default]." )

    parser.add_option("--tree", dest="tree", type="string",
                      help="use tree for estimation [%default]." )

    parser.set_defaults(
        input_format="fasta",
        omega = None,
        codon_frequencies = None,
        paml_method = None,
        optimization_threshold = None,
        seqtype = "codon",
        dump = False,
        clean_data = False,
        min_overlap = 60,
        gap_chars = "-.",
        mask_chars = "nN",
        pairwise = False,
        kappa = None,
        fix_kappa = False,
        fix_omega = False,
        clean_mali = True,
        method = "paml",
        report_step = 1000,
        loglevel = 1,
        xrate_insert_frequencies = False,
        xrate_fix_frequencies = False,
        write = [],
        output_pattern = "%s.eg",
        value_format = "%6.4f",
        fix_rates = None,
        xrate_from_parameters = False,
        xrate_model = "f3x4-four",
        with_rho = False,
        with_counts = False,
        iteration = "all-vs-all",
        remove_stops = False,
        xrate_min_increment = 0.000001,
        replicates = None,
        tree = None,
        )

    (options, args) = E.Start( parser )

    if options.method == "xrate":
        # imports for xrate computation
        from XGram.Generator.Prebuilt import Codons
        from XGram.Model import Annotation
        import XGram.Run
        import Bio.Data.CodonTable

        # paml like estimation using xrate
        if options.codon_frequencies == "uniform":
            options.xrate_fix_frequencies = True
            options.xrate_insert_frequencies = False
        elif options.codon_frequencies == "f3x4":
            options.xrate_fix_frequencies = True
            options.xrate_insert_frequencies = True
    elif options.method == "paml":
        if not options.codon_frequencies:
            options.codon_frequencies = "F3X4"
    
    if options.fix_rates: options.fix_rates = map( float, options.fix_rates.split(",") )

    if options.pairwise or options.replicates:
        ## read sequences, but not as a multiple alignment. This permits multiple names.
        mali = Mali.SequenceCollection()
    else:
        mali = Mali.Mali()

    mali.readFromFile( sys.stdin, format = options.input_format )

    E.info("read multiple alignment")
    
    if mali.getLength() == 0:
        raise "refusing to process empty alignment."

    ################################################################
    ################################################################
    ################################################################
    ## setup methods
    ################################################################

    options.stdout.write("seq1\tseq2\tdN\tdS\tdNdS\tN\tS\tdN_err\tdS_err\tkappa\tlnL\ttau" )
        
    if options.with_rho:
        options.stdout.write("\trN\trS\tt\trN0\trS0\tt0" )

    if options.with_counts:
        options.stdout.write("\t%s" % Genomics.SequencePairInfo().getHeader() )

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

    E.Stop()


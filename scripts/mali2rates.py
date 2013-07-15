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
mali2rates.py - calculate rates for sequences in a multiple alignment
=====================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

calculate evolutionary rates for sequences in a multiple alignment.

Default output is as a list of pairwise distances. 

REV: is NOT the general 12 parameter model, but - like in PAML -
denotes the 9 parameter GTR model.

Usage
-----

Example::

   python mali2rates.py --help

Type::

   python mali2rates.py --help

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
import math
import tempfile
import shutil
import CGAT.Experiment as E
import CGAT.Mali as Mali
import CGAT.Genomics as Genomics
import CGAT.WrapperPhylip as WrapperPhylip
import CGAT.WrapperCodeML as WrapperCodeML
import CGAT.RateEstimation as RateEstimation

##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------
def writeModel( grammar, section, options):
    """write a model to output file."""
    if section in options.write or "all" in options.write:
        outfile = open( options.output_pattern % section, "w" )
        outfile.write( "%s\n" % grammar.getGrammar())
        outfile.close()


def CalculateDistancePOVL( seq1, seq2, gap_chars = (".", "-") ):
    """calculate overlap distance between two sequences."""
    naligned, nunaligned = 0, 0
    for c1, c2 in zip( seq1.upper(), seq2.upper()):
        if c1 in gap_chars or c2 in gap_chars:
            if not (c1 in gap_chars and c2 in gap_chars):
                nunaligned += 1
        else:
            naligned += 1
        
    return 100.0 * float( naligned) / float(nunaligned + naligned), 0

def CalculateDistancePID( seq1, seq2,
                          gap_chars = ("-", "."),
                          mask_chars = ("X"),
                          ):

    nidentical, naligned, nunaligned = 0, 0, 0
    for c1, c2 in zip( seq1.upper(), seq2.upper()):
        if c1 in gap_chars or c2 in gap_chars:
            nunaligned += 1
            continue
        naligned += 1        
        if c1 == c2:
            nidentical += 1

    return 100.0 * float(nidentical) / naligned, 0

def CalculateDistanceJC69( info, do_gamma = False, alpha = None ):
    """return Jukes-Cantor distance.
    """

    p = float(info.mNDifferent) / info.mNAligned
    
    if do_gamma:
        ## not done yet
        distance = 0.75 * alpha * (pow(1 - 4*p/3, -1/ alpha) - 1)
        variance = p*(1 - p)/(pow(1 - 4*p/3, -2/(alpha + 1)) * L)

    else:
        distance = -0.75 * math.log(1.0 - 4.0 * p / 3.0);
        variance = p * ( 1.0 - p) / ( math.pow(1.0 - 4.0 * p / 3, 2.0) * info.mNAligned);
        
    return distance, variance

def CalculateDistanceT92( info ):
    """
    P,Q: transitions, transversions frequencies
    q: G+C content
    
    d = -2q(1 - q)loge(1 - P/[2q(1 - q)] - Q) -[1 -2q(1 -q)]loge(1 - 2Q)/2,(4.18)
    V(d) = [c12P + c32Q - (c1P + c3Q)2]/n,(4.19)
    where c1 = 1/(1 - P/[2q(1 - q)] - Q), c2 = 1/(1 - 2Q), c3 = 2q(1 - q)(c1 - c2) + c2, and q is the G+C content

    Note: result is undefined if
        the number of transversions is >= 0.5
        the G+C content is 0
    
    """
    gc = info.getGCContent()

    # if there are no GC or no AT pairs: result is undefined
    if gc == 0 or gc == 1:
        return -1, -1
    
    wg = 2.0 * gc * (1.0 - gc)

    P = float(info.mNTransitions) / info.mNAligned
    Q = float(info.mNTransversions) / info.mNAligned

    a1 = 1.0 - P / wg - Q;
    if a1 <= 0: return -1, -1
    
    a2 = 1.0 - 2.0 * Q;
    if a2 <= 0: return -1, -1

    # print a1, a2, wg, gc, "p=", P, "q=", Q, str(info)

    distance = -wg * math.log(a1) - 0.5 * (1.0 - wg) * math.log(a2);

    c1 = 1/a1;
    c2 = 1/a2;
    c3 = wg*(c1 - c2) + c2;
    
    variance = (c1 * c1 * P + c3 * c3 * Q - math.pow(c1 * P + c3 * Q, 2.0))/ info.mNAligned;

    return distance, variance

def writePhylipResult( result, options ):
    """write result of phylip run."""
    
    options.stdout.write( "seq1\tseq2\tdist\n" )

    for x in range( len(result.mRowHeaders) - 1 ):
        for y in range( x+1, len(result.mColHeaders) ):
            options.stdout.write( "\t".join( map(str, (
                result.mRowHeaders[x],
                result.mColHeaders[y],
                options.format % result.mMatrix[x,y] ))) + "\n" )

def runDNADIST( mali, pairs, options ):
    """run dnadist."""
    # use phylip for these
    phylip = WrapperPhylip.Phylip()
    phylip.setProgram( "dnadist" )
    phylip.setMali( mali )

    phylip_options = []
    if options.distance == "K80":
        phylip_options += ["D"] * 1
    elif options.distance == "JC69":
        phylip_options += ["D"] * 2
    elif options.distance == "LogDet":
        phylip_options += ["D"] * 3

    phylip_options.append("Y")
    phylip.setOptions( phylip_options )

    if options.dump:
        phylip.setLogLevel(2)
    
    result = phylip.run()

    writePhylipResult( result, options )


def printValue( value, format = "%f" ):
    if value == "na":
        return "na"
    else:
        return format % value

##-------------------------------------------------------------------
##-------------------------------------------------------------------
##-------------------------------------------------------------------
def printPair( pair, mali, map_new2old, options, msg ="" ):
    """print pairs form codeml."""

    ids = mali.getIdentifiers()
    if options.output_format == "list":
        options.stdout.write( "\t".join( (map_new2old[ids[0]],
                                          map_new2old[ids[1]],
                                          options.format % pair.mDistanceMatrix[ids[0]][ids[1]],
                                          options.format % pair.mLogLikelihood,
                                          printValue( pair.mAlpha, options.format),
                                          printValue( pair.mKappa, options.format),
                                          msg ) ) )

    elif options.output_format == "tree":
        options.stdout.write( ">pair%i" % (noutput + 1))
        options.stdout.write( "%s\n" % pair.mTree )

    if options.with_counts:
        info = Genomics.CalculatePairIndices( mali[ids[1]], mali[ids[0]], with_codons = options.is_codons )
        options.stdout.write( "\t%s" % (str(info)) )

    options.stdout.write( "\n" )

    return 1

def runBaseML( mali, pairs, options ):

    baseml = WrapperCodeML.BaseML()

    paml_options = {}

    map_new2old = mali.mapIdentifiers()
    ids = mali.getIdentifiers()

    if options.kappa != None:
        paml_options[ "kappa" ] = str(options.kappa)

    if options.fix_kappa:
        paml_options[ "fix_kappa" ] = "1"

    if options.alpha != None:
       paml_options[ "alpha" ] = str(options.alpha)

    if options.fix_alpha:
        paml_options[ "fix_alpha" ] = "1"

    if options.clean_data:
        paml_options[ "cleandata" ] = options.clean_data

    map_distance2index = {}
    for key,val in baseml.mOptions["model"].items():
        map_distance2index[val] = key
        
    if options.distance.upper() in map_distance2index:
        paml_options[ "model"] = map_distance2index[options.distance]
    else:
        raise "unknown distance for baseml: %s" % options.distance

    if options.filename_tree:
        result = baseml.Run( mali,
                             tree = options.filename_tree,
                             dump = options.dump,
                             test = options.test,
                             options = paml_options )

    elif options.pairwise:
        
        noutput = 0
        ninput = 0
        ntotal = (len(ids) * (len(ids) -1 ) ) / 2

        if options.output_format == "list":
            options.stdout.write( "\t".join( ("seq1", "seq2", "distance", "lnL", "alpha", "kappa", "msg" ) ) )

            if options.with_counts:
                options.stdout.write("\t%s" % Genomics.SequencePairInfo().getHeader() )
            options.stdout.write( "\n" )

        for x,y in pairs:
            m1 = mali.getSequence( ids[x] )            
            ninput += 1
            temp_mali = Mali.Mali()
            m2 = mali.getSequence( ids[y] )

            temp_mali.addSequence( ids[x], m1.mFrom, m1.mTo, m1.mString )
            temp_mali.addSequence( ids[y], m2.mFrom, m2.mTo, m2.mString )

            result = baseml.Run( temp_mali,
                                 tree = "(%s,%s);" % (ids[x], ids[y]),
                                 dump = options.dump,
                                 test = options.test,
                                 options = paml_options )

            if options.loglevel >= 1 and ninput % options.report_step == 0:
                options.stdlog.write( "# pairwise computation: %i/%i -> %i%% in %i seconds.\n" % (ninput, ntotal, 100.0 * ninput / ntotal, time.time() - tstart ) )
                options.stdlog.flush()

            noutput += printPair( result,
                                  temp_mali,
                                  map_new2old,
                                  options )

            options.stdout.flush()
            
    else:
        ## assume that there are only two sequences
        if mali.getLength() == 2:
            id1, id2 = mali.getIdentifiers()
            result = baseml.Run( mali,
                                 tree = "(%s,%s);" % (id1, id2),
                                 dump = options.dump,
                                 test = options.test,
                                 options = paml_options )

        else:
            raise "please supply tree if there are more than two sequences and pairwise mode is not selected."

        if options.output_format == "list":
            all_identifiers = mali.getIdentifiers()
            options.stdout.write( "\t".join( ("seq1", "seq2", "distance", "lnL", "alpha", "kappa" ) ) + "\n" )    
            for x in range( len(all_identifiers) -1 ):
                id_x = all_identifiers[x]
                for y in range(x + 1, len(all_identifiers)):
                    id_y = all_identifiers[y]

                    options.stdout.write( "\t".join( (id_x,
                                                      id_y,
                                                      options.format % result.mDistanceMatrix[id_x][id_y],
                                                      options.format % result.mLogLikelihood,
                                                      options.format % result.mAlpha,
                                                      options.format % result.mKappa ) ) + "\n")

        elif options.output_format == "tree":
            options.stdout.write( "%s\n" % result.mTree )

##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------        
def runXrate( mali, pairs, options ):

    from XGram.Generator.Prebuilt import DNA
    from XGram.Model import Annotation
    import XGram.Run

    xgram = XGram.XGram()
    if options.xrate_min_increment:
        xgram.setMinIncrement( options.xrate_min_increment )

    ninput, noutput, nskipped = 0, 0, 0

    tempdir = tempfile.mkdtemp()
    data = tempdir + "/data"

    if options.distance == "K80":
        model = DNA.buildModel( substitution_model = "k80" )
    elif options.distance == "JC69":
        model = DNA.buildModel( substitution_model = "jc69" )        
    elif options.distance == "REV":
        model = DNA.buildModel( substitution_model = "gtr" )        
    else:
        raise "distance %s not implemented for xrate" % (options.distance)

    writeModel( model, "input", options )
    
    if options.output_format == "list":
        options.stdout.write( "\t".join( ("seq1", "seq2", "distance", "lnL", "alpha", "kappa", "msg" ) ) )

        if options.with_counts:
            options.stdout.write("\t%s" % Genomics.SequencePairInfo().getHeader() )
        options.stdout.write( "\n" )

    for x,y in pairs:
        
        m1 = mali.getSequence( ids[x] )            
        ninput += 1
        temp_mali = Mali.Mali()
        m2 = mali.getSequence( ids[y] )

        temp_mali.addSequence( m1.mId, m1.mFrom, m1.mTo, m1.mString )
        temp_mali.addSequence( m2.mId, m2.mFrom, m2.mTo, m2.mString )

##         if temp_mali.getWidth() < options.min_overlap:
##             if options.loglevel >= 1:
##                 options.stdlog.write("# pair %s-%s: not computed because only %i residues overlap\n" % (mali.getEntry(ids[x]).mId,
##                                                                                                         mali.getEntry(ids[y]).mId,
##                                                                                                         temp_mali.getWidth()) )

##             nskipped += 1
##             continue

        outfile = open(data, "w" )
        temp_mali.writeToFile( outfile, format="stockholm",
                               write_ranges = False,
                               options = ("#=GF NH (%s:1.0)%s;" % tuple(temp_mali.getIdentifiers()),) )
        outfile.close()

        o_alpha, o_kappa = "na", "na"
        o_distance = "na"
        msg = ""

        if options.test_xrate:
            for alpha in (0.1, 0.5, 1.0, 1.5):
                for beta in (0.1, 0.5, 1.0, 1.5):
                    model.mGrammar.setParameter( "alpha", alpha )
                    model.mGrammar.setParameter( "beta", beta )                        
                    result = xgram.train( model, data )
                    trained_model = result.getModel()
                    xalpha, xbeta = \
                            ( trained_model.mGrammar.getParameter( 'alpha' ),
                              trained_model.mGrammar.getParameter( 'beta' ) )
                    # this assumes that the branch length in the input is normalized to 1
                    # this is the normalization constant
                    o_distance = options.format % (2 * xbeta + xalpha)
                    o_kappa = options.format % (xalpha / xbeta)

                    msg = "alpha=%6.4f, beta=%6.4f" % (xalpha, xbeta)

                    options.stdout.write( "\t".join( ("%f" % alpha,
                                                      "%f" % beta,
                                                      o_distance,
                                                      options.format % result.getLogLikelihood(),
                                                      o_alpha,
                                                      o_kappa,
                                                      msg) ) )
                    options.stdout.write("\n" )
            continue
        
        options.stdout.write( "%s\t%s\t" % (m1.mId, m2.mId) )
        
        if options.distance in ( "K80", ):
            result = xgram.train( model, data )
            trained_model = result.getModel()

        elif options.distance in ( "REV", ):
            result = xgram.train( model, data )
            trained_model = result.getModel()
            alpha, beta, gamma, delta, epsilon, theta = \
                ( trained_model.mGrammar.getParameter( 'alpha' ),
                  trained_model.mGrammar.getParameter( 'beta' ),
                  trained_model.mGrammar.getParameter( 'gamma' ),
                  trained_model.mGrammar.getParameter( 'delta' ),
                  trained_model.mGrammar.getParameter( 'epsilon' ),
                  trained_model.mGrammar.getParameter( 'theta' ) )

            pi = trained_model.evaluateTerminalFrequencies( ('A0',) )[('A0',)]
            matrix = trained_model.evaluateRateMatrix( ('A0',) )[('A0',)]
            q, d = RateEstimation.getDistanceGTR( pi, matrix )
            o_distance = options.format % (d)
            o_kappa = ""
            msg = "alpha=%6.4f, beta=%6.4f, gamma=%6.4f, delta=%6.4f, epsilon=%6.4f, theta=%6.4f" % (alpha, beta, gamma, delta, epsilon, theta)

        elif options.distance in ('JC69', ):
            result = xgram.buildTree( model, data )

        if options.distance == "K80":
            alpha, beta = \
                   ( trained_model.mGrammar.getParameter( 'alpha' ),
                     trained_model.mGrammar.getParameter( 'beta' ) )
            # this assumes that the branch length in the input is normalized to 1
            # this is the normalization constant
            o_distance = options.format % (2 * beta + alpha)
            o_kappa = options.format % (alpha / beta)

            msg = "alpha=%6.4f, beta=%6.4f" % (alpha, beta)
            alpha = "na"
            
        elif options.distance == "JC69":

            tree = result.getTree()
            ## multiply distance by tree, as rates are set to 1 and 
            ## thus the matrix is scaled by a factor of 3
            o_distance = options.format % (3.0 * float(re.search("\(\S+:([0-9.]+)\)", tree).groups()[0]))
            o_kappa = "na"
            msg = ""

        writeModel( result.mModel, "trained", options )
    
        options.stdout.write("\t".join( (o_distance,
                                         options.format % result.getLogLikelihood(),
                                         o_alpha,
                                         o_kappa,
                                         msg) ) )

        if options.with_counts:
            info = Genomics.CalculatePairIndices( mali[ids[x]], mali[ids[y]], with_codons = options.is_codons )
            options.stdout.write( "\t%s" % (str(info)) )

        options.stdout.write( "\n" )
    
    shutil.rmtree( tempdir )

if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: mali2rates.py 2781 2009-09-10 11:33:14Z andreas $",
                                    usage = globals()["__doc__"])

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal", "stockholm", "phylip" ),
                      help="input format of multiple alignment"  )

    parser.add_option("-s", "--sites", dest="sites", type="string",
                      help="sites to use [default=%default].", )

    parser.add_option("-f", "--file", dest="filename", type="string",
                      help="filename of multiple alignment (- for stdin) [default=%default].",
                      metavar = "FILE" )

    parser.add_option("-o", "--format", dest="format", type="string",
                      help="format [default=%default].",
                      metavar = "format" )

    parser.add_option("-d", "--distance", dest="distance", type="choice",
                      choices=("PID", "T92", "JC69", "POVL", "F84", "LogDet",
                               "K80", "F81", "HKY85", "TN93", "REV", "UNREST", "REVU",
                               "UNRESTU",
                               "JTT", "PMB", "PAM", "Kimura", "CategoriesModel" ),
                      help="method to use for distance calculation [default=%default]." )

    parser.add_option( "--method", dest="method", type="choice",
                      choices=("phylip", "baseml", "own", "xrate"),
                      help = "program to use for rate calculation.")

    parser.add_option( "--output-format", dest="output_format", type="choice",
                      choices=("list", "tree"),
                      help = "output format.")
                     
    parser.add_option("-m", "--min-sites", dest="min_sites", type="int",
                      help="minimum number of sites for output[default=%default].", )

    parser.add_option("-a", "--alphabet", dest="alphabet", type="choice",
                      choices=("aa", "na", "auto"),
                      help="alphabet to use.", )

    parser.add_option("-t", "--filename-tree", dest="filename_tree", type="string",
                      help="filename with tree information."  )

    parser.add_option( "--set-alpha", dest="alpha", type="float",
                      help="initial alpha value.")

    parser.add_option( "--fix-alpha", dest="fix_alpha", action = "store_true",
                       help="do not estimate alpha.")

    parser.add_option( "--set-kappa", dest="kappa", type="float",
                       help="initial kappa value.")

    parser.add_option( "--fix-kappa", dest="fix_kappa", action = "store_true",
                       help="do not estimate kappa.")

    parser.add_option( "--dump", dest="dump", action="store_true",
                      help="dump output."  )

    parser.add_option( "--test", dest="test", action="store_true",
                      help="test run - does not clean up."  )

    parser.add_option( "--pairwise", dest="pairwise", action="store_true",
                      help="force pairwise comparison.")

    parser.add_option( "--set-clean-data", dest="clean_data", type="choice",
                       choices = ("0", "1"),
                       help="PAML should cleanup data:  0=only gaps within pair are removed, 1=columns in the mali with gaps are removed.")

    parser.add_option( "--with-counts", dest="with_counts", action="store_true",
                      help="output counts of aligned positions, transitions and transversions."  )

    parser.add_option("-w", "--write", dest="write", type="choice", action="append",
                      choices=("input", "trained", "all" ),
                      help="output sections to write for xrate."  )
    
    parser.add_option( "--output-pattern", dest="output_pattern", type="string",
                      help="output pattern for output files." )

    parser.add_option( "--xrate-min-increment", dest="xrate_min_increment", type = float,
                       help="minimum increment to stop iteration in xrate." )

    parser.set_defaults( \
        input_format = "fasta",
        filename_tree = None,
        with_counts = False,
        sites = "d4",
        distance = "T92",
        min_sites = 1,
        filename = "-",
        alphabet="auto",
        format= "%6.4f",
        method="phylip",
        kappa = None,
        fix_kappa = False,
        alpha = None,
        fix_alpha = False,
        dump = False,
        clean_data = None,
        output_format = "list",
        iteration="all-vs-all",
        pairwise=False,
        report_step = 1000,
        output_pattern = "%s.eg",
        write = [],
        test_xrate = False,
        xrate_min_increment = None,
        is_codons = False,
        )
    
    (options, args) = E.Start( parser )

    if options.filename != "-":
        infile = open(options.filename, "r")
    else:
        infile = sys.stdin

    ## read multiple alignment
    if options.pairwise:
        ## read sequences, but not as a multiple alignment. This permits multiple names.
        mali = Mali.SequenceCollection()
        options.iteration = "pairwise"
    else:
        mali = Mali.Mali()

    mali.readFromFile( infile, format = options.input_format )
    
    ids = mali.getIdentifiers()

    if options.alphabet == "auto":
        s = "".join( map( lambda x: x.mString, mali.values())).lower()
        ss = re.sub( "[acgtxn]", "", s)
        if float(len(ss)) < (len(s) * 0.1):
            options.alphabet = "na"
            if mali.getNumColumns() % 3 == 0:
                options.is_codons = True
        else:
            options.alphabet = "aa"

        if options.loglevel >= 1:
            options.stdlog.write("# autodetected alphabet: %s\n" % options.alphabet)
            
    if options.filename != "-":
        infile.close()

    npairs = 0
    nskipped_length = 0
    nskipped_distance = 0

    pairs = []
    if options.iteration == "all-vs-all":
        for x in range(len(ids) - 1):
            for y in range(x+1, len(ids)):
                pairs.append( (x,y) )
    elif options.iteration == "first-vs-all":
        for y in range(1, len(ids)):
            pairs.append( (0,y) )
    elif options.iteration == "pairwise":
        if len(ids) % 2 != 0:
            raise "uneven number of sequences (%i) not compatible with --iteration=pairwise" % len(ids)
        for x in range(0, len(ids), 2):
            pairs.append( (x, x+1) )

    if options.alphabet == "na":

        if options.method=="baseml":
            runBaseML( mali, pairs, options )
        elif options.method == "phylip" and options.distance in ("F84", "K80", "JC69", "LogDet" ):
            runDNADIST( mali, pairs, options )
        elif options.method == "xrate":
            runXrate( mali, pairs, options )
        else:
            if options.is_codons:
                h = Genomics.SequencePairInfoCodons().getHeader()
            else:
                h = Genomics.SequencePairInfo().getHeader()
            options.stdout.write( "seq1\tseq2\tdist\tvar\t%s\n" % ( h ))

            for x, y in pairs:
                id_x = ids[x]
                npairs += 1

                id_y = ids[y]

                info = Genomics.CalculatePairIndices( mali[id_x], mali[id_y], with_codons = options.is_codons )

                if options.distance in ("T92", "JC69"):
                    if options.sites == "d4":
                        seq1, seq2 = Genomics.GetDegenerateSites( mali[id_x], mali[id_y], position=3, degeneracy=4 )

                        if len(seq1) < options.min_sites:
                            nskipped_length += 1
                            continue
                    else:
                        raise "unknown sites %s" % options.sites

                if options.distance == "T92":
                    distance, variance = CalculateDistanceT92( info )
                elif options.distance == "JC69":
                    distance, variance = CalculateDistanceJC69( info )
                elif options.distance == "PID":
                    distance, variance = CalculateDistancePID( mali[id_x], mali[id_y] )
                elif options.distance == "POVL":
                    distance, variance = CalculateDistancePOVL( mali[id_x], mali[id_y] )

                if distance >= 0:
                    options.stdout.write( "\t".join( map(str, (id_x, id_y,
                                                               options.format % distance,
                                                               options.format % variance, info) ) ) + "\n" )
                else:
                    nskipped_distance += 1

    elif options.alphabet == "aa":

        if options.distance in ( "JTT", "PMB", "PAM", "Kimura", "CategoriesModel" ):

            # use phylip for these
            phylip = WrapperPhylip.Phylip()
            phylip.setProgram( "protdist" )
            phylip.setMali( mali )
            
            phylip_options = []
            if options.distance == "PMG":
                phylip_options += ["D"] * 1
            elif options.distance == "PAM":
                phylip_options += ["D"] * 2
            elif options.distance == "Kimura":
                phylip_options += ["D"] * 3
            elif options.distance == "CategoriesModel":
                phylip_options += ["D"] * 4
                
            phylip_options.append("Y")
            phylip.setOptions( phylip_options )
            result = phylip.run()

            writePhylipResult( result, options )

        else:
            options.stdout.write( "id1\tid2\tdist\tvar\n" )

            ## iterate over all pairs of sequences
            for x, y in pairs:
                id_x = ids[x]
                npairs += 1

                id_y = ids[y]

                if options.distance == "PID":
                    distance, variance = CalculateDistancePID( mali[id_x], mali[id_y] )
                elif options.distance == "POVL":
                    ## percentage overlap
                    distance, variance = CalculateDistancePOVL( mali[id_x], mali[id_y] )

                if distance >= 0:
                    options.stdout.write("\t".join( (id_x, id_y,
                                      options.format % distance,
                                      options.format % variance) ) + "\n" ) 
                else:
                    nskipped_distance += 1
        
    if options.loglevel >= 1:
        options.stdlog.write("# nseqs=%i, npairs=%i, nskipped_length=%i, nskipped_distance=%i\n" % ( len(ids), npairs, nskipped_length, nskipped_distance))
    
    E.Stop()

"""
/* Copyright 2005 Emmanuel Paradis <paradis@isem.univ-montp2.fr> */

/* This file is part of the `ape' library for R and related languages. */
/* It is made available under the terms of the GNU General Public */
/* License, version 2, or at your option, any later version, */
/* incorporated herein by reference. */

/* This program is distributed in the hope that it will be */
/* useful, but WITHOUT ANY WARRANTY; without even the implied */
/* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR */
/* PURPOSE.  See the GNU General Public License for more */
/* details. */

/* You should have received a copy of the GNU General Public */
/* License along with this program; if not, write to the Free */
/* Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, */
/* MA 02111-1307, USA */

#include <R.h>

#define both_non_n(a, b) (strcmp(a, "n") && strcmp(b, "n"))

void dist_dna_JC69(char **x, int *n, int *s, double *d, int *pairdel,
		   int *variance, double *var, int *gamma, double *alpha)
{
  int i, j, k, s1, s2, target, Nd, L;
  double p;

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) Nd += 1;
      }
      if (!*pairdel) L = *s;
      p = ((double) Nd/L);
      if (*gamma) d[target] = 0.75 * *alpha * (pow(1 - 4*p/3, -1/ *alpha) - 1);
      else d[target] = -0.75 * log(1 - 4 * p/3);
      if (*variance) {
	if (*gamma) var[target] = p*(1 - p)/(pow(1 - 4*p/3, -2/(*alpha + 1)) * L);
	else var[target] = p*(1 - p)/(pow(1 - 4*p/3, 2)*L);
      }
    }
  }
}

void dist_dna_K80(char **x, int *n, int *s, double *d, int *pairdel,
		   int *variance, double *var, int *gamma, double *alpha)
{
  int i, j, k, s1, s2, target, Nd, Ns, L;
  double P, Q, a1, a2, b, c1, c2, c3;

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P = ((double) Ns/L);
      Q = ((double) (Nd - Ns)/L);
      a1 = 1 - 2 * P - Q;
      a2 = 1 - 2 * Q;
      if (*gamma) {
	b = -1 / *alpha;
	d[target] = *alpha * (pow(a1, b) + 0.5 * pow(a2, b) - 1.5) / 2;
      }
      else d[target] = -0.5 * log(a1 * sqrt(a2));
      if (*variance) {
	if (*gamma) {
	  b = -(1 / *alpha + 1);
	  c1 = pow(a1, b);
	  c2 = pow(a2, b);
	  c3 = (c1 + c2) / 2;
	} else {
	  c1 = 1 / a1;
	  c2 = 1 / a2;
	  c3 = (c1 + c2) / 2;
	}
	var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;
      }
    }
  }
}

void dist_dna_F81(char **x, int *n, int *s, double *d, double *BF,
		  int *pairdel, int *variance, double *var,
		  int *gamma, double *alpha)
{
  int i, j, k, Nd, L, target, s1, s2;
  double E, p;

  E = 1 - BF[0]*BF[0] - BF[1]*BF[1] - BF[2]*BF[2] - BF[3]*BF[3];

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) Nd += 1;
      }
      if (!*pairdel) L = *s;
      p = ((double) Nd/L);
      if (*gamma) d[target] = E * *alpha * (pow(1 - p/E, -1/ *alpha) - 1);
      else d[target] = -E*log(1 - p/E);
      if (*variance) {
	if (*gamma) var[target] = p*(1 - p)/(pow(1 - p/E, -2/(*alpha + 1)) * L);
	else var[target] = p*(1 - p)/(pow(1 - p/E, 2)*L);
      }
    }
  }
}

void dist_dna_K81(char **x, int *n, int *s, double *d, int *pairdel,
		   int *variance, double *var)
{
  int i, j, k, Nd, Nv1, Nv2, L, s1, s2, target;
  double P, Q, R, a1, a2, a3, a, b, c;

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Nv1 = 0;
      Nv2 = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a")) {
	    if (!strcmp(x[s2], "t")) Nv1 += 1;
	    if (!strcmp(x[s2], "c")) Nv2 += 1;
	  } else {
	    if (!strcmp(x[s1], "g")) {
	      if (!strcmp(x[s2], "c")) Nv1 += 1;
	      if (!strcmp(x[s2], "t")) Nv2 += 1;
	    } else {
	      if (!strcmp(x[s1], "c")) {
		if (!strcmp(x[j * *s + k], "g")) Nv1 += 1;
		if (!strcmp(x[j * *s + k], "a")) Nv2 += 1;
	      } else {
		if (!strcmp(x[s1], "t")) {
		  if (!strcmp(x[s2], "a")) Nv1 += 1;
		  if (!strcmp(x[s2], "g")) Nv2 += 1;
		}
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P = ((double) (Nd - Nv1 - Nv2)/L);
      Q = ((double) Nv1/L);
      R = ((double) Nv2/L);
      a1 = 1 - 2*P - 2*Q;
      a2 = 1 - 2*P - 2*R;
      a3 = 1 - 2*Q - 2*R;
      d[target] = -0.25 * log(a1 * a2 * a3);
      if (*variance) {
        a = (1/a1 + 1/a2)/2;
	b = (1/a1 + 1/a3)/2;
	c = (1/a2 + 1/a3)/2;
        var[target] = (a*a*P + b*b*Q + c*c*R - pow(a*P + b*Q + c*R, 2))/2;
      }
    }
  }
}

void dist_dna_F84(char **x, int *n, int *s, double *d, double *BF,
		  int *pairdel, int *variance, double *var)
{
  int i, j, k, Nd, Ns, L, target, s1, s2;
  double P, Q, A, B, C, a, b, t1, t2, t3;

  A = (BF[0] * BF[2])/(BF[0] + BF[2]) + (BF[1] * BF[3])/(BF[1] + BF[3]);
  B = BF[0] * BF[2] + BF[1] * BF[3];
  C = (BF[0] + BF[2]) * (BF[1] + BF[3]);

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P = ((double) Ns/L);
      Q = ((double) (Nd - Ns)/L);
      d[target] = -2*A*log(1 - (P/(2*A) - (A - B)*Q/(2*A*C))) + 2*(A - B - C)*log(1 - Q/(2*C));
      if (*variance) {
        t1 = A*C;
	t2 = C*P/2;
	t3 = (A - B)*Q/2;
        a = t1/(t1 - t2 - t3);
	b = A*(A - B)/(t1 - t2 - t3) - (A - B - C)/(C - Q/2);
	var[target] = (a*a*P + b*b*Q - pow(a*P + b*Q, 2))/2;
      }
    }
  }
}

void dist_dna_T92(char **x, int *n, int *s, double *d, double *BF,
                  int *pairdel, int *variance, double *var)
{
  int i, j, k, Nd, Ns, L, target, s1, s2;
  double P, Q, wg, a1, a2, c1, c2, c3;

  wg = 2 * (BF[1] + BF[2]) * (1 - (BF[1] + BF[2]));

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P = ((double) Ns/L);
      Q = ((double) (Nd - Ns)/L);
      a1 = 1 - P/wg - Q;
      a2 = 1 - 2*Q;
      d[target] = -wg*log(a1) - 0.5*(1 - wg)*log(a2);
      if (*variance) {
        c1 = 1/a1;
        c2 = 1/a2;
        c3 = wg*(c1 - c2) + c2;
        var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;
      }
    }
  }
}

void dist_dna_TN93(char **x, int *n, int *s, double *d, double *BF,
		  int *pairdel, int *variance, double *var,
		  int *gamma, double *alpha)
{
  int i, j, k, Nd, Ns1, Ns2, L, target, s1, s2;
  double P1, P2, Q, A, B, C, gR, gY, k1, k2, k3, k4, w1, w2, w3, c1, c2, c3, c4, b;

  gR = BF[0] + BF[2];
  gY = BF[1] + BF[3];
  k1 = 2 * BF[0] * BF[2] / gR;
  k2 = 2 * BF[1] * BF[3] / gY;
  k3 = 2 * (gR * gY - BF[0]*BF[2]*gY/gR - BF[1]*BF[3]*gR/gY);

  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns1 = 0;
      Ns2 = 0;
      L = 0;
      for (k = 0; k < *s; k++) {
        s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns1 += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns1 += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns2 += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns2 += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L = *s;
      P1 = ((double) Ns1/L);
      P2 = ((double) Ns2/L);
      Q = ((double) (Nd - Ns1 - Ns2)/L);
      w1 = 1 - P1/k1 - Q/(2*gR);
      w2 = 1 - P2/k2 - Q/(2*gY);
      w3 = 1 - Q/(2*gR*gY);
      if (*gamma) {
        k4 = 2*(BF[0]*BF[2] + BF[1]*BF[3] + gR*gY);
        b = -1 / *alpha;
        c1 = pow(w1, b);
        c2 = pow(w2, b);
        c3 = pow(w3, b);
        c4 = k1*c1/(2*gR) + k2*c2/(2*gY) + k3*c3/(2*gR*gY);
        d[target] = *alpha * (k1*pow(w1, b) + k2*pow(w2, b) + k3*pow(w3, b) - k4);
      } else {
        k4 = 2*((BF[0]*BF[0] + BF[2]*BF[2])/(2*gR*gR) + (BF[2]*BF[2] + BF[3]*BF[3])/(2*gY*gY));
        c1 = 1/w1;
        c2 = 1/w2;
        c3 = 1/w3;
        c4 = k1 * c1/(2 * gR) + k2 * c2/(2 * gY) + k4 * c3;
        d[target] = -k1*log(w1) - k2*log(w2) - k3*log(w3);
      }
      if (*variance) {
        var[target] = (c1*c1*P1 + c2*c2*P2 + c4*c4*Q - pow(c1*P1 + c2*P2 + c4*Q, 2))/L;
      }
    }
  }
}

void dist_dna_GG95(char **x, int *n, int *s, double *d,
		   int *pairdel, int *variance, double *var)
{
  int i, j, k, s1, s2, target, *GC, gccount, Nd, Ns, *L, length, tl, npair;
  double *theta, gcprop, *P, pp, *Q, qq, *tstvr, svr, A,
    sum, ma /* mean alpha */, K1, K2;

  GC = &gccount;
  theta = &gcprop;
  L = &length;
  P = &pp;
  Q = &qq;
  tstvr = &svr;

  npair = *n * (*n - 1) / 2;

  GC = (int*)R_alloc(*n, sizeof(int));
  theta = (double*)R_alloc(*n, sizeof(double));
  L = (int*)R_alloc(npair, sizeof(int));
  P = (double*)R_alloc(npair, sizeof(double));
  Q = (double*)R_alloc(npair, sizeof(double));
  tstvr = (double*)R_alloc(npair, sizeof(double));

  /* get the proportion of GC, theta, in each sequence */
  for (i = 0; i < *n; i++) {
    tl = 0;
    GC[i] = 0;
    for (k = 0; k < *s; k++) {
      s1 = i * *s + k;
      if (*pairdel) {
	if (strcmp(x[s1], "n")) tl += 1;
	else continue;
      }
      if (!strcmp(x[s1], "c") || !strcmp(x[s1], "g"))
	GC[i] += 1;
    }
    if (!*pairdel) tl = *s;
    theta[i] = ((double) GC[i]/tl);
  }

  /* get the proportions of transitions and transversions,
     and the estimates of their ratio for each pair; we
     also get the sample size for each pair in L  */
  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      Nd = 0;
      Ns = 0;
      L[target] = 0;
      for (k = 0; k < *s; k++) {
	s1 = i * *s + k;
	s2 = j * *s + k;
	if (*pairdel) {
	  if (both_non_n(x[s1], x[s2])) L[target] += 1;
	  else continue;
	}
	if (strcmp(x[s1], x[s2])) {
	  Nd += 1;
	  if (!strcmp(x[s1], "a") && !strcmp(x[s2], "g"))
	    Ns += 1;
	  else {
	    if (!strcmp(x[s1], "g") && !strcmp(x[s2], "a"))
	      Ns += 1;
	    else {
	      if (!strcmp(x[s1], "c") && !strcmp(x[s2], "t"))
		Ns += 1;
	      else {
		if (!strcmp(x[s1], "t") && !strcmp(x[s2], "c"))
		  Ns += 1;
	      }
	    }
	  }
	}
      }
      if (!*pairdel) L[target] = *s;
      P[target] = ((double) Ns/L[target]);
      Q[target] = ((double) (Nd - Ns)/L[target]);
      A = log(1 - 2*Q[target]);
      tstvr[target] = 2*(log(1 - 2*P[target] - Q[target]) - 0.5*A)/A;
    }
  }

  /* compute the mean alpha (ma) = mean Ts/Tv */
  sum = 0;
  tl = 0;
  for (i = 0; i < npair; i++)
    /* some values of tstvr are -inf if there is no
       transversions observed: we exclude them */
    if (R_FINITE(tstvr[i])) {
      sum += tstvr[i];
      tl += 1;
    }
  ma = sum/tl;

  /* compute the distance for each pair */
  for (i = 0; i < *n - 1; i++) {
    for (j = i + 1; j < *n; j++) {
      target = *n * i - i*(i + 1)/2 + j - i - 1;
      A = 1 - 2*Q[target];
      K1 = 1 + ma*(theta[i]*(1 - theta[i]) + theta[j]*(1 - theta[j]));
      K2 = ma*pow(theta[i] - theta[j], 2)/(ma + 1);
      d[target] = -0.5*K1*log(A) + K2*(1 - pow(A, 0.25*(ma + 1)));
      if (*variance)
	var[target] = pow(K1 + K2*0.5*(ma + 1)*pow(A, 0.25*(ma + 1)), 2)*Q[target]*(1 - Q[target])/(A*A*L[target]);
    }
  }
}

BF: base frequencies: a c g t
void dist_dna(char **x, int *n, int *s, int *model, double *d,
	      double *BF, int *pairdel, int *variance, double *var,
	      int *gamma, double *alpha)
{
  switch (*model) {
  case 1 : dist_dna_JC69(x, n, s, d, pairdel, variance,
			 var, gamma, alpha); break;
  case 2 : dist_dna_K80(x, n, s, d, pairdel, variance,
			var, gamma, alpha); break;
  case 3 : dist_dna_F81(x, n, s, d, BF, pairdel, variance,
			var, gamma, alpha); break;
  case 4 : dist_dna_K81(x, n, s, d, pairdel, variance,
			var); break;
  case 5 : dist_dna_F84(x, n, s, d, BF, pairdel, variance,
			var); break;
  case 6 : dist_dna_T92(x, n, s, d, BF, pairdel, variance,
			var); break;
  case 7 : dist_dna_TN93(x, n, s, d, BF, pairdel, variance,
			 var, gamma, alpha); break;
  case 8 : dist_dna_GG95(x, n, s, d, pairdel,
			 variance, var); break;
  }
}
"""
        

        

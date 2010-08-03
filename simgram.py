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
simgram.py - 
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

   python simgram.py --help

Type::

   python simgram.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os, sys, string, re, tempfile, subprocess, optparse, time, math, shutil, random

from types import *

USAGE = """python simgram.py [OPTIONS] < mali.in > out

Simulate alignments using simgram
"""

import Experiment as Experiment
import Mali
import Genomics
import TreeTools

# imports for xrate computation
from XGram.Generator.Prebuilt import Codons
from XGram.Model import Annotation
import XGram.Run
from XGram.Generator.Prebuilt import DNA

import Bio.Data.CodonTable

class Error(Exception):
    """Base class for exceptions in this module."""
    def __str__(self):
        return str(self.message)
    def _get_message(self, message): return self._message
    def _set_message(self, message): self._message = message
    message = property(_get_message, _set_message)

class ParsingError(Error):
    """Exception raised for errors while parsing

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message, line):
        self.message = message + " at line " + line

class UsageError(Error):
    """Exception raised for errors while starting

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class SimgramResult:
    def __init__(self):
        self.mStockholm = None
        self.mMali = Mali.Mali()

class WrapperSimgram:

    mExecutable = "simgram"
    
    def __init__(self):
        pass

    ##------------------------------------------------------------------------
    def run( self,
             grammar,
             tree = None,
             dump = 0,
             test = False,
             options = {} ):

        self.mTempdir = tempfile.mkdtemp()
        self.mFilenameGrammar = "grammar.eg"
        self.mFilenameTree = "tree.nh"        
        self.mFilenameOutput = None
        self.mWarnings = []
        
        if test:
            print "# temporary directory is %s" % self.mTempdir
            
        outfile = open(self.mTempdir + "/" + self.mFilenameGrammar, "w")
        outfile.write( grammar.getGrammar() )
        outfile.close()

        if tree:

            outfile = open(self.mTempdir + "/" + self.mFilenameTree, "w" )
            
            ## check what kind of tree is given.
            if type(tree) == StringType:
                t = tree.strip()
                if t[0] == "(" and t[-1] in ");":
                    outfile.write("%s\n" % t)
                
                else:
                    nexus = TreeTools.Newick2Nexus( open(tree, "r" ) )                    
                    t = nexus.trees[0]
                    outfile.write("%s\n" % TreeTools.Tree2Newick(t))                    

            outfile.close()

        # use your own random seed. Time won't do, if simgram
        # is called in quick succession.
        # Are there any restrictions on seeds? Ian using an even number.
        statement = "%s -rndseed %i -g %s -t %s" % (self.mExecutable,
                                                    random.randint(0, 4294967296),
                                                    self.mFilenameGrammar,
                                                    self.mFilenameTree )
        
        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = self.mTempdir,
                              close_fds = True)                              

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise UsageError, "Error in running %s \n%s\n%s\nTemporary directory in %s" % (self.mExecutable, err, out, self.mTempdir)

        if dump: 
            print "# stdout output of %s:\n%s\n######################################" % (self.mExecutable, out)

        if not test:
            shutil.rmtree( self.mTempdir )

        return self.parseOutput( out.split("\n") )

    def parseOutput( self, lines ):
        """parse stdout output from simgram."""

        result = SimgramResult()

        result.mStockholm = lines
        result.mMali.readFromFile( lines, format="stockholm" )

        return result

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
def countSites( model ):
    """count number of expected synonymous/nonsynonymous sites in a grammar.
    """

    ## number of synonymous/non-synonymous sites
    n, s = 0.0, 0.0

    xpi = model.evaluateTerminalFrequencies()[('COD0', 'COD1', 'COD2')]
    pi = {}
    for codon, f in xpi.items():
        pi["".join(codon).upper()] = f

    ## translate pi and the matrix to codons
    for key, value in pi.items():
        del pi[key]
        pi["".join(key).upper()] = value
    
    for codon, freq in pi.items():
        
        try:
            degeneracy = Genomics.GetDegeneracy( codon )
        except KeyError:
            continue

        for x in range(1,4):
            d = (degeneracy[x] - 1.0 ) / 3.0
            s += freq * d
            n += freq * (1.0-d)
             
##              if degeneracy[x] > 1:
##                  s += freq
##              else:
##                  n += freq

    assert( float("%5.2f" % (n + s)) == 3.0 )

##     print s / (n+s)
    
##     n = 184.9
##     s = 76.1
##     t = n + s
##     n /= t
##     s /= t
##     print s / (n+s)
    
    return n, s

if __name__ == "__main__":
    
    parser = optparse.OptionParser( version = "%prog version: $Id: simgram.py 2781 2009-09-10 11:33:14Z andreas $" )

    parser.add_option( "-m", "--model", dest="model", type="choice", 
                      choices=( "custom", "sn", "akaksgc", "K80", "f3x4-four" ),
                      help="grammar to use for simulation."  )

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal", "stockholm", "phylip" ),
                      help="input format of multiple alignment"  )

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=("plain", "fasta", "clustal", "stockholm", "phylip" ),
                      help="output format of multiple alignment"  )

    parser.add_option( "--set-omega", dest="omega", type="float",
                      help="set omega (non-synonymous/synonymous rate ratio).")

    parser.add_option( "--set-kappa", dest="kappa", type="float",
                       help="set kappa (transition/transversion ratio).")

    parser.add_option( "--set-ds", dest="ds", type="float",
                       help="set divergence.")

    parser.add_option("-w", "--write", dest="write", type="choice", action="append",
                      choices=("input_fixed", "trained_fixed", "input_variable", "trained_variable", "all" ),
                      help="output sections to write."  )

    parser.add_option( "--output-pattern", dest="output_pattern", type="string",
                       help="output pattern for output files." )

    parser.add_option("--insert-frequencies", dest="insert_frequencies", type="string",
                      help="insert frequencies from a multiple alignment." )

    parser.add_option("--uniform-frequencies", dest="insert_frequencies", action="store_false",
                      help="use uniform codon frequencies." )

    parser.add_option( "--test", dest="test", action="store_true",
                      help="run a test."  )

    parser.add_option( "--dump", dest="dump", action="store_true",
                      help="dump output."  )

    parser.add_option( "--num-replicates", dest="num_replicates", type="int",
                       help="number of replicates to output." )

    parser.add_option( "--length", dest="length", type="int",
                       help="length of simulated alignment. Set to 0 (default) for variable length." )

    parser.add_option( "--remove-stop-codons", dest="remove_stop_codons", action="store_true",
                       help="remove positions in mali with stop-codons." )

    parser.set_defaults(
        input_format="fasta",
        output_format="fasta",        
        model = "akaksgc",
        omega = None,
        kappa = None,
        ds = None,
        dump = False,
        test = False,
        insert_frequencies = None,
        write = [],
        output_pattern = "%s.eg",
        value_format = "%6.4f",
        debug = False,
        fix_rates = None,
        filename_grammar = None,
        num_replicates = 1,
        remove_stop_codons = False,
        length = 0,
        )

    (options, args) = Experiment.Start( parser )

    if options.fix_rates: options.fix_rates = map( float, options.fix_rates.split(",") )

    if options.loglevel >= 1:
        if options.omega:
            o_omega = "%6.4f" % options.omega
        else:
            o_omega = "na"
        if options.kappa:
            o_kappa = "%6.4f" % options.kappa
        else:
            o_kappa = "na"
        if options.ds:
            o_ds = "%6.4f" % options.ds
        else:
            o_ds = "na"
            
        options.stdlog.write("# input parameters: model=%s, ds=%s, omega=%s, kappa=%s\n" % ( options.model, o_ds, o_omega, o_kappa) )

    ## load a grammar
    if options.model in ( "sn" , "akaksgc", "f3x4-four" ):

        if options.model in ("sn", ):
            infile = open(XGram.PATH_DATA + "/sn.eg", "r")
            input_model = XGram.Parser.parseGrammar( infile.readlines() )
            
        elif options.model in ( "akaksgc", ):
            infile = open(XGram.PATH_DATA + "/akaksgc.eg", "r")
            input_model = XGram.Parser.parseGrammar( infile.readlines() )
            
        elif options.model in ( "f3x4-four", ):
            input_model = Codons.buildCodonML(codon_model = options.model,
                                              explicit_extension = True,
                                              fix_kappa = options.kappa == None,
                                              fix_omega = options.omega == None )


        ## set codon usage frequencies
        if options.insert_frequencies:

            mali = Mali.Mali()

            mali.readFromFile( open(options.insert_frequencies, "r"),
                               format = options.input_format )

            if mali.getLength() == 0:
                raise "refusing to process empty alignment."

            frequencies = Codons.getFrequenciesPerCodonPosition( map( lambda x: x.mString, mali.values() ))

            # build a dummy grammar to insert frequencies
            dummy_grammar = XGram.Model.Grammar()
            for x in range(0,3):
                params = []
                for a in ('A', 'C', 'G', 'T'):
                    params.append( ("p%s%i" % (a.lower(), x), frequencies[x][a]) )
                dummy_grammar.addVariable( params )

            input_model.mGrammar.copyParameters( dummy_grammar, 
                                                 ignore_missing = True)

        ## count the number of synonymous and non-synonymous sites
        if options.omega or options.kappa:
            n, s = countSites( input_model )
            ps = s / (n + s )

            if options.omega:
                branchlength = 3 * options.ds * ( ps + options.omega * (1 - ps ) )
            else:
                branchlength = 3 * options.ds 
                
            if options.loglevel >= 1:
                options.stdlog.write("# derived  parameters: n=%6.4f, s=%6.4f, ps=%6.4f, branchlength=%6.4f\n" % (n, s, ps, branchlength ) )

            if options.kappa and options.omega:

                if options.model in ("akaksgc", ) :
                    n /= 3.0
                    s /= 3.0

                    rni = branchlength * n * options.omega * options.kappa
                    rsi = branchlength * s * options.kappa
                    rnv = branchlength * n * options.omega
                    rsv = branchlength * s

                    input_model.mGrammar.setParameter( "Rsi", rsi )
                    input_model.mGrammar.setParameter( "Rni", rni )
                    input_model.mGrammar.setParameter( "Rnv", rnv )
                    input_model.mGrammar.setParameter( "Rsv", rsv )

                    if options.loglevel >= 1:
                        options.stdlog.write("# computed parameters: rsi=%6.4f, rsv=%6.4f, rni=%6.4f, rnv=%6.4f\n" % (rsi, rsv, rni, rnv) )
                
                elif options.model in ("f3x4-four", ):
                    
                    # branchlength = 3 * options.ds

                    rni = 25.0 * branchlength * options.omega * options.kappa 
                    rsi = 25.0 * branchlength * options.kappa 
                    rnv = 25.0 * branchlength * options.omega 
                    rsv = 25.0 * branchlength 

                    input_model.mGrammar.setParameter( "Rsi", rsi )
                    input_model.mGrammar.setParameter( "Rni", rni )
                    input_model.mGrammar.setParameter( "Rnv", rnv )
                    input_model.mGrammar.setParameter( "Rsv", rsv )

                    if options.loglevel >= 1:
                        options.stdlog.write("# computed parameters: rsi=%6.4f, rsv=%6.4f, rni=%6.4f, rnv=%6.4f\n" % (rsi, rsv, rni, rnv) )

            elif options.kappa:

                ## without omega, we have a plain nucleotide model. 
                ## Because we work in codon space,
                ## the branch length needs to be 20 * 20 / 4 * 4 as long = 25
                alpha = 25.0 * branchlength * options.kappa / (1.0 + 2.0 * options.kappa )
                beta =  25.0 * branchlength / (1.0 + 2.0 * options.kappa )        
                
                input_model.mGrammar.setParameter( "kappa", alpha )
                input_model.mGrammar.setParameter( "not_kappa", beta )                
                
            elif options.omega:

                ## without omega, we have a plain nucleotide model. Because we work in codon space,
                ## the branch length needs to be 20 * 20 / 4 * 4 as long = 25
                omega     = 25.0 * options.ds 
                not_omega = 25.0 * options.ds * options.omega

                input_model.mGrammar.setParameter( "omega", omega )
                input_model.mGrammar.setParameter( "not_omega", not_omega )

    elif options.model in ("K80" ):

        if not (options.omega == None or options.omega == 1.0):
            raise "can only accept 1.0 for omega using the kimura model."

        if options.model == "K80":
            input_model = DNA.buildModel( substitution_model = "k80",
                                          explicit_extension = True )
            
        alpha = options.ds * options.kappa / (1.0 + 2.0 * options.kappa )
        beta  = options.ds / (1.0 + 2.0 * options.kappa )        

        if options.loglevel >= 1:
            options.stdlog.write("# computed parameters: alpha=%6.4f, beta=%6.4f\n" % (alpha, beta) )

        input_model.mGrammar.setParameter( "alpha", alpha )
        input_model.mGrammar.setParameter( "beta", beta )

    ## set ext and not_ext to allow for long chains
    input_model.mGrammar.setParameter( "ext", "0.999" )
    input_model.mGrammar.setParameter( "not_ext", "0.001" )        
        
    writeModel( input_model, "input", options )

    simgram = WrapperSimgram()

    noutput = 0
    last_mali = None
    while noutput < options.num_replicates:
        
        result = simgram.run( input_model, tree = "(seq1:1.0)seq2;",
                              dump = options.dump,
                              test = options.test )

        mali = result.mMali

        if options.remove_stop_codons:
            mali.removePattern( lambda x: x.upper() in ("TAG", "TAA", "TGA"),
                                allowed_matches = 0,
                                minimum_matches = 1,
                                frame = 3 )

            width = mali.getWidth()
            mali.truncate( 0, 3 * int(math.floor( width / 3.0 ) ) )
        
        if options.loglevel >= 1:
            options.stdlog.write("# received mali: %i sequences, %i columns.\n" % (mali.getLength(), mali.getWidth()))

        if last_mali:
            for key, value in last_mali.items():
                mali.getEntry(key).mString += value.mString 

        if options.loglevel >= 1:
            options.stdlog.write("# cumulative  mali: %i sequences, %i columns.\n" % (mali.getLength(), mali.getWidth()))            

        output = True
        if options.length:
            if mali.getWidth() > options.length:
                mali.truncate( 0, options.length )
                output = True
            elif mali.getWidth() < options.length:
                output = False
            
        if output:
            noutput += 1
            mali.writeToFile( sys.stdout, format = options.output_format )

        last_mali = mali
        
    Experiment.Stop()


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
codonbias_acai2tsv.py - codonbias analysis
=========================================

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

   python codonbias_acai2tsv.py --help

Type::

   python codonbias_acai2tsv.py --help

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
import random

"""Wrapper for adaptive codon bias program
"""

import CGAT.Experiment as E
import CGAT.Genomics as Genomics

import CGAT.WrapperAdaptiveCAI as WrapperAdaptiveCAI
import CGAT.IOTools as IOTools
import CGAT.CSV as CSV

## order of codon matrix as expected by caijava.
OUTPUT_ORDER_CODON_MATRIX = ( \
    ("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG"),
    ("TAT","TAC","TAA","TAG","TGT","TGC","TGA","TGG"),
    ("CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG"),
    ("CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG"),
    ("ATT","ATC","ATA","ATG","ACT","ACC","ACA","ACG"),
    ("AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG"),
    ("GTT","GTC","GTA","GTG","GCT","GCC","GCA","GCG"),
    ("GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"),
)

# a list of codon preferences (counts)
CODON_PREFERENCES = {
    # taken from Shields et al. (1988), table 2
    'dmelanogaster' : {
    "TTT" : 13,
    "TTC" : 145,
    "TTA" : 3,
    "TTG" : 39,
    "TCT" : 30,
    "TCC" : 127,
    "TCA" : 5,
    "TCG" : 57,
    "TAT" : 32,
    "TAC" : 140,
    "TAA" : 14, # stop
    "TAG" : 1, # stop
    "TGT" : 3,
    "TGC" : 54,
    "TGA" : 0, # stop
    "TGG" : 39,
    "CTT" : 14,
    "CTC" : 43,
    "CTA" : 7,
    "CTG" : 236,
    "CCT" : 21,
    "CCC" : 138,
    "CCA" : 31,
    "CCG" : 14,
    "CAT" : 13,
    "CAC" : 72,
    "CAA" : 12,
    "CAG" : 156,
    "CGT" : 75,
    "CGC" : 97,
    "CGA" : 4,
    "CGG" : 0,
    "ATT" : 61,
    "ATC" : 213,
    "ATA" : 0,
    "ATG" : 114,
    "ACT" : 34,
    "ACC" : 216,
    "ACA" : 5,
    "ACG" : 12,
    "AAT" : 17,
    "AAC" : 185,
    "AAA" : 13,
    "AAG" : 346,
    "AGT" : 1,
    "AGC" : 63,
    "AGA" : 1,
    "AGG" : 9,
    "GTT" : 56,
    "GTC" : 131,
    "GTA" : 9,
    "GTG" : 163,
    "GCT" : 95,
    "GCC" : 299,
    "GCA" : 18,
    "GCG" : 19,
    "GAT" : 126,
    "GAC" : 160,
    "GAA" : 31,
    "GAG" : 313,
    "GGT" : 130,
    "GGC" : 160,
    "GGA" : 83,
    "GGG" : 0,
    }
    }


if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: codonbias_acai2tsv.py 865 2007-01-15 13:44:43Z andreas $")

    parser.add_option("-o", "--input-file-trace", dest="input_filename_trace", type="string",
                      help="input filename for cai.",
                      metavar="FILE" )

    parser.add_option("-e", "--input-file-genes", dest="input_filename_genes", type="string",
                      help="input filename for genes information from cai.",
                      metavar="FILE" )

    parser.add_option("-c", "--input-file-codons", dest="input_filename_codons", type="string",
                      help="input filename for codon usage information.",
                      metavar="FILE" )

    parser.add_option( "--input-file-sequences", dest="input_filename_sequences", type="string",
                      help="input filename with sequences.",
                      metavar="FILE" )

    parser.add_option("-t", "--input-file-subset", dest="input_filename_subset", type="string",
                      help="input filename with subset.",
                      metavar="FILE" )

    parser.add_option("--codon-table-format", dest="codon_table_format", type="choice",
                      choices=("list", "matrix"),
                      help="output options for output codon tables." )

    parser.add_option("--codon-table-type", dest="codon_table_type", type="choice",
                      choices=("counts", "frequencies", "weights", "absolute-frequencies"),
                      help="type of codon table." )

    parser.add_option("-r", "--reference", dest="reference", type="string",
                      help="dump CAI reference weights for species." )

    parser.add_option("-s", "--select", dest="select", type="string",
                      help="fields to select from genes table." )

    parser.add_option("-m", "--map", dest="input_filename_map", type="string",
                      help="filename with mapping information for gene names.",
                      metavar="FILE")

    parser.add_option("-i", "--invert-map", dest="invert_map", action="store_true",
                      help="invert map.")

    parser.add_option("-d", "--dominant-set", dest="dominant_set", type ="float",
                      help="only print out dominant set (# fraction of most biased genes).")

    parser.add_option("--reverse-set", dest="reverse_set", action="store_true",
                      help="print the reverse set, i.e., then non-dominant set.")

    parser.add_option("-u", "--codon-usage", dest="codon_usage", type="string",
                      help="print codon usage for the full/biased set of genes [full|biased].")

    parser.add_option("-w", "--weights", dest="weights", type="string",
                      help="print weights [final-list|final-matrix|random|compute|weights|frequencies|absolute-frequencies].")

    parser.add_option( "--weights-matrix2table", dest="weights_matrix2table", action="store_true",
                      help="convert a weights matrix to a weights table.")

    parser.add_option( "--get-preferred-codons", dest="get_preferred_codons", type="string",
                      help="compute overview of preferred codons.")

    parser.set_defaults( \
        input_filename = "-",
        input_filename_trace = None,
        input_filename_genes = None,
        input_filename_codons = None,
        input_filename_map = None,
        input_filename_subset = None,
        input_filename_sequences = None,
        invert_map = False,
        select = None,
        codon_usage = None,
        weights = None,
        revserse_set = False,
        pseudocounts = 1,
        codon_table_format = "list",
        codon_table_type = "weights",
        weights_matrix2table = False,
        random_size = 1000,
        get_preferred_codons = None,
        dominant_set = 0.0 )

    (options, args) = E.Start( parser )
    if options.select:
        options.select = options.select.split(",")

    outfile = options.stdout

    ###################################################################
    ## convert weights table to a codon table
    if options.weights_matrix2table:
        lines = options.stdin.readlines()
        data = []
        for line in lines:
            if line[0] == "#": continue
            data += list(map(float, line[:-1].split(",")))

        weights = {}
        x = 0
        for cc in OUTPUT_ORDER_CODON_MATRIX:
            for c in cc:
                weights[c] = data[x]
                x += 1

        outfile.write ("CODON\tWEIGHT\n")
        codons = weights.keys()
        codons.sort()
        for codon in codons:
            outfile.write( "%s\t%f\n" % ( codon, weights[codon]) )

        E.Stop()
        sys.exit(1)

    ###################################################################
    map_genes = {}

    if options.input_filename_map:
        data = map( lambda x: x[:-1].split("\t")[:2], filter( lambda x: x[0] != "#", open( options.input_filename_map, "r").readlines()))

        for a, b in data:
            if options.invert_map: a,b = b,a
            map_genes[a] = b

    result = WrapperAdaptiveCAI.AdaptiveCAIResult()

    if options.input_filename_genes:
        gene_file = open(options.input_filename_genes,"r")
    else:
        gene_file = None

    if options.input_filename_codons:
        codon_file = open(options.input_filename_codons,"r")
    else:
        codon_file = None

    if options.input_filename_trace:
        trace_file = open(options.input_filename_trace,"r")
    else:
        trace_file = None

    if options.input_filename_subset:
        l, e = IOTools.ReadList( open(options.input_filename_subset,"r") )
        subset = set(l)
        if options.loglevel >= 1:
            options.stdlog.write( "# read %i entries into subset from %s.\n" % (len(subset), options.input_filename_subset))
    else:
        subset = None
        
    result.Read( gene_file=gene_file, codon_file=codon_file, trace_file = trace_file )

    if gene_file: gene_file.close()
    if codon_file: codon_file.close()
    if trace_file: trace_file.close()    

    if options.reference:
        if options.reference not in CODON_PREFERENCES:
            raise "unknown species %s: possibles species are: %s" % (options.reference, str(CODON_PREFERNCES.keys()))

        weights = Genomics.CalculateCAIWeightsFromCounts( CODON_PREFERENCES[options.reference], options.pseudocounts )
        
        for x in range(len(OUTPUT_ORDER_CODON_MATRIX)):
            outfile.write(",".join( map(lambda z: "%5.3f" % z, [ weights[codon.upper()] for codon in OUTPUT_ORDER_CODON_MATRIX[x]])))
            outfile.write("\n")

    if options.dominant_set and gene_file:
        cai_threshold = result.GetDominantThreshold( options.dominant_set )
    else:
        if options.reverse_set:
            cai_threshold = 1.0
        else:
            cai_threshold = 0.0
        
    if options.select:

        fields = []
        titles = []
        for x in options.select:
            f = re.match("(\S+) (AS|as) (\S+)", x)
            if f:
                fields.append(f.groups()[0].upper())
                titles.append(f.groups()[2])
            else:
                fields.append(x.upper())
                titles.append(x)
                
        outfile.write( "GENENAME\t" + string.join(titles, "\t") + "\n")
        
        for genename,data in result.mGeneInfo.items():
            if genename in map_genes: genename = map_genes[genename]

            if options.reverse_set:
                if data["CAICLASS"] >= cai_threshold: continue                
            else:
                if data["CAICLASS"] < cai_threshold: continue
            
            outfile.write( genename )
            for c in fields:
                outfile.write("\t%s" % str(data[c]))
            outfile.write("\n")

    if options.weights:

        format = options.codon_table_format

        if options.weights in ("compute-counts", "compute-weights", "compute-frequencies"):
            ## compute codon usage weights from a set of sequences
            codons = CODON_PREFERENCES["dmelanogaster"].keys()
            counts = {}
            for x in codons: counts[x] = 0
            
            if options.input_filename_sequences:
                sequences = Genomics.ReadPeptideSequences( open( options.input_filename_sequences, "r"), filter = subset )
                for key, sequence in sequences.items():
                    sequence = re.sub(" ", "", sequence )
                    if len(sequence) % 3 != 0:
                        raise "warning: sequence %s is not multiple of 3" % key
                    for codon in [ sequence[x:x+3] for x in range(0, len(sequence), 3)]:
                        counts[codon.upper()] += 1
                        
            if options.weights == "compute-frequencies":
                weights = Genomics.CalculateCodonFrequenciesFromCounts( counts, options.pseudocounts )
            elif options.weights == "compute-weights":
                weights = Genomics.CalculateCAIWeightsFromCounts( counts, options.pseudocounts )
            else:
                weights = counts
            
        elif options.weights in ("final-list", "final-matrix"):
            
            weights = result.mFinalWeights
            if options.weights == "final-list":
                format = "list"
            else:
                format = "matrix"

        elif options.weights == "random":
            ## get random weights
            codons = CODON_PREFERENCES["dmelanogaster"].keys()
            counts = {}
            for x in codons:
                counts[x] = random.randint( 1, options.random_size )
                
            weights = Genomics.CalculateCAIWeightsFromCounts( counts, options.pseudocounts )
            format = "matrix"

        elif options.weights == "biased":
            ## get biased weights
            codons = Genomics.GetUniformCodonUsage()

            weights = Genomics.CalculateCAIWeightsFromCounts( counts, options.pseudocounts )
            format = "matrix"

        elif options.weights in ( "uniform-weights", "uniform-frequencies"):
            ## get uniform weights
            codons = Genomics.GetUniformCodonUsage()

            if options.weights == ( "uniform-weights" ):
                weights = Genomics.CalculateCAIWeightsFromCounts( counts, options.pseudocounts )
                format = "matrix"
            else:
                weights = codons
                format = "list"

        elif options.weights in ("counts", "frequencies", "absolute-frequencies"):
            ## get weights as frequencies
            ## compute from scratch. In the caijava file, the absolute frequencey f / gene_length is
            ## given. Thus the total number of codons is f * gene_length.
            codons = CODON_PREFERENCES["dmelanogaster"].keys()            
            counts = {}
            for c in codons: counts[c] = 0
            
            for genename, data in result.mGeneInfo.items():
                
                if options.reverse_set:
                    if data["CAICLASS"] >= cai_threshold: continue                
                else:
                    if data["CAICLASS"] < cai_threshold: continue
                    
                l = data["GENELENGTH"]
                for c in codons:
                    counts[c] += int(data[c] * l)

            if options.weights == "frequencies":
                weights = Genomics.CalculateCodonFrequenciesFromCounts( counts, options.pseudocounts )
            elif options.weights == "counts":
                weights = counts
            elif options.weights == "absolute-frequencies":
                ## compute absolute frequencies (with pseudo-counts, but do not normalize per aa)
                weights = {}
                m = sum(counts.values())
                for k,v in counts.items():
                    weights[k] = float(v) / m
                
            format = "list"

        elif options.weights == "subset":

            codons = CODON_PREFERENCES["dmelanogaster"].keys()            
            counts = {}
            for c in codons: counts[c] = 0
            
            for genename, data in result.mGeneInfo.items():

                found = genename in subset
                if (not found and not options.reverse_set) or (found and options.reverse_set):
                    continue
                    
                l = data["GENELENGTH"]
                for c in codons:
                    counts[c] += int(data[c] * l)
                    
            if options.codon_table_type == "frequencies":
                weights = Genomics.CalculateCodonFrequenciesFromCounts( counts, options.pseudocounts )                
            elif options.codon_table_type == "weights":                
                weights = Genomics.CalculateCAIWeightsFromCounts( counts, options.pseudocounts )
            elif options.codon_table_type == "counts":
                weights = counts
            if options.codon_table_type == "absolute-frequencies":
                ## compute absolute frequencies (with pseudo-counts, but do not normalize per aa)
                weights = {}
                m = sum(counts.values())
                for k,v in counts.items():
                    weights[k] = float(v) / m
                
        else:
            raise "unknown weights %s" % options.weights
        
        if format == "list":
            outfile.write ("CODON\tWEIGHT\n")
            codons = weights.keys()
            codons.sort()
            for codon in codons:
                outfile.write( "%s\t%f\n" % ( codon, weights[codon]) )
                
        elif format == "matrix":

            for x in range(len(OUTPUT_ORDER_CODON_MATRIX)):
                outfile.write(",".join( map(lambda z: "%5.3f" % z, [ weights[codon.upper()] for codon in OUTPUT_ORDER_CODON_MATRIX[x]])))
                outfile.write("\n")
        
    if options.codon_usage:
        outfile.write( "CODON\tFREQUENCY\n" )

        if options.codon_usage == "biased":
            usages = result.mCodonUsages[-1]
        elif options.codon_usage == "full":
            usages = result.mCodonUsages[0]
        elif options.codon_usage == "weights":
            usages = WrapperAdaptiveCAI.CalculateWeightsFromUsage(result.mCodonUsages[0])
        else:
            raise "unknown option '%s' for codon-usage." % options.codon_usage

        codons = usages.keys()
        codons.sort()
        for codon in codons:
            outfile.write( "%s\t%f\n" % ( codon, usages[codon]) )

    E.Stop()
    

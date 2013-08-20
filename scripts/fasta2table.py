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
fasta2table.py - analyze sequences for codon bias and other sequence properties
====================================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads a collection of sequences in :term:`fasta` format and computes
various sequence properties (length, information content, codon bias, ...).

Usage
-----

Example::

   python cat in.fasta | fasta2table.py --sections=cpg

In this example we inout a fasta file and compute the sequence composition i.e.
%C, %G, %A, %T as well as dinucleotide (CpG) composition for each sequence in the 
set.


Type::

   python fasta2table.py --help

for command line help.

Documentation
-------------

This script was formerly called :file:`analyze_codonbias_shannon.py`.

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
import math

import CGAT.Experiment as E
import CGAT.Genomics as Genomics
import CGAT.IOTools as IOTools
from CGAT.SequenceProperties import *
import CGAT.FastaIterator as FastaIterator

##---------------------------------------------------------------------------------------------
    
def main( argv = None ):

    parser = E.OptionParser( version = "%prog version: $Id: analyze_codonbias_shannon.py 2864 2010-03-03 10:18:16Z andreas $",
                                    usage = globals()["__doc__"] )

    parser.add_option("-w", "--filename-weights", dest="filename_weights", type="string",
                      help="filename with codon frequencies. Multiple filenames can be separated by comma." )
    
    parser.add_option("-s", "--section", dest="sections", type="string", action="append",
                      help="which sections to output. Possible choices are length|na|aa|cpg|degeneracy|bias|codons|codon-usage|codon-translator|sequence, [%default]" )

    parser.add_option("-t", "--type", dest="seqtype", type="choice",
                      choices=("na", "aa"),
                      help="type of sequence: na=nucleotides, aa=amino acids [%default].")

    parser.add_option("-e", "--regex-identifier", dest="regex_identifier", type="string",
                      help="regular expression to extract identifier from fasta description line." )

    parser.set_defaults(
        filename_weights = "uniform",
        pseudocounts = 1,
        sections = [],
        regex_identifier = "(.+)",
        seqtype = "na",
        )
    
    (options, args) = E.Start( parser, argv = argv )
    options.filename_weights = options.filename_weights.split(",")

    rx = re.compile( options.regex_identifier )

    reference_codons = []
    if options.filename_weights:
        for filename in options.filename_weights:
            if filename == "uniform":
                reference_codons.append( Genomics.GetUniformCodonUsage() )
            else:
                reference_codons.append( IOTools.ReadMap( open(filename, "r"), has_header = True, map_functions=(str, float) ) )

        ## print codon table differences
        options.stdlog.write("# Difference between supplied codon usage preferences.\n")
        for x in range(0, len(reference_codons)):
            for y in range(0, len(reference_codons)):
                if x == y: continue
                # calculate KL distance
                a = reference_codons[x]
                b = reference_codons[y]
                d = 0
                for codon, p in a.items():
                    if Genomics.IsStopCodon( codon ): continue
                    d += b[codon] * math.log( b[codon] / p )

                options.stdlog.write( "# tablediff\t%s\t%s\t%f\n" % (options.filename_weights[x],
                                                                     options.filename_weights[y],
                                                                     d ) )

    iterator = FastaIterator.FastaIterator( options.stdin )

    def getCounter( section ):

        if options.seqtype == "na":
            if section == "length":
                s = SequencePropertiesLength()
            elif section == "sequence":
                s = SequencePropertiesSequence()
            elif section == "hid":
                s = SequencePropertiesHid()
            elif section == "na":
                s = SequencePropertiesNA()
            elif section == "cpg":
                s = SequencePropertiesCpg()
            elif section == "aa":
                s = SequencePropertiesAA()
            elif section == "degeneracy":
                s = SequencePropertiesDegeneracy()
            elif section == "bias":
                s = SequencePropertiesBias( reference_codons )
            elif section == "codons":
                s = SequencePropertiesCodons()
            elif section == "codon-usage":
                s = SequencePropertiesCodonUsage()
            elif section == "codon-translator":
                s = SequencePropertiesCodonTranslator()
            else:
                raise ValueError("unknown section %s" % section)
        elif options.seqtype == "aa":
            if section == "length":
                s = SequencePropertiesLength()
            elif section == "sequence":
                s = SequencePropertiesSequence()
            elif section == "hid":
                s = SequencePropertiesHid()
            elif section == "aa":
                s = SequencePropertiesAminoAcids()
            else:
                raise ValueError("unknown section %s" % section)
        return s

    ## setup totals
    totals = {}
    for section in options.sections:
        totals[section] = getCounter( section )
        
    options.stdout.write ("id" )
    for section in options.sections:
        options.stdout.write( "\t" + "\t".join(totals[section].getHeaders()) )
        
    options.stdout.write("\n")
    options.stdout.flush()

    for cur_record in iterator:
        
        sequence = re.sub(" ", "", cur_record.sequence).upper()

        if len(sequence) == 0:
            raise ValueError( "empty sequence %s" % cur_record.title )

        id = rx.search( cur_record.title ).groups()[0]

        options.stdout.write("%s" % id)
        options.stdout.flush()

        for section in options.sections:
            s = getCounter( section )
            s.loadSequence( sequence )
            totals[section].addProperties( s )

            options.stdout.write( "\t" + "\t".join(s.getFields()) )
            
        options.stdout.write("\n")
        
    options.stdout.write( "total" )
    for section in options.sections:
        options.stdout.write( "\t" + "\t".join(totals[section].getFields()))
    options.stdout.write("\n")
    
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

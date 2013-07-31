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
codonbias_weights2tsv.py - 
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

   python codonbias_weights2tsv.py --help

Type::

   python codonbias_weights2tsv.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import string
import os
import optparse

import CGAT.Experiment as E
import CGAT.CSV as CSV
import CGAT.Genomics as Genomics

def WriteHeader( options ):
    if options.is_frequencies:
        cat = "frequency"
    else:
        cat = "weight"

    options.stdout.write("""#
# Changes between preferred codons. The output
# is the preferred codon together with the %s of the unpreferred codon. I.e.:
# speciesA        speciesB      R       CGT     0.646000        CGC     0.587000
# means:
# CGT is preferred in speciesA and has a %s of 0.64 in speciesB.
# CGC is preferred in speciesB and has a %s of 0.59 in speciesA.
species1\tspecies2\taa\tcodon1\tweight1\tcodon2\tweight2\tdifference
""" % (cat, cat, cat))

def WriteOutput( output, options ):
    """sort output."""

    if options.global_sort:
        output.sort()

    for o in output:
        options.stdout.write( "%s\t%s\t%s\t%s\t%f\t%s\t%f\t%5.2f\t%s\n" % o[1])

def WriteChanges( genome1, genome2, changed, options ):
    """write changed preference codons
    """

    output = []
    
    for aa, changes in changed.items():
        ## change will be from -> to
        changes.sort()

        if len(changes) != 2:
            raise "how is more than a change possible?"

        is_preferred1, pref1, codon1 = changes[0]
        is_preferred2, pref2, codon2 = changes[1]

        # check, if thought out correctly
        if is_preferred1 != False and is_preferred2 != True:
            raise "preferred/unpreferred mixed up."

        percent_difference = abs(pref2-pref1) / (pref1 + pref2) * 200.0

        code = "*" * (int(percent_difference) / 10)

        if options.sort == "percent-difference":
            output.append( (-percent_difference, (genome1, genome2, aa, codon1, pref2, codon2, pref1, percent_difference, code ) ) )
        elif options.sort == "aa":
            output.append( (aa, (genome1, genome2, aa, codon1, pref2, codon2, pref1, percent_difference, code ) ) )

        output.sort()
        
    return output

def WriteOverviewWeights( fields, table, options ):

    output = []
    WriteHeader(options)
    for x in range(1, len(fields)-1):
        for y in range(x+1, len(fields)):
            changed = {}

            for c in table:
                codon = c[0]
                w1 = c[x]
                w2 = c[y]
                t1 = w1 == 1.0 and w2 != 1.0
                t2 = w1 != 1.0 and w2 == 1.0
                if t1 or t2:
                    aa = Genomics.MapCodon2AA( codon )
                    if aa not in changed: changed[aa] = []
                    if t1:
                        changed[aa].append( (t1, w2, codon) )
                    else:
                        changed[aa].append( (t1, w1, codon) )

            output += WriteChanged( fields[x], fields[y], changes, options )            

    WriteOutput( output, options)
    
def WriteOverviewFrequencies( fields, table, options ):

    WriteHeader(options)
    output = []
    
    for x in range(1, len(fields)-1):
        for y in range(x+1, len(fields)):
            frequencies = {}

            ## collect frequencies per amino acid
            for c in table:
                codon = c[0]
                f1 = c[x]
                f2 = c[y]
                aa = Genomics.MapCodon2AA( codon )
                if aa not in frequencies: frequencies[aa] = []
                    
                frequencies[aa].append( (codon, f1, f2) )

            changed = {}
            
            ## sort for both genomes, and check if preference has changed
            for aa, codons in frequencies.items():
                codons.sort( lambda x, y: cmp(x[1], y[1]) )
                pref_codon1 = codons[-1]
                codons.sort( lambda x, y: cmp(x[2], y[2]) )                
                pref_codon2 = codons[-1]

                if pref_codon1 == pref_codon2:
                    continue
                else:
                    changed[aa] = [ (True, pref_codon1[2], pref_codon1[0]),
                                    (False, pref_codon2[1], pref_codon2[0]) ]
                    
            output += WriteChanges( fields[x], fields[y], changed, options )            

    WriteOutput( output, options )
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: codonbias_weights2tsv.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option( "--methods", dest="methods", type="string",
                      help="methods to apply.")

    parser.add_option( "--is-frequencies", dest="is_frequencies", action="store_true",
                      help="data is frequencies (default: weights).")

    parser.add_option( "-s", "--sort", dest="sort", type="choice",
                       choices=("percent-difference", "aa"),
                       help="sort order of output table.")

    parser.add_option( "-g", "--global-sort", dest="global_sort", action="store_true",
                       help="globally sort results (otherwise: by species pair).")

    parser.set_defaults( \
       methods = "",
       is_frequencies = False,
       sort = "percent-difference",
       global_sort= False,
       )

    (options, args) = E.Start( parser )
    if options.methods:
        options.methods = options.methods.split(",")

    fields, table = CSV.ReadTable(sys.stdin)

    ## convert weights to floats
    table = CSV.getConvertedTable( table, range( 1, len(fields) ) )

    for method in options.methods:

        if method == "overview":
            if options.is_frequencies:
                WriteOverviewFrequencies( fields, table, options )
            else:
                WriteOverviewWeights( fields, table, options )            
        

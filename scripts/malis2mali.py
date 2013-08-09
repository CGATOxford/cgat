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
malis2mali.py - concatenate multiple alignments
===============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

concatenate multiple alignments into a single multiple alignment. 

This script requires the following input:

1. Multiple alignments 

The multiple alignments can be given as separate files (option --pattern-mali contains
a '%%s') or in a single file.

2. A list of components

This file maps sequence identifiers in the multiple alignment file(s)
to multiple alignment identifiers. The format for this file is a tab separated table 
containing the following fields::

    1       sequence        a sequence identifier
    2       input_id        The id under which the alignment is found. This value is substituted
                            for the '%%s' in --pattern-mali.
    3       component_id    The component_id or multiple alignment identifier.
                            This field is optional. If no third column is specified, the input_id
                            is used.

Methods:

filter-variants
   only keep variant columns in the multiple alignment. Set ``--width`` to 3 to test for variant
   codons. Gaps will be ignored.

Usage
-----

Example::

   python malis2mali.py --help

Type::

   python malis2mali.py --help

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


import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Mali as Mali

## import shared helper functions from malis2malis
from malis2malis import *

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: malis2mali.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    addOptions( parser )

    parser.add_option( "--filename-coordinates", dest="filename_coordinates", type="string",
                      help="filename of coordinates that constitute the multiple alignment."  )

    parser.add_option( "--filename-identifiers", dest="filename_identifiers", type="string",
                      help="filename with list of identifiers to use."  )

    parser.add_option("-x", "--pattern-identifier", dest="pattern_identifier", type="string",
                      help="pattern to extract identifier from a sequence header."  )

    parser.add_option("-w", "--width", dest="width", type="int",
                      help="width of an alignment column (choose 3 for codon alignments) [default=%default]."  )

    parser.add_option("-m", "--method", dest="methods", type="choice",
                      choices = ("filter-variants", ),
                      help = "methods to apply" )

    parser.add_option("-p", "--parameters", dest="parameters", type="string",
                      help="parameter stack for methods that require one."  )

    parser.add_option("--mask-acgtn", dest="mask_actgn", action="store_true",
                      help="mask. Anything not [ACGTN] will be N."  )

    parser.set_defaults(
        pattern_identifier="(^\S+)",
        methods = [],
        parameters = "",
        filename_identifiers = None,
        filename_coordinates = None,
        mask_acgtn = False,
        )

    (options, args) = E.Start( parser )

    options.parameters = options.parameters.split(",")    

    if not options.pattern_mali:
        raise ValueError("Please specifiy a pattern to find the malis using --pattern-mali")

    ####################################################################
    ####################################################################
    ####################################################################
    ## Read components
    ####################################################################
    map_seq_id2component, map_component2seq_id, map_component2input_id = \
        readComponents( options )

    ####################################################################
    ####################################################################
    ####################################################################
    ## Read regions to mask
    ####################################################################
    map_component2masks = readMasks( options, map_component2input_id )

    ####################################################################
    ####################################################################
    ####################################################################
    ## Read regions to extract
    ####################################################################
    map_component2extracts = readExtracts( options, map_component2input_id )

    ####################################################################
    ####################################################################
    ####################################################################
    ## read identifiers
    ####################################################################
    if options.filename_identifiers:
        identifiers, nerrors = IOTools.ReadList( open(options.filename_identifiers, "r") )
        identifiers_set = set(identifiers)
    else:
        identifiers = None
        identifiers_set = None
        
    ####################################################################
    ####################################################################
    ####################################################################
    ## Prepare for run
    ####################################################################

    rx_identifier = re.compile(options.pattern_identifier)

    ## build list of concatenated malis
    sequences = {}
    if identifiers:
        for id in identifiers_set:
            sequences[id] = []
    else:
        identifiers_set = set()
        for seq_id in map_seq_id2component.keys():
            id = rx_identifier.search(seq_id).groups()[0]
            sequences[id] = []
            identifiers_set.add( id )
        identifiers = list(identifiers_set)
        identifiers.sort()
    
    component_ids = map_component2seq_id.keys()
    component_ids.sort()

    if options.test:
        component_ids = component_ids[:options.test]
        
    ####################################################################
    ####################################################################
    ####################################################################
    ## Build list of components to output.
    ####################################################################
    component_ids, map_sample2reference = selectComponents( component_ids, 
                                                            map_component2seq_id,
                                                            map_component2input_id,
                                                            None,
                                                            options )

    nskipped = 0
    new_component_ids = []

    for component_id in component_ids:

        try:
            mali = getMali( component_id, 
                            map_component2seq_id, 
                            map_component2input_id,
                            None,
                            options )
        except OSError, msg:
            E.warn("could not find mali %s: %s" % (component_id, msg))
            nskipped += 1
            continue

        ###############################################################
        ###############################################################
        ###############################################################
        ## check if all identifiers in component are present in mali
        ## and build a temporary alignment with all of those found
        component_set = set(map_component2seq_id[component_id])
        if len(component_set.difference( set(mali.getIdentifiers()))) != 0:
            nskipped += 1
            continue

        found = {}
        is_double = None
        temp_mali = Mali.Mali()
        temp_mali.setName( str(component_id) )

        for seq_id in map_component2seq_id[component_id]:
            id = rx_identifier.search(seq_id).groups()[0]
            if id not in identifiers_set:
                continue
            if id in found:

                if options.skip_doubles:
                    if options.loglevel >= 1:
                        options.stdlog.write("# component %s: removed double entry %s\n" % (component_id, seq_id ))
                    continue
                else:
                    is_double = id
                    break
                
            if options.output_format == "codeml":
                if len(mali[seq_id]) % 3 != 0:
                    raise "length of sequence %s is not a multiple of 3: %i" % (seq_id, len(mali[seq_id]))

            ## change identifier to id
            found[id] = True
            entry = mali.getEntry( seq_id ) 
            temp_mali.addSequence( id, entry.mFrom, entry.mTo, entry.mString )

        if is_double:
            nskipped += 1
            if options.loglevel >= 1:
                options.stdlog.write("# component %s: skipped because it contains double entry %s\n" % (component_id, is_double ))
            continue
        
        if set(found.keys()) != identifiers_set:
            nskipped += 1
            if options.loglevel >= 1:
                options.stdlog.write("# component %s: skipped because incomplete: %s\n" % (component_id, str(found.keys())))
            continue

        ###############################################################
        ###############################################################
        ###############################################################
        ## mask the temporary alignment
        maskAlignment( temp_mali, 
                       map_component2masks, 
                       map_component2extracts, 
                       map_sample2reference,
                       options )
        
        for id, o in temp_mali.items():
            if options.mask_acgtn:
                s = re.sub( "[^ACGTNactgn]", "N", o.mString)
            else:
                s = o.mString
            sequences[id].append( s )

        new_component_ids.append(component_id)

        ## if we only sample, stop if you have reached
        ## the desired number
        if options.sample and len(new_component_ids) == options.sample:
            break
        
    component_ids = new_component_ids
    
    nnucleotides = reduce( lambda x, y: x + y, map(len, sequences[identifiers[0]] ) )
    nseqs = len(sequences)
    ngenes = len(sequences)

    if options.output_format == "fasta":
        
        for id in identifiers:
            options.stdout.write( ">%s\n" % id )
            for seq in sequences[id]:
                options.stdout.write( "%s\n" % seq )                
            
    elif options.output_format == "codeml":

        ## Codeml output with G option. The format is
        ## "nseqs nnucleotides"
        ## "G ngenes genelengths in codons"

        options.stdout.write( "%i %i G\n" % (nseqs, nnucleotides))
        options.stdout.write( "G %i %s\n" % (ngenes, " ".join(map(lambda x: str(len(x) / 3), sequences[identifiers[0]]) ) ) )

        for id in identifiers:
            options.stdout.write( "%s\n" % id )
            for seq in sequences[id]:
                options.stdout.write( "%s\n" % seq )                

    elif options.output_format == "phylip":

        options.stdout.write( "%i %i\n" % (nseqs, nnucleotides))

        for id in identifiers:
            options.stdout.write( "%s\n" % id )
            for seq in sequences[id]:
                options.stdout.write( "%s\n" % seq )                
                
    if options.filename_coordinates:
            
        outfile = open(options.filename_coordinates, "w")
        outfile.write("component\tlength\tposition\n" )
        p = 0
        lsequences = map(lambda x: len(x), sequences[identifiers[0]])

        if len(lsequences) != len(component_ids):
            ll = []
            for id in identifiers:
                ll.append( len(map(lambda x: len(x), sequences[id]) ) )
            raise "%s: %s, lsequences=%i, ncomponents=%i" % (identifiers[0], str(ll), len(lsequences), len(component_ids))
        
        assert( len(lsequences) == len(component_ids) )
        
        for c in range(len(component_ids)):
            
            outfile.write( "%s\t%i\t%i\n" % (component_ids[c], lsequences[c], p) )
            p += lsequences[c]
            
        outfile.close()

    if options.loglevel >= 1:
        options.stdlog.write("# nsample=%i, nskipped=%i\n" % (len(component_ids), nskipped))

    E.Stop()
    

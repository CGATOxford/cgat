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
pdb2sequence.py - 
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

   python pdb2sequence.py --help

Type::

   python pdb2sequence.py --help

for command line help.

Documentation
-------------

Code
----

'''
USAGE="""extracts the sequence from a PDB file

python pdb2sequence.py [OPTIONS] < in.pdb 

Options:
-f, --format    output format (plain|fasta)
-s, --seqres    take sequence from seqres fields [default]
-a, --atoms     take sequence from atom coordinates
-c, --chain     only take specified chain
"""

import sys, string, os, getopt, re

from Pairsdb import *
import alignlib
import Experiment, Dots, Fasta

param_format    = "plain"
param_seqres    = 1
param_atoms     = None
param_chain     = None

param_filename_sequence = None

## translation tables for amino acids
MAP_MODELLER = { 'ALA' : 'A',
                 'ARG' : 'R',
                 'ASN' : 'N',
                 'ASP' : 'D',
                 'CYS' : 'C',
                 'CYX' : 'X',
                 'GLN' : 'Q',
                 'GLU' : 'E',
                 'GLY' : 'G',
                 'HIS' : 'H',
                 'HID' : 'H',     # histidine, delta H
                 'HIE' : 'H',     # histidine, epsilon H
                 'HIP' : 'H',     # histidine, protonated
                 'HSD' : 'H',     # histidine, delta H
                 'HSE' : 'H',     # histidine, epsilon H
                 'HSP' : 'H',     # histidine, protonated
                 'ILE' : 'I',
                 'LEU' : 'L',
                 'LYS' : 'K',
                 'MET' : 'M',
                 'NHE' : 'X',     # amine ending group
                 'NME' : 'X',     # N-methylamine ending group
                 'PHE' : 'F',
                 'PRO' : 'P',
                 'SER' : 'S',
                 'THR' : 'T',
                 'TRP' : 'W',
                 'TYR' : 'Y',
                 'VAL' : 'V',
                 'ACE' : 'X',     # acetyl beginning group
                 'NME' : 'X',
                 'UNK' : 'X',
                 'KCX' : 'X',
                 'MSE' : 'M',     # Selenomethionine
                 'KCX' : 'K',     # LYSINE NZ-CARBOXYLIC ACID, skip!
                 }

MAP_IUPAC = { 'ALA' : 'A',
              'ARG' : 'R',
              'ASN' : 'N',
              'ASP' : 'D',
              'CYS' : 'C',
              'GLN' : 'Q',
              'GLU' : 'E',
              'GLY' : 'G',
              'HIS' : 'H',
              'ILE' : 'I',
              'LEU' : 'L',
              'LYS' : 'K',
              'MET' : 'M',
              'PHE' : 'F',
              'PRO' : 'P',
              'SER' : 'S',
              'THR' : 'T',
              'TRP' : 'W',
              'TYR' : 'Y',
              'VAL' : 'V',
              }


param_map_amino_acids = MAP_MODELLER

def PrintAlignedSequences( sequence1, sequence2, chain = None, format="modeller" ):

    ## align sequences by identity

    seq_row = alignlib.makeSequence( sequence1 )
    seq_col = alignlib.makeSequence( sequence2 )
    alignator = alignlib.makeAlignatorFullDP( -0.0, -0.0 )
    map_row2col = alignlib.makeAlignataVector()
    alignator.Align( seq_row, seq_col, map_row2col )

    lines = string.split(alignlib.writePairAlignment( seq_row, seq_col, map_row2col ), "\n")

    if format == "modeller":
        
        first_res, sequence, last_res = string.split( lines[0], "\t" )
        
        print ">P1;structure"  
        print "structureX: %s : %s : %s : %s : %s : : : : " % ("structure", first_res, "" , last_res, "" )
        print "%s*" % sequence

        first_res, sequence, last_res = string.split( lines[1], "\t" )
        
        print ">P1;sequence"
        print "sequence:%s : %s : %s : %s : %s : : : : " % ("sequence" , first_res, "", last_res, "")
        print "%s*" % sequence
    else:
        print lines

def PrintSequence( sequence, chain = None, format="plain", alignment_sequence = None ):        

    if alignment_sequence:
        PrintAlignedSequences( sequence, alignment_sequence, chain, format )
    else:
        if format == "plain":
            print sequence
        elif format == "fasta":
            print ">%s %i\n%s" % (chain, len(sequence), sequence)

if __name__ == '__main__':
     
    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "V:f:c:ail:",
                                      ["Verbose=", "format=", "chain=", "atoms",
                                       "iupac", "align="])
    except getopt.error, msg:
        print USAGE
        sys.exit(2)


    for o,a in optlist:
        if o in ("-f", "--format"):
            param_format = a
        elif o == "--format":
            param_seqres = 1
            param_atoms = None
        elif o in ("-c", "--chain"):
            param_chain = a
        elif o in ("-a", "--atoms"):
            param_atoms = 1
            param_seqres = None
        elif o in ("-i", "--iupac"):            
            param_map_amino_acids = MAP_IUPAC
        elif o in ("-l", "-align"):
            infile = open(a, "r")
            infile.readline()
            alignment_sequence = ""
            for line in infile:
                alignment_sequence += re.sub( "[^ACDEFGHIKLMNPQRSTVWYZ]", "", line )
            
    ## you can also submit pdb_id-chain:
    if param_chain and len(param_chain) > 1:
        if string.find(param_chain, "-") != -1:
            param_chain = param_chain[-1]
        else:
            param_chain = None
            
    sequence = None
    last_chain = None
    last_residue_number = None
    
    for line in sys.stdin:

        if param_seqres and line[:6] == "SEQRES":
            
            chain = line[11]
            if chain != last_chain:
                if sequence:
                    if not param_chain or param_chain == last_chain:
                        PrintSequence( sequence, last_chain, param_format, alignment_sequence )
                sequence = ""
                
            last_chain = chain
            
            for residue in string.split(string.strip(line[19:-1]), " "):
                if not param_map_amino_acids.has_key( residue ):
                    print "# unknown residue %s" % residue
                else:
                    c = param_map_amino_acids[residue]
                    if c: sequence += param_map_amino_acids[residue]
                    
        if param_atoms and line[:6] in ("ATOM  ", "HETATM"):

            chain = line[21]
            residue_number = line[22:26]
            residue = line[17:20]
            
            if chain != last_chain:
                if sequence:
                    if not param_chain or param_chain == last_chain:
                        PrintSequence( sequence, last_chain, sequence )
                sequence = ""
                
            last_chain = chain

            if residue_number != last_residue_number:
                if param_map_amino_acids.has_key( residue ):
                    c = param_map_amino_acids[residue]
                    if c: sequence += param_map_amino_acids[residue]

            last_residue_number = residue_number
            
    if sequence:
        if not param_chain or param_chain == last_chain:        
            PrintSequence( sequence, last_chain, param_format, alignment_sequence )


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
MaliIO.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import string, sys, re
try: import alignlib_lite
except ImportError: pass

##---------------------------------------------------------------------------------------
def writeFasta( outfile,
                mali,
                line_width = 60,
                identifiers = None,
                skip_first = 0,
                gap_char = None
                ):
    
    """print multiple alignment in fasta format.
    """
    
    mali_width = mali.getWidth()
    mali_length = mali.getLength()
    
    for x in range(skip_first, mali_width):
        if identifiers:
            outfile.write(">" + str(identifiers[x]) + "\n")
        else:
            outfile.write(">%i\n" % x)
        if not gap_char:
            outfile.write(mali.getRow(x).getString() + "\n" )
        else:
            outfile.write(string.replace(mali.getRow(x).getString(), "-", gap_char) + "\n" )            

##---------------------------------------------------------------------------------------
def WriteMODELLER( outfile,
                    mali,
                    line_width = 60,
                    identifiers = None,
                    skip_first = 0,
                    gap_char = None
                    ):
    """print multiple alignment in PIR (MODELLER) format.
example:

>P1;5fd1
structureX:5fd1:1    : :106  : :ferredoxin:Azotobacter vinelandii: 1.90: 0.19
AFVVTDNCIKCKYTDCVEVCPVDCFYEGPNFLVIHPDECIDCALCEPECPAQAIFSEDEVPEDMQEFIQLNAELA
EVWPNITEKKDPLPDAEDWDGVKGKLQHLER*
    """
    
    mali_width = mali.getWidth()
    mali_length = mali.getLength()
    
    for x in range(skip_first, mali_width):
        
        if identifiers:
            outfile.write(">" + str(identifiers[x]) + "\n")
        else:
            outfile.write(">%i\n" % x)
            
        if not gap_char:
            outfile.write(mali.getRow(x).getString() + "\n" )
        else:
            outfile.write(string.replace(mali.getRow(x).getString(), "-", gap_char) + "\n" )            
    
              
##---------------------------------------------------------------------------------------
def writeClustalW( outfile, mali, line_width = 60, identifiers = None ):
    """print alignment in ClustalW format.
    (have to still look up the format, dsc can read it.)
    """

    mali_width = mali.getWidth()
    mali_length = mali.getLength()

    lines = []
    for x in range(0, mali_width):
        lines.append(mali.getRow(x).getString()) 

    outfile.write( "CLUSTAL W(1.60) multiple sequence alignment\n" )
    outfile.write( "\n" * 5 )

    for i in range(0, mali_length, line_width):
        i_to = i + line_width

        for x in range( 0, mali_width ):
            if identifiers:
                outfile.write("%-16s" % identifiers[x])
            else:
                outfile.write("%-16i" % x)

            outfile.write(lines[x][i:i_to] + "\n")
            
        outfile.write("\n" * 2)

##---------------------------------------------------------------------------------------
def readPicasso( infile ):
    """read alignment in the non-defined picasso format.
    """

    mali = alignlib_lite.py_makeMultipleAlignment()

    while 1:
        line = infile.readline()
        if not line: break

        x = re.search( "\d+\s+([A-Z\-\.]*)\s+\d+", line)
        if x:
            s = x.groups()[0]
            a = alignlib_lite.py_makeAlignatumFromString(s)
            a.thisown = 0
            mali.addAlignatum( a )
            
    return mali

##---------------------------------------------------------------------------------------
def compressAlignment( alignment, gap_character = "-", ignore_beginning = 0 ):
    """compress an alignment string.
    Lower-case characters at the beginning are ignored if so wished.
    --xxabBCDEfgHI
    becomes:
    -6+4

    This was necessary for radar output (e.g., 46497)
    """

    # subsitute all initial lower case characters with a gap
    if ignore_beginning:
        x = re.search("[A-Z]", alignment)
        if x:
            alignment = gap_character * x.start() + alignment[x.start():] 
        
    d = 0
    if alignment[0] == gap_character:
        gap = 1
    else:
        gap = 0
        
    result = ""
    
    for char in alignment:
        if char == gap_character and not gap:
            result += "+%i" % d
            d = 1
            gap = 1
            continue

        if char != gap_character and gap:
            result += "-%i" % d
            d = 1
            gap = 0
            continue

        d = d + 1
        
    if gap:
        result = result + "-%i" % d
    else:
        result = result + "+%i" % d
        
    return result

##---------------------------------------------------------------------------------------
def readFasta( infile, pattern_identifier="\S+" ):
    """read alignment in fasta format.
    """

    mali = {}
    id = None
    fragments = []
    identifiers = []
    for line in infile:
        if line[0] == "#": continue
        if line[0] == ">":
            if id: mali[id] = re.sub( "\s", "", string.join( fragments, ""))
            id = re.search( "^(%s)" % pattern_identifier, line[1:-1]).group(0)
            identifiers.append(id)
            fragments = []
            continue
        fragments.append( line[:-1] )

    if id: mali[id] = re.sub( "\s", "", string.join( fragments, ""))
    return convertGaps( mali ), identifiers

##---------------------------------------------------------------------------------------
def readPlain( infile ):
    
    mali = {}
    for line in infile:
        if line[0] == "#": continue
        data = line[:-1].split("\t")
        identifiers.append( data[3] )
        mali[data[3]] = data[2]
        
    return convertGaps( mali ), identifiers

##---------------------------------------------------------------------------------------
def readMali( infile, format = "fasta" ):

    if format == "fasta": return readFasta( infile )
    elif format == "plain": return readPlain( infile )
    else:
        raise ValueError, "unknown format %s" % format

##---------------------------------------------------------------------------------------
def convertGaps( mali, old_gap=".", new_gap = "-"):
    """convert gaps characters in mali.
    """
    for id in mali.keys():
        mali[id] = string.replace( mali[id], old_gap, new_gap,  )
    return mali

##---------------------------------------------------------------------------------------
def removeGappedColumns( mali, gap_char = "-" ):
    """remove all gapped columns in mali."""

    if len(mali) == 0: return mali
    
    keys = mali.keys()
    lmali = len(mali[keys[0]])
    wmali = len(keys)

    gaps_per_column = [0] * lmali

    for m in mali.values():
        for x in range(lmali):
            if m[x] == gap_char: gaps_per_column[x] += 1

    new_mali = {}
    
    for k in keys:
        s = []
        m = mali[k]
        for x in range(lmali):
            if gaps_per_column[x] < wmali:
                s += m[x]
        new_mali[k] = string.join(s, "")
        
    return new_mali

##---------------------------------------------------------------------------------------
def getSubset( mali, identifiers, not_in_set = False ):
    """return subset of mali which only contains identifiers."""

    new_mali = {}

    if not_in_set:
        for k,v in mali.items():
            if k not in identifiers: new_mali[k] = v
    else:
        for k,v in mali.items():
            if k in identifiers: new_mali[k] = v
        
        
    return new_mali

##---------------------------------------------------------------------------------------
def getFrameColumnsForMaster( mali, master, gap_char = "-" ):
    """get columns in frame according to master.
    """
    columns = []
    sequence = mali[master]
    x = 0
    while x < len(sequence):
        c = 0
        codon = []
        while x < len(sequence) and c < 3:
            if  sequence[x] != gap_char: #  and sequence[x] in string.uppercase:
                codon.append( x )
                c += 1
            x += 1
            
        if len(codon) == 3: columns.append( codon )

    return columns

##---------------------------------------------------------------------------------------
def getFrameColumnsForMasterPattern( mali, identifiers, master_pattern, gap_char = "-"):
    """get columns in frame for all masters matching the pattern.
    """
    columns = []
    for id in identifiers:
        if re.search( master_pattern, id ):
            columns += getFrameColumnsForMaster( mali, id )

    if len(columns) == 0:
        columns += getFrameColumnsForMaster( mali, identifiers[0] )

    # sort all columns by tuple. The "shortest" codon will be first (1,2,3) before (1,2,100)
    columns.sort()

    # select codons
    frame_columns = []
    last_codon = columns[0]
    for codon in columns[1:]:
        # skip identical codons
        if codon == last_codon: continue

        # take first (shortest) codon in case of identical first residue
        if codon[0] == last_codon[0]: continue

        # if not overlapping, keep
        if codon[0] > last_codon[2]:
            frame_columns.append( last_codon )

        # if overlapping, but out of register: skip
        last_codon = codon

    frame_columns.append( last_codon )
    
    return frame_columns

##---------------------------------------------------------------------------------------
def getMapFromMali( seq1, seq2, gap_char = "-" ):
    """build map of positions between mali."""
    xpos = 0
    ypos = 0

    map_a2b = alignlib_lite.py_makeAlignataVector()
    # build map between genomic sequences:
    for p in range(len(seq1)):

        if     seq1[p] != gap_char and \
               seq2[p] != gap_char and \
               seq1[p] in string.uppercase and \
               seq2[p] in string.uppercase:
            map_a2b.addPairExplicit( xpos + 1, ypos + 1, 0)
            
        if seq1[p] != gap_char:
            xpos += 1
        if seq2[p] != gap_char:
            ypos += 1
    return map_a2b


##---------------------------------------------------------------------------------------
def getCodonSequence( sequence, frame_columns, gap_char = "-", remove_stops = True ):
    """return a pruned sequence given frame columns.

    everything not in frame is deleted, only complete codons are kept.
    """
    
    fragments = []
    nstops, ncodons, naligned = 0, 0, 0
    for a,b,c in frame_columns:
        codon = sequence[a] + sequence[b] + sequence[c]

        codon_is_aligned = False
        codon_is_ok = True
        for x in codon:
            ## a codon will be masked, if it either
            ## 1. contains a gap character
            ## 2. is an unaligned character (lowercase)
            residue_is_unaligned = (x == gap_char) or \
                                   (x in string.lowercase)
            codon_is_aligned = codon_is_aligned or not residue_is_unaligned
            codon_is_ok = codon_is_ok and not residue_is_unaligned

        if codon_is_aligned: naligned += 1

        if codon_is_ok:
            ncodons += 1
            if string.upper(codon) in ("TAG", "TAA", "TGA"):
                if remove_stops:
                    fragments.append( gap_char * 3 )
                else:
                    fragments.append( codon )                        
                nstops += 1
            else:
                fragments.append( codon )                                            
        else:
            fragments.append( gap_char * 3 )

    return string.join(fragments, ""), naligned, ncodons, nstops

##---------------------------------------------------------------------------------------
def getPercentIdentity( seq1, seq2, gap_char = "-" ):
    """get number of identical residues between seq1 and seq2."""

    ntotal = 0
    nidentical = 0
    for a in range(len(seq1)):
        if seq1[a] != gap_char and seq2[a] != gap_char:
            ntotal += 1
            if seq1[a] == seq2[a]:
                nidentical += 1

    if ntotal == 0: return 0.0
    return float(nidentical) / ntotal
            
    
if __name__ == '__main__':

    c = pairsdblib.ConnectionPicasso()
    c.Connect()

    import pairsdblib
    mali = pairsdblib.makeMultipleAlignmentNeighbours( c, 54, "pairsdb_90x90")
    
    writeClustalW( sys.stdout, mali )

    c.Disconnect()

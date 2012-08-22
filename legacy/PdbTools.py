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
PdbTools.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, math, numpy

from ftplib import FTP

import Scientific.IO.PDB
import alignlib

CMD_SUPERIMPOSE="u3b"
CMD_ROTATE="rotate_pdb"

## PDB_SERVER = "ftp.ebi.ac.uk"
## PDB_ACCESS = "/pub/databases/msd/pdb_uncompressed/pdb%s.ent"
## USER = None
## PASSWORD = None
PDB_SERVER = "mozart.ebi.ac.uk"
PDB_ACCESS = "/ebi/ftp/pub/databases/msd/pdb_uncompressed/pdb%s.ent"
USER=""
PASSWORD=""
RASMOL_CMD = "rasmol"

PALETTE = ( "black", 
            "blue",             
            "cyan",
            "green",
            "greenblue",
            "magenta",
            "orange",
            "purple",  
            "red",
            "redorange",
            "violet",
            "white",
            "yellow")

AMINOACIDS = { 'ALA' : 'A',
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
               'HSD' : 'X',     # histidine, delta H
               'HSE' : 'X',     # histidine, epsilon H
               'HSP' : 'X',     # histidine, protonated
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
               'MSE' : 'X',     # Selenomethionine
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
    

def GetPdbFile( pdb_id, store_filename ):

    ftp = FTP(PDB_SERVER)       # connect to host, default port
    ftp.login()                 # user anonymous, passwd user@hostname
    
    # retrieve and save
    ftp.retrbinary('RETR ' + PDB_ACCESS % string.lower(pdb_id), open(store_filename, "wb").write)

def GetPdbFileLine( pdb_id ):

    ftp = FTP(PDB_SERVER)       # connect to host, default port
    ftp.login(USER,PASSWORD)                 # user anonymous, passwd user@hostname

    result = []
    # retrieve and save
    ftp.retrlines('RETR ' + PDB_ACCESS % string.lower(pdb_id), result.append)
    return result

class RasmolView:

    def __init__ (self, filename, pipe):
        
        self.mPdbFileName = filename       # filename to view

        self.mOutfile = pipe

        self.mOutfile.write("load %s\n" % self.mPdbFileName )
        self.mOutfile.write("""
        set background white
        renumber
        wireframe off
        ribbons off
        color grey
        color structure
        backbone on
        """)

    def WriteScript( self ):
        pass
        
    def RestrictChain(self, chain):
        self.mOutfile.write("restrict *:%s\n" % chain)

    def ColorRange(self, first_res, last_res, chain, color):
        self.mOutfile.write("select %i-%i:%s\ncolor %s\n" % (first_res, last_res, chain, color))

    def ColorResidue( self, residue, chain, color):
        self.mOutfile.write("select %i:%s\ncolor %s\n" % (residue, chain, color))

    def ColorSet( self, name, color):
        self.mOutfile.write("select %s\ncolor %s\n" % (name, color))

    def DefineSet( self, residues, name):
        
        self.mOutfile.write("define " + name + " ")
        step_size = 20
        i = 1
        while i < len(residues):
            self.mOutfile.write(string.join( map(str,residues[i:i+step_size]), ","))
            i += step_size
            if i < len(residues):
                self.mOutfile.write( ",\n")
        self.mOutfile.write("\n")
            
    def highlightResidues(self, first_res, last_res, chain, colour = None, shape = None):
        
        self.mOutfile.write("select %i-%i:%s\n" % (first_res, last_res, chain))
        if colour:
            self.mOutfile.write( "color %s\n" % (colour))
        if shape:
            if shape == "ball-and-stick":
                self.mOutfile.write( "spacefill %i\n" % (200) )
            elif shape == "cpk" or shape == "spacefill" :
                self.mOutfile.write( "spacefill\n" )
            
    def SelectSet(self, name ):
        self.mOutfile.write("select %s\n" % name)

    def Command( self, command ):
        self.mOutfile.write(command + "\n")
                
##print OUT "restrict *:$chain\n";
##print OUT "select all\ncolor grey\n";
##select $last-$next:$chain\ncolor $color[$index]\n";


class RasmolViewInline( RasmolView ):
    
    def __init__ (self, pdb_lines, pipe):
        
        self.mPdbLines    = pdb_lines

        self.mOutfile = pipe

        self.mOutfile.write("load inline" )
        self.mOutfile.write("""
        set background white
        set ambient 60
        select all
        renumber
        wireframe off
        ribbons off
        backbone off
        cartoon on
        color grey
        """)

    def WriteScript(self):

        self.mOutfile.write("""
        select all
        exit\n""")
        
        self.mOutfile.write( string.join(self.mPdbLines, "") + "\n")
        
        RasmolView.WriteScript( self )
        

def buildMapPdb2Sequence( sequence, filename_pdb, options, pdb_chain = ""):
    """build a map for residue numbers in pdb file to residue numbers on
    a sequence.

    returns the following maps:

    map_structure2seq: mapping of residue numbers between structure and
        sequence. These are mappings that will work if you "renumber" the
        structure.
        
    map_pdb2seq, map_seq2pdb: mapping according to residue numbers in pdb file.
    """

    if not os.path.exists( filename_pdb ):
        return None, None
    
    structure = Scientific.IO.PDB.Structure( filename_pdb )
    
    map_pdb2seq = {}
    map_seq2pdb = {}
    
    for chain in structure.peptide_chains:

        if chain.chain_id == pdb_chain:
            
            ## align pdb sequence to sequence
            map_structure2seq = alignlib.makeAlignataVector()
            alignator = alignlib.makeFullDP( -10.0, -2.0 )

            ## build sequence of pdb file
            structure = ""
            
            for residue in chain.sequence():
                structure += AMINOACIDS[residue]

            ## align reference sequence to sequence of pdb file
            row = alignlib.makeSequence( structure )
            col = alignlib.makeSequence( sequence )
            alignator.Align(row, col, map_structure2seq)

            if options.loglevel >= 3:
                options.stdlog.write( "structure: %s\n" % structure )                
                options.stdlog.write( "sequence : %s\n" % sequence )
                options.stdlog.write( "alignment of structure to sequence:\n" )
                options.stdlog.write( alignlib.writePairAlignment( row, col, map_structure2seq ) + "\n" )
                
            # print alignlib.writeAlignataTable(map_structure2seq)

            residue_number = 0
            
            for residue in chain.residues:

                residue_number += 1
                
                mapped_residue = map_structure2seq.mapRowToCol(residue_number)
                
                if not mapped_residue:
                    if options.loglevel >= 3:
                        options.stdlog.write( "# skipped residue %s=%s %i\n" % (str(residue.number), residue.name, residue_number))
                    continue

                r = str(residue.number)
                map_pdb2seq[r] = mapped_residue
                map_seq2pdb[mapped_residue] = r
                
            return map_structure2seq, map_pdb2seq, map_seq2pdb, residue_number-1, str(chain.residues[0].number), str(chain.residues[-1].number), structure

def TranslateAminoAcid( aa ):
    return AMINOACIDS[aa]

def GetPdbCoordinates( filename,
                       select_atom = ("CA",),
                       select_chain = None,
                       renumber = None,
                       only_coordinates = None):
    """read a pdb file and return coordinates of selected atoms
    """

    if not os.path.exists(filename):
        raise "pdb file %s does not exist" % filename

    if filename[-3:] == ".gz":
        lines = os.popen("gunzip < %s" % filename).readlines()
    else:
        lines = open(filename,"r").readlines()

    result = []

    current_number = 1
    
    for line in lines:
        if line[:6] not in ("ATOM  ", "HETATM"): continue

        chain   = line[21]
        number  = line[22:26]
        aa      = line[17:20]
        atom    = string.strip(line[13:17])
        
        x,y,z = map(string.atof, (line[30:38], line[38:46], line[46:54]))

        if select_chain and chain not in select_chain: continue
        if select_atom and atom not in select_atom: continue
        
        if renumber:
            number = current_number
            current_number += 1

        if AMINOACIDS.has_key(aa):
            aminoacid = AMINOACIDS[aa]
        else:
            sys.stderr.write( "# error in PdbCoordinates: aminoacid %s not known\n" % aa )
            continue

        if only_coordinates:
            result.append( (x, y, z) )
        else:
            result.append( (number, aminoacid, x, y, z) )            
        
    return result
        

def ConvertSequence2StructuralAlignment( src1, src2, source=None, format="plain", check_residues = 1):
    """calculate a structural alignment from two pdb files.
    """

    ca1 = GetPdbCoordinates( src1, renumber = 1)

    if len(ca1) == 0:
        raise "no coordinates found for %s" % src1

    ca2 = GetPdbCoordinates( src2, renumber = 1 )

    if len(ca2) == 0:
        raise "no coordinates found for %s" % src2

    if string.lower(format) not in ("plain",):
        raise "unknown alignment format %s" % format

    if source:
        lines = open(source, "r").readlines()
    else:
        lines = sys.stdin.readlines()

    ## replace gap characters
    lines = map(lambda x: re.sub( "\s", "", string.replace(x, ".", "-")), lines)
    if not lines:
        raise ValueError, "alignment is empty"

    lali = len(lines[0])

    current1 = 0
    current2 = 0

    index1 = 0
    index2 = 0

    output = []

    alignment = []

    for x in range(0, lali):

        res1 = lines[0][x]
        res2 = lines[1][x]

        if res1 != "-": current1+=1
        if res2 != "-": current2+=1

        try:
            while (ca1[index1][0] < current1): index1 += 1
            while (ca2[index2][0] < current2): index2 += 1                    
        except IndexError:
            break

        if res1 == "-" or res2 == "-":
            continue

        (i1, aa1, x1, y1, z1) = ca1[index1]
        (i2, aa2, x2, y2, z2) = ca2[index2]        

        if check_residues:
            if aa1 != res1:
                sys.stderr.write("# mismatch in 1:%s at residue alignment %i(%s) -> structure %i(%s)\n" %\
                                 (source, current1, res1, index1, aa1))
            if aa2 != res2:
                sys.stderr.write("# mismatch in 2:%s at residue %i(%s) -> %i(%s)\n" %\
                                 (source, current2, res2, index2, aa2))

        alignment.append( (x1, y1, z1, x2, y2, z2, 1) )

    return alignment


def GetSuperImposition( alignment ):
    """get superimposition based on aligned residue positions.
    """
    (infile, outfile) = os.popen2( CMD_SUPERIMPOSE )

    ## write input data
    infile.write( "%i\n" % len(alignment))
    for data in alignment:
        infile.write( "%f\t%f\t%f\t%f\t%f\t%f\t%i\n" % data)
    infile.close()

    ## get output data
    lines = outfile.readlines()

    rms = string.atof(re.split("\s+", lines[0] )[1])
    u = map(string.atof, re.split("\s+", lines[1])[1:-1] )
    t = map(string.atof, re.split("\s+", lines[2])[1:-1] )

    return (math.sqrt(rms / float(len(alignment))),
            numpy.array(t, numpy.float),
            numpy.reshape(numpy.array( u, numpy.float, 2), (3,3)))

def RotateCoordinates( coordinates, t, u ):
    """rotates a set of coordinates.
    """

    rotated_coordinates = []

    for x, y, z in coordinates:
        v = numpy.dot( numpy.array( (x, y, z), numpy.float), u) + t
        rotated_coordinates.append( (v[0], v[1], v[2]) )

    return rotated_coordinates
    

def GetAlignmentBetweenCorrespondingAtoms( coordinates1, coordinates2, cutoff ):
    """returns a list of atom positions, which are close to each other.

    This is done via a dynamic programming step. First all versus all comparison
    between atom positions is done. Only those positions are kept below cutoff.
    """

    dots = alignlib.makeAlignataMatrixRow()
    for i in range(len(coordinates1)):
        x1,y1,z1 = coordinates1[i]
        for j in range(len(coordinates2)):
            x2,y2,z2 = coordinates2[j]
            d = math.sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
            if d <= cutoff:
                dots.addPairExplicit(i+1, j+1, 1)
                
    seq1 = alignlib.makeSequence ("A" * len(coordinates1))
    seq2 = alignlib.makeSequence ("A" * len(coordinates2))    

    if dots.getLength() <= 3:
        return None
    
    dottor = alignlib.makeAlignatorDummy( dots )
    alignator = alignlib.makeAlignatorDotsSquared( 0, 0, dottor)
    map_a2b = alignlib.makeAlignataVector()
    
    alignator.Align( seq1, seq2, map_a2b)

    return map_a2b

##-----------------------------------------------------------------------
def IterativeSuperImposition( src1, src2,
                              initial_alignment, format = "plain",
                              cutoff = 4.0, max_iterations = 10 ):
    """superimposes two structures.
    """

    alignment = ConvertSequence2StructuralAlignment( src1, src2, initial_alignment, format)

    if len(alignment) == 0:
        raise "alignment is empty"

    rms = 10000

    iteration = 0

    ## get list of CA coordinates for structure1
    ca1 = GetPdbCoordinates( src1, renumber = 1, only_coordinates = 1)

    ## get list of CA coordinates for structure2 and rotate them
    ca2 = GetPdbCoordinates( src2, renumber = 1, only_coordinates = 1)

    last_alignment = None
    
    while alignment != last_alignment:
        
        # print "# aligning in iteration %i with %i pairs" % (iteration, len(alignment))
        
        rmsd, translation, rotation = GetSuperImposition( alignment )
    
##         print "# iteration=%i, atoms=%i, rmsd=%f, rotation=%s, translation=%s" % (iteration,
##                                                                                   len(alignment),
##                                                                                   rmsd,
##                                                                                   str(rotation),
##                                                                                   str(translation))

        ca1_rotated = RotateCoordinates( ca1, translation, rotation)
        
        last_alignment = alignment
        
        map_a2b = GetAlignmentBetweenCorrespondingAtoms( ca1_rotated, ca2, cutoff)

        if not map_a2b:
            raise ValueError, "no corresponding atoms within %f Angstroms found in iteration %i" % (cutoff, iteration)
            
        alignment = []
    
        for row in range( map_a2b.getRowFrom(), map_a2b.getRowTo()):
            col = map_a2b.mapRowToCol(row)
            if col:
                alignment.append( ca1[row-1] + ca2[col-1] + (1,))

        iteration += 1

        if iteration >= max_iterations:
            break

    return rmsd, len(alignment), translation, rotation, iteration


################################################################################
def MergePdbFiles( filenames,
                   rotations = None,
                   target = None,
                   chains = None,
                   do_gunzip = None):
    """Merge several pdb files.

    Rotates them if so desired.
    """
    
    current_chain = 0

    map_old2new = {}

    current_id = ord("A")
    index = 0
    coord = []

    for index in range(len(filenames)):
        
        coord.append("MODEL            %i\n" % (index))

        if not os.path.exists( filenames[index] ):
            raise "pdb file for %s not found" % pdb_id
        
        if do_gunzip:
            statement = "gunzip < %s" % filenames[index]
        else:
            statement = "cat %s" % filenames[index]

        if rotations:
            if rotations[index]:
                t = rotations[index][0]
                u = numpy.transpose(rotations[index][1])
                statement += " | %s %s %s" % (CMD_ROTATE,
                                           string.join(map(str, t), " "),
                                           string.join(map(str, numpy.reshape(u, (9,))), " "))

        lines = os.popen(statement).readlines()

        for line in lines:
            if line[:4] == "ATOM":
                chain = string.upper(line[21])
                if chains and chain in chains[index]: continue

                key = str(index) + "-" + chain
                
                if not map_old2new.has_key(key):
                    map_old2new[key] = chr(current_id)
                    current_id += 1

                new_chain = map_old2new[key]
                
                coord.append(line[:21] + new_chain + line[22:])

        coord.append("ENDMDL\n")

    result=["HEADER    PROTEIN STRUCTURE ALIGNMENT\n"]
        
    for key in map_old2new.keys():
        result.append("COMPND    (%s) %s\n" % (map_old2new[key], key))
        
    result += coord 

    if target:
        file = open(target,"w")
    else:
        file = sys.stdout

    file.write(string.join(result, "") + "\n")

    if target:
        file.close()
        




def RotatePdbFile( input_filename, output_filename,
                   translation, rotation,
                   do_gunzip = None):


    if do_gunzip:
        statement = "gunzip < %s" % input_filename
    else:
        statement = "cat %s" % output_filename

    u = numpy.transpose(rotation)
    statement += " | %s %s %s > %s" % (CMD_ROTATE,
                                       string.join(map(str, translation), " "),
                                       string.join(map(str, numpy.reshape(u, (9,))), " "),
                                       output_filename)

    return os.system(statement)
                

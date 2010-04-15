"""renumber a pdb file

python renumber_pdb_file.py [OPTIONS] < in > out

-s, --sequence: renumber according to sequence (fasta format)
-c, --compress: compress gaps due to missing residues
-g, --replace_glycine:  replace modified amino acids with glycine
-a, --chain: only take a specified chain
"""
import sys, os, getopt, re, string

import PdbTools

param_sequence = None
param_loglevel = 0
param_compress = None
param_glycine  = None
param_chain    = None

#  FORMAT(6A1,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,1X,I3)
#      Column             Content
#       1-6    'ATOM' or 'HETATM'
#       7-11   Atom serial number (may have gaps)
#      13-16   Atom name, in IUPAC standard format
#        17    Alternate location indicator indicated by A, B or C
#      18-20   Residue name, in IUPAC standard format
#      23-26   Residue sequence number (ordered as below)
#        27    Code for insertions of residues (i.e. 66A & 66B)
#      31-38   X coordinate
#      39-46   Y coordinate
#      47-54   Z coordinate
#      55-60   Occupancy
#      61-66   Temperature factor
#      68-70   Footnote number 

AMINOACIDS_TO_REPLACE = {
    'KCX' : 1,
    }

if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "s:V:c:ga:",
                                      ["sequence=","Verbose=", "compress",
                                       "replace_glycine", "chain="])
                                      
        
    except getopt.error, msg:
        print USAGE
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-s", "--sequence" ):
            param_sequence = a
        elif o in ("-V", "--Verbose"):
            param_loglevel = string.atoi(a)
        elif o in ("-g", "--replace_glycine"):
            param_glycine = 1
        elif o in ("-a", "--chain"):
            if len(a) > 1:
                if string.find(a, "-"):
                    a = string.split(a, "-")[1]
                elif len(a) == 5:
                    a = a[4]
            param_chain = a
            
    if param_sequence:
        param_sequence = re.sub("\s+", "", string.join(open(param_sequence, "r").readlines()[1:], ""))
        
    last_chain = None
    last_number = None
    current_number = 1
    current_atom_number = 1
    
    while 1:
        line = sys.stdin.readline()
        if not line: break

        if line[:6] not in ("ATOM  ", "HETATM"):
            if param_loglevel == 0:
                print line[:-1]
            continue

        header = line[:6]
        aminoacid = line[17:20]
        number=line[22:26]
        chain =line[21]
        atomtype = string.strip(line[13:16])
            
##         if param_loglevel >= 1:
##             print "%s %s %i %s" % (aa, number, current_number, param_sequence[current_number-1])
        
        if (last_chain != chain):
            current_number = 1
        else:
            if last_number != number:
                
                delta = string.atoi(number) - string.atoi(last_number)
                    
                if param_loglevel >= 2:
                    if delta > 1:
                        print "# gap between %s:%s and %s:%s %s" % (last_number, last_chain, number, chain, aminoacid )
                        
                if param_compress:
                    delta = 1

                if param_sequence:
                    aa = PdbTools.TranslateAminoAcid( aminoacid )        
                    if aa != param_sequence[current_number - 1]:
                        continue
                current_number += delta
                    
        last_number = number
        last_chain  = chain
        
        new_line = None

        if param_glycine:
            if aminoacid in AMINOACIDS_TO_REPLACE:
                if atomtype in ['N', 'CA', 'C', 'O']:
                    header = "ATOM  "
                    aminoacid = "GLY"
                else:
                    continue

                    
        new_line = "%-6s" % header + "%5i " % current_atom_number + " %-4s" % atomtype + "%-3s" % aminoacid + " %s" % chain +\
                   "%4i" % current_number + line[26:-1]
        current_atom_number += 1

        if param_chain and chain != param_chain:
            continue
        
        print new_line
        
    








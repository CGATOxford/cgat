USAGE = """filter out a list of residues from a
pdb file.
"""


import sys, re, os, string, getopt

import alignlib

param_master = 0
param_format = None
param_multiple_alignment = None

GAPCHARS=(".", "-")

if __name__ == '__main__':
     
    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "V:m:f:",
                                      ["Verbose=", "master=", 
                                       "format="])
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)


    for o,a in optlist:
        if o in ("-f", "--format"):
            param_format = a
        elif o in ("-i", "--iterations"):
            param_max_iterations = string.atoi(a)
        

    if len(args) != 1:
        print "please specify file with residue number."
        print USAGE
        sys.exit(1)

    residues = map(lambda x: string.split(x[:-1], "\t")[0], filter(lambda x: not re.match("#",x), open(args[0], "r").readlines()))
    
    residues_map = {}
    for r in residues:
        residues_map[r] = 1

    for line in sys.stdin:
        
        if line[:6] not in ("ATOM  ", "HETATM"): continue

        chain   = line[21]
        number  = string.strip(line[22:26])
        aa      = line[17:20]
        atom    = string.strip(line[13:17])

        if not residues_map.has_key( number ):
            continue

        print line[:-1]

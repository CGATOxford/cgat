"""append several pdb files.

python merge_pdb_files.py [OPTIONS] pdb_id1 pdb_id2 ...

pdb_id: four letter code. If a chain is given, only this
chain is taken, otherwise everything is output.

The chains are renamed in order of occurence.
Any other fields other than ATOM are removed.

OPTIONS:
-g, --gunzip:   files need to be decompressed
-r, --rotate:   translation vector + rotation matrix, separted by ,
"""
import sys
import os
import getopt
import re
import string

param_pattern_pdb = "/data/databases/pdb/pdb%s.ent.gz"
param_cmd_rotate  = "/data/bin/rotate_pdb"
param_gunzip = 0
param_rotate = None
param_target = None

def MergePdbFiles( pdb_ids,
                   filenames,
                   rotations = None,
                   target = None,
                   do_gunzip = None):

    current_chain = 0

    map_old2new = {}

    current_id = ord("A")
    index = 0
    coord = []

    for index in range(len(pdb_ids)):
        
        coord.append("MODEL            %i\n" % (index))

        pdb_id = pdb_ids[index]
        
        pdb = pdb_id[0:4]
        
        if len(pdb_id) > 4:
            pdb_chain = string.upper(pdb_id[4])
        else:
            pdb_chain = None

        if not os.path.exists( filenames[index] ):
            raise "pdb file for %s not found" % pdb_id
        
        if do_gunzip:
            statement = "gunzip < %s" % filenames[index]
        else:
            statement = "cat %s" % filenames[index]

        if rotations:
            if rotations[index]:
                statement += " | %s %s" % (param_cmd_rotate, rotations[index])

        lines = os.popen(statement).readlines()

        for line in lines:
            if line[:4] == "ATOM":
                chain = string.upper(line[21])
                if pdb_chain and pdb_chain != chain: continue

                key = pdb_id + chain
                
                if not map_old2new.has_key(key):
                    map_old2new[key] = chr(current_id)
                    current_id += 1

                new_chain = map_old2new[key]
                
                coord.append(line[:21] + new_chain + line[22:])

        coord.append("ENDMDL\n")

    result=["HEADER    PROTEIN STRUCTURE ALIGNMENT"]
        
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
        
                
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "r:t:gf:",
                                      ["rotate=", "target=", "gunzip", "filenames="])
        
    except getopt.error, msg:
        print USAGE
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-g", "--gunzip" ):
            param_gunzip = 1
        elif o in ( "-t", "--target" ):
            param_target = a
        elif o in ("-r", "--rotate"):
            param_rotate = a
        elif o in ("-f", "--filenames"):
            param_filenames = string.split(a, ",")

    pdb_ids = []
    pdb_files = []

    for a in args:
        pdb_ids.append(a)
        pdb_files.append(param_pattern_pdb % a[:4] )

    if param_rotate:
        lines = map(string.strip, string.split(param_rotate, "@"))
        rotations = []
        for line in lines:
            data = re.split("\s*,\s*", line[:-1])
            ## empty separator will give ''
            if len(data) == 1:
                rotations.append(None)
            else:
                if len(data) == 12:
                    rotations.append(string.join(data, " "))
                else:
                    raise "incomplete rotation matrix!"


        if len(rotations) != len(pdb_ids):
            raise "number of pdbids and rotations differ!"
    else:
        rotations = None

    if param_filenames:
        pdb_files = param_filenames

    MergePdbFiles( pdb_ids, pdb_files,
                   rotations,
                   target = param_target,
                   do_gunzip = param_gunzip)
    
    




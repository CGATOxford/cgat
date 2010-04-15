USAGE="""superimposes two pdb structures

python pdb2superimpose.py [OPTIONS] filename1 filename2 alignment

Options:
-i, --iterations #       do iterative superimposition, number specifies
                        maximum number of iterations
-c, --cutoff            proximity cutoff for iterative superimposition
                         

"""

import sys, string, os, getopt

import Experiment, PdbTools, numpy

param_loglevel   = 1
param_format     = "plain"
param_cutoff     = 3.0
param_max_iterations = 0
param_write_pdb  = None
param_write_rotated = None

if __name__ == '__main__':
     
    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "V:i:c:p:r:",
                                      ["Verbose=", "cutoff", "iterations", "write_pdb",
                                       "write_rotated"])
                                       
    except getopt.error, msg:
        print USAGE
        sys.exit(2)


    for o,a in optlist:
        if o in ("-f", "--format"):
            param_format = a
        elif o in ("-i", "--iterations"):
            param_max_iterations = string.atoi(a)
        elif o in ("-c", "--cutff"):
            param_cutoff = string.atof( cutoff )
        elif o in ("-p", "--write_pdb"):
            param_write_pdb = a
        elif o in ("-p", "--write_rotated"):
            param_write_rotated = a

    if len(args) < 3:
        print USAGE
        print "not enough arguments specified"
        sys.exit(1)

    param_filename_pdb1 = args[0]
    param_filename_pdb2 = args[1]    
    param_filename_alignment = args[2]

    print Experiment.GetHeader()
    print Experiment.GetParams()    

    try:
        rmsd, natoms, translation, rotation, iterations = PdbTools.IterativeSuperImposition( param_filename_pdb1,
                                                                                             param_filename_pdb2,
                                                                                             param_filename_alignment,
                                                                                             max_iterations = param_max_iterations,
                                                                                             cutoff = param_cutoff)
    except ValueError, msg:
        print "# ERROR: %s" % msg
        print "rmsd\tna" 
        print "natoms\tna" 
        print "iterations\tna" 
        sys.exit(1)

    print "rmsd\t%5.2f" % rmsd
    print "natoms\t%i" % natoms
    print "iterations\t%i" % iterations
    print translation
    print rotation
    
    if param_write_pdb:
        PdbTools.MergePdbFiles( (param_filename_pdb1, param_filename_pdb2),
                                rotations=( ( translation, rotation ), None),
                                target = param_write_pdb )

    if param_write_rotated:
        PdbTools.RotatePdbFile( param_filename_pdb1,
                                param_write_rotated,
                                translation,
                                rotation)                                






import os, sys, string, re, getopt, math

USAGE="""python %s [OPTIONS] < mali > filtered

Extract sequences from a multiple alignment

Version = $Id: mali_extract.py 2782 2009-09-10 11:40:29Z andreas $

Options:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-s, --subset=                   subset of ids to select
-c, --components                filename with components to be analyses separately in the multiple alignment
""" % sys.argv[0]

param_long_options=["verbose=", "help", "subset=", "components=" ]

param_short_options="v:ho:s:c::"

param_loglevel = 1

param_gap_char = "-"
param_mask_char = "x"

param_subset = None
param_filename_components = None

import Experiment
import MaliIO
            
##------------------------------------------------------------
if __name__ == '__main__':

    try:
        optlist, args = getopt.getopt(sys.argv[1:], param_short_options, param_long_options)
    except getopt.error, msg:
        print USAGE, msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-v", "--verbose" ):
            param_loglevel = int(a)
        elif o in ( "-h", "--help" ):
            print USAGE
            sys.exit(0)
        elif o in ("-s", "--subset"):
            param_subset = a
        elif o == ("-c", "--components"):
            param_filename_components = a

    if param_loglevel >= 1:
        print Experiment.GetHeader()
        print Experiment.GetParams()
            
    ## 1. read multiple alignment in fasta format
    all_mali, all_identifiers = MaliIO.readFasta( sys.stdin )

    if len(all_identifiers) == 0:
        raise "alignment is empty."

    if param_loglevel >= 1:
        print "# read mali with %i entries." % len(all_identifiers)

    if param_filename_components:

        infile = open(param_filename_components, "r")
        components = {}
        for line in infile:
            if line[0] == "#": continue
            if line[0] == ">": continue
            a, b = line[:-1].split("\t")[:2]
            if b not in components:
                components[b] = []
            components[b].append(a)
        
        if param_loglevel >= 1:
            print "# read %i components." % len(components)

    if param_subset:
        components = { 'all' : string.split( param_subset, "," ) }

    for key, identifiers in components.items():
        # 1. remove gaps in multiple alignment
        mali = MaliIO.removeGappedColumns( MaliIO.getSubset( all_mali, identifiers),
                                           param_gap_char )

        for i in identifiers:
            print ">%s\n%s\n" % (i, mali[i] )

    if param_loglevel >= 1:        
        print Experiment.GetFooter()
    
    

#!/sw/arch/bin/python
####
####
##
## Project PairsDBTools
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: pdb_domains2molscript.py 2782 2009-09-10 11:40:29Z andreas $
##
##
####
####


#----------------------------------------------------------------

USAGE="""
python pdb_domains2molscript.py -s|--script=rasmol.script -p|--pdb=file.pdb < molscript.in > molscript.out

adds to a molscript script the domain definitions in
a rasmol script.

1. run view_pdb_domains and save as in.script
3. view script and save pdb and molscript file
        write molscript in.mol
        save in.pdb
4. run command above:
python pdb_domains2molscript.py --script=in.script --pdb=in.pdb < in.mol > out.mol

"""
import sys, string, re, getopt

param_rasmol_script = None
param_pdb_filename  = None

if __name__ == '__main__':

    
    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "s:p:",
                                      ["script=",
                                       "pdb=",])
        
    except getopt.error, msg:
        print USAGE
        print msg
        sys.exit(2)

    for o,a in optlist:
        if o in ( "-s", "--script" ):
            param_rasmol_script = a
        elif o in ( "-p", "--pdb" ):
            param_pdb_filename = a


    if not param_rasmol_script:
        print USAGE
        print "please define --script="
        sys.exit(1)

    if not param_pdb_filename:
        print USAGE
        print "please define --pdb="
        sys.exit(1)

    
        
    infile = open(param_rasmol_script, "r")

    domains = []
    pdb_lines = []
    while 1:
        line = infile.readline()
        
        if not line:
            break

        x =  re.search("echo Domain (\d+) \(did=(\S+), (\d+)-(\d+):(\S+);.*= color (\S+)", line)
        
        if x:
            domains.append( x.groups() )
            
    infile.close()

    while 1:
        line = sys.stdin.readline()
        if not line: break

        if re.search("read mol \"inline\";", line):
            print '   read mol %s;' % param_pdb_filename
            print '   set colourparts on;'
            for domain in domains:
                (n, did, first_res, last_res, chain, color) = domain
                print "set residuecolour from %s%s to %s%s %s;" % (chain, first_res, chain, last_res, color)
            continue
        
        if not re.match("end_plot", line):
            print line,
            continue

        print line
        break
    
    
    
            











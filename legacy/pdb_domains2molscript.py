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
pdb_domains2molscript.py - 
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

   python pdb_domains2molscript.py --help

Type::

   python pdb_domains2molscript.py --help

for command line help.

Documentation
-------------

Code
----

'''
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
    
    
    
            











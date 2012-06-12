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
WrapperDialign.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile, subprocess

"""Wrapper for Dialign.
"""

import Genomics

class Dialign:

    mOptions = ""
    mExecutable = "dialign -stdo "
    mEnvironment = "DIALIGN2_DIR=/net/cpp-group/legacy/bin/dialign2_dir; export DIALIGN2_DIR;"
    
    def __init__( self, options ):
        self.mOptions = options

    def Align( self, s1, s2, result ):

        result.clear()
        
        handle_tmpfile, filename_tmpfile = tempfile.mkstemp()
        os.write( handle_tmpfile, ">s1\n%s\n" % (s1))
        os.write( handle_tmpfile, ">s2\n%s\n" % (s2))
        os.close( handle_tmpfile )

        statement = string.join( ( "(", self.mEnvironment, 
                                   self.mExecutable, self.mOptions, filename_tmpfile, ")" ),
                                 " ")

        p = subprocess.Popen( statement , 
                              shell=True, 
                              stdin=subprocess.PIPE, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              close_fds=True)

        (file_stdout, file_stdin, file_stderr) = (p.stdin, p.stdout, p.stderr)
        file_stdin.close()
        lines = file_stdout.readlines()
        lines_stderr = file_stderr.readlines()
        exit_code = file_stdout.close()
        file_stderr.close()
        
        if exit_code:
            raise "Error while executing statement %s" % statement

        r = None
        for x in range(len(lines)):
            if re.search( "Alignment \(FASTA format\):", lines[x]):
                r = Genomics.ParseFasta2Hash( lines[x+2:])
                break

        if not r: return None
        
        a1 = r['s1']
        a2 = r['s2']

        x1 = 1
        x2 = 1
        for pos in range(len(a1)):
            if a1[pos] in string.uppercase and a2[pos] in string.uppercase:
                result.addPairExplicit( x1, x2, 0 )
                x1 += 1
                x2 += 1
                continue
            
            if a1[pos] != "-": x1 += 1
            if a2[pos] != "-": x2 += 1            


        os.remove( filename_tmpfile )
        
        return result
                          
            
        
                                 
        
        
        
        
    
        

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
WrapperDBA.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile, subprocess

"""Wrapper for DBA.
"""

import Genomics

class DBA:

    mOptions = ""
    mExecutable = "dba -quiet -nomatchn -align "
    mEnvironment = "WISECONFIGDIR=/net/cpp-group/src/wise2.2.0/wisecfg; export WISECONFIGDIR;"
    
    def __init__( self, options = ""):
        self.mOptions = options

    def Align( self, s1, s2, result ):

        result.clear()
        
        handle_tmpfile1, filename_tmpfile1 = tempfile.mkstemp()
        handle_tmpfile2, filename_tmpfile2 = tempfile.mkstemp()
        os.write( handle_tmpfile1, ">s1\n%s\n" % (s1))
        os.close( handle_tmpfile1 )
        os.write( handle_tmpfile2, ">s2\n%s\n" % (s2))
        os.close( handle_tmpfile2 )

        statement = string.join( ( "(",
                                   self.mEnvironment, 
                                   self.mExecutable,
                                   self.mOptions,
                                   filename_tmpfile1,
                                   filename_tmpfile2,
                                   ")" ),
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

        for line in lines:
            data = re.split("\s+", line[:-1] )
            if len(data) != 6: continue
            x1 = int(data[1]) 
            x2 = int(data[4])
            result.addPair( x1, x2, 0 )

        os.remove( filename_tmpfile1 )
        os.remove( filename_tmpfile2 )
        
        return result
                          
            
        
                                 
        
        
        
        
    
        

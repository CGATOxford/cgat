################################################################################
#   Gene prediction pipeline 
#
#   $Id: WrapperGblocks.py 1799 2008-03-28 11:44:19Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
import os, sys, string, re, tempfile, popen2

"""Wrapper for Gblocks
"""

import Genomics
import alignlib

class Gblocks:

    mOptions = ""
    mExecutable = "Gblocks %s -t d"
    mEnvironment = ""
    
    def __init__( self, options = ""):
        self.mOptions = options

    def GetBlocks( self, s1, s2 ):
        """the strings have to be already aligned!!!"""
        
        handle_tmpfile, filename_tmpfile = tempfile.mkstemp()
        os.write( handle_tmpfile, ">s1\n%s\n" % (s1))
        os.write( handle_tmpfile, ">s2\n%s\n" % (s2))
        os.close( handle_tmpfile )

        statement = string.join( ( "(", self.mEnvironment, 
                                   self.mExecutable % filename_tmpfile,
                                   self.mOptions, ")" ),
                                 " ")

        (file_stdout, file_stdin, file_stderr) = popen2.popen3( statement )
        file_stdin.close()
        lines = file_stdout.readlines()
        lines_stderr = file_stderr.readlines()
        exit_code = file_stdout.close()
        file_stderr.close()
        if exit_code:
            raise "Error while executing statement %s" % statement

        if not os.path.exists( filename_tmpfile + "-gb"):
            os.remove( filename_tmpfile )            
            return "", ""
        
        lines = open( filename_tmpfile + "-gb").readlines()
        r = Genomics.ParseFasta2Hash( lines)

        if not r: return "", ""

        os.remove( filename_tmpfile )
        os.remove( filename_tmpfile + "-gb" )
        os.remove( filename_tmpfile + "-gb.htm")        
        
        return r['s1'], r['s2']


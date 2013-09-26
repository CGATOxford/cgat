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
WrapperClustal.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile

try:
    from Bio.Clustalw import MultipleAlignCL
    from Bio.Clustalw import do_alignment
except ImportError:
    # code needs updating for newer biopython versions
    pass
"""Wrapper for Clustal based on Biopython
"""

import Genomics

class Clustal:

    mOptions = ""
    
    def __init__( self, options = ""):
        self.mOptions = options

    def Align( self, s1, s2, result ):

        result.clear()

        handle_tmpfile1, filename_tmpfile1 = tempfile.mkstemp()
        handle_tmpfile2, filename_tmpfile2 = tempfile.mkstemp()        
        os.write( handle_tmpfile1, ">s1\n%s\n" % (s1))
        os.write( handle_tmpfile1, ">s2\n%s\n" % (s2))
        os.close( handle_tmpfile1 )
        os.close( handle_tmpfile2 )
        cline = MultipleAlignCL( filename_tmpfile1 )
        cline.set_output( filename_tmpfile2 )
                  
        align = do_alignment(cline)
                  
        seqs = align.get_all_seqs()
        if len(seqs) != 2: return result
        
        a1 = seqs[0].seq.tostring()
        a2 = seqs[1].seq.tostring()

        x1 = 0
        x2 = 0
        for pos in range(len(a1)):
            
            if a1[pos] not in "Nn-" and a2[pos] not in "Nn-":
                result.addPair( x1, x2, 0 )
                x1 += 1
                x2 += 1
                continue
            
            if a1[pos] != "-": x1 += 1
            if a2[pos] != "-": x2 += 1            

            
        os.remove( filename_tmpfile1 )
        os.remove( filename_tmpfile2 )
            
        return result


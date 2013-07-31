################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
"""
Fasta.py - Methods for dealing with fasta files.
================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

"""

import string, os, sys, re

class Fasta:

    def __init__ (self,file):
        self.mFile = file
        self.mLastLine = " "
        
    def FetchOne( self ):
        """returns a tuple (description, sequence).

        Returns (None, None) if no more data is there.
        """
        while self.mLastLine != None:
            if self.mLastLine[0] == ">": break
            self.mLastLine = self.mFile.readline()
        else:
            return (None, None)
        
        description = self.mLastLine[1:-1]
        sequence = ""
        while 1:
            line = self.mFile.readline()
            if not line: break
            if line[0] == ">":
                self.mLastLine = line
                break

            sequence += line[:-1]

        return (description, re.sub("\s", "", sequence))
            
            
    

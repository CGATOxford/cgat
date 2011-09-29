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
pipeline_vitaminD_motifs.py - 
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

   python pipeline_vitaminD_motifs.py --help

Type::

   python pipeline_vitaminD_motifs.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys, re, os, tempfile, collections

import Experiment as E
import Pipeline as P
import sqlite3

PARAMS = P.getParameters()

def filterMotifsFromMEME( infile, outfile, selected ):
    '''select motifs from a MEME file and save into outfile
    '''

    outs = open(outfile, "w" )
    if len(selected) == 0:
        outs.close()
        return

    keep = True

    for line in open(infile, "r"):
        if line.startswith( "MOTIF" ):
            motif = re.match( "MOTIF\s+(\d+)", line ).groups()[0]
            if motif in selected:
                keep = True
            else:
                keep = False
        if keep: outs.write( line )
    outs.close()

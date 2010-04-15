################################################################################
#   Gene prediction pipeline 
#
#   $Id: rotate_table.py 2782 2009-09-10 11:40:29Z andreas $
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
import os, sys, string, re, tempfile, subprocess, optparse, time, math

"""rotate a table
"""

import Experiment

parser = optparse.OptionParser( version = "%prog version: $Id: rotate_table.py 2782 2009-09-10 11:40:29Z andreas $")

if __name__ == "__main__":

    (options, args) = Experiment.Start( parser )
        
    data = []
    
    data = map( lambda x: x[:-1].split("\t"), filter( lambda x: x[0] != "#", sys.stdin.readlines()))
    
    for x in zip(*data):
        print "\t".join( x )

    Experiment.Stop()
        
            

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
'''
cgat.py - frontend for cgat scripts
===================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

To get help for a specific command, type::

    cgat <command> --help

'''

import os
import sys
import glob
import imp

def main():

    argv = sys.argv

    path = os.path.abspath( os.path.dirname(__file__) )

    if argv[1] == "--help" or argv[1] == "-h":
        print(globals()["__doc__"])
        print("The list of available commands is:\n" )
        print( "\n".join( sorted([os.path.basename(x)[:-3] for x in glob.glob( os.path.join( path, "*.py") )]) ) )
        return

    command = argv[1]

    (file, pathname, description ) = imp.find_module( command, [path,] )
    module = imp.load_module( command, file, pathname, description)
    module.main( sys.argv )
    
if __name__ == "__main__":
    sys.exit( main() )

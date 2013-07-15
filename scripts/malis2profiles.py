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
malis2profiles.py - build profiles from malis
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

convert a set of plain alignments to profiles.

Usage
-----

Example::

   python malis2profiles.py --help

Type::

   python malis2profiles.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import optparse
import math
import time
import random
import types

USAGE="""python %s [OPTIONS]



""" % sys.argv[0]

import CGAT.Experiment as E
import CGAT.Mali as Mali

##------------------------------------------------------------
if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: malis2profiles.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.set_defaults(
        )

    (options, args) = E.Start( parser )
    
    mali = Mali.SequenceCollection()
    last_id = None
    ninput, noutput, nskipped = 0, 0, 0

    for line in sys.stdin:
        if line[0] == "#": continue

        start, ali, end, id = line[:-1].split("\t")
        ninput += 1
        if id != last_id:
            if last_id:
                mali.setName(last_id)
                mali.writeToFile( sys.stdout, format="profile" )
                noutput += 1
            mali = Mali.SequenceCollection()
            last_id = id

        mali.addSequence( id, start, end, ali )
    
    if last_id:
        mali.setName(last_id)
        mali.writeToFile( sys.stdout, format="profile" )
        noutput += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nskipped=%i.\n" % (ninput, noutput, nskipped))
        
    E.Stop()
    

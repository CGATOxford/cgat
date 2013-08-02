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
convert_time2seconds.py -
=============================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

convert output of the time command into seconds:

input: for example something like

1h33m47.32s

Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

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

import CGAT.Experiment as E

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    for line in sys.stdin:

        if line[0] == "#": continue

        data = string.split(line[:-1], "\t")

        new_fields = []

        for field in data:

            x = re.match("^(\d+)h(\d+)m([\d.]+)s", field)
            if not x:
                x = re.match("^(\d+)m([\d.]+)s", field)
                if x:
                    minutes, seconds = x.groups()
                else:
                    new_fields.append(field)
                    continue
                hours = "0"
            else:
                hours, minutes, seconds = x.groups()


            new_fields.append( "%.2f" % (string.atof(hours) * 3600.0 + string.atof(minutes) * 60.0 + string.atof(seconds)))

        print string.join(new_fields, "\t")

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )



    



    


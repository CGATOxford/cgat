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
nofarm.py - stand-in for ``farm.py``
====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script is a stand-in for :doc:`farm`. Instead of running
splitting the input and running jobs on the cluster, the command
is executed as is.

This script accepts all the options of ``farm.py``, but ignores
most of them.

Usage
-----

Example::

   python nofarm.py --help

Type::

   python nofarm.py --help

for command line help.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import string
import optparse
import time
import glob
import subprocess
import tempfile
import shutil

import farm

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import threadpool

##--------------------------------------------------------------------
def main():

    parser = farm.getOptionParser()

    (options, args) = E.Start( parser, 
                                        add_cluster_options = True )


    cmd = args[0]
    if len(args) > 1:
        cmd += " '" + "' '".join(args[1:]) + "'" 

    cmd = re.sub( "%DIR%", "", cmd )
    retcode = subprocess.call( cmd,
                               shell = True,
                               stdin = sys.stdin,
                               stdout = sys.stdout,
                               cwd = os.getcwd(),
                               close_fds = True)                              
    E.Stop()

if __name__ == '__main__':
    main()

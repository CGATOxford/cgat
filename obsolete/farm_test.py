##########################################################################
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
##########################################################################
'''
farm_test.py - tests for :ref:`farm.py`
=======================================

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

   python farm_test.py --help

Type::

   python farm_test.py --help

for command line help.

Documentation
-------------

Code
----

'''
import random
import os
import time


if __name__ == "__main__":
    queue = "test.q"

    datafile = "tmp.data"
    scriptfile = "tmp.scriptfile"
    resultfile = "tmp.out"
    logfile = "tmp.log"

    num_jobs = 100
    num_parts = 20
    max_wait = 4

    outfile = open(datafile, "w")
    outfile.write("id\twait\n")
    for x in range(0, num_parts):
        outfile.write("%i\t%i\n" % (x, random.randint(0, max_wait)))
    outfile.close()

    outfile = open(scriptfile, "w")
    outfile.write("""import sys,time
    logfile=open("my.log", "w")
    print"id\twaited"

    for line in sys.stdin:
       if line.startswith("id"): continue
       id, wait = map(int, line[:-1].split("\t"))
       time.sleep(wait)
       print "%i\twaited %i" % (id, wait)
       logfile.write("logging for %i\\n"% id)

    logfile.close()
    """)

    outfile.close()

    cmd = "farm.py --cluster-priority=-10 --cluster-queue=%s --cluster-num-jobs=%i --split-at-column=1 --input-header --output-header --log=%s python %s --log=my.log < %s > %s" %\
        (queue, num_jobs, logfile, scriptfile, datafile, resultfile)
    print cmd
    os.system(cmd)

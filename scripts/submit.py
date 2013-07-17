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
submit.py - run a collection of qsub scripts on the cluster
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script will submit a collection of qsub scripts on the cluster.
The script will only submit those jobs for which a result is not 
present.

Usage
-----

Example::

   python submit.py *.qsub

Type::

   python submit.py --help

for command line help.

Documentation
-------------

Code
----

'''
USAGE="""python submit.py [OPTIONS] 

submit a list of jobs to the cluster.
"""

import os
import sys
import re
import string
import optparse
import time
import glob
import subprocess
import stat

import CGAT.Experiment as E
import CGAT.IOTools as IOTools

def checkPythonRuns( filename ):
    """returns true if a python run is complete."""
    last_line = IOTools.getLastLine( filename )
    return re.match( "# job finished", last_line )
    
def isNewer( a, b ):
    """return true if file a is newer than file b."""

    ## get times of most recent access
    at = os.stat( a )[7]
    bt = os.stat( b )[7]

    return at > bt

if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: submit.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"] )

    parser.add_option( "-g", "--glob", dest="glob_pattern", type="string" ,
                       help="glob pattern to use for collecting cluster jobs descriptions [%default]")

    parser.add_option( "-i", "--input-pattern", dest="input_pattern", type="string" ,
                       help="regular expression to extract job id from filename [%default].")

    parser.add_option( "-o", "--output-pattern", dest="output_pattern", type="string" ,
                       help="string to convert a job id to a filename [%default].")


    parser.set_defaults( glob_pattern = "job*.qsub",
                         input_pattern = "(\S+).qsub",
                         output_pattern = "%s.stdout",
                         remove_old = True,
                         force = False,
                         check_completeness = "python",
                         )

    (options, args) = E.Start( parser, 
                                        add_pipe_options = True )

    if args:
        filenames = args
    elif options.glob_pattern:
        filenames = glob.glob( options.glob_pattern )
        
    ninput, nrun, nskipped, nerrors = 0, 0, 0, 0
    ndeleted = 0

    if options.check_completeness == "python":
        isComplete = checkPythonRuns

    ##############################################################
    ##############################################################
    ##############################################################
    ## decide what to do
    ##############################################################
    jobs = []
    files_to_delete = []

    for filename in filenames:
        
        ninput += 1
        try:
            job_name = re.search( options.input_pattern, filename ).groups()[0]
        except AttributeError:
            options.stderr.write( "# could not extract invariant job name from %s\n" % filename )
            nerrors += 1
            continue
        
        result_filename = options.output_pattern % job_name
        
        do = False
        status = "up-to-date"

        if options.force:
            status = "force"
            do = True
            
        if not do:
            if os.path.exists( result_filename ):
                if isNewer( filename, result_filename ):
                    status = "newer"
                    do = True
                    if options.remove_old:
                        files_to_delete.append( result_filename )
                if not do and not isComplete( result_filename ):
                    status = "incomplete"
                    do = True
                    if options.remove_old:
                        files_to_delete.append( result_filename )
            else:
                status = "missing"
                do = True

        if options.loglevel >= 1:
            options.stdlog.write("# %s->%s (%s)\n" % (filename, result_filename, status))
                
        if not do:
            nskipped += 1
            continue

        jobs.append( filename )

    ##############################################################
    ##############################################################
    ##############################################################
    ## delete old files
    ##############################################################
    for filename in files_to_delete:
        if os.path.exists( filename ):
            os.remove( filename )
            ndeleted += 1

    ##############################################################
    ##############################################################
    ##############################################################
    ## start jobs
    ##############################################################
    for filename in jobs:

        cmd = "qsub %s" % filename
        try:
            retcode = subprocess.call(cmd, shell=True)
            if retcode != 0:
                if options.loglevel >= 1:
                    options.stdlog.write("# ERROR: failed to execute %s\n" % cmd)
                nerrors += 1
                continue
        except OSError, e:
            if options.loglevel >= 1:
                options.stdlog.write("# ERROR: failed to execute %s with msg %s\n" % (cmd, e))
        nrun += 1

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, nrun=%i, nskipped=%i, ndeleted=%i, nerrors=%i\n" %\
                                 (ninput, nrun, nskipped, ndeleted, nerrors ) )

    E.Stop()

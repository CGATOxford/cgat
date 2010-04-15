#!/bin/env python
####
####
##
## Copyright (C) 2007 Andreas Heger All rights reserved
##
## Author: Andreas Heger <andreas.heger@dpag.ox.ac.uk>
##
## $Id: nofarm.py 2782 2009-09-10 11:40:29Z andreas $
##
##
####
####

USAGE="""farm.py [OPTIONS] cmd [ARGS] < stdin > stdout

execute a cmd on the cluster.

The input on stdin is split for embarrasingly parallel jobs.
The --split-at-.. options describe how standard input is to
be split. A temporary directory is created in the current
directory. This directory has be visible on the cluster nodes
and accessible under the same name.

The output is written to stdout. Results are returned in the
same order as they are submitted.

On error, error messages are echoed and nothing is returned.
The temporary directory is not deleted to allow manual recovery.

Examples:

The following command will split the file "go" at the first column 
and execute the command perl -p -e "s/GO/gaga/".

   cat go | farm.py --split-at-colum=1 perl -p -e "s/GO/gaga/"

The following command will split a fasta file at each entry and
compute an approximate sequence length:

   cat genome.fasta | farm.py --split-at-regex="^>(\S+)" "wc -c"

The following command will split a fasta file at every 10 sequences

   cat genome.fasta | farm.py --split-at-regex="^>(\S+)" --chunksize=10 "wc -c"

TODO: 

implement continuation of jobs
implement better error messages
"""

import os, sys, re, string, optparse, time, glob, subprocess, tempfile, shutil

import farm

import Experiment
import IOTools
import threadpool

##--------------------------------------------------------------------
def main():

    parser = farm.getOptionParser()

    (options, args) = Experiment.Start( parser, 
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
    Experiment.Stop()

if __name__ == '__main__':
    main()

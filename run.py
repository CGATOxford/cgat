#!/bin/env python
####
####
##
## Copyright (C) 2007 Andreas Heger All rights reserved
##
## Author: Andreas Heger <andreas.heger@dpag.ox.ac.uk>
##
## $Id: run.py 2782 2009-09-10 11:40:29Z andreas $
##
##
####
####

USAGE="""python run.py [OPTIONS] 

generic wrapper around a command-line. 

The wrapper records where a command was executed and how long
it took.
"""

import os, sys, re, string, optparse, time, glob, subprocess

import Experiment

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: run.py 2782 2009-09-10 11:40:29Z andreas $", 
                                    usage=USAGE )

    ## stop parsing options at the first argument
    parser.disable_interspersed_args()

    (options, args) = Experiment.Start( parser, 
                                        add_pipe_options = True )
    
    if len(args) > 0:
        
        cmd = args[0]
        if len(args) > 1:
            cmd += " '" + "' '".join(args[1:]) + "'" 
    
        s = subprocess.Popen( cmd,
                              shell = True,
                              cwd = os.getcwd(),
                              close_fds = True)                              

        (out, err) = s.communicate()
        returncode = s.returncode
    else:
        returncode = 0

    Experiment.Stop()

    sys.exit( returncode )

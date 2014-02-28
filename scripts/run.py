'''
run.py - wrapper around command line
====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script is a generic wrapper for tools with a command line interface.
The wrapper records execution times and logs successfull completion
of the command.

Usage
-----

Example::

   python run.py ls

Type::

   python run.py --help

for command line help.

Command line options
--------------------

'''
USAGE = """python run.py [OPTIONS] 

generic wrapper around a command-line. 

The wrapper records where a command was executed and how long
it took.
"""

import os
import sys
import re
import string
import optparse
import time
import glob
import subprocess

import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv == None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id: run.py 2782 2009-09-10 11:40:29Z andreas $",
                            usage=globals()["__doc__"])

    # stop parsing options at the first argument
    parser.disable_interspersed_args()

    (options, args) = E.Start(parser,
                              add_pipe_options=True)

    if len(args) > 0:

        cmd = args[0]
        if len(args) > 1:
            cmd += " '" + "' '".join(args[1:]) + "'"

        s = subprocess.Popen(cmd,
                             shell=True,
                             cwd=os.getcwd(),
                             close_fds=True)

        (out, err) = s.communicate()
        returncode = s.returncode
    else:
        returncode = 0

    E.Stop()

    sys.exit(returncode)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

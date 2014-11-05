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
gpipe/setup.py - setup Makefile pipelines
===================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script will setup a Makefile pipeline.

Usage
-----

Example::

   python gpipe/setup.py --help

Type::

   python gpipe/setup.py --help

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
import tempfile
import optparse
import time

import CGAT.Experiment as E


def AddOptions(outfile,  infile, source_directory, options):
    """add options from infile to outfile."""

    keep = False
    for line in infile:
        if re.search("Section parameters: start", line):
            keep = True
            continue
        elif re.search("Section parameters: end", line):
            keep = False
            continue
        elif re.match("include", line):
            other_file = re.match("include\s+(\S+)", line).groups()[0]
            other_file = re.sub(
                "\$\(DIR_SCRIPTS_GENEPREDICTION\)", source_directory + "/", other_file)
            AddOptions(
                outfile, open(other_file, "r"), source_directory, options)
        elif re.match("PARAM_PROJECT_NAME", line):
            if options.project:
                line = "PARAM_PROJECT_NAME=%s\n" % options.project
        elif re.match("DIR_SCRIPTS_GENEPREDICTION", line):
            line = "DIR_SCRIPTS_GENEPREDICTION=%s/\n" % source_directory
        elif re.match("DIR_ROOT", line):
            if options.root_dir:
                line = "DIR_ROOT=%s/\n" % os.path.abspath(options.root_dir)
            else:
                line = "DIR_ROOT=%s/\n" % os.path.abspath(options.destination)
        elif re.match("DIR_TMP_SHARED_LOCAL", line):
            if options.temporary_local:
                line = "DIR_SCRIPTS_TMP_SHARED_LOCAL=%s\n" % options.temporary_local
        elif re.match("DIR_TMP_SHARED_REMOTE", line):
            if options.temporary_remote:
                line = "DIR_SCRIPTS_TMP_SHARED_REMOTE=%s\n" % options.temporary_remote
        if keep:
            outfile.write(re.sub("\?=", "=", line))

    outfile.write("\n")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/setup.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-f", "--force-output", dest="force", action="store_true",
                      help="force overwrite of existing Makefile.")
    parser.add_option("-d", "--dest", dest="destination", type="string",
                      help="destination directory.")
    parser.add_option("-p", "--project", dest="project", type="string",
                      help="name of the project.")
    parser.add_option("-t", "--tempdir", dest="temporary", type="string",
                      help="temporary directory, used if local/remote are not set.")
    parser.add_option("-l", "--tempdir-local", dest="temporary_local", type="string",
                      help="pathname to local temporary directory.")
    parser.add_option("-r", "--tempdir-remote", dest="temporary_remote", type="string",
                      help="pathname to remote temporary directory.")
    parser.add_option("-o", "--rootdir", dest="root_dir", type="string",
                      help="root directory.")

    parser.add_option("-y", "--python-libraries", dest="python_libraries", type="string",
                      help="pathname to directory with python libraries.",
                      action="append")

    parser.add_option("-i", "--ld-libraries", dest="ld_libraries", type="string",
                      help="pathname to directory with libaries (LD_LIBRARY_PATH).",
                      action="append")

    parser.add_option("-b", "--binaries", dest="binaries", type="string",
                      help="pathname with binaries.",
                      action="append")

    parser.add_option("-c", "--include-pattern", dest="includes", type="string",
                      help="pathnames of Makefile to be included BEFORE the parameter section of the Makefile.",
                      action="append")

    parser.add_option("-m", "--method", dest="method", type="string",
                      help="""method to install [gpipe|gpipe_export|blast|introns|orthology|structures|synteny|
                      orthology_pairwise|orthology_pairwise_multiple|orthology_multiple|orthology_malis|
                      codonbias|codonbias_multiple|analyze_codonbias_duplications].""" )

    parser.set_defaults(
        destination=".",
        force=False,
        project=None,
        temporary="/net/cpp-group/gpipe/tmp",
        temporary_local=None,
        temporary_remote=None,
        method=None,
        includes=[],
        root_dir=None,
    )

    (options, args) = E.Start(parser)

    if not options.method:
        raise ValueError("please specify a method.")

    if len(args) >= 1 and not options.destination:
        options.destination = args[0]
        del args[0]

    if not options.destination:
        raise ValueError(
            "please specify a destination directory for installing the scripts.")

    if not options.temporary_remote:
        options.temporary_remote = options.temporary_local

    source_directory = os.path.realpath(
        os.path.abspath(os.path.dirname(__file__)))
    target_directory = options.destination

    # create working directories
    if not os.path.exists(target_directory):
        os.makedirs(target_directory)
    if options.temporary_local:
        if not os.path.exists(options.temporary_local):
            os.makedirs(options.temporary_local)
        if not os.path.exists(options.temporary_remote):
            os.makedirs(options.temporary_remote)

    target_makefile = target_directory + "/Makefile"
    target_script = target_directory + "/setup"

    if os.path.exists(target_makefile) and not options.force:
        print "Makefile %s already exists. Use -f to overwrite" % target_makefile
        sys.exit(1)

    # Create makefile
    source_makefile = source_directory + \
        "/makefiles/Makefile.%s" % options.method
    if not os.path.exists(source_makefile):
        raise ValueError("unknown method %s: %s not found" %
                         (options.method, source_makefile))

    # 1. write parameters
    keep = 0
    outfile = open(target_makefile, "w")
    outfile.write("############################################\n")
    outfile.write("# created from %s\n#\n" % source_makefile)
    outfile.write("# created at %s\n#\n" %
                  time.asctime(time.localtime(time.time())))
    outfile.write("# python %s\n#\n" % " ".join(sys.argv))
    outfile.write("############################################\n")

    AddOptions(outfile,
               open(source_makefile, "r"),
               source_directory,
               options)

    # add include makefiles
    if options.includes:
        outfile.write("############################################\n")
        outfile.write("# default parameter includes\n")
        outfile.write("# These overwrite all previous parameters.\n")

        for filename in options.includes:
            outfile.write("include %s\n\n" % os.path.abspath(filename))

    # 2. write user specified variables
    outfile.write("############################################\n")
    outfile.write("# user specified options\n\n")
    if len(args) > 0:
        for arg in args:
            outfile.write(arg + "\n\n")

    # 3. write paths
    outfile.write("############################################\n")
    outfile.write("# include master makefile\n\n")
    outfile.write("include %s\n\n" % source_makefile)

    outfile.write("############################################\n")
    outfile.close()

    print "Setup of module '%s' complete.\n" % (options.method)
    print "Working directory is %s.\n" % (os.path.realpath(target_directory))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

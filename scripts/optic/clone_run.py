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
optic/clone_run.py - 
======================================================

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

   python optic/clone_run.py --help

Type::

   python optic/clone_run.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import subprocess
import CGAT.Experiment as E

""" program $Id: optic/clone_run.py 2781 2009-09-10 11:33:14Z andreas $

clone a prediction run.

1. copy all tables from new schema into old schema.

"""


def Run(statement, stdout=sys.stdout, stderr=sys.stderr):

    s = subprocess.Popen(statement,
                         shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         close_fds=True)

    (out, err) = s.communicate()

    if s.returncode != 0:
        raise "Error %s while executing %s " % (err, statement)

    stdout.write(out)
    stdout.write(err)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: optic/clone_run.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-o", "--old-schema", dest="old_schema", type="string",
                      help="old_schema.")

    parser.add_option("-n", "--new-schema", dest="new_schema", type="string",
                      help="new_schema.")

    parser.add_option("-O", "--old-directory", dest="old_directory", type="string",
                      help="old directory.")

    parser.add_option("-N", "--new-directory", dest="new_directory", type="string",
                      help="new directory.")

    parser.set_defaults(
        old_schema=None,
        new_schema=None,
        old_directory=None,
        new_directory=None,
    )

    files = "peptides2genes", "peptides.fasta", "genome*.fasta", "reference.exons"

    (options, args) = E.Start(parser,
                              add_database_options=True,
                              add_pipe_options=True)

    if not options.old_schema:
        raise "please specify an old schema to copy data from."

    if not options.old_schema:
        raise "please specify a new schema to copy data to."

    if not options.old_directory:
        options.old_directory = options.old_schema

    if not os.path.exists(options.old_directory):
        raise "old directory %s does not exist." % options.old_directory

    if not options.new_directory:
        options.new_directory = options.new_schema

    host, db = options.psql_connection.split(":")

    print "######### creating directory %s" % (options.new_directory)

    Run("mkdir %s" % options.new_directory)

    print "######### checking out code into %s" % (options.new_directory)

    Run("svn co svn://fgu202/andreas/gpipe/trunk %s/src" %
        options.new_directory)

    print "######### creating makefile from %s to %s" % (options.old_directory, options.new_directory)

    statement = 'perl -p -e "s/%s/%s/g; s/%s/%s/g" < %s/Makefile > %s/Makefile' %\
                (options.old_schema, options.new_schema,
                 options.old_directory, options.new_directory,
                 options.old_directory, options.new_directory)

    Run(statement, stdout=options.stdout, stderr=options.stderr)

    print "########## setting up links from %s to %s" % (options.old_schema, options.new_schema)

    for f in files:
        Run("cd %s; ln -fs ../%s/%s .; cd .." %
            (options.new_directory, options.old_directory, f))

    print "########### running make prepare"

    # this will create temporary directories

    Run("make -C %s prepare" %
        options.new_directory, stdout=options.stdout, stderr=options.stderr)

    E.Stop()

    print "########### running make all"

    Run("make -C %s -t all" %
        options.new_directory, stdout=options.stdout, stderr=options.stderr)

    E.Stop()

    print "########## cloning schema %s to %s" % (options.old_schema, options.new_schema)

    statement = 'psql -h %s %s -c "DROP SCHEMA %s CASCADE"' %\
                (host, db, options.new_schema)

    Run(statement, stdout=options.stdout, stderr=options.stderr)

    statement = '/usr/bin/pg_dump -h %s -v -n %s %s | perl -p -e "s/%s/%s/g" | psql -h %s %s' %\
                (host, options.old_schema, db,
                 options.old_schema, options.new_schema,
                 host, db)

    Run(statement, stdout=options.stdout, stderr=options.stderr)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

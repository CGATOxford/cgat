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
gpipe/gbrowser_delete_features.py - 
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

   python gpipe/gbrowser_delete_features.py --help

Type::

   python gpipe/gbrowser_delete_features.py --help

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
import random
import optparse
import CGAT.Experiment as E
import MySQLdb
import _mysql

USAGE = """python %s.py [OPTIONS] value

delete a section from a gbrowser database.
""" % (sys.argv[0])


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/gbrowser_delete_features.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-d", "--dry-run", dest="dry_run", action="store_true",
                      help="dry run, do not execute any commands.")
    parser.add_option("-f", "--source-type", dest="source_type", type="string",
                      help="source type to remove.")

    parser.set_defaults(
        dry_run=False,
        source_type=None,
    )

    (options, args) = E.Start(parser, add_database_options=True)

    dbhandle = MySQLdb.Connect(host=options.host,
                               user=options.user,
                               passwd=options.password,
                               port=options.port)

    if not options.database:
        raise "please supply a database."

    dbhandle.cursor().execute("USE %s" % options.database)

    cc = dbhandle.cursor()
    cc.execute("SHOW TABLES")
    tables = map(lambda x: x[0], cc.fetchall())
    cc.close()

    if "ftype" not in tables:
        raise "table ftype not found in %s - is it a gbrowser database?" % (
            self.database)

    if options.source_type:

        source_types = options.source_type.split(",")

        if len(source_types) == 1 and "%" in source_types[0]:
            statement = "SELECT ftypeid FROM ftype WHERE fsource LIKE '%s'" % source_types[
                0]
        else:
            statement = "SELECT ftypeid FROM ftype WHERE fsource IN ('%s')" % "','".join(
                source_types)

        if options.dry_run:
            print statement

        cc = dbhandle.cursor()
        cc.execute(statement)
        typeids = map(lambda x: str(x[0]), cc.fetchall())
        cc.close()

        if not typeids:
            raise "no features of type %s found" % (options.source_type)

        # find annotations
        statement = "SELECT DISTINCTROW a.fid FROM fattribute_to_feature AS a, fdata AS b WHERE b.ftypeid IN (%s) AND a.fid = b.fid" % ",".join(
            typeids)
        if options.dry_run:
            print statement

        cc = dbhandle.cursor()
        cc.execute(statement)
        fids = map(lambda x: str(x[0]), cc.fetchall())
        cc.close()

        # find gene identifiers (gid) in fdata for group
        statement = "SELECT DISTINCTROW b.gid FROM fdata AS b WHERE b.ftypeid IN (%s)" % ",".join(
            typeids)
        if options.dry_run:
            print statement

        cc = dbhandle.cursor()
        cc.execute(statement)
        gids = map(lambda x: str(x[0]), cc.fetchall())
        cc.close()

        options.stdout.write(
            "deleting %i genes from %s.\n" % (len(gids), options.database))
        options.stdout.flush()

        statement = "DELETE FROM fgroup WHERE gid IN (%s)" % ",".join(gids)
        if options.dry_run:
            print statement
        else:
            dbhandle.cursor().execute(statement)

        options.stdout.write(
            "deleting attributes for %i features from %s.\n" % (len(fids), options.database))
        options.stdout.flush()

        statement = "DELETE FROM fattribute_to_feature WHERE fid IN (%s)" % ",".join(
            fids)
        if options.dry_run:
            print statement
        else:
            dbhandle.cursor().execute(statement)

        options.stdout.write(
            "deleting features %s from %s.\n" % (",".join(typeids), options.database))
        options.stdout.flush()

        statement = "DELETE FROM fdata WHERE ftypeid IN (%s)" % ",".join(
            typeids)
        if options.dry_run:
            print statement
        else:
            dbhandle.cursor().execute(statement)

        options.stdout.write("deleting types %s from %s.\n" %
                             (",".join(typeids), options.database))
        options.stdout.flush()

        statement = "DELETE FROM ftype WHERE ftypeid IN (%s)" % ",".join(
            typeids)
        if options.dry_run:
            print statement
        else:
            dbhandle.cursor().execute(statement)

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

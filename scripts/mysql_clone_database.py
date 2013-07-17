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
mysql_clone_database.py - clone a mysql database
================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script clones a mysql database.

Usage
-----

Example::

   python mysql_clone_database.py src dest

will clone the database ``src`` to ``dest``. 

Type::

   python mysql_clone_database.py --help

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
import MySQLdb
import _mysql

USAGE="""python clone_database.py src dest

Clone a database by copying tables from one
database to another.
"""

def copyTables(dbhandle, dest, src, dry_run = False ):
    """copy tables from dest into src."""
    
    cc = dbhandle.cursor()
    cc.execute("SHOW tables FROM %s" % src)
    src_tables = map(lambda x: x[0], cc.fetchall())        
    cc.close()

    cc = dbhandle.cursor()
    cc.execute("SHOW tables FROM %s" % dest)
    dest_tables = set(map(lambda x: x[0], cc.fetchall()))
    cc.close()

    cc = dbhandle.cursor()
    cc.execute( "USE %s" % dest )
        
    for table in src_tables:

        if table in dest_tables:
            delete_statement = "DROP TABLE %s.%s" % (dest, table)
            if dry_run:
                print delete_statement
            else:
                cc.execute( delete_statement )
                
        cc.execute( "SHOW CREATE TABLE %s.%s" % (src, table) )
        create_statement = cc.fetchone()[1]
        copy_statement = "INSERT INTO %s.%s SELECT * FROM %s.%s" % (dest, table, src, table )
        if dry_run:
            print create_statement
            print copy_statement
        else:
            cc.execute( create_statement )
            cc.execute( copy_statement )
    
    cc.close()
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: mysql_clone_database.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-c", "--create", dest="create", action="store_true",
                      help="create target database if it does not exist."  )

    parser.add_option("-d", "--dry-run", dest="dry_run", action="store_true",
                      help="dry run, do not execute any commands."  )


    parser.set_defaults (
        connection = "db",
        create = False,
        dry_run = False,
    )

    (options, args) = E.Start( parser, add_mysql_options = True )

    if len(args) != 2:
        print USAGE
        sys.exit(1)
        
    src, dest = args

    if src == dest:
        if options.loglevel >= 1:
            options.stdlog.write("# src and dest identical - no copying performed")            
        
        sys.exit(0)
        E.Stop()

    dbhandle = MySQLdb.Connect( host = options.host,
                                user = options.user,
                                passwd = options.password,
                                port = options.port )
                                               
    cc = dbhandle.cursor()

    cc.execute("SHOW databases")
    databases = map(lambda x: x[0], cc.fetchall())
    cc.close()
    
    databases = map(lambda x: x[0], cc.fetchall())

    if src not in databases:
        raise "unknown database %s" % src

    if dest not in databases:
        if options.create:
            if options.loglevel >= 1:
                options.stdlog.write("# creating database %s" % src)
                
            if not options.dry_run:
                cc = dbhandle.cursor()
                cc.execute("CREATE DATABASE %s" % dest)
                cc.close()

    copyTables( dbhandle, dest, src, options.dry_run )
        
    E.Stop()

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
mysql_copy_tables.py - 
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

   python mysql_copy_tables.py --help

Type::

   python mysql_copy_tables.py --help

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

# name, database 
data = (
    ( "Dmelanogaster", "dmel_vs_dmel4" ),
    ( "Dyakuba", "dyak_vs_dmel7" ),
    ( "Dpseudoobscura", "dpse_vs_dmel8" ),
    ( "Dvirilis", "dvir_vs_dmel6" ),
    ( "Dmojavanesis", "dmoj_vs_dmel6" ),
    ( "Dananassae", "dana_vs_dmel6" ),
    ( "Dgrimshawi", "dgri_vs_dmel5" ),
    ( "Dsimulans", "dsim_vs_dmel6" ),
    ( "Derecta", "dere_vs_dmel6" ),
    ( "Dsechellia", "dsec_vs_dmel3" ),
    ( "Dpersimilis", "dper_vs_dmel3" ),
    ( "Dwillistoni", "dwil_vs_dmel1" ),    
    )

#data = ( ( "Monodelphis", "mono2" ), )
#data = ( ( "Monodelphis", "mono_p1" ), )
#data = ( ( "Dsimulans", "dsim_vs_dmel6"), )

USAGE="""python gpipe/gbrowser_clone_devel.py

Clone the devel installation into the target installation.

The following cloing operations will be done:

%s
""" % (str(data))


class X:
    def __init__(self):
        pass

    def execute(self, msg):
        print msg

    def fetchall(self):
        return ( ("test",), )
    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: mysql_copy_tables.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-t", "--target-dir", dest="target", type="string",
                      help="target directory.", metavar = "FILE"  )

    parser.set_defaults (
        connection = "db",
        source= "/var/www/conf/gbrowse-devel.conf",
        target= "/var/www/conf/gbrowse.conf",
        source_database = "gbrowser_devel_%s",
        target_database = "gbrowser_%s",        
    )

    (options, args) = E.Start( parser, add_mysql_options = True )

    dbhandle = MySQLdb.Connect( host = options.host,
                                user = options.user,
                                passwd = options.password,
                                port = options.port )
                                               
    cc = dbhandle.cursor()
    cc.execute("SHOW databases")
    databases = map(lambda x: x[0], cc.fetchall())
    cc.close()
    
    databases = map(lambda x: x[0], cc.fetchall())

    for name, schema in data:

        old_name = options.source_database % schema
        new_name = options.target_database % schema

        if old_name in databases:
            print "cloning: %s -> %s" % (old_name, new_name)

        cc = dbhandle.cursor()
        cc.execute("SHOW tables FROM %s" % old_name)
        tables = map(lambda x: x[0], cc.fetchall())        
        cc.close()

        ## I want to do this without file system access and
        ## without delay
        ## 1. delete target database
        ## 2. recreate target database
        ## 3. rename tables to target database
        ## 4. recreate tables in source database
        ## 5. copy data to source database
        cc = dbhandle.cursor()
        if new_name in databases:
            cc.execute( "DROP DATABASE %s" % new_name )

        cc.execute( "CREATE DATABASE %s" % new_name )
        
        cc = dbhandle.cursor()
        for table in tables:
            cc.execute( "ALTER TABLE %s.%s RENAME AS %s.%s" % ( old_name, table, new_name, table))

        print "copying has finished."
        cc.close()

        cc = dbhandle.cursor()
        cc.execute( "USE %s" % old_name )
        
        for table in tables:
            print "recreating table %s in %s" % (table, old_name)
            cc.execute( "SHOW CREATE TABLE %s.%s" % (new_name, table) )
            create_statement = cc.fetchone()[1]
            cc.execute( create_statement )
            cc.execute( "INSERT INTO %s.%s SELECT * FROM %s.%s" % (old_name, table, new_name, table ))
            
        cc.close()
        
    E.Stop()

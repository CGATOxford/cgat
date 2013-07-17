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
psql_clone_database.py - clone a psql schema
============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script will clone a psql schema

Usage
-----

Example::

   python psql_clone_database.py src dest

will clone the schema ``src`` to ``dest``. 

Type::

   python psql_clone_database.py --help

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
import time
import optparse
import math
import tempfile
import copy
import subprocess

import pgdb

import CGAT.Experiment as E

##---------------------------------------------------------------------------------------------
def DbExecute( dbhandle, statement ):

    cc = dbhandle.cursor()
    cc.execute(statement)
    result = cc.fetchall()
    cc.close()
    return result

##---------------------------------------------------------------------------------------------
def DbDo( dbhandle, statement ):

    cc = dbhandle.cursor()
    cc.execute(statement)
    dbhandle.commit()
    cc.close()

##------------------------------------------------------------------------------------------
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: psql_clone_database.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.set_defaults( \
        start = None,
        fields = [] )

    (options, args) = E.Start( parser, add_psql_options = True )
    
    dbhandle = pgdb.connect( options.psql_connection )

    if len(args) != 2:
        raise ValueError("please supply src and dest." )
    
    src, dest = args

    sql_get_schemas = "SELECT DISTINCT schemaname FROM pg_catalog.pg_stat_all_tables"
    sql_create_schema = "CREATE SCHEMA %s"

    cc = dbhandle.cursor()

    cc.execute( sql_get_schemas )
    databases = map(lambda x: x[0], cc.fetchall())
    cc.close()

    if src not in databases:
        raise ValueError("unknown schema %s" % src )

    if dest in databases:
        raise ValueError( "schema %s already exists" % (dest ) )

    host, database = options.psql_connection.split(":")

    cmd = 'pg_dump --ignore-version --schema-only --host=%(host)s --schema=%(src)s | sed "s/%(src)s/%(dest)s/g" | psql --host=%(host)s %(database)s' % (locals() ) 
    p = subprocess.Popen( cmd,
                          shell = True )
    sts = os.waitpid(p.pid, 0)

    cmd = 'pg_dump --ignore-version --data-only --host=%(host)s --schema=%(src)s | psql --host=%(host)s %(database)s' % (locals() ) 
    p = subprocess.Popen( cmd,
                          shell = True )
    sts = os.waitpid(p.pid, 0)

    E.Stop()

    

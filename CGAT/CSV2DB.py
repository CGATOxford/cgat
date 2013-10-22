#! /bin/env python
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
CSV2DB.py - utilities for uploading a table to database
=======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

create a table from a csv separated file and load data into it.

This module supports backends for postgres and sqlite3. Column types are
auto-detected.

.. todo::

   Use file import where appropriate to speed up loading. Currently, this is
   not always the case.

Usage
-----

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
import tempfile
import subprocess
import types

import Experiment as E
import csv
import CSV as CSV
import sqlite3

def executewait( dbhandle, statement, error, retry = False, wait=5):
    '''execute sql statement.

    Retry on error, if retry is True.
    Returns a cursor object.
    '''

    cc = dbhandle.cursor()    
    i = 20
    while i>0:
        try:
            cc.execute( statement )
            return cc
        except sqlite3.OperationalError as e:
            msg = e.message
            E.warn("import failed: msg=%s, statement=\n  %s" % (msg, statement ) )
        # TODO: check for database locked msg
            if not retry:
                raise error, msg
            if not re.search("locked", str(msg)):
                raise error, msg
            time.sleep(wait)
            i -= 1
            continue
        break
    raise sqlite.OperationalError("Database locked and too many retries")

def quoteRow( row, take, 
              map_column2type, 
              missing_values, 
              null = "NULL",
              string_value = "%s" ):
    """return a dictionary with properly quoted values."""

    # set empty values for int/float to NULL
    d = {}
    for t in take:
        v = row[t]
        if v == "": 
            d[t] = null
        elif v in missing_values: 
            d[t] = null
        elif map_column2type[t] in (types.IntType, types.FloatType):
            d[t] = str(row[t])
        else:
            d[t] = string_value % row[t]
            
    return d

def quoteTableName( name, quote_char = "_", backend="sqlite" ):
    if backend == "sqlite":
        # no special characters. Column names can not start with a number.
        if name[0] in "0123456789": name = "_" + name
        return re.sub( "[-(),\[\].:]", "_", name )
    elif backend in ("mysql","pg"):
        if name[0] in "0123456789": name = "_" + name
        return re.sub( "[-(),\[\]]:", "_", name )

def createTable( dbhandle, error, options, rows = None, headers = None,
                 first_column = None,
                 existing_tables = []):

    ## create table by guessing column types from data type.
    if rows:
        map_column2type, ignored, max_values = CSV.GetMapColumn2Type( rows,
                                                                      ignore_empty = options.ignore_empty,
                                                                      get_max_values = True )
        if ignored:
            E.info( "ignored columns: %s" % str(ignored) )

        headers = map_column2type.keys()
        headers.sort()

    elif headers:
        map_column2type = dict( zip( headers, [None,] * len(headers) ) )
        ignored = 0

    columns_to_ignore = set( [ x.lower() for x in options.ignore_columns] )
    columns_to_rename = dict( [x.lower().split(":") for x in options.rename_columns ] )

    take = []
    ## associate headers to field names
    columns = []
    present = {}
    for header_index, h in enumerate(headers):
        hh = h
        if options.lowercase:
            hh = string.lower(h)

        if hh in columns_to_ignore: continue
    
        if hh in present:
            if options.ignore_duplicates:
                continue
            else:
                raise ValueError("duplicate column %s" % hh)

        present[hh] = 1
        take.append( h )

        if map_column2type[h] == types.IntType:
            max_value = max_values[h]
            if h > 2147483647:
                t = "BIGINT DEFAULT '0'"
            elif h > 32767:
                t = "INTEGER DEFAULT '0'"
            else:
                t = "SMALLINT DEFAULT '0'"

        elif map_column2type[h] == types.FloatType:
            t = "FLOAT DEFAULT '0'"
        else:
            t = "TEXT"

        # remove special characters from column names
        if hh == "": 
            if first_column != None and header_index == 0:
                hh = first_column
            else:
                raise ValueError("column '%s' without header " % h )
        hh = columns_to_rename.get( hh, hh )
        hh = re.sub( '''['"]''', "", hh)
        hh = re.sub( "[,;.:\-\+/ ()%?]", "_", hh)
        if hh[0] in "0123456789": hh = "_" + hh
        columns.append( "%s %s" % (hh, t))

    # delete old table if it exists
    while 1:
        try:
            cc = dbhandle.cursor()
            cc.execute("DROP TABLE IF EXISTS '%s'" % options.tablename)
            dbhandle.commit()
            cc.close()
            E.info("existing table %s deleted" % options.tablename )
        except sqlite3.OperationalError, msg:
            E.warn( msg )
            time.sleep(5)
            continue
        except error, msg:
            E.warn( "could not delete existing table %s: %s" % (options.tablename, str(msg) ))
            dbhandle.rollback()
            if not options.retry:
                raise error, msg
            elif options.tablename in existing_tables:
                # table exists, but drop did not work (e.g. database lock)
                time.sleep(5)
                continue
            else:
                # table might not have existed
                break
        break
    
    ## create new table
    statement = "CREATE TABLE %s ( %s );" % (options.tablename, ", ".join( columns))

    E.debug( "table create:\n# %s" % (statement ) )
        
    while 1:
        try:
            cc = dbhandle.cursor()
            cc.execute(statement)
            cc.close()
            dbhandle.commit()
        except error, msg:
            E.warn("table creation failed: msg=%s, statement=\n  %s" % (msg, statement ) )
            # TODO: check for database locked msg
            if not options.retry:
                raise error, msg
            if not re.search("locked", str(msg)):
                raise error( "%s: %s" % (msg, statement))
            time.sleep(5)
            continue
        break

    E.info("table %s created successfully." % options.tablename )
    
    return take, map_column2type, ignored

def run( infile, options ):

    options.tablename = quoteTableName( options.tablename, backend = options.backend )

    if options.map:
        m = {}
        for x in options.map:
            f,t = x.split(":")
            m[f] = t
        options.map = m
    else:
        options.map = {}
        
    existing_tables = None
    
    if options.backend == "pg":
        import pgdb
        dbhandle = pgdb.connect( options.psql_connection )
        error = pgdb.DatabaseError
        options.null = "NULL"
        options.string_value = "'%s'"
        if options.insert_quick:
            raise ValueError("quick import not implemented.")

    elif options.backend == "sqlite":
        import sqlite3
        dbhandle = sqlite3.connect( options.database )
        try:
            os.chmod( options.database, 0664 )
        except OSError, msg:
            E.warn("could not change permissions of database: %s" % msg )

        # Avoid the following error:
        # sqlite3.ProgrammingError: You must not use 8-bit bytestrings unless you use a text_factory that can interpret 8-bit bytestrings (like text_factory = str). It is highly recommended that you instead just switch your application to Unicode strings
        # Note: might be better to make csv2db unicode aware.
        dbhandle.text_factory = str

        error = sqlite3.OperationalError
        options.insert_many = True  # False
        options.null = None # "NULL" 
        options.string_value = "%s" # "'%s'"

        statement = "SELECT name FROM sqlite_master WHERE type='table'"
        cc = executewait( dbhandle, statement, error, options.retry )
        existing_tables = set( [ x[0] for x in cc ] )
        cc.close()

        quick_import_statement = "sqlite3 -header -csv -separator '\t' %s '.import %%s %s'" % (options.database, options.tablename)

    if options.header != None:
        options.header = [x.strip() for x in options.header.split(",")]

    if options.utf:
        reader = CSV.UnicodeDictReader( infile, 
                             dialect=options.dialect, 
                             fieldnames = options.header )
    else:
        reader = CSV.DictReader( infile, 
                                 dialect=options.dialect, 
                                 fieldnames = options.header )

    if options.replace_header:
        reader.next()

    E.info("reading %i columns to guess column types" % options.guess_size)

    rows = []
    for row in reader:
        if None in row: 
            raise ValueError( "undefined columns in input file at row: %s" % row )
        
        try:
            rows.append( CSV.ConvertDictionary( row , map=options.map ))
        except TypeError, msg:
            E.warn( "incomplete line? Type error in conversion: '%s' with data: %s" % (msg, str(row) ) )
        except ValueError, msg:
            E.warn( "incomplete line? Type error in conversion: '%s' with data: %s" % (msg, str(row) ) )

        if len(rows) >= options.guess_size:
            break

    E.info("read %i rows for type guessing" % len(rows) )
    E.info("creating table" )

    if len(rows) == 0:
        if options.allow_empty:
            if not reader.fieldnames:
                E.warn( "no data - no table created")
            else:
                # create empty table and exit
                take, map_column2type, ignored = createTable( dbhandle, error, headers = reader.fieldnames,
                                                              options = options,
                                                              first_column = options.first_column,
                                                              existing_tables = existing_tables )
                E.info( "empty table created" )
            return
        else:
            raise ValueError( "empty table" )
    else:
        take, map_column2type, ignored = createTable( dbhandle, error,
                                                      rows = rows,
                                                      options = options,
                                                      first_column = options.first_column,
                                                      existing_tables = existing_tables )


    def row_iter( rows, reader):
        for row in rows: 
            yield quoteRow( row, take, map_column2type, 
                            options.missing_values, 
                            null = options.null,
                            string_value = options.string_value )
        for data in reader: 
            yield quoteRow( CSV.ConvertDictionary( data , map=options.map ), take, map_column2type, 
                            options.missing_values, 
                            null = options.null,
                            string_value = options.string_value )

    ninput = 0

    E.info("inserting data")

    if options.insert_quick:
        E.info( "using quick insert" )

        outfile, filename = tempfile.mkstemp()

        E.debug("dumping data into %s" % filename )

        for d in row_iter( rows, reader ):

            ninput += 1
            os.write( outfile, "\t".join( [ str(d[x]) for x in take ] ) + "\n" )

            if options.loglevel >= 1 and ninput % options.report_step == 0:
                options.stdlog.write( "# iteration %i\n" % ninput )
                options.stdlog.flush()

        os.close( outfile )

        statement = quick_import_statement % filename
        E.debug( statement )

        # infinite loop possible
        while 1:
            
            retcode = E.run( statement, cwd = os.getcwd(), close_fds = True)
        
            if retcode != 0:
                E.warn("import error using statement: %s" % statement)

                if not options.retry:
                    raise ValueError( "import error using statement: %s" % statement )
            
                time.sleep(5)
                continue

            break
        
        os.remove( filename )
        
        # there is no way to .insert NULL values into sqlite. The only solution is
        # to update all colums.
        for column in take:
            executewait( dbhandle, 
                         "UPDATE %s SET %s = NULL WHERE %s = 'None'" % (options.tablename, column, column),
                         error, 
                         options.retry)

    elif options.insert_many:
        data = []
        for d in row_iter( rows, reader ):
            ninput += 1

            data.append( [d[x] for x in take] )

            if ninput % options.report_step == 0:
                E.info( "iteration %i" % ninput )
                
        statement = "INSERT INTO %s VALUES (%s)" % (options.tablename, ",".join( "?" * len(take ) ) )

        E.info( "inserting %i rows" % len(data) )
        E.debug( "multiple insert:\n# %s" % statement )

        while 1:
            try:
                dbhandle.executemany( statement, data )
            except error, msg:
                E.warn("import failed: msg=%s, statement=\n  %s" % (msg, statement ) )
                # TODO: check for database locked msg
                if not options.retry:
                    raise error, msg
                if not re.search("locked", str(msg)):
                    raise error, msg
                time.sleep(5)
                continue
            break


    else:
        ## insert line by line (could not figure out how to do bulk loading with subprocess and COPY FROM STDIN)
        statement = "INSERT INTO %s VALUES (%%(%s)s)" % (options.tablename,
                                                         ')s, %('.join( take ))

        # output data used for guessing:
        for d in row_iter( rows, reader):

            ninput += 1

            E.debug( "single insert:\n# %s" % (statement % d) )
            cc = executewait( dbhandle, statement % d, error, options.retry )
            cc.close()
            
            if options.loglevel >= 1 and ninput % options.report_step == 0:
                E.info("iteration %i" % ninput )

    E.info("building indices")
    nindex = 0
    for index in options.indices:
        
        nindex += 1
        try:
            statement = "CREATE INDEX %s_index%i ON %s (%s)" % ( options.tablename, nindex, options.tablename, index )
            cc = executewait( dbhandle, statement, error, options.retry )
            cc.close()
            E.info( "added index on column %s" % (index) )
        except error, msg: 
            E.info( "adding index on column %s failed: %s" % (index, msg) )

    statement = "SELECT COUNT(*) FROM %s" % ( options.tablename )
    cc = executewait( dbhandle, statement, error, options.retry )
    result = cc.fetchone()
    cc.close()

    noutput = result[0]
    
    E.info("ninput=%i, noutput=%i, nskipped_columns=%i" % (ninput, noutput, len(ignored)) )

    dbhandle.commit()

def buildParser( ):

    parser = E.OptionParser( version = "%prog version: $Id: csv2db.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "--dialect", dest="dialect", type="string",
                      help="csv dialect to use [default=%default]." )

    parser.add_option("-m", "--map", dest="map", type="string", action="append",
                      help="explicit mapping function for columns The format is column:type (e.g.: length:int) [default=%default]." )

    parser.add_option("-t", "--table", dest="tablename", type="string",
                      help="table name for all backends [default=%default]." )

    parser.add_option("-d", "--database", dest="database", type="string",
                      help="database name for sqlite3 [default=%default]." )

    parser.add_option("-H", "--header", dest="header", type="string",
                      help="',' separated list of header for files without header [default=%default]." )

    parser.add_option( "--replace-header", dest="replace_header", action="store_true",
                       help="replace header with --header instead of adding it [default=%default]." )

    parser.add_option("-l", "--lowercase", dest="lowercase", action="store_true",
                      help="force lower case column names [default=%default]." )

    parser.add_option("-u", "--ignore-duplicates", dest="ignore_duplicates", action="store_true",
                      help="ignore columns with duplicate names [default=%default]." )

    parser.add_option("-s", "--ignore-same", dest="ignore_same", action="store_true",
                      help="ignore columns with identical values [default=%default]." )

    parser.add_option("--ignore-column", dest="ignore_columns", type="string", action="append",
                      help="ignore columns [default=%default]." )

    parser.add_option("--rename-column", dest="rename_columns", type="string", action="append",
                      help="rename columns [default=%default]." )

    parser.add_option("--first-column", dest="first_column", type="string", 
                      help="name of first column - permits loading CSV table where the first "
                      "column name is the empty string [default=%default]." )
    
    parser.add_option("-e", "--ignore-empty", dest="ignore_empty", action="store_true",
                      help="ignore columns which are all empty [default=%default]." )

    parser.add_option("-q", "--quick", dest="insert_quick", action="store_true",
                      help="try quick file based import - needs to be supported by the backend [default=%default]." )

    parser.add_option("-b", "--backend", dest="backend", type="choice",
                      choices=("pg", "sqlite"),
                      help="database backend to choose [default=%default]." )

    parser.add_option("-i", "--index", dest="indices", type="string", action="append",
                      help="create an index for the named column [default=%default]." )

    parser.add_option("-a", "--allow-empty", dest="allow_empty", action="store_true",
                      help="allow empty table [default=%default]." )

    parser.add_option( "--retry", dest="retry", action="store_true",
                      help="retry if an SQL statement fails - warning: THIS MIGHT CAUSE DEADLOCKS [default=%default]." )

    parser.add_option( "-z", "--from-zipped", dest="from_zipped", action="store_true",
                       help="input is zipped.")

    parser.add_option( "--utf8", dest="utf", action="store_true",
                       help="standard in is encoded as UTF8 rather than local default, WARNING: does not strip comment lines yet [default=%default]")

    parser.set_defaults(
        map = [],
        dialect = "excel-tab",
        database = "csvdb",
        lowercase = False,
        tablename = "csv",
        from_zipped = False,
        ignore_duplicates= False,
        ignore_identical = False,
        ignore_empty = False,
        insert_many = False,
        ignore_columns = [],
        rename_columns = [],
        header = None,
        replace_header = False,
        guess_size = 1000,
        report_step = 10000,
        backend="sqlite",
        indices = [],
        missing_values = ("na", "NA", ),
        insert_quick = False,
        allow_empty = False,
        retry = False,
        utr = False
        )
    
    return parser

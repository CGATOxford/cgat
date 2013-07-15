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
check_db.py - do a sanity check for a collection of linked tables
=================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

field-in-table:  count number of entries of field in tables
fields-in-table: count number on non-null values in all fields of a table
tables-on-field: count tables

Usage
-----

Example::

   python check_db.py --help

Type::

   python check_db.py --help

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
import tempfile
import subprocess
import types
import numpy
import CGAT.Experiment as E

def output( counts, row_names, col_names, options ):

    outfile = open(options.output_filename_pattern % "counts", "w" )
    outfile.write ("\t%s\n" % "\t".join( col_names ) )
    for x, head in enumerate(row_names):
        outfile.write("%s\t%s\n" % (head, "\t".join( [ "%i" % i for i in counts[x,] ] )) )
    outfile.close()

    percents = numpy.zeros( (len(row_names), len(col_names)), numpy.float )
    for x in range(len(row_names)):
        m = max(counts[x,:])
        if m != 0:
            for y in range(len(col_names)):
                percents[x,y] = 100.0 * float(counts[x,y]) / m

    def __output( filename, matrix, pattern ):

        outfile = open(filename, "w" )
        outfile.write ("head\t%s\n" % "\t".join( col_names ) )
        for x, head in enumerate(row_names):
            outfile.write("%s\t%s\n" % (head, "\t".join( [ pattern % i for i in matrix[x,] ] )) )
        outfile.close()

    __output(options.output_filename_pattern % "counts", counts, "%i" )
    __output(options.output_filename_pattern % "percents", percents, "%5.2f" )


def doColumnInTables( dbhandle, error, options ):

    if options.backend == "pg":
        pass
    elif options.backend == "sqlite":
        cc = dbhandle.cursor()
        cc.execute("select tbl_name, sql from sqlite_master where type='table'")
        # tables returned are unicode, convert to string and filter
        # for presence of the join column 
        tag = " %s " % options.column
        tables = [ str(x[0]) for x in cc.fetchall() if tag in x[1] ]
        cc.close()

    E.debug( "tablenames from sql=%s" % str(tables) )

    if options.exclude_pattern:
        rx = re.compile( options.exclude_pattern )
        tables = [ x for x in tables if not re.search( rx, x ) ]
        E.debug( "tablenames after filtering=%s" % str(tables) )

    reg_head = re.compile(options.regex_head)
    reg_tail = re.compile(options.regex_tail)
    heads,tails = set(), set()
    for x in tables:
        try: heads.add( reg_head.search(x).groups()[0] )
        except AttributeError: pass
        try: tails.add( reg_tail.search(x).groups()[0] )
        except AttributeError: pass

    heads, tails = list(heads), list(tails)
    heads.sort()
    tails.sort()
    counts = numpy.zeros( (len(heads), len(tails)), numpy.int )

    map_head2index = dict( [(x[1],x[0]) for x in enumerate(heads)  ] )
    map_tail2index = dict( [(x[1],x[0]) for x in enumerate(tails)  ] )

    E.debug( "map_head2index=%s" % str(map_head2index) )
    E.debug( "map_tail2index=%s" % str(map_tail2index) )

    statement = "SELECT COUNT( DISTINCT %s) FROM %%s" % options.column

    E.info( "counting %s tables " % len(tables) )

    for table in tables:

        try: x = map_head2index[ reg_head.search(table).groups()[0] ]
        except AttributeError: continue
        try: y = map_tail2index[ reg_tail.search(table).groups()[0] ]
        except AttributeError: continue

        cc = dbhandle.cursor()
        cc.execute( statement % table )
        count = cc.fetchone()[0]
        cc.close()

        E.debug( "%s: %i: %s" % (table, count, statement % table) )
        counts[x,y]=count

    return counts, heads, tails

def doColumnsInTable( dbhandle, error, options ):

    rx = re.compile( options.table )

    if options.backend == "pg":
        pass
    elif options.backend == "sqlite":
        cc = dbhandle.cursor()
        cc.execute("select tbl_name FROM sqlite_master where type='table'")
        # tables returned are unicode, convert to string and filter
        # for presence of the join column 
        tables = [ str(x[0]) for x in cc.fetchall() if rx.search( x[0] ) ]
        cc.close()

        columns = set()
        for table in tables:
            cc = dbhandle.cursor()
            columns.update( set( [str(x[1]) for x in cc.execute("PRAGMA table_info( %s )" % table).fetchall()] ) )
            cc.close()

    tables, columns = list(tables), list(columns)
    tables.sort()
    columns.sort()

    E.debug( "tables=%s" % str(tables) )
    E.debug( "columns=%s" % str(columns) )

    counts = numpy.zeros( (len(tables), len(columns)), numpy.int )

    statement = "SELECT COUNT( %s) FROM %s WHERE %s IS NOT NULL" 

    E.info( "counting %i columns in %i tables " % (len(columns), len(tables) ))
    
    for x, table in enumerate(tables):

        for y, column in enumerate(columns):
            cc = dbhandle.cursor()
            try:
                cc.execute( statement % (column,table,column) )
                count = cc.fetchone()[0]
                cc.close()
            except error, msg:
                continue
            counts[x,y]=count

    new_tables = []
    rx = re.compile( options.regex_head )
    for table in tables:
        try:
            new_tables.append( rx.search( table ).groups()[0])
        except ValueError:
            new_tables.append( table )
            
    tables = new_tables

    return counts, tables, columns 

def main():

    parser = E.OptionParser( version = "%prog version: $Id: check_db.py 2782 2009-09-10 11:40:29Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-d", "--database", dest="database", type="string",
                      help="database name [default=%default]." )

    parser.add_option("-t", "--table", dest="table", type="string",
                      help="table name [default=%default]." )

    parser.add_option("-c", "--column", dest="column", type="string",
                      help="column name [default=%default]." )

    parser.add_option( "--regex-head", dest="regex_head", type="string",
                      help="regular expression to extract the head of a table, which will be the rows [default=%default]." )

    parser.add_option( "--regex-tail", dest="regex_tail", type="string",
                      help="regular expression to extract the tail of a table, which be the columns [default=%default]." )

    parser.add_option("-b", "--backend", dest="backend", type="choice",
                      choices=("pg", "sqlite"),
                      help="database backend to choose [default=%default]." )

    parser.add_option( "-p", "--output-filename-pattern", dest="output_filename_pattern", type="string" ,
                       help="pattern for output filenames [%default].")

    parser.add_option( "-e", "--exclude", dest="exclude_pattern", type="string" ,
                       help="pattern to exclude tables from the analysis [%default].")

    parser.add_option( "-m", "--method", dest="method", type="choice",
                      choices=("column-in-tables", "columns-in-table", "tables-on-column" ),
                      help="analysis to perform [default=%default]." )

    parser.set_defaults(
        database = "csvdb",
        backend="sqlite",
        column = None,
        table = None,
        regex_head = "^([^_]+)_",
        regex_tail = "_(.*)",
        output_filename_pattern = "%s",
        method = "column-in-tables",
        exclude_pattern = None,
        )

    (options, args) = E.Start( parser, add_psql_options = True )

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
        error = sqlite3.OperationalError
        options.insert_many = True  # False
        options.null = None # "NULL" 
        options.string_value = "%s" # "'%s'"

    if options.method == "column-in-tables":
        if not options.column:
            raise ValueError( "please supply a column name or pattern" )

        counts, row_names, column_names = doColumnInTables( dbhandle, error, options )

    elif options.method == "columns-in-table":
        if not options.table:
            raise ValueError( "please supply a table name or pattern" )
        counts, row_names, column_names = doColumnsInTable( dbhandle, error, options )

    output( counts, row_names, column_names, options )

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

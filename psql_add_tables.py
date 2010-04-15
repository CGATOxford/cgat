
################################################################################
#   Gene prediction pipeline 
#
#   $Id: psql_add_tables.py 2781 2009-09-10 11:33:14Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
import os, sys, string, re, getopt, time, sets, optparse, math, tempfile, copy

import pgdb

import Experiment

USAGE="""python psql_add_tables.py src dest field"""


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

    parser = optparse.OptionParser( version = "%prog version: $Id: psql_add_tables.py 2781 2009-09-10 11:33:14Z andreas $", usage=USAGE)

    parser.add_option( "-s", "--start", dest="start", type="int",
                       help="index to start.")

    parser.add_option( "-f", "--fields", dest="fields", type="string",
                       action="append",
                       help="field with index.")

    parser.set_defaults( \
        start = None,
        fields = [] )

    (options, args) = Experiment.Start( parser, add_psql_options = True )
    
    dbhandle = pgdb.connect( options.psql_connection )

    if len(args) != 2:
        raise "please supply src and dest."
    
    src_table, dest_table = args

    if len(options.fields) == 0:
        raise "please supply a field."

    statement = "SELECT %s FROM %s" % (options.fields[0], dest_table )

    indices = map(lambda x: int(x[0]), DbExecute( dbhandle, statement ))
    
    indices.sort()
    max_index = indices[-1]

    if options.start:
        if options.start < max_index:
            raise "start %i is less than maximum index %i" % (options.start, max_index)
        increment = options.start
    else:
        increment = max_index

    for field in options.fields:
        statement = "UPDATE %s SET %s=CAST(%s AS INT)+%i" % (src_table, field, field, increment )
        options.stdlog.write("# incrementing field %s in %s by %i\n" % (field, src_table, increment ) )
        DbDo( dbhandle, statement )

    statement = "INSERT INTO %s SELECT * FROM %s" % (dest_table, src_table)
    options.stdlog.write("# adding table %s to %s\n" % (src_table, dest_table) )
    DbDo( dbhandle, statement )
    
    Experiment.Stop()

    

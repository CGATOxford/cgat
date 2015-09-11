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
'''Database.py - Database utility functions
===========================================

This module contains convenience functions to work with a relational
database.


Reference
---------

'''
import time
import re


def executewait(dbhandle, statement, error=Exception, regex_error="locked",
                retries=-1, wait=5):
    '''repeatedly execute an SQL statement until it succeeds.


    Arguments
    ---------
    dbhandle : object
        A DB-API conform database handle.
    statement : string
        SQL statement to execute.
    error : string
        Exception to catch and examine for error messages.
    regex_error : string
        Any error message matching `regex_error` will be ignored,
        otherwise the procedure exists.
    retries : int
        Number of retries. If set to negative number, retry indefinitely.
        If set to 0, there will be only one attempt.
    wait : int
        Number of seconds to way between retries.

    Returns
    -------
    A cursor object

    '''
    cc = dbhandle.cursor()

    while 1:
        try:
            cc.execute(statement)
        except error, msg:
            if retries == 0:
                raise
            if not re.search("locked", str(msg)):
                raise
            time.sleep(wait)
            retries -= 1
            continue
        break
    return cc


def getColumnNames(dbhandle, table):
    """return column names of a table from a database.
    """

    cc = executewait(dbhandle, "SELECT * FROM %s LIMIT 1" % table)
    return tuple([x[0] for x in cc.description])


def getTables(dbhandle):
    """get list of tables in an sqlite database"""
    cc = executewait(
        dbhandle, """select name from sqlite_master where type='table'""")
    return tuple([x[0] for x in cc])


def toTSV(dbhandle, outfile, statement, remove_none=True):
    '''execute statement and save as tsv file
    to disk.

    If *remove_none* is true, empty/NULL values will be output as
    empty values.

    '''
    cc = dbhandle.cursor()
    cc.execute(statement)
    outfile.write("\t".join([x[0] for x in cc.description]) + "\n")

    def _str(x):
        if x is None:
            return ""
        else:
            return str(x)

    if remove_none:
        f = _str
    else:
        f = str

    outfile.write("\n".join(
        ["\t".join(map(f, x)) for x in cc]))


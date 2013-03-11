################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: PipelineGO.py 2877 2010-03-27 17:42:26Z andreas $
#
#   Copyright (C) 2013 Steve Sansom
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
"""
=============================================================
PipelineUtilities.py - helper functions for CGAT pipelines
=============================================================

:Author: Steve Sansom
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

To make excuting database commands and files IO as trivial and (TODO: robust) as possible possible. Before changing the default behaviour of any of these functions please discuss!. Better to add new functions and then deprecate.

Usage
-----

Type::

   python <script_name>.py --help

for command line help.

Code
----

"""

import sqlite3
import Pipeline as P
import rpy2.robjects as R

try:
    PARAMS = P.getParameters()
except IOError:
    pass

#######################################################################
#################### Database Queries #################################
#######################################################################

def execute(queries, database=PARAMS["database"], attach=False):
    '''Execute a list of statements sequentially'''

    dbhandle = sqlite3.connect( database )
    cc = dbhandle.cursor()

    if attach:
        for attach_statement in attach:
            cc.execute(attach_statement)

    for statement in queries: cc.execute(statement)
    cc.close()


def fetch(query, database=PARAMS["database"], attach=False):
    '''Fetch all query results and return'''
    dbhandle = sqlite3.connect( database )
    cc = dbhandle.cursor()
    if attach:
        for attach_statement in attach:
            cc.execute(attach_statement)
    sqlresult = cc.execute(query).fetchall()
    cc.close()
    return sqlresult


def fetch_with_names(query, database=PARAMS["database"], attach=False):
    '''Fetch query results and returns them as an array of row arrays, 
       in which the first entry is an array of the field names'''

    dbhandle = sqlite3.connect( database )
    cc = dbhandle.cursor()
    if attach:
        for attach_statement in attach:
            cc.execute(attach_statement)
    sqlresult = cc.execute(query).fetchall()
    data=[]
    # http://stackoverflow.com/questions/4147707/
    # python-mysqldb-sqlite-result-as-dictionary
    field_names = [ d[0] for d in cc.description ]
    data.append( [ name for name in field_names ])
    for record in sqlresult:
        line = [ field for field in record ]
        data.append(line)

    cc.close()
    return data


def write(outfile, lines, header=False):
    ''' expects [[[line1-field1],[line1-field2 ] ],... ]'''
    handle = open(outfile,"w")

    if header:
        handle.write("\t".join([str(title) for title in header])+"\n")

    for line in lines:
        handle.write("\t".join([str(field) for field in line])+"\n")

    handle.close()


###############################################################################
########################### Biomart Access ####################################
###############################################################################

def biomart_iterator( attributes,
                      host,
                      biomart, 
                      dataset,
                      filters=None,
                      values=None
                      ):
    ''' Modified from pipeline biomart... '''

    r = R.r
    r.library("biomaRt")
    mart = r.useMart( biomart,
                      dataset=dataset,
                      host=host,
                      path="/biomart/martservice",
                      archive=False )

    if filters != None: filter_names = R.StrVector(filters)
    else: filter_names = ""

    if values != None: filter_values = values
    else: filter_values = ""

    result = r.getBM( attributes=R.StrVector(attributes), 
                      filters=filter_names,
                      values=filter_values,
                      mart=mart )

    # result is a dataframe.
    # rx returns a dataframe.
    # rx()[0] returns a vector
    for data in zip( *[ result.rx(x)[0] for x in attributes] ):
        yield dict( zip(attributes, data) )

###############################################################################
########################### File Utitilies ####################################
###############################################################################

def tabToDict(filename, key=None, value=None,sep="\t"):
    ''' make a dictionary from a text file on the 
        specified column names '''

    count = 0
    result = {}
    valueidx, keyidx = False, False

    with open(filename,"r") as fh:
        for line in fh:
            if line.startswith("#"): continue
            if count==0:
                fieldn= 0
                for field in line.split(sep):
                    if field == key:
                        keyidx = fieldn
                    if field == value:
                        valueidx = fieldn
                    fieldn += 1
                if not keyidx: raise ValueError("key name not found in header")
                if not valueidx: raise ValueError("value name not found in header")
            else:
                fields = line.split(sep)
                result[fields[keyidx]] = fields[valueidx]
            count +=1

    return(result)

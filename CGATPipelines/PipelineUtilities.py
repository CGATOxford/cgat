"""=============================================================
PipelineUtilities.py - helper functions for CGAT pipelines
=============================================================

:Author: Steve Sansom
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

To make excuting database commands and files IO as trivial and (TODO:
robust) as possible possible. Before changing the default behaviour of
any of these functions please discuss!. Better to add new functions
and then deprecate.

.. note::

   This module needs refactoring as it implements
   some functionality found elsewhere such as the
   Database commands and text IO commands.

Code
----


"""
import sqlite3
try:
    import CGAT.Pipeline as P
except AttributeError, OSError:
    pass
import CGAT.IOTools as IOTools
import pickle
from pandas import DataFrame

try:
    PARAMS = P.getParameters()
except:
    PARAMS = {}


# ################### Database Queries ##############
def db_execute(cc, statements):
    '''excute a statement or statements against a cursor'''

    if type(statements) not in (list, tuple):
        statements = [statements]

    for statement in statements:
        cc.execute(statement)


def execute(queries, database=PARAMS.get("database", ""), attach=False):
    '''Execute a statement or a  list of statements (sequentially)'''

    dbhandle = sqlite3.connect(database)
    cc = dbhandle.cursor()

    if attach:
        db_execute(cc, attach)

    db_execute(cc, statement)
    cc.close()


def fetch(query, database=PARAMS.get("database", ""), attach=False):
    '''Fetch all query results and return'''

    try:
        cc = database.cursor()
    except AttributeError:

        dbhandle = sqlite3.connect(database)
        cc = dbhandle.cursor()

    if attach:
        db_execute(cc, attach)

    sqlresult = cc.execute(query).fetchall()
    cc.close()
    return sqlresult


def fetch_with_names(query,
                     database=PARAMS.get("database", ""),
                     attach=False):
    '''Fetch query results and returns them as an array of row arrays, in
       which the first entry is an array of the field names

    '''

    try:
        cc = database.cursor()
    except AttributeError:
        dbhandle = sqlite3.connect(database)
        cc = dbhandle.cursor()

    if attach:
        db_execute(cc, attach)

    sqlresult = cc.execute(query).fetchall()
    data = []
    # http://stackoverflow.com/questions/4147707/
    # python-mysqldb-sqlite-result-as-dictionary
    field_names = [d[0] for d in cc.description]
    data.append([name for name in field_names])
    for record in sqlresult:
        line = [field for field in record]
        data.append(line)

    cc.close()
    return data


def fetch_DataFrame(query,
                    database=PARAMS.get("database", ""),
                    attach=False):
    '''Fetch query results and returns them as a pandas dataframe'''

    try:
        cc = database.cursor()
    except AttributeError:
        dbhandle = sqlite3.connect(database)
        cc = dbhandle.cursor()

    if attach:
        if type(attach) is str:
            db_execute(cc, attach)
        elif isinstance(attach, (tuple, list)):
            for attach_statement in attach:
                db_execute(cc, attach_statement)

    sqlresult = cc.execute(query).fetchall()
    cc.close()

    # see http://pandas.pydata.org/pandas-docs/dev/generated/
    # pandas.DataFrame.from_records.html#pandas.DataFrame.from_records
    # this method is design to handle sql_records with proper type
    # conversion

    field_names = [d[0] for d in cc.description]
    pandas_DataFrame = DataFrame.from_records(
        sqlresult,
        columns=field_names)
    return pandas_DataFrame


def write_DataFrame(dataframe,
                    tablename,
                    database=PARAMS.get("database", ""),
                    index=False,
                    if_exists='replace'):
    '''write a pandas dataframe to an sqlite db, index on given columns
       index columns given as a string or list eg. "gene_id" or
       ["gene_id", "start"]

    '''

    dbhandle = sqlite3.connect(database)

    dataframe.to_sql(tablename,
                     con=dbhandle,
                     flavor='sqlite',
                     if_exists=if_exists)

    def indexStat(tablename, column):
        istat = ('create index %(tablename)s_%(column)s '
                 'on %(tablename)s(%(column)s)') % locals()
        return istat

    if index:

        cc = dbhandle.cursor()

        if type(index) is str:
            istat = indexStat(tablename, index)
            print istat
            db_execute(cc, istat)
        elif isinstance(index, (tuple, list)):
            for column in index:
                istat = indexStat(tablename, column)
                db_execute(cc, istat)

        cc.close()


def write(outfile, lines, header=False):
    ''' expects [[[line1-field1],[line1-field2 ] ],... ]'''
    handle = IOTools.openFile(outfile, "w")

    if header:
        handle.write("\t".join([str(title) for title in header]) + "\n")

    for line in lines:
        handle.write("\t".join([str(field) for field in line]) + "\n")

    handle.close()

# ########################## Biomart Access


def biomart_iterator(attributes,
                     host,
                     biomart,
                     dataset,
                     filters=None,
                     values=None
                     ):
    ''' Modified from pipeline biomart... '''

    # Use PipelineBiomant.biomart_iterator
    raise NotImplementedError(
        "deprecated, use PipelineBiomart.biomart_iterator")

# ########################## File Utitilies


def txtToDict(filename, key=None, sep="\t"):
    ''' make a dictionary from a text file keyed
        on the specified column '''

    # Please see function in IOTools: readDict()
    count = 0
    result = {}
    valueidx, keyidx = False, False
    field_names = []

    with open(filename, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if count == 0:
                fieldn = 0
                for rawfield in line.split(sep):
                    field = rawfield.strip()
                    if field == key:
                        keyidx = fieldn
                    field_names.append(field)
                    fieldn += 1

                if not keyidx:
                    raise ValueError("key name not found in header")
                # if not valueidx:
                #   raise ValueError(
                #     "value name not found in header")
            else:
                fields = [x.strip() for x in line.split(sep)]
                fieldn = 0
                thiskey = fields[keyidx]
                result[thiskey] = {}
                for field in fields:
                    if fieldn == keyidx:
                        pass
                    else:
                        colkey = field_names[fieldn]
                        result[thiskey][colkey] = field
                    fieldn += 1
            count += 1

    return(result)

# ############################## Objects


def save(file_name, obj):
    '''dump a python object to a file using pickle'''
    with open(file_name, "wb") as pkl_file:
        pickle.dump(obj, pkl_file)
    return


def load(file_name):
    '''retrieve a pickled python object from a file'''
    with open(file_name, "r") as pkl_file:
        data = pickle.load(pkl_file)
    return data

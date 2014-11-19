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

Command line options
--------------------

'''
import os
import sys
import subprocess

import psycopg2

import CGAT.Experiment as E

# ------------------------------------------------------------------------


def DbExecute(dbhandle, statement):

    cc = dbhandle.cursor()
    cc.execute(statement)
    result = cc.fetchall()
    cc.close()
    return result

# ------------------------------------------------------------------------


def DbDo(dbhandle, statement):

    cc = dbhandle.cursor()
    cc.execute(statement)
    dbhandle.commit()
    cc.close()

# ------------------------------------------------------------------------


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id",
        usage=globals()["__doc__"])

    parser.set_defaults(
        start=None,
        fields=[])

    (options, args) = E.Start(parser, add_psql_options=True)

    dbhandle = psycopg2.connect(options.psql_connection)

    if len(args) != 2:
        raise ValueError("please supply src and dest.")

    src, dest = args

    sql_get_schemas = "SELECT DISTINCT schemaname FROM pg_catalog.pg_stat_all_tables"
    sql_create_schema = "CREATE SCHEMA %s"

    cc = dbhandle.cursor()

    cc.execute(sql_get_schemas)
    databases = map(lambda x: x[0], cc.fetchall())
    cc.close()

    if src not in databases:
        raise ValueError("unknown schema %s" % src)

    if dest in databases:
        raise ValueError("schema %s already exists" % (dest))

    host, database = options.psql_connection.split(":")

    cmd = 'pg_dump --ignore-version --schema-only --host=%(host)s --schema=%(src)s | sed "s/%(src)s/%(dest)s/g" | psql --host=%(host)s %(database)s' % (locals())
    p = subprocess.Popen(cmd,
                         shell=True)
    sts = os.waitpid(p.pid, 0)

    cmd = 'pg_dump --ignore-version --data-only --host=%(host)s --schema=%(src)s | psql --host=%(host)s %(database)s' % (
        locals())
    p = subprocess.Popen(cmd,
                         shell=True)
    sts = os.waitpid(p.pid, 0)

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

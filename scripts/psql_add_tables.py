'''
psql_add_tables.py - merge two tables
=====================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script will merge two tables by adding one to the other.
An index ``field`` is automatically incremented such that there
is no overlap between the two tables.

The column ``field`` will not necessary be continuous.

Usage
-----

Example::

   python psql_add_tables.py src dest field

will add table ``src`` to table ``dest`` using ``field``
as an incremental counter. The column ``field`` has to
present in both tables.

Type::

   python psql_add_tables.py --help

for command line help.

Command line options
--------------------

'''
import sys
import psycopg2

import CGAT.Experiment as E

USAGE = """python psql_add_tables.py src dest field"""


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
        version="%prog version: $Id: psql_add_tables.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-s", "--start", dest="start", type="int",
                      help="index to start.")

    parser.add_option("-f", "--fields", dest="fields", type="string",
                      action="append",
                      help="field with index.")

    parser.set_defaults(
        start=None,
        fields=[])

    (options, args) = E.Start(parser, add_psql_options=True)

    dbhandle = psycopg2.connect(options.psql_connection)

    if len(args) != 2:
        raise "please supply src and dest."

    src_table, dest_table = args

    if len(options.fields) == 0:
        raise "please supply a field."

    statement = "SELECT %s FROM %s" % (options.fields[0], dest_table)

    indices = map(lambda x: int(x[0]), DbExecute(dbhandle, statement))

    indices.sort()
    max_index = indices[-1]

    if options.start:
        if options.start < max_index:
            raise "start %i is less than maximum index %i" % (
                options.start, max_index)
        increment = options.start
    else:
        increment = max_index

    for field in options.fields:
        statement = "UPDATE %s SET %s=CAST(%s AS INT)+%i" % (src_table,
                                                             field, field, increment)
        options.stdlog.write(
            "# incrementing field %s in %s by %i\n" % (field, src_table, increment))
        DbDo(dbhandle, statement)

    statement = "INSERT INTO %s SELECT * FROM %s" % (dest_table, src_table)
    options.stdlog.write("# adding table %s to %s\n" % (src_table, dest_table))
    DbDo(dbhandle, statement)

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

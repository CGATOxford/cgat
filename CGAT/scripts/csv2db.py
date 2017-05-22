'''csv2db.py - upload table to database
====================================

:Tags: Python

Purpose
-------

create a table from a csv separated file and load data into it.

This module supports backends for postgres and sqlite3. Column types are
auto-detected.

Read a table from stdin and create an sqlite3 database. By default,
the database will reside in a file called csvdb and in a table csv.

.. todo::

   Use file import where appropriate to speed up loading. Currently, this is
   not always the case.

Usage
-----

Example::

   python csv2db.py -b sqlite < stdin 

Type::

   python csv2db.py --help

for command line help.

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.CSV2DB as CSV2DB

import csv

csv.field_size_limit(sys.maxsize)


def main(argv=sys.argv):

    parser = CSV2DB.buildParser()

    (options, args) = E.Start(parser, argv=argv,
                              add_database_options=True)

    if options.from_zipped:
        import gzip
        infile = gzip.GzipFile(fileobj=options.stdin, mode='r')

    else:
        infile = options.stdin

    CSV2DB.run(infile, options)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

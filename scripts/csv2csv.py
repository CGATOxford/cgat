'''
csv2csv.py - operate on tables
==============================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

operate on tables.

Usage
-----

Example::

   python csv2csv.py --help

Type::

   python csv2csv.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import string
import re
import getopt
import time
import optparse
import math
import tempfile


"""do some table magic.

--sort: sort table 


"""
import CGAT.Experiment as E
import csv

parser = E.OptionParser( version = "%prog version: $Id: csv2csv.py 2782 2009-09-10 11:40:29Z andreas $")

def ConvertDictionary( d ):
    """tries to convert values in a dictionary.
    """

    rx_int = re.compile("^[+-]*[0-9]+$")
    rx_float = re.compile("^[+-]*[0-9.+-eE]+$")
    for k,v in d.items():
        if rx_int.match( v ):
            d[k] = int(v)
        elif rx_float.match( v ):
            d[k] = float(v)

    return d
    
if __name__ == "__main__":

    parser.add_option( "-s", "--sort", dest="sort", type="string" ,
                       help="fields to take (in sorted order).")

    (options, args) = E.Start( parser, add_csv_options  = True)

    reader = csv.DictReader( sys.stdin, dialect=options.csv_dialect )

    if options.sort:
        fields = options.sort.split(",")
    else:
        fields = None
        
    writer = csv.DictWriter( sys.stdout,
                             fields,
                             dialect=options.csv_dialect,
                             lineterminator = options.csv_lineterminator,
                             extrasaction = 'ignore' )

    print "\t".join(fields)
    
    for row in reader:
        row = ConvertDictionary( row )
        writer.writerow(row)
        
    E.Stop()

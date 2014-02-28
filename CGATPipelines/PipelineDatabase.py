'''
PipelineDatabase.py - utility functions for working with a database
===================================================================

'''

import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections

import CGAT.Experiment as E
import CGAT.Pipeline as P
import CGAT.Stats as Stats

try:
    PARAMS = P.getParameters()
except IOError:
    pass


def importFromIterator(
        outfile,
        tablename,
        iterator,
        columns=None,
        indices=None):
    '''import data in *iterator* into *tablename* via temporary file.

    '''

    tmpfile = P.getTempFile(".")

    if columns:
        keys, values = zip(*columns.items())
        tmpfile.write("\t".join(values) + "\n")

    for row in iterator:
        if not columns:
            keys = row[0].keys()
            values = keys
            columns = keys
            tmpfile.write("\t".join(values) + "\n")

        tmpfile.write("\t".join(str(row[x]) for x in keys) + "\n")

    tmpfile.close()

    if indices:
        indices = " ".join("--index=%s" % x for x in indices)
    else:
        indices = ""

    tmpfilename = tmpfile.name

    statement = '''
       python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                     --table=%(tablename)s
                     %(indices)s
        < %(tmpfilename)s > %(outfile)s
    '''

    P.run()

    os.unlink(tmpfilename)

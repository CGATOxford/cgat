'''
cgat_pep8_check_code_quality.py - check PEP8 conformance of CGAT Code
=====================================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script runs pep8.py on the CGAT code collection and outputs
summary statistics of code quality onto stdout.

Usage
-----

To use, simply run the script from the root directory of the
CGAT code collection::

   python cgat_pep8_check_code_quality.py

Type::

   python cgat_pep8_check_code_quality.py --help

for command line help.

Command line options
--------------------


'''

import subprocess
import collections
import re
import sys
import CGAT.Experiment as E

DATA = collections.namedtuple('DATA', 'count code description')

expressions = (
    ('tests', 'tests/*.py'),
    ('scripts', 'scripts/*.py'),
    ('optic', 'scripts/optic/*.py'),
    ('gpipe', 'scripts/gpipe/*.py'),
    ('CGAT', 'CGAT/*.py'),
    ('CGATPipelines', 'CGATPipelines/*.py'),
    ('trackers', 'CGATPipelines/pipeline_docs/*/trackers/*.py'))


def runPep8(expr):
    '''run pep8 on files matching the glob expression
    *expr*
    '''

    args = ['/ifs/devel/andreas/python/bin/pep8',
            '--statistics',
            '--quiet',
            expr]

    args = ' '.join(args)
    # pep8 returns unix error if there are errors
    try:
        output = subprocess.check_output(args, shell=True)
    except subprocess.CalledProcessError, msg:
        output = msg.output

    # count python files
    nchecked = len([x for x in output.split('\n') if x.endswith('.py')])

    # parse data
    # Output is '120 W123 Description'
    lines = [x for x in output.split('\n') if not x.endswith('.py')]
    data = [DATA(*re.split(' +', x, 2)) for x in lines if x]

    return nchecked, data


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    rows = []
    labels = {}
    for label, expr in expressions:
        nchecked, data = runPep8(expr)
        rows.append((label, nchecked, data))
        labels.update(dict([(x.code, x.description) for x in data]))

    # build table
    #
    # each row is data set and each column is a Warning/Error type
    # with some additional columns such as total and n.

    # build dictionary mapping error codes to columns
    # consistently across samples
    map_code2column = dict([(y, x + 3) for x, y in enumerate(labels.keys())])

    # build first row containing the column labels
    results = [['code', 'n', 'total'] + labels.keys()]

    # build array with column totals
    column_totals = [0] * (len(map_code2column) + 3)
    for label, nchecked, data in rows:
        row = [label, nchecked, 0] + [0] * len(map_code2column)
        column_totals[1] += nchecked
        for x in data:
            c = map_code2column[x.code]
            row[c] = x.count
            row[2] += int(x.count)
            column_totals[2] += int(x.count)
            column_totals[c] += int(x.count)

        results.append(row)
    # add column totals
    column_totals[0] = 'total'
    results.append(column_totals)

    # add descriptions as last row
    results.append(['description',
                    'number of files checked',
                    'total errors/warnings in set'] + labels.values())

    # output transposed table
    outfile = sys.stdout
    for row in zip(*results):
        outfile.write('%s\n' % ('\t'.join(map(str, row))))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

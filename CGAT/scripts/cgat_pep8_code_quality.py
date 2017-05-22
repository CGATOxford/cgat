'''
cgat_pep8_check_code_quality.py - check PEP8 conformance of CGAT Code
=====================================================================

:Author:
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

import collections
import sys
import CGAT.Experiment as E
import CGAT.Style

DATA = collections.namedtuple('DATA', 'count code description')

expressions = (
    ('tests', 'tests/*.py'),
    ('scripts', 'scripts/*.py'),
    ('optic', 'scripts/optic/*.py'),
    ('gpipe', 'scripts/gpipe/*.py'),
    ('CGAT', 'CGAT/*.py'),
    ('CGATPipelines', 'CGATPipelines/*.py'),
    ('trackers', 'CGATPipelines/pipeline_docs/*/trackers/*.py'))


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
        nchecked, data = CGAT.Style.runPep8(expr)
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
    results = [['code', 'n', 'total'] + list(labels.keys())]

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
                    'total errors/warnings in set'] + list(labels.values()))

    # output transposed table
    outfile = sys.stdout
    for row in zip(*results):
        outfile.write('%s\n' % ('\t'.join(map(str, row))))

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

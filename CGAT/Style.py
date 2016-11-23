'''Style.py
============

utility functions for checking coding style in
CGAT scripts.

'''

import subprocess
import collections
import re

DATA = collections.namedtuple('DATA', 'count code description')


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
    except subprocess.CalledProcessError as msg:
        output = msg.output

    # count python files
    nchecked = len([x for x in output.split('\n') if x.endswith('.py')])

    # parse data
    # Output is '120 W123 Description'
    lines = [x for x in output.split('\n') if not x.endswith('.py')]
    data = [DATA(*re.split(' +', x, 2)) for x in lines if x]

    return nchecked, data

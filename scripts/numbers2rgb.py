'''
numbers2rgb.py - map numbers to RGB values
==========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Reads numbers from stdin and maps to random RGB values.

Usage
-----

Example::

   python numbers2rgb.py --help

Type::

   python numbers2rgb.py --help

for command line help.

Command line options
--------------------

'''
import sys
import random

import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: numbers2rgb.py 2782 2009-09-10 11:40:29Z andreas $")

    parser.add_option("-m", "--method", dest="method", type="string",
                      help="method to use.")

    parser.set_defaults(
        method="random_rgb",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    numbers = map(lambda x: int(x[:-1].split("\t")[0]),
                  filter(lambda x: x[0] != "#", sys.stdin.readlines()))

    if options.method == "random_rgb":
        f = lambda x: "%i,%i,%i" % (random.randint(0, 256),
                                    random.randint(0, 256),
                                    random.randint(0, 256))

    map_number2output = {}
    for x in numbers:
        if x not in map_number2output:
            map_number2output[x] = f(x)

    for x, y in map_number2output.items():
        options.stdout.write("%i\t%s\n" % (x, y))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

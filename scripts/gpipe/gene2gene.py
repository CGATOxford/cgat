##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
##########################################################################
'''
gpipe/gene2gene.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python gpipe/gene2gene.py --help

Type::

   python gpipe/gene2gene.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import CGAT.Experiment as E

USAGE = """python %s [OPTIONS] < gene_list > graph

print list of all transcripts within a gene.
""" % sys.argv[0]

# add links between genes


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/gene2gene.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-q", "--restrict-quality", dest="restrict_quality", type="string",
                      help="restrict genes to given quality codes.")

    parser.set_defaults(separator="|",
                        restrict_quality=None)

    options, args = E.Start(parser)

    if options.restrict_quality:
        options.restrict_quality = set(options.restrict_quality.split(","))

    ninput, noutput, nskipped, nerrors = 0, 0, 0, 0

    def print_lines(lines):
        global noutput
        if not lines:
            return
        for x in range(len(lines) - 1):
            for y in range(x + 1, len(lines)):
                options.stdout.write(options.separator.join(
                    lines[x]) + "\t" + options.separator.join(lines[y]) + "\t0\n")
                noutput += 1

    transcripts = []

    for line in sys.stdin:
        try:
            schema, prediction_id, gene_id, quality = line[
                :-1].split(options.separator)
        except ValueError:
            nerrors += 1
            if options.loglevel >= 1:
                options.stdlog.write("# PARSING ERROR in line %s" % line)
            continue
        transcripts.append((schema, prediction_id, gene_id, quality))

    transcripts.sort(lambda x, y: cmp((x[0], x[2]), (y[0], y[2])))

    last_gene_id = None
    last_schema = None
    lines = []

    ninput = len(transcripts)
    for schema, prediction_id, gene_id, quality in transcripts:
        if last_gene_id != gene_id or last_schema != schema:
            print_lines(lines)
            lines = []
            last_gene_id = gene_id
            last_schema = schema

        if options.restrict_quality and quality not in options.restrict_quality:
            nskipped += 1
            continue

        lines.append((schema, prediction_id, gene_id, quality))

    print_lines(lines)

    E.info("ninput=%i, noutput=%i, nskipped=%i, nerrors=%i" %
           (ninput, noutput, nskipped, nerrors))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

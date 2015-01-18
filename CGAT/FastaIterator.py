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
"""
FastaIterator.py - iterate over fasta files
===========================================

The difference to the biopython iterator is that this one 
skips over comment lines starting with "#".

Code
----

"""
import subprocess
import os


class FastaRecord:

    def __init__(self, title, sequence):

        self.title = title
        self.sequence = sequence


def iterate(infile, comment="#"):
    '''iterate over fasta data in infile

    Lines before the first fasta record are
    ignored (starting with '>') as well as
    lines starting with the comment character.

    yields FastaRecord for each fasta file.
    '''

    h = infile.readline()[:-1]

    if not h:
        raise StopIteration

    # skip everything until first fasta entry starts
    while h[0] != ">":
        h = infile.readline()[:-1]
        if not h:
            raise StopIteration
        continue

    h = h[1:]
    seq = []

    for line in infile:

        if line.startswith(comment):
            continue

        if line.startswith('>'):
            yield FastaRecord(h, ''.join(seq))

            h = line[1:-1]
            seq = []
            continue

        seq.append(line[:-1])

    yield FastaRecord(h, ''.join(seq))


class FastaIterator:

    '''a iterator of :term:`fasta` formatted files.
    '''

    def __init__(self, f, *args, **kwargs):
        self.mIterator = iterate(f)

    def __iter__(self):
        return self

    def next(self):
        return self.mIterator.next()

# ------------------------------------------------------------


def iterate_together(*args):
    """iterate synchronously over one or more fasta files.

    The iteration finishes once any of the files is exhausted.

    yield output tuples of sequences."""
    iterators = [FastaIterator(open(x, "r")) for x in args]

    while 1:
        yield [x.next() for x in iterators]


# ------------------------------------------------------------
def count(filename):
    '''count number of sequences in fasta file.'''
    if filename.endswith(".gz"):
        statement = "zcat %s | grep -c '>'" % filename
    else:
        statement = "cat %s | grep -c '>'" % filename

    if not os.path.exists(filename):
        raise OSError("file '%s' does not exist" % filename)

    # grep returns error if no match is found
    try:
        return subprocess.check_output(statement, shell=True)
    except subprocess.CalledProcessError:
        return 0

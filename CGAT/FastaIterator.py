################################################################################
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
#################################################################################
"""
FastaIterator.py - iterate over fasta files
===========================================

The difference to the biopython iterator is that this one 
skips over comment lines starting with "#".

Code
----

"""
import subprocess, os

class FastaRecord:

    def __init__(self, title, sequence ):

        self.title = title
        self.sequence = sequence

def iterate( infile ):

    h = infile.readline()[:-1]

    if not h: raise StopIteration
    
    while h[0] != ">":
        h = infile.readline()[:-1]
        if not h: raise StopIteration
        continue

    h = h[1:]
    
    seq = []
    
    for line in infile:

        if line[0] == "#": continue
        
        line = line[:-1] # remove newline
        if not line: continue

        if line[0] == '>':
            yield FastaRecord(h,''.join(seq))
            
            h = line[1:]
            seq = []
            continue
        
        seq.append(line)
        
    yield FastaRecord(h,''.join(seq))

class FastaIterator:
    '''a iterator of :term:`fasta` formatted files.
    '''

    def __init__(self, f, *args, **kwargs):
        self.mIterator = iterate(f)
        
    def __iter__(self):
        return self

    def next(self):
        return self.mIterator.next()

##------------------------------------------------------------
def iterate_together( *args ):
    """iterate synchronously over one or more fasta files.

    The iteration finishes once any of the files is exhausted.

    yield output tuples of sequences."""
    iterators = [ FastaIterator( open(x, "r") ) for x in args ]

    while 1:
        yield [ x.next() for x in iterators ]


##------------------------------------------------------------
def count( filename ):
    '''count number of sequences in fasta file.'''
    if filename.endswith(".gz"):
        statement = "zcat %s | grep -c '>'" % filename
    else:
        statement = "cat %s | grep -c '>'" % filename

    if not os.path.exists( filename ):
        raise OSError( "file '%s' does not exist" % filename )

    # grep returns error if no match is found
    try:
        return subprocess.check_output( statement, shell = True )
    except subprocess.CalledProcessError:
        return 0

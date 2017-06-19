"""
GFF3 - Classes, functions and iterators for working with GFF3 files
-------------------------------------------------------------------

This module mostly inherits from the :mod:`GTF` and replaces selected
functionality to allow working with :term:`GFF3` formatted files.

"""

import collections
import itertools
import pysam
from CGAT import GTF as GTF


def flat_file_iterator(infile):
    ''' simple iterator that iterators over lines in a field
    and yeilds GFF3 Entry objects '''

    return pysam.tabix_iterator(infile, parser=pysam.asGFF3())


def chrom_iterator(gff3_iterator):
    ''' takes a an iterator and returns an iterator over iterators,
    with a new instance every time a new chromosome is found '''

    seen = set()

    for chrom_set in itertools.groupby(gff3_iterator, lambda x: x.contig):

        if chrom_set[0] in seen:
            raise ValueError("File doesn't appear to be chromosome sorted")

        seen.add(chrom_set[0])

        yield chrom_set[1]

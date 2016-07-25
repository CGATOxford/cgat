"""
GFF3 - Classes, functions and iterators for working with GFF3 files
-------------------------------------------------------------------

This module mostly inherits from the :mod:`GTF` and replaces selected
functionality to allow working with :term:`GFF3` formatted files.

"""

import collections
import itertools
from CGAT import GTF as GTF


class Entry(GTF.Entry):
    """representation of a :term:`GFF3` formatted entry.

    This class inherits from :class:`GTF.Entry`, but changes
    the parsing to reflect GFF3.
    """

    def parseInfo(self, attributes, line=None):
        ''' Parse the attributes line of an entry,
        line parameter provided purely for backwards compatability'''

        # remove comments
        attributes = attributes.split("#")[0]
        # separate into fields
        fields = [x.strip() for x in attributes.split(";")]
        self.attributes = collections.OrderedDict()

        for f in fields:

            d = [x.strip() for x in f.split("=")]

            n, v = d[0], d[1]
            if len(d) > 2:
                v = d[1:]

            if v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            else:
                # try to convert to a value
                try:
                    v = float(v)
                    v = int(v)
                except ValueError:
                    pass
                except TypeError:
                    pass

            # The reversed Parent attribute can contain multiple, "," sepearted
            # values
            if n == "Parent":
                v = v.split(",")

            self.attributes[n] = v

    def getAttributeField(self):
        ''' return the attributes field as a ; delimied field '''

        aa = []
        for k, v in list(self.attributes.items()):

            if k == "Parent":
                v = ",".join(v)

            aa.append('%s=%s' % (k, str(v)))

        return ";".join(aa)


def flat_file_iterator(infile):
    ''' simple iterator that iterators over lines in a field
    and yeilds GFF3 Entry objects '''

    for line in infile:

        if line.startswith("##FASTA"):
            raise StopIteration

        if line.startswith("#"):
            continue

        if len(line.strip()) == 0:
            continue

        gff_entry = Entry()
        gff_entry.read(line)
        yield gff_entry


def iterator_from_gff(gff_iterator):
    ''' to make this slot in with other gtf using scripts,
    allow copying of an entry into gff3 format. Acts via str,
    probably not the most efficient way to do things'''

    for gff in gff_iterator:

        result = Entry()
        result.read(str(gff))
        yield result


def chrom_iterator(gff3_iterator):
    ''' takes a an iterator and returns an iterator over iterators,
    with a new instance every time a new chromosome is found '''

    seen = set()

    for chrom_set in itertools.groupby(gff3_iterator, lambda x: x.contig):

        if chrom_set[0] in seen:
            raise ValueError("File doesn't appear to be chromosome sorted")

        seen.add(chrom_set[0])

        yield chrom_set[1]

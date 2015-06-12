'''
analyze_sequences.py - compute some sequence properties
=======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

compute some sequence properties.

Usage
-----

Example::

   python analyze_sequences.py --help

Type::

   python analyze_sequences.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re

import CGAT.Experiment as E
import CGAT.Genomics as Genomics

import CGAT.FastaIterator as FastaIterator
import Bio.Alphabet.IUPAC


class SequenceProperties:

    def __init__(self):
        pass

    def Load(self, in_sequence):
        """load sequence properties from a sequence."""

        # uppercase all letters
        sequence = in_sequence.upper()

        self.mNCodons = len(sequence) / 3

        self.mNStopCodons = 0

        # setup counting arrays
        # counts of amino acids
        self.mCountsAA = {}

        for x in Bio.Alphabet.IUPAC.extended_protein.letters:
            self.mCountsAA[x] = 0

        # nucleotide counts for each position (is not a sum of the counts
        # per degenerate site, as the codon might be intelligible, e.g. GNN).
        self.mCounts = [{'A': 0, 'C': 0, 'G': 0, 'T': 0, 'X': 0, 'N': 0},
                        {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'X': 0, 'N': 0},
                        {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'X': 0, 'N': 0}]

        # nucleotide counts for each position per degeneracy
        self.mCountsDegeneracy = []

        self.mLength = len(sequence)

        for x in (0, 1, 2):
            xx = []
            for y in range(5):
                yy = {}
                for z in Bio.Alphabet.IUPAC.extended_dna.letters:
                    yy[z] = 0
                xx.append(yy)
            self.mCountsDegeneracy.append(xx)

        for codon in [sequence[x:x + 3] for x in range(0, len(sequence), 3)]:

            for x in (0, 1, 2):
                self.mCounts[x][codon[x]] += 1

            if Genomics.IsStopCodon(codon):
                self.mNStopCodons += 1
                continue

            try:
                aa, deg1, deg2, deg3 = Genomics.GetDegeneracy(codon)
                degrees = (deg1, deg2, deg3)
                for x in range(len(degrees)):
                    self.mCountsDegeneracy[x][degrees[x]][codon[x]] += 1
                self.mCountsAA[aa] += 1

            except KeyError:
                pass

        self.Update()

    def Update(self):
        """update fields from counts."""

        self.mNGC = 0
        self.mNSites1D, self.mNSites2D, self.mNSites3D, self.mNSites4D = 0, 0, 0, 0

        for x in (0, 1, 2):
            self.mNGC += self.mCounts[x]['C'] + self.mCounts[x]['G']
            self.mNSites1D += sum(self.mCountsDegeneracy[x][1].values())
            self.mNSites2D += sum(self.mCountsDegeneracy[x][2].values())
            self.mNSites3D += sum(self.mCountsDegeneracy[x][3].values())
            self.mNSites4D += sum(self.mCountsDegeneracy[x][4].values())

        self.mNGC3 = self.mCounts[2]['C'] + self.mCounts[2]['G']
        self.mN4DGC3 = self.mCountsDegeneracy[2][4][
            'C'] + self.mCountsDegeneracy[2][4]['G']

    def __str__(self):

        if self.mLength > 0:
            pngc = float(self.mNGC) / float(self.mLength)
        else:
            pngc = 0.0

        if self.mNCodons > 0:
            pngc3 = float(self.mNGC3) / float(self.mNCodons)
        else:
            pngc3 = 0.0

        if self.mNSites4D > 0:
            pn4dgc3 = float(self.mN4DGC3) / float(self.mNSites4D)
        else:
            pn4dgc3 = 0.0

        return "\t".join(("%i" % self.mLength,
                          "%i" % self.mNCodons,
                          "%i" % self.mNStopCodons,
                          "%i" % self.mNSites1D,
                          "%i" % self.mNSites2D,
                          "%i" % self.mNSites3D,
                          "%i" % self.mNSites4D,
                          "%i" % self.mNGC,
                          "%i" % self.mNGC3,
                          "%i" % self.mN4DGC3,
                          "%5.4f" % pngc,
                          "%5.4f" % pngc3,
                          "%5.4f" % pn4dgc3,
                          ))

    def GetHeader(self):
        return "\t".join(("length",
                          "ncodons",
                          "nstops",
                          "nsites1d",
                          "nsites2d",
                          "nsites3d",
                          "nsites4d",
                          "ngc",
                          "ngc3",
                          "n4dgc3",
                          "pgc",
                          "pgc3",
                          "p4dgc3",
                          ))


def main(argv=None):

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-s", "--species", dest="species", type="string",
                      help="species to use.")

    parser.add_option("-p", "--column-prefix", dest="prefix", type="string",
                      help="prefix to use for temporary files.")

    parser.set_defaults(
        filename_map=None,
    )

    (options, args) = E.Start(parser, add_database_options=True)

    iterator = FastaIterator.FastaIterator(options.stdin)

    print "id\t" + SequenceProperties().GetHeader()

    while 1:
        cur_record = iterator.next()

        if cur_record is None:
            break

        sequence = re.sub(" ", "", cur_record.sequence)

        if len(sequence) % 3:
            raise "sequence %s is not multiples of 3: length=%i!" % (
                cur_record.title, len(sequence))

        s = SequenceProperties()
        s.Load(sequence)

        print cur_record.title + "\t" + str(s)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

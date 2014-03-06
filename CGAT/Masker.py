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
Masker.py - Wrapper for sequence masking tools
==============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os
import subprocess
import tempfile
import string
import re
import random

from CGAT import Experiment as E
from CGAT import Genomics as Genomics
from CGAT import FastaIterator as FastaIterator
import cStringIO as StringIO

# Class for calling masking programs.


class Masker:

    """a masker preserves gaps, but it does not preserve
    whitespace characters."""

    mExecutable = None
    mOptions = ""
    mHasPeptideMasking = False
    mHasNucleicAcidMasking = False

    # set to true if masker outputs softmasked sequence
    soft_mask = False

    def __init__(self):
        pass

    def getAlphabet(self, sequence):
        """get sequence type (aa,na,codons)."""

        s1 = re.sub("[acgtxn\-.]", "", sequence.lower())
        s2 = re.sub("[-.]", "", sequence.lower())
        if float(len(s1)) < (len(s2) * 0.1):
            alphabet = "na"
            if len(sequence) % 3 == 0:
                alphabet = "codons"
        else:
            alphabet = "aa"

        return alphabet

    def __call__(self, sequence):
        """mask a sequence."""

        sequence = re.sub("\s", "", sequence)

        a = self.getAlphabet(sequence)

        seq = list(sequence)

        if len(seq) < 5:
            # do not mask empty/short sequences
            pass

        elif a == "aa" and self.mHasPeptideMasking:

            c = 0
            m = self.maskSequence(sequence)
            if self.soft_mask:
                m = re.sub("[a-z]", "x", m)
            for p, m in zip(sequence, m):
                if m in "Xx":
                    if p.isupper():
                        seq[c] = "X"
                    else:
                        seq[c] = "x"
                c += 1

        elif a == "codons" and self.mHasPeptideMasking:

            peptide_sequence = Genomics.TranslateDNA2Protein(sequence)
            masked_sequence = self.maskSequence(peptide_sequence)
            if self.soft_mask:
                masked_sequence = re.sub("[a-z]", "x", masked_sequence)

            c = 0
            for p, m in zip(peptide_sequence, masked_sequence):
                if m in "Xx":
                    if p.isupper():
                        seq[c:c + 3] = ["N"] * 3
                    else:
                        seq[c:c + 3] = ["n"] * 3
                c += 3

        elif a in ("na", "codons") and self.mHasNucleicAcidMasking:
            masked_sequence = self.maskSequence(sequence)
            if self.soft_mask:
                masked_sequence = re.sub("[a-z]", "N", masked_sequence)
            return masked_sequence
        else:
            raise ValueError(
                "masking of sequence type %s not implemented." % a)

        return "".join(seq)

    def maskSequence(self, peptide_sequence):
        """mask peptide sequence
        """

        Masker.__init__(self)

        outfile, filename_peptide = tempfile.mkstemp()
        os.write(outfile, ">test\n%s\n" % (peptide_sequence))
        os.close(outfile)

        infile = filename_peptide
        statement = self.mCommand % locals()

        E.debug("statement: %s" % statement)

        s = subprocess.Popen(statement,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             close_fds=True)

        (out, err) = s.communicate()
        if s.returncode != 0:
            raise RuntimeError(
                "Error in running %s \n%s\nTemporary directory" %
                (statement, err))

        os.remove(filename_peptide)

        masked_sequence = re.sub(
            "\s", "", string.join(out.split("\n")[1:], ""))
        return masked_sequence

    def maskSequences(self, sequences):
        '''mask a collection of sequences.'''

        outfile, infile = tempfile.mkstemp()

        for x, s in enumerate(sequences):
            os.write(outfile, ">%i\n%s\n" % (x, s))

        os.close(outfile)

        statement = self.mCommand % locals()

        E.debug("statement: %s" % statement)

        s = subprocess.Popen(statement,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             close_fds=True)

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise RuntimeError(
                "Error in running %s \n%s\nTemporary directory" %
                (statement, err))

        result = [
            x.sequence for x in FastaIterator.iterate(StringIO.StringIO(out))]

        os.remove(infile)

        return result


class MaskerBias (Masker):

    mCommand = "biasdb.pl %(infile)s"
    mHasPeptideMasking = True


class MaskerSeg (Masker):
    # mCommand = "seg %(infile)s 12 2.2 2.5 -x"
    mCommand = ("segmasker -in %(infile)s -window 12 -locut 2.2 "
                "-hicut 2.5 -outfmt fasta")
    mHasPeptideMasking = True
    soft_mask = True


class MaskerDustMasker(Masker):

    '''use dustmasker. masked chars are returned as
    lower case characters.'''

    mCommand = "dustmasker -outfmt fasta -in %(infile)s"
    mHasNucleicAcidMasking = True


class MaskerRandom (Masker):
    """randomly mask a proportion of positions in a sequence
    in multiple alignment."""

    def __init__(self, proportion=10, *args, **kwargs):
        Masker.__init__(self, *args, **kwargs)
        self.mProportion = proportion / 100.0

    def __call__(self, sequence):
        """mask a sequence."""

        sequence = re.sub("\s", "", sequence)

        a = self.getAlphabet(sequence)

        if a == "codons":
            frame = 3
        else:
            frame = 1

        positions = [
            x for x in range(0, len(sequence), frame) if sequence[x] != "-"]
        to_mask = random.sample(
            positions, int(len(positions) * self.mProportion))
        print int(len(positions) * self.mProportion)

        s = list(sequence)
        print to_mask
        for x in to_mask:
            for y in range(x, x + frame):
                s[x] == "x"

        return "".join(s)

if __name__ == "__main__":
    x = MaskerRandom()
    print x("AAA AAA AAA AAA --- AAA AAA AAA AAA")

    x = MaskerDustMasker()
    print x.maskSequences(("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                           "GGGGGGGGGG", ))


def maskSequences(sequences, masker=None):
    '''return a list of masked sequence.

    *masker* can be one of
        dust/dustmasker * run dustmasker on sequences
        softmask        * use softmask to hardmask sequences
    '''

    if masker in ("dust", "dustmasker"):
        masker_object = MaskerDustMasker()
    else:
        masker_object = None

    if masker == "softmask":
        # the genome sequence is repeat soft-masked
        masked_seq = sequences
    elif masker in ("dust", "dustmasker"):
        # run dust
        masked_seq = masker_object.maskSequences(
            [x.upper() for x in sequences])
    elif masker is None:
        masked_seq = [x.upper() for x in sequences]
    else:
        raise ValueError("unknown masker %s" % masker)

    # hard mask softmasked characters
    masked_seq = [re.sub("[a-z]", "N", x) for x in masked_seq]

    return masked_seq

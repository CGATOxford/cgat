##########################################################################
#   Gene prediction pipeline
#
#   $Id: snp2counts_test.py 2855 2010-02-10 09:59:58Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
"""unit testing module for the Tree.py class."""

import sys
import os
import shutil
import optparse
import random
import math
import unittest
import tempfile
import snp2counts

import CGAT.GTF as GTF
import CGAT.Genomics as Genomics
import CGAT.IndexedFasta as IndexedFasta


class getCDSPositionTestPos(unittest.TestCase):

    def setUp(self):

        self.mExons = []

        self.mSplitCodonsNext = {}
        self.mSplitCodonsPrev = {}

        self.mSpliceSize = 4
        self.mExonSize = 100
        self.mIntronSize = 900
        self.strand = "+"
        self.mNExons = 9
        self.mOffset = 1000
        length = 0
        self.frame = 0
        self.mIncrement = self.mIntronSize + self.mExonSize

        seq = list("123" * int((self.mNExons * self.mExonSize) / 3))

        exon_id = 0

        start = self.mOffset
        for x in range(self.mNExons):

            e = GTF.Entry()
            e.contig, e.strand, e.gene_id, e.transcript_id = "chr1", "+", "gene1", "trans1"
            e.start, e.end = start, start + self.mExonSize
            e.frame = (3 - (length % 3)) % 3
            length += e.end - e.start
            self.mExons.append(e)
            if e.frame != 0:
                for y in range(0, e.frame):
                    self.mSplitCodonsPrev[start + y] = start - self.mIntronSize
                for y in range(0, 3 - e.frame):
                    self.mSplitCodonsNext[
                        start - self.mIntronSize - y - 1] = start

            exon_id += 1
            if exon_id < self.mNExons:
                p = exon_id * self.mExonSize + self.mIntronSize * (exon_id - 1)
                seq[p:p] = list("AG")
                seq[p:p] = list("T" * (self.mIntronSize - 4))
                seq[p:p] = list("GT")

            start += self.mIncrement
            # print str(e)
        # print self.mSplitCodonsNext
        # print self.mSplitCodonsPrev
        seq[0:0] = "C" * self.mOffset
        seq.append("G" * self.mOffset)
        tmpfile = tempfile.NamedTemporaryFile()
        tmpfile.close()

        seq = "".join(seq)
        self.mSequence = seq
        self.contigSize = len(seq)
        IndexedFasta.createDatabase(tmpfile.name, iter([("chr1", seq), ]))
        self.mFasta = IndexedFasta.IndexedFasta(tmpfile.name)

    def tearDown(self):
        os.unlink(self.mFasta.getDatabaseName())
        os.unlink(self.mFasta.getDatabaseName()[:-len(".fasta")] + ".idx")

    def toRange(self, x, y):
        '''convert snp to positive strand base.'''
        if self.strand == "+":
            return x, y
        else:
            return self.contigSize - y, self.contigSize - x

    def testCodingSNPs(self):

        length = 0
        framed_length = (3 - self.frame) % 3
        phase = (3 - self.frame) % 3

        if self.strand == "+":
            motif = "123"
        else:
            motif = "321"

        for x in range(self.mOffset, self.contigSize - self.mOffset, self.mIncrement):
            for y in range(0, self.mExonSize):
                base = x + y
                rangex, rangey = self.toRange(base, base + 1)

                result = snp2counts.getCDSPosition(self.mExons, rangex, rangey,
                                                   lcontig=self.contigSize,
                                                   fasta=self.mFasta)

                self.assertEqual(result.strand, self.strand)
                self.assertEqual(result.cds_start, length)
                self.assertEqual(result.cds_end, length + 1)
                self.assertEqual(result.cds_phase, phase)
                self.assertEqual(result.intron_start, None)
                self.assertEqual(result.intron_end, None)
                self.assertEqual(len(result.cds_seq), 3)
                # print x, y, base, str(result)
                if self.frame == 0:
                    self.assertEqual(result.cds_seq, motif)

                self.assertEqual(result.cds_seq_start, framed_length % 3)
                self.assertEqual(result.cds_seq_end, (framed_length % 3) + 1)

                self.assertEqual(result.nc_seq, None)
                self.assertEqual(result.nc_start, None)
                self.assertEqual(result.nc_end, None)

                if base in self.mSplitCodonsPrev:
                    self.assertEqual(
                        result.prev_exon_end, self.mSplitCodonsPrev[base])
                else:
                    self.assertEqual(result.prev_exon_end, None)
                if base in self.mSplitCodonsNext:
                    self.assertEqual(
                        result.next_exon_start, self.mSplitCodonsNext[base])
                else:
                    self.assertEqual(result.next_exon_start, None)

                length += 1
                framed_length += 1
                phase += 1
                if phase >= 3:
                    phase = 0

    def testIntronsSNPs(self):

        length = 0
        t = 0
        exon_id = 0
        for x in range(self.mOffset, self.contigSize - self.mIncrement - self.mOffset, self.mIncrement):

            # exons
            for y in range(0, self.mExonSize):
                base = x + y
                base_x, base_y = self.toRange(base, base + 1)
                result = snp2counts.getCDSPosition(
                    self.mExons, base_x, base_y, lcontig=self.contigSize, fasta=self.mFasta)
                self.assertEqual(result.strand, self.strand)
                self.assertEqual(result.cds_start, t)
                self.assertEqual(result.cds_end, t + 1)
                self.assertEqual(result.intron_start, None)
                self.assertEqual(result.intron_end, None)
                self.assertEqual(len(result.cds_seq) % 3, 0)
                self.assertEqual(result.nc_seq, None)
                self.assertEqual(result.nc_start, None)
                self.assertEqual(result.nc_end, None)
                self.assertEqual(result.exon_id, exon_id)
                self.assertEqual(result.intron_id, None)
                t += 1

            exon_id += 1

            # introns
            for y in range(self.mExonSize, self.mExonSize + self.mIntronSize):
                base = x + y
                base_x, base_y = self.toRange(base, base + 1)
                result = snp2counts.getCDSPosition(
                    self.mExons, base_x, base_y, lcontig=self.contigSize, fasta=self.mFasta)
                self.assertEqual(result.strand, self.strand)
                self.assertEqual(result.cds_start, None)
                self.assertEqual(result.cds_end, None)
                self.assertEqual(result.cds_phase, None)
                self.assertEqual(result.intron_start, x + self.mExonSize)
                self.assertEqual(
                    result.intron_end, x + self.mIntronSize + self.mExonSize)
                self.assertEqual(result.cds_seq, None)
                self.assertEqual(result.cds_seq_start, None)
                self.assertEqual(result.cds_seq_end, None)
                self.assertEqual(len(result.nc_seq), 1)
                self.assert_(result.nc_seq not in "abc")
                self.assertEqual(result.nc_start, base)
                self.assertEqual(result.nc_end, base + 1)
                self.assertEqual(result.exon_id, exon_id)
                self.assertEqual(result.intron_id, exon_id - 1)

    def testIndels(self):
        '''test with segments of size 5'''
        size = 5

        length = 0
        framed_length = (3 - self.frame) % 3
        phase = (3 - self.frame) % 3

        if self.strand == "+":
            motif = "123"
        else:
            motif = "321"

        for x in range(self.mOffset, self.contigSize - self.mIncrement - self.mOffset, self.mIncrement):

            for y in range(-2 * size, self.mExonSize + 2 * size):
                base = x + y
                if base < self.mOffset:
                    continue

                base_x, base_y = self.toRange(base, base + size)
                result = snp2counts.getCDSPosition(
                    self.mExons, base_x, base_y, lcontig=self.contigSize, fasta=self.mFasta)

                if -size < y < self.mExonSize:
                    # overlap with coding sequence
                    self.assertEqual(len(result.cds_seq) % 3, 0)
                    self.assertEqual(result.cds_start, length)
                    if y < 0:
                        self.assertEqual(result.cds_end, length + size + y)
                    else:
                        self.assertEqual(
                            result.cds_end, length + min(size, self.mExonSize - y))

                    self.assertEqual(result.cds_phase, phase)

                    self.assertEqual(result.strand, self.strand)
                    ncodons = int(
                        math.ceil((result.cds_phase + result.cds_end - result.cds_start) / 3.0))
                    if self.frame == 0:
                        self.assertEqual(result.cds_seq, motif * ncodons)
                    self.assertEqual(result.cds_seq_start, framed_length % 3)
                    self.assertEqual(
                        result.cds_seq_end, framed_length % 3 + min(size, size + y, self.mExonSize - y))
                    if result.nc_end != None:
                        self.assertEqual(
                            result.cds_end - result.cds_start + (result.nc_end - result.nc_start), size)
                        self.assertEqual(
                            len(result.nc_seq), (result.nc_end - result.nc_start))
                else:
                    self.assertEqual(result.cds_start, None)
                    self.assertEqual(result.cds_end, None)
                    self.assertEqual(result.cds_phase, None)

                if y > self.mExonSize - size:
                    self.assertEqual(result.intron_start, x + self.mExonSize)
                    self.assertEqual(
                        result.intron_end, x + self.mIntronSize + self.mExonSize)
                elif y < 0:
                    self.assertEqual(result.intron_start, x - self.mIntronSize)
                    self.assertEqual(result.intron_end, x)

                if 0 <= y < self.mExonSize:
                    length += 1
                    framed_length += 1
                    phase += 1
                    if phase >= 3:
                        phase = 0


class getCDSPositionTestNeg(getCDSPositionTestPos):

    def setUp(self):
        getCDSPositionTestPos.setUp(self)

        for x in self.mExons:
            x.start, x.end = self.contigSize - x.end, self.contigSize - x.start
            x.strand = "-"
            # frame remains

        self.mExons.reverse()
        self.strand = "-"


class getCDSPositionTestWithStartingFrame2(getCDSPositionTestPos):

    '''test with a transcript not starting at frame 0, but at frame 2.'''

    def setUp(self):

        getCDSPositionTestPos.setUp(self)

        self.mSplitCodonsNext = {}
        self.mSplitCodonsPrev = {}

        start = self.mOffset
        l = 1
        for exon_id, e in enumerate(self.mExons):
            e.frame = (3 - l % 3) % 3
            l += e.end - e.start
            if e.frame != 0:
                if exon_id > 0:
                    for y in range(0, e.frame):
                        self.mSplitCodonsPrev[
                            start + y] = start - self.mIntronSize
                if exon_id < self.mNExons - 1:
                    for y in range(0, 3 - e.frame):
                        self.mSplitCodonsNext[
                            start - self.mIntronSize - y - 1] = start
            start += self.mIncrement

        self.frame = self.mExons[0].frame
#         for e in self.mExons:
#             print str(e)
#         print self.mSplitCodonsPrev
#         print self.mSplitCodonsNext


class iterateOverFrames(unittest.TestCase):

    def setUp(self):
        self.seq = list("AAA" * 20)
        self.length = len(self.seq)
        self.ncodons = self.length / 3

    def merge(self, result):

        n = []
        last = result[0]
        for this in result[1:]:
            if last[0] == this[0]:
                last[-1] = this[-1]
            else:
                n.append(tuple(last))
                last = this
        n.append(tuple(last))
        return n

    def testDeletion(self):
        '''test single deletion.'''

        for l in range(1, 7):
            for x in range(0, len(self.seq)):
                s = list(self.seq)
                todelete = min(l, self.length - x)
                for y in range(x, x + todelete):
                    s[y] = ""
                ncodons = self.ncodons - todelete // 3
                i = list(snp2counts.iterateOverFrames(s))

                codon_start = (x // 3) * 3
                codon_end = min(self.length, x + l + (3 - (x + l) % 3) % 3)

                result = []
                if codon_start > 0:
                    result.append([True, 0, codon_start])

                if todelete % 3 == 0:
                    if x % 3 != 0:
                        result.append([False, codon_start, codon_end])
                        if codon_end < self.length:
                            result.append([True, codon_end, self.length])
                    else:
                        result.append([True, codon_start, self.length])
                else:
                    o = codon_start
                    if todelete > 3 and x % 3 == 0:
                        o = codon_start + (todelete // 3) * 3
                        result.append([True, codon_start, o])

                    result.append([False, o, codon_end])
                    result.append([False, codon_end, self.length])

                result = self.merge(result)

                self.assertEqual(i, result)

    def testInsertion(self):
        '''test single insertion.'''

        for l in range(1, 7):
            for x in range(len(self.seq)):
                s = list(self.seq)
                s[x] = "A" * l + s[x]
                i = list(snp2counts.iterateOverFrames(s))
                result = []

                codon_start = (x // 3) * 3

                if codon_start > 0:
                    result.append([True, 0, codon_start])

                if l % 3 == 0:
                    result.append([True, 0, self.length])
                else:
                    result.append([False, codon_start, self.length])

                result = self.merge(result)

                self.assertEqual(i, result)


class countEffectsOnTranscript(unittest.TestCase):

    '''test countEffectsOnTranscript'''

    def setUp(self):
        self.seq = list("AAA" * 20)
        self.length = len(self.seq)
        self.ncodons = self.length / 3

    def testEmpty(self):
        r = snp2counts.countEffectsOnTranscript(self.seq, self.seq)
        self.assertEqual(r.ninserted_bases, 0)
        self.assertEqual(r.ninserted_codons, 0)
        self.assertEqual(r.ndeleted_bases, 0)
        self.assertEqual(r.ndeleted_codons, 0)
        self.assertEqual(r.noframe_codons, 0)
        self.assertEqual(r.nwrong_frames, 0)
        self.assertEqual(r.ncorrected_frames, 0)
        self.assertEqual(r.first_stop, self.ncodons)
        self.assertEqual(r.nstop_codons, 0)
        self.assertEqual(r.nstops, 0)
        self.assertEqual(r.nunaffected_codons, self.ncodons)
        self.assertEqual(r.nsynonymous_codons, 0)
        self.assertEqual(r.nnonsynonymous_codons, 0)

    def testInsertion(self):
        '''test single insertion.'''

        for l in range(1, 7):
            for x in range(len(self.seq)):
                s = list(self.seq)
                s[x] = "A" * l + s[x]

                r = snp2counts.countEffectsOnTranscript(s, self.seq)
                # print s, str(r)
                self.assertEqual(r.ninserted_bases, l)
                if l % 3 == 0:
                    self.assertEqual(r.ninserted_codons, l // 3)
                    self.assertEqual(r.nwrong_frames, 0)

                self.assertEqual(r.ndeleted_bases, 0)
                self.assertEqual(r.ndeleted_codons, 0)
                unaffected = x // 3
                if l % 3 == 0:
                    self.assertEqual(r.noframe_codons, 0)
                    self.assertEqual(r.ncorrected_frames, 0)
                    self.assertEqual(r.nunaffected_codons, 20)
                    self.assertEqual(r.nsynonymous_codons, 0)
                else:
                    self.assertEqual(r.noframe_codons, self.ncodons - x / 3)
                    self.assertEqual(r.ncorrected_frames, 0)
                    self.assertEqual(r.nunaffected_codons, unaffected)

                self.assertEqual(r.first_stop, (self.length + l) // 3)
                self.assertEqual(r.nnonsynonymous_codons, 0)

    def testDeletion(self):
        '''test single deletion.'''

        for l in range(1, 7):
            for x in range(0, len(self.seq)):
                s = list(self.seq)
                todelete = min(l, self.length - x)
                for y in range(x, x + todelete):
                    s[y] = ""
                ncodons = self.ncodons - todelete // 3
                r = snp2counts.countEffectsOnTranscript(s, self.seq)
                # print s, str(r)

                self.assert_(r.ndeleted_codons + r.nunaffected_codons + r.nincomplete_codons +
                             r.nnonsynonymous_codons + r.nsynonymous_codons + r.nstop_codons <=
                             self.ncodons)

                self.assertEqual(r.ninserted_bases, 0)
                self.assertEqual(r.ninserted_codons, 0)
                self.assertEqual(r.ndeleted_bases, todelete)

                codon_start = (x // 3) * 3
                codon_end = x + l + (3 - (x + l) % 3) % 3
                affected_codons = (codon_end - codon_start) // 3
                deletion_codon_start = x + (3 - (x % 3)) % 3
                deletion_codon_end = min(self.length, ((x + l) // 3) * 3)

                # subtract fully deleted codons
                deleted_codons = max(
                    0, (deletion_codon_end - deletion_codon_start) // 3)
                self.assertEqual(r.ndeleted_codons, deleted_codons)

                inframe = x // 3

                # delete in-frame, multiple of 3
                if x % 3 == 0 and todelete % 3 == 0:
                    self.assertEqual(r.noframe_codons, 0)
                    self.assertEqual(r.nwrong_frames, 0)
                    self.assertEqual(r.ncorrected_frames, 0)
                    self.assertEqual(r.nsynonymous_codons, 0)

                # delete out-of-frame, multiple of 3
                elif x % 3 != 0 and todelete % 3 == 0:
                    self.assertEqual(
                        r.noframe_codons, affected_codons - deleted_codons)
                    self.assertEqual(r.nwrong_frames, 1)
                    self.assertEqual(r.ncorrected_frames, 1)

                # delete, but not multiple of 3
                else:
                    self.assertEqual(
                        r.noframe_codons, self.ncodons - inframe - deleted_codons)
                    self.assertEqual(r.nwrong_frames, 1)
                    self.assertEqual(r.ncorrected_frames, 0)
#                    self.assertEqual( r.nsynonymous_codons,
# self.ncodons - r.nincomplete_codons - inframe - deleted_codons)

                self.assertEqual(r.first_stop, (self.length - todelete) // 3)
                # self.assertEqual( r.nunaffected_codons, self.ncodons - (int(math.ceil( (x + todelete) / 3.0)) - x // 3) )
                self.assertEqual(r.nnonsynonymous_codons, 0)

    def testFrameCorrectionAdacent(self):
        '''test frame correction within a codon for
        two adjacent bases.

        Strictly speaking this should not happen as these
        would be called as substitutions.
        '''
        return

        for l in range(1, 7):
            for x in range(len(self.seq) - l):
                s = list(self.seq)
                todelete = l
                toinsert = l
                for y in range(x, x + todelete):
                    s[y] = ""
                s[x + todelete] = "A" * toinsert + s[x + todelete]
                ncodons = self.ncodons
                # print l,x,todelete, toinsert
                # print s
                r = snp2counts.countEffectsOnTranscript(s, self.seq)
                # print str(r)

                self.assert_(r.ndeleted_codons + r.nunaffected_codons + r.nincomplete_codons +
                             r.nnonsynonymous_codons + r.nsynonymous_codons + r.nstop_codons <=
                             self.ncodons)

                self.assertEqual(r.ninserted_codons, 0)

                if (x + todelete) % 3 != 0:
                    self.assertEqual(r.ninserted_bases, 0)
                    self.assertEqual(r.ndeleted_bases, 0)
                    self.assertEqual(r.noframe_codons, 0)
                else:
                    self.assertEqual(r.ninserted_bases, toinsert)
                    self.assertEqual(r.ndeleted_bases, todelete)
                    self.assertEqual(r.noframe_codons, 2)

                if x % 3 == 0 and todelete % 3 == 0:
                    self.assertEqual(r.ndeleted_codons, todelete / 3)
                else:
                    self.assertEqual(r.ndeleted_codons, 0)

                self.assert_(r.noframe_codons <= self.ncodons)

                self.assertEqual(r.nwrong_frames, 1)
                self.assertEqual(r.ncorrected_frames, 1)
                self.assertEqual(r.ntruncated_codons_stop, 0)
                # self.assertEqual( r.nunaffected_codons, self.ncodons )
                self.assertEqual(r.nsynonymous_codons, 0)
                self.assertEqual(r.nnonsynonymous_codons, 0)

    def testFrameCorrection(self):
        '''test frame correction within a codon for
        two adjacent bases.

        Strictly speaking this should not happen as these
        would be called as substitutions.
        '''
        return
        for l in range(1, 7):
            for offset in range(1, 5):
                for x in range(len(self.seq) - (l + offset)):
                    s = list(self.seq)
                    todelete = l
                    toinsert = l
                    for y in range(x, x + todelete):
                        s[y] = ""

                    codon_start = (x // 3) * 3
                    codon_end = ((x + offset + todelete) // 3 + 1) * 3

                    s[x + todelete + offset] = "A" * toinsert + s[x + todelete]
                    ncodons = self.ncodons
                    # print "l=",l,"x=",x,"offest=",offset,"del=",todelete,
                    # "ins=",toinsert, "start=",codon_start, "end=", codon_end

                    naffected_codons = (codon_end - codon_start) // 3
                    if todelete % 3 == 0 and (x + todelete) // 3 != (x + todelete + offset) // 3:
                        # affected codons reduced, if offset includes full
                        # codons
                        naffected_codons -= (x + todelete +
                                             offset) // 3 - (x + todelete) // 3
                        # if offset > 3:
                        #    naffected_codons -= 1

                    # print s
                    r = snp2counts.countEffectsOnTranscript(s, self.seq)
                    # print str(r)

                    self.assertEqual(r.ninserted_codons, l // 3)

                    self.assertEqual(r.ninserted_bases, toinsert)
                    self.assertEqual(r.ndeleted_bases, todelete)

                    if l + offset <= 2 and x % 3 == 0 or (x % 3 == 0 and l % 3 == 0):
                        # within codon correction
                        self.assertEqual(r.noframe_codons, 0)
                        self.assertEqual(r.nwrong_frames, 0)
                        self.assertEqual(r.ncorrected_frames, 0)

                    else:
                        # between codon correction
                        self.assertEqual(r.ninserted_bases, toinsert)
                        self.assertEqual(r.ndeleted_bases, todelete)
                        self.assertEqual(r.noframe_codons, naffected_codons)
                        self.assertEqual(r.nwrong_frames, 1)
                        self.assertEqual(r.ncorrected_frames, 1)

                    if x % 3 == 0 and todelete % 3 == 0:
                        self.assertEqual(r.ndeleted_codons, todelete / 3)
                    else:
                        self.assertEqual(r.ndeleted_codons, 0)

                    self.assert_(r.noframe_codons <= self.ncodons)

                    self.assertEqual(r.first_stop, 0)

                    # self.assertEqual( r.nunaffected_codons, self.ncodons )
                    self.assertEqual(r.nsynonymous_codons, 0)
                    self.assertEqual(r.nnonsynonymous_codons, 0)

    def testStop(self):
        '''test one stop codon.'''
        for x in range(len(self.seq)):
            s = list(self.seq)
            s[x] = "T"
            r = snp2counts.countEffectsOnTranscript(s, self.seq)
            # print s, str(r)
            self.assertEqual(r.ninserted_bases, 0)
            self.assertEqual(r.ninserted_codons, 0)
            self.assertEqual(r.ndeleted_bases, 0)
            self.assertEqual(r.ndeleted_codons, 0)
            self.assertEqual(r.noframe_codons, 0)
            self.assertEqual(r.nwrong_frames, 0)
            self.assertEqual(r.ncorrected_frames, 0)
            self.assertEqual(r.nsynonymous_codons, 0)

            if x % 3 == 0:
                self.assertEqual(r.nstops, 1)
                self.assertEqual(r.first_stop, x // 3)
                # ignore last incomplete codon
                if x < self.length - 3:
                    self.assertEqual(r.nstop_codons, 1)
                self.assertEqual(r.nnonsynonymous_codons, 0)
            else:
                self.assertEqual(r.nstops, 0)
                self.assertEqual(r.first_stop, self.ncodons)
                self.assertEqual(r.nstop_codons, 0)
                self.assertEqual(r.nnonsynonymous_codons, 1)

    def testMutation(self):
        '''test synonymous/nonsynonymous mutation.'''

        for x in range(len(self.seq)):
            s = list(self.seq)
            # aaa = K, aag = N
            s[x] = "G"
            r = snp2counts.countEffectsOnTranscript(s, self.seq)
            self.assertEqual(r.ninserted_bases, 0)
            self.assertEqual(r.ninserted_codons, 0)
            self.assertEqual(r.ndeleted_bases, 0)
            self.assertEqual(r.ndeleted_codons, 0)
            self.assertEqual(r.noframe_codons, 0)
            self.assertEqual(r.nwrong_frames, 0)
            self.assertEqual(r.ncorrected_frames, 0)
            self.assertEqual(r.nstops, 0)
            self.assertEqual(r.first_stop, self.ncodons)
            self.assertEqual(r.nstop_codons, 0)

            if x % 3 == 2:
                self.assertEqual(r.nsynonymous_codons, 1)
                self.assertEqual(r.nnonsynonymous_codons, 0)
            else:
                self.assertEqual(r.nsynonymous_codons, 0)
                self.assertEqual(r.nnonsynonymous_codons, 1)

    def testStopDouble(self):
        '''test two stop codons.'''
        ref = list(self.seq)
        ref[-3] = "T"
        for x in range(len(ref) - 3):
            s = list(ref)
            s[x] = "T"
            r = snp2counts.countEffectsOnTranscript(s, ref)
            self.assertEqual(r.ninserted_bases, 0)
            self.assertEqual(r.ninserted_codons, 0)
            self.assertEqual(r.ndeleted_bases, 0)
            self.assertEqual(r.ndeleted_codons, 0)
            self.assertEqual(r.noframe_codons, 0)
            self.assertEqual(r.nwrong_frames, 0)
            self.assertEqual(r.ncorrected_frames, 0)

            if x % 3 == 0:
                self.assertEqual(r.nstops, 2)
                self.assertEqual(r.first_stop, x // 3)
                # ignore last incomplete codon
                self.assertEqual(r.nstop_codons, 1)
            else:
                self.assertEqual(r.nstops, 1)
                self.assertEqual(r.first_stop, self.ncodons - 1)
                self.assertEqual(r.nstop_codons, 0)


if __name__ == "__main__":
    unittest.main()

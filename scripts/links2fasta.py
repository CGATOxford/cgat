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
'''
links2fasta.py - convert links into alignments
==============================================

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

   python links2fasta.py --help

Type::

   python links2fasta.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import getopt
import tempfile
import time
import optparse

import CGAT.Experiment as E
import CGAT.BlastAlignments as BlastAlignments
import alignlib
import CGAT.Genomics as Genomics
import CGAT.Exons as Exons
import CGAT.FastaIterator as FastaIterator

class Map:
    def __init__(self):
        pass

    def read( self, line ):
        try:
            ( self.mToken,
              self.mOldFrom, self.mOldTo, self.mOldAli,
              self.mNewFrom, self.mNewTo, self.mNewAli,
              self.mOldLength, self.mNewLength) = line[:-1].split("\t")
        except ValueError:
            raise ValueError("parsing error in line\n%s" % line)

        (self.mOldFrom, self.mOldTo, self.mNewFrom, self.mNewTo) = \
                        map(int, (self.mOldFrom, self.mOldTo, self.mNewFrom, self.mNewTo))
        self.mMapOld2New = None
        
    def expand( self ):
        if not self.mMapOld2New:
            self.mMapOld2New = alignlib.makeAlignmentVector()
        
            alignlib.AlignmentFormatEmissions( 
                self.mOldFrom, self.mOldAli,
                self.mNewFrom, self.mNewAli).copy( self.mMapOld2New )
            
    def clear( self ):
        if self.mMapOld2New:
            self.mMapOld2New.clear()
        self.mMapOld2New = None

    def __str__(self):
        return string.join(map(str, (self.mToken,
                                     self.mOldFrom, self.mOldTo, self.mOldAli,
                                     self.mNewFrom, self.mNewTo, self.mNewAli,
                                     self.mOldLength, self.mNewLength)), "\t")
            
def ScaleAlignment( alignment, factor ):
    """scale alignment string."""

    data = re.split("[+-]", alignment[1:])
    
    data = map( lambda x: int(x) * factor, data )
    signs = [ "+", "-" ] * (1 + len(data) / 2)
    
    if alignment[0] == "+":
        del signs[-1]
    else:
        del signs[0]

    s = map( lambda x,y: "%s%i" % (x,y), signs, data)
    return string.join(s, "")

##-------------------------------------------------------------------------
def Write( map_row2col, row_seq, col_seq, link,
           no_gaps = False, no_identical = False,
           min_length = 0,
           suffix1="", suffix2="",
           outfile = None,
           pair_filter = None,
           format = "fasta" ):
    """write alignment based on map_row2col."""
        
    status = None

    filter_status = "new"
    
    if map_row2col.getLength() == 0:
        status = "empty"

    if not status:

        f = alignlib.AlignmentFormatExplicit( map_row2col, row_seq, col_seq )
        
        row_from = map_row2col.getRowFrom()
        row_to = map_row2col.getRowTo()
        col_from = map_row2col.getColFrom()
        col_to = map_row2col.getColTo()
        row_ali, col_ali = f.mRowAlignment, f.mColAlignment

    if not status:
        if no_gaps:
            # remove gaps from fasta
            r = []
            c = []
            for x in range(len(row_ali)):
                if row_ali[x] != "-" and col_ali[x] != "-":
                    r.append( row_ali[x] )
                    c.append( col_ali[x] )
            row_ali = string.join(r, "")
            col_ali = string.join(c, "")            

    if not status and len(row_ali) < min_length:
        status = "length"

    if not status and no_identical:
            if row_ali == col_ali:
                status = "identical"

    if not status:

        if pair_filter:
            id = "%s-%s" % (link.mQueryToken, link.mSbjctToken)
            if id in pair_filter:
                h = Genomics.GetHID( row_ali + ";" + col_ali ) 
                if h in pair_filter[id]:
                    filter_status = "old"
        
        translation1 = Genomics.TranslateDNA2Protein( row_ali )
        translation2 = Genomics.TranslateDNA2Protein( col_ali )        

        if "X" in translation1 or "x" in translation2:
            status = "stops"
        else:
            status = "success"

        if filter_status == "new":
            if format == "fasta":
                print ">%s%s %s %s\n%s\n>%s%s %s %s\n%s" % (link.mQueryToken, suffix1, row_from, row_to, row_ali, 
                                                            link.mSbjctToken, suffix2, col_from, col_to, col_ali )
            elif format == "dummy":
                pass
            else:
                raise ValueError("unknown format")

    if outfile:
        outfile.write( "%s%s\t%s%s\t%s\t%i\t%s\n" % (link.mQueryToken, suffix1, link.mSbjctToken, suffix2,
                                                     status, map_row2col.getLength(), filter_status ) )

    return status

def GetAdjustedBoundaries( id, exons ):
    """return codon adjusted exon boundaries."""

    f, t = exons[id].mPeptideFrom, exons[id].mPeptideTo

    f += exons[id].frame
    
    if id < len(exons) -1:
        next_frame = exons[id+1].frame
    else:
        next_frame = 0
        
    if next_frame:
        t -= 3 - next_frame

    return f, t

if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: links2fasta.py 2446 2009-01-27 16:32:35Z andreas $", usage = globals()["__doc__"] )

    parser.add_option( "-s", "--sequences", dest="filename_sequences", type="string",
                       help="peptide sequence [Default=%default]" )

    parser.add_option( "-f", "--format", dest="format", type="string",
                       help="output format [Default=%default]" )

    parser.add_option( "-e", "--expand",  dest="expand", action="store_true",
                       help="expand positions from peptide to nucleotide alignment [Default=%default]")

    parser.add_option( "-m", "--map",  dest="filename_map", type="string",
                       help="map alignments [Default=%default]")
    
    parser.add_option( "-c", "--codons",  dest="require_codons", action="store_true",
                       help="require codons [Default=%default]")

    parser.add_option( "--one-based-coordinates",  dest="one_based_coordinates", action="store_true",
                       help="expect one-based coordinates. The default are zero based coordinates [Default=%default].")

    parser.add_option( "--no-identical",  dest="no_identical", action="store_true",
                       help="do not output identical pairs [Default=%default]" )

    parser.add_option( "-g", "--no-gaps",  dest="no_gaps", action="store_true",
                       help="remove all gaps from aligned sequences [Default=%default]")

    parser.add_option( "-x", "--exons",  dest="filename_exons", type="string",
                       help="filename with exon boundaries [Default=%default]")
    
    parser.add_option( "-o", "--outfile",  dest="filename_outfile", type="string",
                       help="filename to save links [Default=%default]")

    parser.add_option( "--min-length",  dest="min_length", type="int",
                       help="minimum length of alignment [Default=%default]")

    parser.add_option( "--filter",  dest="filename_filter", type="string",
                       help="given a set of previous alignments, only write new pairs [Default=%default].")

    parser.set_defaults(
        filename_sequences = None,
        filename_exons = None,
        filename_map = None,
        filename_outfile = None,
        no_gaps = False,
        format = "fasta",
        expand = False,
        require_codons = False,
        no_identical = False,
        min_length = 0,
        report_step = 100,
        one_based_coordinates = False,
        filename_filter = None)

    (options, args) = E.Start( parser, add_mysql_options = True )

    t0 = time.time()
    if options.filename_sequences:
        sequences = Genomics.ReadPeptideSequences( open(options.filename_sequences, "r") )
    else:
        sequences = {}

    if options.loglevel >= 1:
        options.stdlog.write( "# read %i sequences\n" % len(sequences) )
        sys.stdout.flush()

    if options.filename_exons:
        exons = Exons.ReadExonBoundaries( open(options.filename_exons, "r") )
    else:
        exons = {}

    if options.loglevel >= 1:
        options.stdlog.write( "# read %i exons\n" % len(exons) )
        sys.stdout.flush()

    if options.filename_map:
        map_old2new = {}
        for line in open(options.filename_map, "r"):
            if line[0] == "#": continue
            m = Map()
            m.read( line )
            map_old2new[m.mToken] = m
    else:
        map_old2new = {}

    if options.loglevel >= 1:
        options.stdlog.write( "# read %i maps\n" % len(map_old2new) )
        sys.stdout.flush()

    if options.filename_filter:
        if options.loglevel >= 1:        
            options.stdlog.write( "# reading filtering information.\n" )
            sys.stdout.flush()
            
        map_pair2hids = {}

        if os.path.exists( options.filename_filter ):
            
            infile = open(options.filename_filter, "r")

            iterator = FastaIterator.FastaIterator( infile )

            while 1:
                cur_record = iterator.next()
                if cur_record is None: break

                record1 = cur_record

                cur_record = iterator.next()
                if cur_record is None: break

                record2 = cur_record

                identifier1 = re.match("(\S+)", record1.title).groups()[0]
                identifier2 = re.match("(\S+)", record2.title).groups()[0]

                id = "%s-%s" % (identifier1, identifier2)
                s = Genomics.GetHID(record1.sequence + ";" + record2.sequence)

                if id not in map_pair2hids: map_pair2hids[id] = []

                map_pair2hids[id].append( s )

            infile.close()
            
        if options.loglevel >= 1:        
            options.stdlog.write( "# read filtering information for %i pairs.\n" % len(map_pair2hids) )
            sys.stdout.flush()
    else:
        map_pair2hids = None
        
    if options.loglevel >= 1:
        options.stdlog.write( "# finished input in %i seconds.\n" % (time.time() - t0))

    if options.filename_outfile:
        outfile = open(options.filename_outfile, "w")
    else:
        outfile = None
        
    map_row2col = alignlib.makeAlignmentVector()
    tmp1_map_row2col = alignlib.makeAlignmentVector()
    counts = {}

    iterations = 0

    t1 = time.time()
    ninput, nskipped, noutput = 0, 0, 0

    for link in BlastAlignments.iterator_links( sys.stdin ):

        iterations += 1
        ninput += 1

        if options.loglevel >= 1:
            if (iterations % options.report_step == 0):
                options.stdlog.write( "# iterations: %i in %i seconds.\n" % (iterations, time.time() - t1) )
                sys.stdout.flush()
                
        if link.mQueryToken not in sequences or \
           link.mSbjctToken not in sequences:
            nskipped += 1
            continue

        if options.loglevel >= 3:
            options.stdlog.write( "# read link %s\n" %  str(link) )
            
        row_seq = alignlib.makeSequence( sequences[link.mQueryToken] )
        col_seq = alignlib.makeSequence( sequences[link.mSbjctToken] )

        if options.one_based_coordinates:
            link.mQueryFrom -= 1
            link.mSbjctFrom -= 1

        if options.expand:
            link.mQueryFrom = link.mQueryFrom * 3 
            link.mSbjctFrom = link.mSbjctFrom * 3
            link.mQueryAli = ScaleAlignment( link.mQueryAli, 3 )
            link.mSbjctAli = ScaleAlignment( link.mSbjctAli, 3 )            
            
        map_row2col.clear()

        alignlib.AlignmentFormatEmissions(
            link.mQueryFrom, link.mQueryAli,
            link.mSbjctFrom, link.mSbjctAli ).copy(  map_row2col )
        
        if link.mQueryToken in map_old2new:
            tmp1_map_row2col.clear()
            map_old2new[link.mQueryToken].expand()
            if options.loglevel >= 3:
                options.stdlog.write( "# combining in row with %s\n" %\
                                      str(alignlib.AlignmentFormatEmissions(map_old2new[link.mQueryToken].mMapOld2New ) ))

            alignlib.combineAlignment( tmp1_map_row2col,
                                      map_old2new[link.mQueryToken].mMapOld2New,
                                      map_row2col,
                                      alignlib.RR )
            map_old2new[link.mQueryToken].clear()
            alignlib.copyAlignment( map_row2col, tmp1_map_row2col )

        if link.mSbjctToken in map_old2new:
            tmp1_map_row2col.clear()
            map_old2new[link.mSbjctToken].expand()            
            if options.loglevel >= 3:
                options.stdlog.write( "# combining in col with %s\n" %\
                                      str(alignlib.AlignmentFormatEmissions(map_old2new[link.mSbjctToken].mMapOld2New ) ))

            alignlib.combineAlignment( tmp1_map_row2col,
                                       map_row2col,
                                       map_old2new[link.mSbjctToken].mMapOld2New,
                                       alignlib.CR )
            map_old2new[link.mSbjctToken].clear()
            alignlib.copyAlignment( map_row2col, tmp1_map_row2col )

        dr = row_seq.getLength() - map_row2col.getRowTo() 
        dc = col_seq.getLength() - map_row2col.getColTo() 
        if dr < 0 or dc < 0:
            raise ValueError("out of bounds alignment: %s-%s: alignment out of bounds. row=%i col=%i ali=%s" %\
                                          (link.mQueryToken,
                                           link.mSbjctToken,
                                           row_seq.getLength(),
                                           col_seq.getLength(),
                                           str(alignlib.AlignmentFormatEmissions(map_row2col))))
            

        if options.loglevel >= 2:
            options.stdlog.write( str( alignlib.AlignmentFormatExplicit( map_row2col, 
                                                                         row_seq, 
                                                                         col_seq )) + "\n" )
        ## check for incomplete codons
        if options.require_codons:

            naligned = map_row2col.getNumAligned()
            
            # turned off, while fixing alignlib
            if naligned % 3 != 0:
                options.stdlog.write( "# %s\n" % str(map_row2col) )
                options.stdlog.write( "# %s\n" % str(link) )
                options.stdlog.write( "# %s\n" % str(map_old2new[link.mQueryToken]) )
                options.stdlog.write( "# %s\n" % str(map_old2new[link.mSbjctToken]) )
                options.stdlog.write( "#\n%s\n" % alignlib.AlignmentFormatExplicit( map_row2col, 
                                                                                    row_seq,
                                                                                    col_seq ) )

                raise ValueError("incomplete codons %i in pair %s - %s" % (naligned, link.mQueryToken, link.mSbjctToken))

        ## if so desired, write on a per exon level:
        if exons:
            if link.mQueryToken not in exons:
                raise IndexError("%s not found in exons" % (link.mQueryToken))
            if link.mSbjctToken not in exons:
                raise IndexError("%s not found in exons" % (link.mSbjctToken))
            exons1 = exons[link.mQueryToken]
            exons2 = exons[link.mSbjctToken]

            ## Get overlapping segments
            segments = Exons.MatchExons( map_row2col, exons1, exons2 )
            
            for a,b in segments:
                tmp1_map_row2col.clear()

                # make sure you got codon boundaries. Note that frameshifts
                # in previous exons will cause the codons to start at positions
                # different from mod 3. The problem is that I don't know where
                # the frameshifts occur exactly. The exon boundaries are given
                # with respect to the cds, which include the frame shifts.
                # Unfortunately, phase information seems to be incomplete in the input files.

                from1, to1 = GetAdjustedBoundaries( a, exons1 )
                from2, to2 = GetAdjustedBoundaries( b, exons2 )

                alignlib.copyAlignment( tmp1_map_row2col, map_row2col,
                                       from1+1, to1, from2+1, to2 )
                
                mode = Write( tmp1_map_row2col, row_seq, col_seq, link,
                              no_gaps = options.no_gaps,
                              no_identical = options.no_identical,
                              min_length = options.min_length,
                              suffix1="_%s" % str(a),
                              suffix2="_%s" % str(b),
                              outfile = outfile,
                              pair_filter = map_pair2hid,
                              format = options.format )

                if mode not in counts: counts[mode] = 0
                counts[mode] += 1

        else:
            mode = Write( map_row2col, row_seq, col_seq, link,
                          min_length = options.min_length,                          
                          no_gaps = options.no_gaps,
                          no_identical = options.no_identical,
                          outfile = outfile,
                          pair_filter = map_pair2hids,
                          format = options.format )
            
            if mode not in counts: counts[mode] = 0
            counts[mode] += 1

        noutput += 1
        
    if outfile: outfile.close()
    
    if options.loglevel >= 1:
        options.stdlog.write("# %s\n" % ", ".join( map( lambda x,y: "%s=%i" % (x,y), counts.keys(), counts.values() ) ))
        options.stdlog.write("# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped) )

    E.Stop()

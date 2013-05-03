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
profile_vs_profile.py - compare two profile libraries
=====================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

compare two profile libraries.

Usage
-----

Example::

   python profile_vs_profile.py --help

Type::

   python profile_vs_profile.py --help

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
import tempfile
import subprocess
import optparse
import time
import math
import shutil

#--------------------------------------------------------
#--------------------------------------------------------
#--------------------------------------------------------
# import of user libraries
#--------------------------------------------------------
import CGAT.Experiment as Experiment
import alignlib
from ProfileLibrary import ProfileLibrary
from ProfileLibraryCompass import ProfileLibraryCompass

def getKeys( plib, start = None, end = None ):
    """get keys of profiles to compare."""
    k = plib.keys()
    k.sort()
    if not start: start = 0
    if not end: end = len(k)
    return k[max(0,start):min(end,len(k))], start, end

class CompassResult:
    def __init__(self):
        pass

class AlignatorCompass:

    mAligner = "compass_db1Xdb2"
    mReferenceLength = 1000000

    def __init__(self):
        self.mTempdir = tempfile.mkdtemp()
        self.mFilenameQuery = self.mTempdir + "/query"
        self.mFilenameSbjct = self.mTempdir + "/sbjct"
        self.mFilenameQueryLength = self.mFilenameQuery + ".len"
        self.mFilenameSbjctLength = self.mFilenameSbjct + ".len"
        
        outfile = open( self.mFilenameQueryLength, "w" )
        outfile.write( "%i\n" % self.mReferenceLength )
        outfile.close()

        outfile = open( self.mFilenameSbjctLength, "w" )
        outfile.write( "%i\n" % self.mReferenceLength )
        outfile.close()
        
    def __del__( self ):
        
        # shutil.rmtree( self.mTempdir ) did not work
        for f in (self.mFilenameQuery, self.mFilenameSbjct, 
                  self.mFilenameQueryLength, self.mFilenameSbjctLength ):
            if os.path.exists(f):
                os.remove(f)
        os.rmdir( self.mTempdir )

    def writeProfile( self, filename, profile, name = None ):
        
        if name:
            old_name = profile.getName()
            profile.setName( name )
        outfile = open( filename, "w" )
        profile.save( outfile )
        outfile.close()

        if name: profile.setName( old_name )

    def align( self, query, sbjct, map_querysbjct ):
        """align query and sbjct profile. 

        Result is stored in map_query2sbjct. In addition,
        a method specific result object is returned.
        """
        
        self.writeProfile( self.mFilenameQuery, query, "query" )
        self.writeProfile( self.mFilenameSbjct, sbjct, "sbjct" )

        statement = "%s -i %s -j %s" % (self.mAligner, self.mFilenameQuery, self.mFilenameSbjct )
        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = self.mTempdir,
                              close_fds = True)                              

        (out, err) = s.communicate()
        
        if s.returncode != 0:
            raise "Error in running %s \n%s\n%s\nTemporary directory in %s" % (statement, err, out, self.mTempdir)

        return self.parseResult( out, err, map_query2sbjct )

    def addBlocks( self, 
                   query_index, query_ali, 
                   sbjct_index, sbjct_ali,
                   map_query2sbjct ):
        """ parse alignment. From the COMPASS website:

CAPITAL letters: residues at positions aligned by COMPASS, i.e. at input alignment positions 
                 with gap content < threshold of gap fraction (see above); 
lower-case letters: residues at positions not used by COMPASS, i.e. at input alignment positions 
                    with gap content >= threshold of gap fraction (see above);
'-' : gaps retained from original alignments at positions aligned by COMPASS, i.e. at positions 
      with gap content < threshold;
'.' : gaps retained from original alignments at positions not used by COMPASS, i.e. at positions 
      with gap content >= threshold;
'=' : gaps introduced by COMPASS in profile-profile alignment;
'~' : gaps introduced by COMPASS against positions that are not used in the construction of 
      profile-profile alignment (positions with gap content >= threshold);
      """
        
        gap_chars = "=~"

        for x in range( 0, len(query_ali) ):

            # skip over gaps introduced by compass
            if query_ali[x] in gap_chars:
                sbjct_index += 1
                continue
            elif sbjct_ali[x] in gap_chars:
                query_index += 1
                continue

            is_unaligned = False
            # deal with unaligned positions - these can be matched up
            if query_ali[x] in string.lowercase:
                query_index += 1
                is_unaligned = True

            if sbjct_ali[x] in string.lowercase:
                sbjct_index += 1
                is_unaligned = True

            if is_unaligned: continue

            map_query2sbjct.addPair( query_index, sbjct_index )
            query_index += 1
            sbjct_index += 1
            
        return query_index, sbjct_index

    def parseResult( self, out, err, map_query2sbjct ):
        """parse result from compass."""
        
        result = CompassResult()
        map_query2sbjct.clear()

        lines = out.split("\n")

        result.mQuery, result.mSbjct = re.match( "Ali1:\s+(\S+)\s+Ali2:\s+(\S+)", lines[0]).groups()
        result.mQueryLength, result.mQueryLengthFiltered, result.mSbjctLength, result.mSbjctLengthFiltered = \
            map( int, re.match("length1=(\d+)\s+filtered_length1=(\d+)\s+length2=(\d+)\s+filtered_length2=(\d+)", lines[2] ).groups() )
        result.mQueryNSeqs, result.mQueryNEffective, result.mSbjctNSeqs, result.mSbjctNEffective = \
            map( float, re.match("Nseqs1=(\S+)\s+Neff1=(\S+)\s+Nseqs2=(\d+)\s+Neff2=(\S+)", lines[3] ).groups() )
        result.score, result.mEvalue = \
            map( float, re.match("Smith-Waterman score = (\S+)\s+Evalue = (\S+)", lines[4]).groups() )
        
        x = 6

        d, query_index, query_ali = re.split("\s+", lines[x] )
        d, sbjct_index, sbjct_ali = re.split("\s+", lines[x+2] )
        query_index, sbjct_index = self.addBlocks( int(query_index) - 1, query_ali, 
                                                   int(sbjct_index) - 1, sbjct_ali,
                                                   map_query2sbjct )

        for x in range( 11, len(lines), 5):
            d, query_ali = re.split("\s+", lines[x] )
            d, sbjct_ali = re.split("\s+", lines[x+2] )
            
            query_index, sbjct_index = self.addBlocks( query_index, query_ali, 
                                                       sbjct_index, sbjct_ali,
                                                       map_query2sbjct )

        map_query2sbjct.setScore( result.score )

#--------------------------------------------------------
#--------------------------------------------------------
#--------------------------------------------------------
# main part of script
#--------------------------------------------------------
if __name__ == "__main__":

    #--------------------------------------------------------
    # command line parsing options
    parser = optparse.OptionParser( version = "%prog version: $Id: profile_vs_profile.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-q", "--query", dest="query", type="string",
                      help="query profile library." )

    parser.add_option("-s", "--sbjct", dest="sbjct", type="string",
                      help="sbjct profile library." )

    parser.add_option("-e", "--self-compare", dest="self_compare", action="store_true",
                      help="self-comparison. Only compare one direction." )

    parser.add_option( "--query-start", dest="query_start", type="int",
                      help="start at xth entry of query." )

    parser.add_option( "--query-end", dest="query_end", type="int",
                      help="stop at xth entry of query." )

    parser.add_option( "--sbjct-start", dest="sbjct_start", type="int",
                      help="start at xth entry of sbjct." )

    parser.add_option( "--sbjct-end", dest="sbjct_end", type="int",
                      help="stop at xth entry of sbjct." )

    parser.add_option( "--filename-pairs", dest="filename_pairs", type="string",
                       help="align a list of pairs." )

    parser.add_option( "--iterative-min-score", dest="iterative_min_score", type="float",
                       help="score threshold for iterative alignment." )

    parser.add_option( "--alignment-mode", dest="alignment_mode", type="choice",
                       choices=("iterative-profile", "iterative-sequence", "compass"),
                       help="alignment mode." )

    parser.set_defaults( query = None,
                         sbjct = None,
                         query_start = None,
                         query_end = None,
                         sbjct_start = None,
                         sbjct_end = None,
                         report_step = 100,
                         filename_pairs= None,
                         iterative_min_score = 40.0,
                         alignment_mode = "iterative-profile",
                         )

    (options, args) = Experiment.Start( parser )

    #--------------------------------------------------------
    # main part of script
    
    if not options.query:
        print USAGE
        raise "please supply a query."
    
    if options.self_compare:
        options.sbjct = options.query
        if options.sbjct_end and options.query_start and \
                options.sbjct_end < options.query_start:
            if options.loglevel >= 1:
                options.stdlog.write( "# subsections to compare are out of range for self comparison." )
            Experiment.Stop()
            sys.exit(0)

        ## adjust sbjct start to upper diagonal
        if options.query_start and options.sbjct_start:
            options.sbjct_start = max( options.query_start, options.sbjct_start )
    else:
        if not options.sbjct:
            print USAGE
            raise "please supply both a query and a sbjct."
    
    if options.alignment_mode == "compass":
        plib_query = ProfileLibraryCompass( options.query, "r" )
        plib_sbjct = ProfileLibraryCompass( options.sbjct, "r" )
    else:
        plib_query = ProfileLibrary( options.query, "r" )
        plib_sbjct = ProfileLibrary( options.sbjct, "r" )

    if options.alignment_mode == "iterative-profile":
        alignator1 = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_LOCAL, -10.0, -2.0 )
        alignator = alignlib.makeAlignatorIterative( alignator1, options.iterative_min_score )

    elif options.alignment_mode == "iterative-sequence":
        class AlignatorSequence:
            def __init__(self):
                self.mAlignator1 = alignlib.makeAlignatorDPFull( alignlib.ALIGNMENT_LOCAL, -10.0, -2.0 )
                self.mAlignator = alignlib.makeAlignatorIterative( self.mAlignator1, options.iterative_min_score )

            def align(self, query, sbjct, map_query2sbjct):
                xrow = alignlib.makeSequence(query.asString())
                xcol = alignlib.makeSequence(sbjct.asString())
                self.mAlignator.align( xrow, xcol, map_query2sbjct)
                
        alignator = AlignatorSequence()
    elif options.alignment_mode == "compass":
        alignator = AlignatorCompass()
    else:
        raise "unknown alignment mode %s" % options.alignment_mode

    map_query2sbjct = alignlib.makeAlignmentVector()

    def __align( query_profile, sbjct_profile ):
        """align two profiles and output the result."""
        
        alignator.align( query_profile, sbjct_profile, map_query2sbjct )
        
        blocks = alignlib.AlignedBlocks( map_query2sbjct )
        
        if options.loglevel >= 3:
            options.stdlog.write( str(map_query2sbjct) )

        if map_query2sbjct.getLength() > 0:
            options.stdout.write("%s\t%s\t%i\t%s\n" % (
                    query, sbjct, map_query2sbjct.getScore(), str(blocks) ) )
            return 1

        return 0

    t_start = time.time()        
    def __report( noutput, ntotal ):

        global t_start
        if options.loglevel >= 1 and noutput % options.report_step == 0:
            t = time.time() - t_start
            options.stdlog.write( "# alignment: %5i (%5.2f)%%, query=%s, sbjct=%s, t=%i, <t>=%5.2fs, etf=%5.2fs, %5.2fh, et=%5.2fh\n" % \
                                      (noutput, 100.0 * noutput / ntotal,
                                       query, sbjct,
                                       t, 
                                       float(t)/noutput,
                                       float(t)/noutput * (ntotal-noutput),
                                       float(t)/noutput * (ntotal-noutput) / 3600,
                                       float(t)/noutput * ntotal / 3600) )
            options.stdlog.flush()
            options.stdout.flush()


    noutput = 0
    nempty = 0
    npairs = 0

    if options.filename_pairs:
        
        pairs = []
        infile = open( options.filename_pairs, "r" )
        for line in infile:
            if line[0] == "#": continue
            query, sbjct = line[:-1].split("\t")[:2]
            pairs.append( (query, sbjct) )
        infile.close()

        ntotal = len(pairs)

        if options.loglevel >= 1:
            options.stdlog.write( "# work: alignments=%i\n" % ( ntotal ) )
            options.stdlog.flush()
            
        last_query, last_sbjct = None, None
        for query, sbjct in pairs:

            if query != last_query:
                query_profile = plib_query.getProfile( query )
                last_query = query
            if sbjct != last_sbjct:
                sbjct_profile = plib_query.getProfile( sbjct )
                last_sbjct = sbjct
                
            npairs += 1
            if __align( query_profile, sbjct_profile ):
                noutput += 1
            else:
                nempty += 1
            __report( npairs, ntotal )

    else:

        query_keys, query_start, query_end = getKeys( plib_query, options.query_start, options.query_end )
        sbjct_keys, sbjct_start, sbjct_end = getKeys( plib_sbjct, options.sbjct_start, options.sbjct_end )

        ntotal = len(query_keys) * len(sbjct_keys)

        ## subtract half-diagonal for self-comparisons. If query_end is smaller than
        ## sbjct_start, the full square is computed
        if options.self_compare:
            d = max( query_end - sbjct_start, 0 ) 
            ntotal -= d * d / 2

        if options.loglevel >= 1:
            options.stdlog.write( "# work: queries=%i, sbjcts=%i, alignments=%i\n" % (len(query_keys), len(sbjct_keys), ntotal ) )
            options.stdlog.flush()

        for query in query_keys:

            query_profile = plib_query.getProfile( query )

            for sbjct in sbjct_keys:

                if options.self_compare and query > sbjct: continue

                sbjct_profile = plib_sbjct.getProfile( sbjct )
                
                npairs += 1
                if __align( query_profile, sbjct_profile ):
                    noutput += 1
                else:
                    nempty += 1

                __report( npairs, ntotal )
                
                break
            break

    if options.loglevel >= 1:
        t = time.time() - t_start
        options.stdlog.write( "# alignment: %5i (%5.2f)%%, t=%is, t=%ih\n" %\
                                  (noutput, 100.0 * noutput / ntotal,
                                   t, t / 3600.0 ) )


    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nempty=%i\n" % (ntotal, noutput, nempty) )

    #--------------------------------------------------------
    # general cleaning up
    Experiment.Stop()

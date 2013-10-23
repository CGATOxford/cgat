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
ProfileLibrary.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile, subprocess, optparse, time, math

#--------------------------------------------------------
#--------------------------------------------------------
#--------------------------------------------------------
# usage message
#--------------------------------------------------------
USAGE="""python %s [OPTIONS]

work with library of profiles.

Actions include: 
   * create: create a profile library
   * verify: verify a profile library
   * merge:  merge several profile libraries
   * split:  split a profile library into smaller sections

This is version $Id: ProfileLibrary.py 2781 2009-09-10 11:33:14Z andreas $.
""" % sys.argv[0]

#--------------------------------------------------------
#--------------------------------------------------------
#--------------------------------------------------------
# import of user libraries
#--------------------------------------------------------
import Experiment
import Mali
try: import alignlib_lite
except ImportError: pass

class ProfileLibrary:

    
    mSuffixDatabase = ".pdb"
    mSuffixIndex = ".pix"

    def __init__(self, name, mode = "r" ):
        self.mName = name

        self.mFilenameProfiles = self.mName + self.mSuffixDatabase
        self.mFilenameIndex = self.mName + self.mSuffixIndex
        self.mIndex = {}
        self.mWeightor = None
        self.mRegularizor = None

        self.mOutfileDatabase = None
        self.mOutfileIndex = None

        if mode == "r":
            self.__loadIndex()
        elif mode == "w":
            if os.path.exists( self.mFilenameProfiles):
                raise "profile database %s already exists." % self.mFilenameProfiles
            self.mOutfileDatabase = open( self.mFilenameProfiles, "wb" )
            self.mOutfileIndex = open( self.mFilenameIndex, "wb" )
        else:
            raise "unknown mode '%s'" % mode
        
    def __getitem__(self, key):
        return self.getProfile( key )

    def __len__(self):
        return len(self.mIndex )

    def __iter__(self):
        for name in self.mIndex.keys():
            yield (name, self[name])
            
    def keys(self):
        return list(self.iterkeys())

    def iterkeys( self ):
        for name in self.mIndex.iterkeys():
            yield name

    def iteritems( self ):
        for name in self.mIndex.keys():
            yield (name, self[name])

    def __del__(self):

        if self.mOutfileDatabase:
            self.mOutfileDatabase.close()
        if self.mOutfileIndex:
            self.mOutfileIndex.close()

    def setOptions( self, options ):
        """set options - access to command line options."""
        self.mLogLevel = options.loglevel
        self.mStdOut = options.stdout
        self.mStdLog = options.stdlog
        self.mStdErr = options.stderr
        
    def setWeightor( self, weightor ):
        """set the sequence weightor to use for profile creation."""
        if weightor == None or weightor == "none":
            self.mWeightor = alignlib_lite.py_makeNoWeightor() 
        elif weightor == "Henikoff" :
            self.mWeightor = alignlib_lite.py_makeWeightorHenikoff()
        elif weightor == "HenikoffKimmen":
            self.mWeightor = alignlib_lite.py_makeWeightorHenikoffKimmen()

    def __loadIndex( self ):

        if not os.path.exists( self.mFilenameProfiles):
            raise "profile database %s could not be found." % self.mFilenameProfiles

        if not os.path.exists( self.mFilenameIndex):
            raise "index %s could not be found." % self.mFilenameIndex
        
        infile = open( self.mFilenameIndex, "r" )
        self.mIndex = {}

        for line in infile:
            if line[0] == "#": continue
            name, first_pos, last_pos = line[:-1].split("\t")
            self.mIndex[name] = (int(first_pos), int(last_pos) )
            
        self.mInfileDatabase = open( self.mFilenameProfiles, "rb" )

    def appendProfile( self, name, profile ):
        """append a profile to this library."""

        start = self.mOutfileDatabase.tell()
        profile.save( self.mOutfileDatabase )
        self.mOutfileIndex.write( "%s\t%s\t%s\n" % (name, 
                                                    str(start),
                                                    str(self.mOutfileDatabase.tell()) ))

    def getProfile( self, name ):
        """append a profile to this library."""

        if name not in self.mIndex: raise KeyError

        self.mInfileDatabase.seek( self.mIndex[name][0] )
        return alignlib_lite.py_loadAlignandum( self.mInfileDatabase )

    def create( self, infile ):
        """create profile library from file."""

        self.mOutfileDatabase = open( self.mFilenameProfiles, "wb" )
        outfile_index = open( self.mFilenameIndex, "w" )

        ninput, noutput = 0, 0

        while mali.readFromFile( sys.stdin, format="profile" ):

            ninput += 1

            m = Mali.convertMali2Alignlib( mali )
            p = alignlib_lite.py_makeProfile( m, weightor = self.mWeightor )
            p.prepare()

            self.appendProfile( mali.getName(), p )
            
            noutput += 1

        return ninput, noutput

    def verify( self, infile ):
        """verify data in database against original data."""

        if not self.mIndex: self.__loadIndex()
        
        ninput, nfound, nnotfound, ndifferent = 0,0,0,0
        while mali.readFromFile( sys.stdin, format="profile" ):

            ninput += 1
            m = Mali.convertMali2Alignlib( mali )
            p1 = alignlib_lite.py_makeProfile( m )
            p1.prepare()
            
            p2 = self.getProfile( mali.getName() )

            if p1.getLength() != p2.getLength() or \
                    str(p1) != str(p2):
                ndifferent += 1
                continue
            
            nfound += 1

        return ninput, nfound, nnotfound, ndifferent

#--------------------------------------------------------
#--------------------------------------------------------
#--------------------------------------------------------
# main part of script
#--------------------------------------------------------
if __name__ == "__main__":

    #--------------------------------------------------------
    # command line parsing options
    parser = E.OptionParser( version = "%prog version: $Id: ProfileLibrary.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-p", "--prefix", dest="prefix", type="string",
                      help="prefix to use for the profile library." )

    parser.add_option("-a", "--action", dest="action", type="choice",
                      choices=("create", "verify", "merge", "split", "stats" ),
                      help="action to undertake." )

    parser.add_option("-w", "--weightor", dest="weightor", type="choice",
                      choices=("none", "Henikoff", "HenikoffKimmen" ),
                      help="sequence weightor to choose." )

    parser.set_defaults( prefix = "profiles", 
                         action = "create",
                         weightor = None,
                         )

    (options, args) = Experiment.Start( parser )

    #--------------------------------------------------------
    # main part of script
    mali = Mali.Mali()

    if options.action in ("create", "merge" ):
        mode = "w"
    else:
        mode = "r"

    plib = ProfileLibrary( options.prefix, mode )
    plib.setWeightor( options.weightor )

    if options.action == "verify":
        ninput, nfound, nnotfound, ndifferent = plib.verify(sys.stdin)
        if options.loglevel >= 1:
            options.stdlog.write( "# verify: ninput=%i, nfound=%i, nnotfound=%i, ndifferent=%i\n" % (ninput, nfound, nnotfound, ndifferent))

    elif options.action == "create":
        # create a new profile library
        ninput, noutput = plib.create( sys.stdin )
        if options.loglevel >= 1:
            options.stdlog.write( "# ninput=%i, noutput=%i\n" % (ninput, noutput) )

    elif options.action == "merge":
        for library in args:
            if options.loglevel >= 1:
                options.stdlog.write("# adding library %s\n" % library )
                options.stdlog.flush()

            other_lib = ProfileLibrary( library )

            for name, profile in other_lib.iteritems():
                plib.appendProfile( name, profile )

    elif options.action == "stats":
        options.stdout.write("profiles\t%i\n" % len(plib) )
        
    #--------------------------------------------------------
    # general cleaning up
    Experiment.Stop()

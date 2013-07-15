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
ProfileLibraryCompass.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import os, sys, string, re, tempfile, subprocess, optparse, time, math, shutil

#--------------------------------------------------------
#--------------------------------------------------------
#--------------------------------------------------------
# usage message
#--------------------------------------------------------
USAGE="""python %s [OPTIONS]

work with library of compass profiles.

Actions include: 
   * create: create a profile library
   * verify: verify a profile library
   * merge:  merge several profile libraries
   * split:  split a profile library into smaller sections

This is version $Id: ProfileLibraryCompass.py 2781 2009-09-10 11:33:14Z andreas $.
""" % sys.argv[0]

#--------------------------------------------------------
#--------------------------------------------------------
#--------------------------------------------------------
# import of user libraries
#--------------------------------------------------------
import Experiment
import Mali
from ProfileLibrary import ProfileLibrary

class CompassProfile:

    def __init__(self, profile ):
        if profile[0] != "#":
            raise "malformatted profile."
        
        pos = profile.find("\n")
        self.mName = profile[2:pos]
        self.mProfile = profile[pos+1:]
            
    def save( self, outfile ):
        outfile.write( "# %s\n" % self.mName )
        outfile.write( self.mProfile )
        
    def setName( self, name ):
        self.mName = name

    def getName( self ):
        return self.mName

class ProfileLibraryCompass( ProfileLibrary ):

    mIndexer = "mk_compass_db"
    
    mSuffixDatabase = ".cdb"
    mSuffixIndex = ".cix"

    def __init__(self, name, mode = "r" ):

        ProfileLibrary.__init__( self, name, mode )

        self.mMaxFilesPerChunk = 2000

    def getProfile( self, name ):
        """get a profile from this library.

        This function returns a compass profile as a string.
        """

        if name not in self.mIndex: raise KeyError

        first_pos, last_pos = self.mIndex[name]
        self.mInfileDatabase.seek( first_pos )
        return CompassProfile(self.mInfileDatabase.read( last_pos - first_pos ) )

    def __create( self, tempdir, nchunk ):
        """create chunk of a profile library from file.
        This is necessary due to a 2Gb file size limit in the compass
        executable that I downloaded.
        """
        
        # run mk_compass_db
        fn = "%s_%i" % (os.path.abspath(self.mFilenameProfiles ), nchunk)
        statement = "%s -i files.list -o %s" % (self.mIndexer, fn )
        s = subprocess.Popen( statement,
                              shell = True,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = tempdir,
                              close_fds = True)                              

        (out, err) = s.communicate()
        
        if s.returncode != 0:
            raise "Error in running %s \n%s\n%s\nTemporary directory in %s" % (statement, err, out, tempdir)

        ## clean up
        shutil.rmtree( tempdir )
        
        return fn

    def create( self, infile ):
        """create profile library from file."""
        
        ninput, noutput = 0, 0

        outfile_components = None

        # write multiple alignment as individual files
        nchunk = 0
        chunks = []

        while mali.readFromFile( infile, format="profile" ):

            ninput += 1
            
            if not outfile_components:
                tempdir = tempfile.mkdtemp()
                outfile_components = open( tempdir + "/files.list", "w" )
                
            fn = tempdir + "/%s" % mali.getName()
            outfile = open( fn, "w" )
            mali.writeToFile( outfile, format="plain-fasta" )
            outfile.close()
            outfile_components.write( "%s\n" % mali.getName() )
            noutput += 1            

            if ninput % self.mMaxFilesPerChunk == 0:
                outfile_components.close()
                if self.mLogLevel >= 2:
                    self.mStdLog.write("# building profiles for chunk %i at %i\n" % (nchunk, ninput) )
                    self.mStdLog.flush()

                chunks.append( self.__create( tempdir, nchunk ) )
                nchunk += 1
                outfile_components = None

        outfile_components.close()
        chunks.append( self.__create( tempdir, nchunk ) )
        
        # combine into a single file and remove temporary files
        if self.mLogLevel >= 2:
            self.mStdLog.write("# concatenating %i chunks\n" % (len(chunks) ))
            self.mStdLog.flush()

        outfile = open( os.path.abspath( self.mFilenameProfiles), "w" )
        for filename in chunks:
            infile = open(filename, "r" )
            for line in infile:
                outfile.write( line )
            infile.close()
            os.remove( filename )
        outfile.close()
        
        # index the profile database
        if self.mLogLevel >= 2:
            self.mStdLog.write("# indexing the database\n" )
            self.mStdLog.flush()

        self.index()

        return ninput, noutput

    def index( self ):
        """build index for database."""

        outfile_index = open( self.mFilenameIndex, "w" )
        infile_database = open( self.mFilenameProfiles, "r" )
        
        pos = infile_database.tell()
        last_name = None
        skip = False
        while 1:

            pos = infile_database.tell()
            line = infile_database.readline()
            if not line: break

            if line[0] == "#": 
                if skip: 
                    skip = False
                    continue
                
                if last_name:
                    outfile_index.write( "%s\t%i\t%i\n" % (last_name, first_pos, pos ) )
                last_name = line[2:-1]
                skip = True
                first_pos = pos

        if last_name:
            pos = infile_database.tell()
            outfile_index.write( "%s\t%i\t%i\n" % (last_name, first_pos, pos ) )
            
        outfile_index.close()
        infile_database.close()
        
    def verify( self, infile ):
        """verify data in database against original data."""

        if not self.mIndex: self.__loadIndex()
        
        ninput, nfound, nnotfound, ndifferent = 0,0,0,0
        while mali.readFromFile( sys.stdin, format="profile" ):

            ninput += 1
            
            raise "not implemented"
            
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
    parser = E.OptionParser( version = "%prog version: $Id: ProfileLibraryCompass.py 2781 2009-09-10 11:33:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option("-p", "--prefix", dest="prefix", type="string",
                      help="prefix to use for the profile library." )

    parser.add_option("-a", "--action", dest="action", type="choice",
                      choices=("create", "verify", "merge", "split", "stats", "index" ),
                      help="action to undertake." )

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

    plib = ProfileLibraryCompass( options.prefix, mode )
    plib.setOptions( options )

    if options.action == "verify":
        ninput, nfound, nnotfound, ndifferent = plib.verify(sys.stdin)
        if options.loglevel >= 1:
            options.stdlog.write( "# verify: ninput=%i, nfound=%i, nnotfound=%i, ndifferent=%i\n" % (ninput, nfound, nnotfound, ndifferent))

    elif options.action == "create":
        # create a new profile library
        ninput, noutput = plib.create( sys.stdin )
        if options.loglevel >= 1:
            options.stdlog.write( "# ninput=%i, noutput=%i\n" % (ninput, noutput) )

    elif options.action == "index":
        plib.index()

    elif options.action == "merge":
        for library in args:
            if options.loglevel >= 1:
                options.stdlog.write("# adding library %s\n" % library )
                options.stdlog.flush()

            other_lib = ProfileLibraryCompass( library )

            for name, profile in other_lib.iteritems():
                plib.appendProfile( name, profile )

    elif options.action == "stats":
        options.stdout.write("profiles\t%i\n" % len(plib) )
        
    #--------------------------------------------------------
    # general cleaning up
    Experiment.Stop()

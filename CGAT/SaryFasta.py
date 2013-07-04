################################################################################
#
#   Gene prediction pipeline 
#
#   $Id: SaryFasta.py 2784 2009-09-10 11:41:14Z andreas $
#
#   Copyright (C) 2007 Andreas Heger
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
#
#################################################################################
"""
SaryFasta.py - index fasta files by suffix array
===================================================

Subroutines for working on I/O of large genomic files.

Index a fasta file to retrieve sequences by suffix-array fragment
search.

python SaryFasta.py [options] name [ files ]
"""

import os, sys, array, string, re, types, optparse, time, struct, hashlib, base64, shutil, subprocess
import glob, random

##------------------------------------------------------------
def getHID ( sequence ):
    """returns a hash identifier for a sequence.
    """
    
    # do the encryption
    h = hashlib.md5(sequence).digest()
    
    # map to printable letters: hid has length 22, so the padded '=' are
    # truncated. You have to add them, if you ever want to decode,
    # but who would do such a thing :=)

    r = base64.encodestring(h)[0:22]
    
    # finally substitute some characters:
    # '/' for '_', so we have legal file names
    # '[' for '+' and ']' for '=' for internet-applications
    
    hid = string.replace(r  , '/', '_')
    hid = string.replace(hid, '+', '[')
    hid = string.replace(hid, '=', ']')
                                                
    return hid

##------------------------------------------------------------
def createDatabase( db, filenames,
                    buf_size=400000000,
                    force = False,
                    regex_identifier = None):
    """index files in filenames to create database.

    buf_size: buffer size for a sary chunk.

    Two new files are created - db.fasta and db_name.idx
    
    regex_identifier: pattern to extract identifier from description line.
    If None, the part until the first white-space character is used.
    """

    index_name = db + ".idx"
    db_name = db + ".dat"
    if db in filenames:
        raise ValueError( "database (%s) is part of input set." % db)

    if os.path.exists( db_name ) and not force:
        raise ValueError( "database %s already exists." % db)

    if os.path.exists( index_name ) and not force:
        raise ValueError( "database index %s already exists." % index_name )
    
    outfile_index = open( index_name, "w" )

    os.mkdir( db_name )
    files_to_index = []
    fn = db_name + "/part%i" % len(files_to_index)
    outfile_src = open( fn, "w" )
    files_to_index.append(fn)
    
    if type(filenames) == types.StringType:
        filenames = [filenames]

    identifiers = {}
    lsequence = 0
    identifier_pos, sequence_pos = 0, 0
    noutput = 0
    
    for filename in filenames:

        if filename[:-3] == ".gz":
            infile = gzip.open(filename, "r")
        else:
            infile = open( filename, "r")

        fragments = []
        first = True
        
        for line in infile:

            if line[0] == "#":  continue
            
            if line[0] == ">" :
                
                if not first:
                    
                    f = "".join(fragments)
                    
                    if len(f) + noutput > buf_size:
                        outfile_src.close()
                        fn = db_name + "/part%i" % len(files_to_index)                        
                        outfile_src = open( fn, "w" )
                        files_to_index.append(fn)
                        noutput = 0

                    outfile_src.write( f + "\n" )
                    outfile_index.write("%s\t%s\n" % (getHID(f), identifier) )

                    noutput += len(f)
                    fragments = []

                first = False
                
                if regex_identifier:
                    try:
                        identifier = re.search(regex_identifier, line[1:-1]).groups()[0]
                    except AttributeError:
                        raise "could not parse identifer from line %s" % line[1:-1]
                else:
                    identifier = re.split("\s", line[1:-1])[0]
                    
                ## check for duplicate identifiers
                if identifier in identifiers:
                    raise ValueError, "%s occurs more than once in %s and %s: line=%s" %\
                          (identifier, identifiers[identifier], filename, line[1:-1])
                identifiers[identifier] = filename
                
            else:
                s = re.sub( "\s", "", line.strip() )
                fragments.append( s )
        
        f = "".join(fragments)
        outfile_src.write( f + "\n" )
        outfile_index.write("%s\t%s\n" % (getHID(f), identifier) )

        infile.close()
                    
    ## build indicies
    for filename in files_to_index:
        
        try:
            retcode = subprocess.Popen( "mksary --quiet %s" % (filename), shell=True)
            
            if retcode < 0:
                print >>sys.stderr, "Child was terminated by signal", -retcode
            else:
                print >>sys.stderr, "Child returned", retcode

        except OSError, e:
            print >>sys.stderr, "Execution failed:", e

class SaryFasta:

    def __init__( self, dbname ):

        self.mNameIndex = dbname + ".idx"
        if not os.path.exists( self.mNameIndex ):
            raise "index file %s does not exist." % (self.mNameIndex)

        self.mNameData = dbname + ".dat"
        if not os.path.exists( self.mNameData ):
            raise "data dir %s does not exist." % (self.mNameData)
        
        self.mDbname = dbname
        self.mIsLoaded = False
        self.mFilenames = map( lambda x: x[:-4], glob.glob( self.mNameData + "/part*.ary"))
        
        if len(self.mFilenames) == 0:
            raise "empty database."

    def __getitem__(self, pattern ):
        """return full length sequence."""
        return self.search( pattern )
        
    def __loadIndex( self ):
        """load complete index into memory."""

        self.mIndex = {}

        for line in open(self.mNameIndex, "r"):

            hid, identifier = line[:-1].split("\t")
            
            self.mIndex[hid] = identifier

        self.mIsLoaded = True
            
    def getDatabaseName( self ):
        """returns the name of the database."""
        return self.mDbname
    
    def search( self, pattern ):
        """search indices with a pattern.
        """

        if not self.mIsLoaded: self.__loadIndex()
        
        result = []
        for filename in self.mFilenames:

            cmd = "sary '%s' %s" % (pattern, filename)
            for line in subprocess.Popen( cmd, shell=True, stdout=subprocess.PIPE).stdout:
                hid = getHID( line[:-1] )
                if hid not in self.mIndex:
                    raise "inconsistency: can not find identifier for hid %s in index.\n" % (hid)
                
                result.append( self.mIndex[hid] )

        return result

## Test function for benchmarking purposes
def benchmarkRandomFragment( fasta, size ):
    """returns a random fragment of size."""

    contig, strand, start, end = fasta.getRandomCoordinates( size )
    s = fasta.getSequence( contig, strand, start, end )
    return s

def verify( reference, fasta, 
            num_iterations, fragment_size,
            stdout = sys.stdout, quiet = False ):
    """verify two databases.

    Get segment from fasta and check for presence in fasta2.
    """
    if not quiet:
        options.stdout.write("verifying %s with %s using %i random segments per sequence of length %i\n" %\
                             (fasta.getDatabaseName(),
                              reference,
                              num_iterations,
                              fragment_size ))
        options.stdout.flush()
    nerrors = 0
    f = []
    for line in open(reference,"r"):
        if line[0] == ">":
            for x in range(num_iterations):
                s = "".join(f)
                if len(s) > fragment_size:
                    c = random.randint( 0, len(s) - fragment_size)
                    pattern = s[c:c+fragment_size]
                    found = set(fasta.search( pattern ))
                    if id not in found:
                        nerrors += 1
                        if not quiet:
                            options.stdout.write("can't find pattern '%s' for %s in %s\n" %\
                                                     (pattern, id, str(found)))
            f = []
            id = re.search(">(\S+)", line[:-1]).groups()[0]
        else:
            f.append( re.sub("\s", "", line[:-1]) )
    return nerrors

if __name__ == "__main__":

    import Experiment

    parser = E.OptionParser( version = "%prog version: $Id: SaryFasta.py 2784 2009-09-10 11:41:14Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "-e", "--extract", dest="extract", type="string",
                       help="extract region ( for testing purposes. Format is contig:strand:from:to." )

    parser.add_option( "-c", "--compression", dest="compression", type="choice",
                       choices=("lzo", "zlib", "gzip", "dictzip", "bzip2", "debug"),
                       help="compress database." )

    parser.add_option( "--random-access-points", dest="random_access_points", type="int",
                       help="save random access points every # number of nucleotides." )
    

    parser.add_option( "-f", "--input-format", dest="input_format", type="choice",
                       choices=("one-forward-open", "zero-both-open" ),
                       help="coordinate format of input." )

    parser.add_option( "-s", "--synonyms", dest="synonyms", type="string",
                       help="list of synonyms, comma separated with =, for example, chr1=chr1b" )

    parser.add_option( "-b", "--benchmark", dest="benchmark", action="store_true",
                       help="benchmark read access." )
    
    parser.add_option( "--benchmark-num-iterations", dest="benchmark_num_iterations", type="int",
                       help="number of iterations for benchmark [%DEFAULT%]." )

    parser.add_option( "--benchmark-fragment-size", dest="benchmark_fragment_size", type="int",
                       help="benchmark: fragment size [%DEFAULT%]." )

    parser.add_option( "--verify", dest="verify", type="string",
                       help="verify against other database.")

    parser.add_option( "-a", "--clean-sequence", dest="clean_sequence", action="store_true",
                       help="remove X/x from DNA sequences - they cause errors in exonerate." )

    parser.add_option( "--regex-identifier", dest="regex_identifier", type="string",
                       help="regular expression for extracting the identifier from fasta description line." )

    parser.set_defaults(
        extract = None,
        input_format = "zero-both-open",
        benchmark_fragment_size = 1000,
        benchmark_num_iterations = 1000000,
        benchmark = False,
        compression = None,
        random_access_points = 0,
        synonyms = None,
        verify = None,
        verify_num_iterations = 1,
        verify_fragment_size = 20,
        clean_sequence = False,
        regex_identifier = None)
    
    (options, args) = Experiment.Start( parser )

    if options.extract:
        fasta = SaryFasta( args[0] )
        converter = getConverter( options.input_format )
        
        contig, strand, start, end = options.extract.split(":")
        start, end = map( int, (start, end) )
        sequence = fasta.getSequence( contig, strand,
                                      start, end,
                                      converter = converter )
        options.stdout.write( ">%s\n%s\n" % \
                              ( options.extract, sequence ) )
    elif options.benchmark:
        import timeit
        timer = timeit.Timer( stmt="benchmarkRandomFragment( fasta = fasta, size = %i)" % (options.benchmark_fragment_size),
                              setup="""from __main__ import benchmarkRandomFragment,SaryFasta\nfasta=SaryFasta( "%s" )""" % (args[0] ) )

        t = timer.timeit( number = options.benchmark_num_iterations )
        options.stdout.write("iter\tsize\ttime\n" )
        options.stdout.write("%i\t%i\t%i\n" % (options.benchmark_num_iterations, options.benchmark_fragment_size, t ) )
    elif options.verify:
        fasta = SaryFasta( args[0] ) 
        nerrors1 = verify( options.verify, 
                           fasta,
                           options.verify_num_iterations,
                           options.verify_fragment_size,
                           stdout=options.stdout )
        options.stdout.write("errors=%i\n" % (nerrors1) )        
    else:
        if options.loglevel >= 1:
            options.stdlog.write("# creating database %s\n" % args[0])            
            options.stdlog.write("# indexing the following files: \n# %s\n" %\
                                 (" \n# ".join( args[1:] ) ))
            options.stdlog.flush()

        if len(args) < 2:
            print USAGE
            sys.exit(1)
            
        createDatabase( args[0], args[1:], 
                        regex_identifier = options.regex_identifier )
    
    Experiment.Stop()

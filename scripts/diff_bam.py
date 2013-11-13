'''
diff_bam.py - compare multiple bam files against each other
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS BAM Comparison

Purpose
-------

Compare RNASeq reads in multiple BAM files against each other.

.. note::
    BAM files need to be sorted by read name. samtools sort
    does NOT work as it uses a custom comparison function
    (strnum_cmp) that is incompatible with the standard 
    lexicographical order in python. See the example below
    on how to get sorted files.

This script is for validation purposes. It might take a while
for large BAM files.

Usage
-----

If you have two sorted :term:`sam` or :term:`bam` formatted
files, type::

   python diff_bam.py a.bam b.bam > out

If they are not sorted, you can use samtools sort to do an
inplace sort::

   python diff_bam.py <( samtools view -h a.bam | hsort 0 -k1,1) 
                  <( samtools view -h b.bam | hsort 0 -k1,1) 

The ``-h`` option outputs the header, and the hsort command sorts without 
disturbing the header. 

Type::

   python diff_bam.py --help

for command line help.

Documentation
-------------

For read counts to be correct the NH flag to be set correctly.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import collections
import itertools
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam
import CGAT.GTF as GTF

class multiway_groupby(object):
    # [k for k, g in groupby('AAAABBBCCDAABBB')] --> A B C D A B
    # [list(g) for k, g in groupby('AAAABBBCCD')] --> AAAA BBB CC D
    def __init__(self, iterables, key=None):
        # set up iterators
        self.it = [ itertools.groupby(iterable, key) for iterable in iterables ]

        # todo: catch StopIterator for empty iterables
        self.current = [ next(self.it[x]) for x in xrange(len(iterables)) ]

        # decide on target key
        self.targetkey = min( [x[0] for x in self.current if x[0] != None ] )

    def __iter__(self):
        return self

    def next(self):

        # check if all iterators exhausted
        if not any( self.it ): raise StopIteration

        # yield all that correspond to target key and update
        result = []
        targetkey = self.targetkey
        updated = 0

        for x, y in enumerate( self.current ):
            key, val = y
            if key and key == targetkey: 
                # save result (instantiate to prevent overwriting)
                # in next method
                result.append( list(val) )
                # advance
                try:
                    self.current[x] = next(self.it[x])
                except StopIteration:
                    self.it[x] = None
                    self.current[x] = (None,None)
                updated += 1
            else: 
                # return empty result
                result.append( [] )

        assert updated > 0, "no updates - infinite loop"

        # decide which is target key
        try: 
            self.targetkey = min( [x[0] for x in self.current if x[0] != None ] )
        except ValueError:
            # if all are None, sequence is empty
            self.targetkey = None
        
        return targetkey, result

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option( "-e", "--headers", dest="headers", type="string",
                       help = "',' separated list of labels used as headers. "
                              " Should correspond in order to command line arguments [%default]" )

    parser.set_defaults(
        headers = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv, add_output_options = True )

    read_length = 50
    max_distance = read_length

    if len(args) < 2:
        raise ValueError( "please specify at least two BAM files" )

    infiles = []
    for arg in args: infiles.append( pysam.Samfile( arg, 'rb' ) )

    if options.headers:
        headers = options.headers.split(",")
        if len(headers) != len(args):
            raise ValueError("number of headers and files differrent")
    else:
        headers = [ "file%i" % x for x in range( 1, len(infiles) + 1) ]

    options.stdout.write( "read\tnlocations\tnmatched\t%s\t%s\n" % \
                              ("\t".join( ["%s_nh" % x for x in headers] ),
                               "\t".join( ["%s_loc" % x for x in headers] ) ) )

    ninput, noutput = 0, 0
    for readname, result in multiway_groupby( infiles, key = lambda x: x.qname ):
        ninput += 1
        nmatches = 0
        locations = set()

        # build location directory
        for r in result: 
            for rr in r: 
                if rr.is_unmapped: continue
                locations.add( (rr.tid, rr.pos) )

        # permit a certain fuzzyness in locations
        # as some mappers clip and others don't.
        locations = list(locations)
        locations.sort()
        pos2loc = {}
        last_tid, last_pos, nlocations = None, 0, -1
        for tid, pos in locations:
            if tid != last_tid or pos - last_pos > max_distance:
                nlocations += 1
            pos2loc[(tid,pos)] = nlocations
            last_tid, last_pos = tid, pos
        nlocations += 1

        # code locations
        codes, nh = [], []
        for r in result:
            c = []
            for rr in r:
                if rr.is_unmapped: continue
                c.append( pos2loc[ (rr.tid, rr.pos) ] )
            codes.append( ",".join( map(str, sorted(c)) ) )
            nh.append( len(c) )
            if len(c): nmatches += 1

        noutput += 1
        options.stdout.write( "%s\t%i\t%i\t%s\t%s\n" % (readname, 
                                                        nlocations,
                                                        nmatches, 
                                                        "\t".join([ "%i" % x for x in nh ]),
                                                        "\t".join( codes ) ) )

    E.info( "ninput=%i, noutput=%i" % (ninput,noutput) )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

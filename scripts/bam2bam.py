'''
bam2bam.py - modify bam files
=============================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Genomics NGS BAM Manipulation

Purpose
-------

This script reads a :term:`bam` formatted file from stdin and outputs a
modified :term:`bam` formatted file on stdout.

.. note::
   You need to redirect logging information to a file or turn it off
   via -v 0 in order to get a valid sam/bam file.

Usage
-----

For example::

   python bam2bam.py --filter=mapped < in.bam > out.bam

will remove all unmapped reads from the bam-file. To remove
unmapped reads from multiple bam-files, try::

   python bam2bam.py --inplace --filter=mapped *.bam

Type::

   python bam2bam.py --help

for command line help.

Documentation
-------------

The script implements the following methods:

``set-nh``
   set the NH flag. Some tools (bowtie, bwa) do not set the NH flag. 
   If set, this option will set the NH flag (for mapped reads).
   This option requires the bam/sam file to be sorted by read name.

``unset-unmapped_mapq`` 
   some tools set the mapping quality of unmapped reads. This
   causes a violation in the Picard tools.

``remove-better``
   remove alignments that are worse than alignments present in
   an alternative :term:`bam` formatted file (``--filter-bam``). 

``filter``
   remove alignments based on a variety of flags.

``strip`` 
   remove the sequence and/or quality scores from all reads in
   a bam-file. If ``match``, sequence and quality are stripped 
   for alignments without mismatches. The latter requires the
   ``NM`` tag. If not present, nothing will be stripped. Note 
   that stripping is not reversible if the read names are not
   unique.

``set-sequence``
   set the sequence and quality scores in the bam file to some 
   dummy sequence. Necessary for some tools that can not work
   with bam-files without sequence.

``unstrip``
   add sequence and quality scores back to a bam file. The
   argument to this option is a :term:`fastq` formatted file
   with the sequences and quality scores to insert.

By default, the script works from stdin and outputs to stdout.
If the ``--inplace option`` is given, the script will modify
bam-files given as command line arguments in-place - a temporary
copy will be written to :file:`/tmp` and once completed copied
over the original file. The script will automatically re-index
the modified files.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import collections
import itertools
import tempfile
import shutil

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import pysam

try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _bam2bam
except ImportError:
    import CGAT._bam2bam as _bam2bam

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                             usage = globals()["__doc__"] )

    parser.add_option( "--set-nh", dest="set_nh", action="store_true",
                       help = "sets the NH flag. The file needs to be sorted by readname [%default]" )

    parser.add_option( "--unset-unmapped-mapq", dest="unset_unmapped_mapq", action="store_true",
                       help = "sets the mapping quality of unmapped reads to 0 [%default]" )

    parser.add_option( "--set-sequence", dest="set_sequence", action="store_true",
                       help = "sets the sequence to 'A's (a valid base) and the quality to 'F's"
                       ",which is defined in all fastq scoring schemes "
                       "[%default]" )

    parser.add_option( "--strip", dest="strip", type="choice",
                       choices = ("sequence", "quality", "match" ),
                       help = "remove parts of the bam-file. Note that stripping the sequence will "
                       "also strip the quality values [%default]" )

    parser.add_option( "--unstrip", dest="unstrip", action="store_true",
                       help = "add sequence and quality into bam file [%default]" )

    parser.add_option( "--filter", dest="filter", action="append", type="choice",
                       choices=('NM', 'CM', 'mapped', 'unique', "non-unique" ),
                       help = "filter bam file. The option denotes the property that is  "
                       "used to determine better match [%default]" )

    parser.add_option( "--reference-bam", dest="reference_bam", type="string",
                       help = "bam-file to filter with [%default]" )
    
    parser.add_option( "--sam", dest="output_sam", action="store_true",
                       help = "output in sam format [%default]" )

    parser.add_option( "--inplace", dest="inplace", action="store_true",
                       help = "modify bam files in-place. Bam files need to be given "
                       "as arguments. Temporary bam files are written to /tmp [%default]" )

    parser.add_option( "--fastq1", "-1", dest="fastq_pair1", type="string",
                       help = "fastq file with read information for first in pair or unpaired [%default]" )

    parser.add_option( "--fastq2", "-2", dest="fastq_pair2", type="string",
                       help = "fastq file with read information for second in pair [%default]" )


    parser.set_defaults(
        filter = [],
        set_nh = False,
        unset_unmapped_mapq = False,
        output_sam = False,
        reference_bam = None,
        strip = None,
        unstrip = None,
        set_sequence = False,
        inplace = False,
        fastq_pair1 = None,
        fastq_pair2 = None,
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if options.inplace:
        bamfiles = args
        if len(bamfiles) == 0:
            raise ValueError("please one or more bam-files as command line arguments" )

        if "-" in bamfiles:
            raise ValueError("can not read from stdin if ``--inplace`` is selected" )
    else:
        bamfiles = ["-"]

    for bamfile in bamfiles:

        E.info('processing %s' % bamfile )

        # reading bam from stdin does not work with only the "r" tag
        pysam_in = pysam.Samfile( bamfile, "rb" )

        if bamfile == "-":
            if options.output_sam:
                pysam_out = pysam.Samfile( "-", "wh", template = pysam_in )
            else:
                pysam_out = pysam.Samfile( "-", "wb", template = pysam_in )
        else:
            tmpfile = tempfile.NamedTemporaryFile( delete=False, prefix = "ctmp" )
            tmpfile.close()

            E.debug("writing temporary bam-file to %s" % tmpfile.name )
            pysam_out = pysam.Samfile( tmpfile.name, "wb", template = pysam_in )

        if options.filter:

            remove_mismatches, colour_mismatches = False, False

            if "NM" in options.filter:
                remove_mismatches = True

            elif "CM" in options.filter:
                remove_mismatches = True
                colour_mismatches = True

            if remove_mismatches:
                if not options.reference_bam:
                    raise ValueError( "requiring reference bam file for removing by mismatches")

                pysam_ref = pysam.Samfile( options.reference_bam, "rb" )
            else:
                pysam_ref = None

            # filter and flags are the opposite way around
            c = _bam2bam.filter_bam( pysam_in, pysam_out, pysam_ref,
                                     remove_nonunique = "unique" in options.filter,
                                     remove_unique = "non-unique" in options.filter,
                                     remove_contigs = None,
                                     remove_unmapped = "mapped" in options.filter,
                                     remove_mismatches = remove_mismatches,
                                     colour_mismatches = colour_mismatches )

            options.stdlog.write( "category\tcounts\n%s\n" % c.asTable() )

        else:
            # set up the modifying iterators
            it = pysam_in.fetch( until_eof = True )

            if options.unset_unmapped_mapq:
                def unset_unmapped_mapq( i ):
                    for read in i:
                        if read.is_unmapped:
                            read.mapq = 0
                        yield read
                it = unset_unmapped_mapq( it )

            if options.set_nh and False:
                def set_nh( i ):

                    for key, reads in itertools.groupby( i, lambda x: x.qname ):
                        l = list(reads)
                        nh = len(l)
                        for read in l:
                            if not read.is_unmapped:
                                t = dict( read.tags )
                                t['NH'] = nh
                                read.tags = list(t.iteritems())
                            yield read
                it = set_nh( it )

            if options.set_sequence:
                def set_sequence( i ):
                    for read in i:
                        # can't get at length of unmapped reads
                        if read.is_unmapped: 
                            read.seq = "A"
                            read.qual = "F"
                        else:
                            read.seq = "A" * read.inferred_length
                            read.qual = "F" * read.inferred_length

                        yield read
                it = set_sequence( it )

            if options.strip != None:
                def strip_sequence( i ):
                    for read in i:
                        read.seq = None
                        yield read

                def strip_quality( i ):
                    for read in i:
                        read.qual = None
                        yield read

                def strip_match( i ):
                    for read in i:
                        try:
                            nm = read.opt( 'NM' )
                        except KeyError:
                            nm = 1
                        if nm == 0:
                            read.seq = None
                        yield read

                if options.strip == "sequence":
                    it = strip_sequence( it )
                elif options.strip == "quality":
                    it = strip_quality( it )
                elif options.strip == "match":
                    it = strip_match( it )

            if options.unstrip:
                def buildReadDictionary( filename ):
                    if not os.path.exists( filename ):
                        raise OSError( "file not found: %s" % filename )
                    fastqfile = pysam.Fastqfile( filename )
                    fastq2sequence = {}
                    for x in fastqfile:
                        if x.name in fastq2sequence: 
                            raise ValueError("read %s duplicate - can not unstrip" % x.name )

                        fastq2sequence[x.name] = (x.sequence, x.quality )
                    return fastq2sequence
            
                if not options.fastq_pair1:
                    raise ValueError("please supply fastq file(s) for unstripping")
                fastq2sequence1 = buildReadDictionary( options.fastq_pair1 )
                if options.fastq_pair2:
                    fastq2sequence2 = buildReadDictionary( options.fastq_pair2 )
                    
                def unstrip_unpaired( i ):
                    for read in i:
                        read.seq, read.qual = fastq2sequence1[read.qname]
                        yield read
                
                def unstrip_pair( i ):
                    for read in i:
                        if read.is_read1:
                            read.seq, read.qual = fastq2sequence1[read.qname]
                        else: 
                            read.seq, read.qual = fastq2sequence2[read.qname]
                        yield read

                if options.fastq_pair2:
                    it = unstrip_pair( it )
                else:
                    it = unstrip_unpaired( it )

            if options.set_nh:
                it = _bam2bam.SetNH( it )

            # read and output
            for read in it: 
                pysam_out.write( read )

            pysam_in.close()
            pysam_out.close()

        if options.inplace:
            # set date and file permissions according to original
            # Note: currently it will not update user and group.
            original = os.stat(bamfile)
            os.utime( tmpfile.name, (original.st_atime, original.st_mtime) )
            os.chmod( tmpfile.name, original.st_mode ) 
            # move new file over original copy
            shutil.move( tmpfile.name, bamfile )
            # re-index
            pysam.index( bamfile )

    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )


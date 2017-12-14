'''bam2bam.py - modify bam files
=============================

:Tags: Genomics NGS BAM Manipulation

Purpose
-------

This script reads a :term:`bam` formatted file from stdin, performs an
action (see methods below) then outputs a modified :term:`bam`
formatted file on stdout.

.. note::
   You need to redirect logging information to a file (via -L) or turn it off
   via -v 0 in order to get a valid sam/bam file.

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

``filter``

   remove alignments based on a variety of flags. The filtering method
   is determined by the ``--filter-method`` option. These may be
   ``unique``, ``non-unique``, ``mapped``, ``NM`` or ``CM``.  If
   ``unique`` is set, only uniquely mapping reads will be output. If
   ``non-unique`` is set then only multi-mapping reads will be
   output. This method first checks for the NH flag - if set, a unique
   match should have at most NH=1 hits.  If not set, the method checks
   for BWA flags. Currently it checks if X0 is set (X0=Number of best
   hits found by BWA).  If ``mapped`` is given, unmapped reads will be
   removed. If ``NM`` or ``CM`` is set, the alignment of reads in two
   sam files (input and reference) is compared and only reads with a
   lower number of mismatches in the input compared to the reference
   sam file will be kept. If ``CM`` is set, the colourspace mismatch
   tag (for ABI Solid reads) will be used to count differences to the
   reference sam file. By default, the ``NM`` (number of mismatches)
   tag is used. The tag that is used needs to present in both input
   sam file and the reference sam file. If ``unique`` is given this
   wil NOT remove any unmapped reads.  This can be achieved by
   providing the ``filter`` option twice, once each with ``mapped``
   and ``unique``.

   .. note::

      The filter methods can't currently combined with any of
      the other methods - this is work in progress.

``strip-sequence``

   remove the sequence from all reads in a bam-file. Note that
   stripping the sequence will also remove the quality scores.
   Stripping is not reversible if the read names are not unique.

``strip-quality``

   remove the quality scores from all reads in a bam-file.
   Stripping is not reversible if the read names are not unique.

``set-sequence``

   set the sequence and quality scores in the bam file to some dummy
   values ('A' for sequence, 'F' for quality which is a valid score in
   most fastq encodings. Necessary for some tools that can not work
   with bam-files without sequence.

``unstrip``

   add sequence and quality scores back to a bam file. Requires a
   :term:`fastq` formatted file with the sequences and quality scores
   to insert.

``unset-unmapped-mapq``

   sets the mapping quality of unmapped reads to 0.

``keep-first-base``

   keep only the first base of reads so that read counting tools will
   only consider the first base in the counts

``downsample-single``

   generates a downsampled :term:`bam` file by randomly subsampling
   reads from a single ended :term:`bam` file. The downsmpling
   retains multimapping reads. The use of this requires downsampling
   parameter to be set and optionally randomseed.

``downsample-paired``

   generates a downsampled :term:`bam` file by randomly subsampling
   reads from a paired ended :term:`bam` file. The downsampling
   retains multimapping reads. The use of this requires downsampling
   parameter to be set and optionally randomseed.


By default, the script works from stdin and outputs to stdout.
If the ``--inplace option`` is given, the script will modify
bam-files given as command line arguments in-place - a temporary
copy will be written to :file:`/tmp` and once completed copied
over the original file. The script will automatically re-index
the modified files.

Usage
-----

For example::

   cgat bam2bam --method=filter --filter-method=mapped < in.bam > out.bam

will remove all unmapped reads from the bam-file.

To remove unmapped reads from multiple bam-files, try::

   cgat bam2bam --inplace --method=filter --filter-method=mapped *.bam

example for running downsample::

   cgat bam2bam --method=downsample-paired --downsample=30000
        --randomseed=1 -L out.log < Paired.bam > out.bam

Type::

   cgat bam2bam --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import tempfile
import shutil
import random
import itertools
import pysam
import CGAT.Experiment as E
import CGAT.IOTools as IOTools


try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _bam2bam
except ImportError:
    import CGAT.scripts._bam2bam as _bam2bam


class SetNH:

    def __init__(self, iter):
        self.iter = itertools.groupby(iter, lambda x: x.qname)
        self.stack = []

    def __iter__(self):
        return self

    def __next__(self): 
        """python version of next().
        """
        
        while 1:
            if self.stack:
                return self.stack.pop(0)
            else:
                key, x = next(self.iter)
                self.stack = list(x)
                nh = len(self.stack)
                for read in self.stack:
                    if not read.is_unmapped:
                        # deal with paired end reads counted
                        # as multi-mapping
                        if read.is_proper_pair and nh > 1:
                            nh -= 1
                        read.set_tag("NH", nh)


class SubsetBam(object):

    ''' base class for performing downsampling on single and
    paired bam file

    A dictionary of the read names is made and then an
    array matching the numbers of reads in the dictionary will be
    created. This array contains 1s that match the number of reads
    to be downsampled to and 0s to fill out the rest of the array.

    The read is kept if the value in the dictionary matches 1. The
    read is yielded and iterated over to produce the output bam.

    The script will handle multimapping reads.
    '''

    def __init__(self, pysam_in, downsample, paired_end=None,
                 single_end=None, random_seed=None):

        self.pysam_in1, self.pysam_in2 = itertools.tee(pysam_in)
        self.downsample = downsample
        self.paired_end = paired_end
        self.single_end = single_end
        self.random_seed = random_seed

    def list_of_reads(self, paired=None):

        '''
        This will create a dictionary of uniqe reads in the bam
        '''
        read_list = []

        for read in self.pysam_in1:
            if paired is True:
                if read.is_proper_pair:
                    read_list.append(read.qname)
            else:
                read_list.append(read.qname)

        return sorted(set(read_list))

    def downsample_paired(self):

        '''
        This function will downsample a paired bam file.
        It will retain multimapping reads if they have not been
        pre-filtered
        '''
        if self.random_seed is not None:
            random.seed(self.random_seed)

        collect_list = self.list_of_reads(paired=True)
        read_list = random.sample(collect_list, self.downsample)

        if self.downsample == len(collect_list):
            E.warn('''The downsample reads is equal to the
            number of unique reads''')
            for read in self.pysam_in2:
                yield read

        # yield read if it is in read_list, if the read is multimapped
        # then all multimapping reads will be yielded
        else:
            E.warn('''Multimaping reads have been detected and these will
            be output to the final bam file''')
            for read in self.pysam_in2:
                if read.qname in read_list:
                    yield read

    def downsample_single(self):

        '''
        This function will downsample a single bam file.
        It will retain multimapping reads if not pre-filtered
        '''
        if self.random_seed is not None:
            random.seed(self.random_seed)

        collect_list = self.list_of_reads(paired=False)
        read_list = random.sample(collect_list, self.downsample)

        if self.downsample == len(collect_list):
            E.warn('''The downsample reads is equal to the
            number of unique reads''')
            for read in self.pysam_in2:
                yield read

        # yield read if it is in read_list, if the reads is multimapped
        # then all multimapping reads will be yielded
        else:
            E.warn('''Multimaping reads have been detected and these will
            be output to the final bam file''')
            for read in self.pysam_in2:
                if read.qname in read_list:
                    yield read


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--methods", dest="methods", type="choice",
                      action="append",
                      choices=("filter",
                               "keep-first-base",
                               "set-nh",
                               "set-sequence",
                               "strip-sequence",
                               "strip-quality",
                               "unstrip",
                               "unset-unmapped-mapq",
                               "downsample-single",
                               "downsample-paired"),
                      help="methods to apply [%default]")

    parser.add_option("--strip-method", dest="strip_method", type="choice",
                      choices=("all", "match"),
                      help="define which sequences/qualities to strip. "
                      "match means that stripping only applies to entries "
                      "without mismatches (requires NM tag to be present). "
                      "[%default]")

    parser.add_option("--filter-method", dest="filter_methods",
                      action="append", type="choice",
                      choices=('NM', 'CM', 'mapped', 'unique', "non-unique"),
                      help="filter method to apply to remove alignments "
                      "from a bam file. Multiple methods can be supplied "
                      "[%default]")

    parser.add_option("--reference-bam-file", dest="reference_bam",
                      type="string",
                      help="bam-file to filter with [%default]")

    parser.add_option("--force-output", dest="force", action="store_true",
                      help="force processing. Some methods such "
                      "as strip/unstrip will stop processing if "
                      "they think it not necessary "
                      "[%default]")

    parser.add_option("--output-sam", dest="output_sam", action="store_true",
                      help="output in sam format [%default]")

    parser.add_option("--inplace", dest="inplace", action="store_true",
                      help="modify bam files in-place. Bam files need "
                      "to be given "
                      "as arguments. Temporary bam files are written "
                      "to /tmp [%default]")

    parser.add_option(
        "--first-fastq-file", "-1", dest="fastq_pair1", type="string",
        help="fastq file with read information for first "
        "in pair or unpaired. Used for unstripping sequence "
        "and quality scores [%default]")

    parser.add_option(
        "--second-fastq-file", "-2", dest="fastq_pair2", type="string",
        help="fastq file with read information for second "
        "in pair. Used for unstripping sequence "
        "and quality scores  [%default]")

    parser.add_option(
        "--downsample", dest="downsample",
        type="int",
        help="Number of reads to downsample to")

    parser.set_defaults(
        methods=[],
        output_sam=False,
        reference_bam=None,
        filter_methods=[],
        strip_method="all",
        force=False,
        inplace=False,
        fastq_pair1=None,
        fastq_pair2=None,
        downsample=None,
        random_seed=None
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)
    # random.seed(options.random_seed)
    bamfiles = []

    if options.stdin != sys.stdin:
        from_stdin = True
        bamfiles.append(options.stdin.name)
    else:
        from_stdin = False
        
    if options.inplace:
        bamfiles.extend(args)
        if len(bamfiles) == 0:
            raise ValueError(
                "please one or more bam-files as command line arguments")

        if "-" in bamfiles:
            raise ValueError(
                "can not read from stdin if ``--inplace`` is selected")

    if len(bamfiles) == 0:
        bamfiles = ["-"]

    to_stdout = False

    for bamfile in bamfiles:

        E.info('processing %s' % bamfile)

        if os.path.islink(bamfile):
            E.warn('ignoring link %s' % bamfile)
            continue

        if IOTools.isEmpty(bamfile):
            E.warn('ignoring empty file %s' % bamfile)
            continue

        # reading bam from stdin does not work with only the "r" tag
        pysam_in = pysam.AlignmentFile(bamfile, "rb")
        if bamfile == "-" or (from_stdin and bamfile == options.stdin.name):
            to_stdout = True
            if options.output_sam:
                pysam_out = pysam.AlignmentFile("-", "wh", template=pysam_in)
            else:
                pysam_out = pysam.AlignmentFile("-", "wb", template=pysam_in)
        else:
            if IOTools.isEmpty(bamfile):
                E.warn('skipping empty file %s' % bamfile)
                continue
            tmpfile = tempfile.NamedTemporaryFile(delete=False, prefix="ctmp")
            tmpfile.close()

            E.debug("writing temporary bam-file to %s" % tmpfile.name)
            pysam_out = pysam.AlignmentFile(tmpfile.name,
                                            "wb",
                                            template=pysam_in)

        if "filter" in options.methods:

            remove_mismatches, colour_mismatches = False, False

            if "NM" in options.filter_methods:
                remove_mismatches = True

            elif "CM" in options.filter_methods:
                remove_mismatches = True
                colour_mismatches = True

            if remove_mismatches:
                if not options.reference_bam:
                    raise ValueError(
                        "requiring reference bam file for removing by "
                        "mismatches")

                pysam_ref = pysam.AlignmentFile(options.reference_bam, "rb")
            else:
                pysam_ref = None

            # filter and flags are the opposite way around
            c = _bam2bam.filter_bam(
                pysam_in, pysam_out, pysam_ref,
                remove_nonunique="unique" in options.filter_methods,
                remove_unique="non-unique" in options.filter_methods,
                remove_contigs=None,
                remove_unmapped="mapped" in options.filter_methods,
                remove_mismatches=remove_mismatches,
                colour_mismatches=colour_mismatches)

            if pysam_ref:
                pysam_ref.close()

            # do not write to stdlog in the middle of a SAM/BAM stdout stream.
            if options.stdlog != options.stdout:
                E.info("category\tcounts\n%s\n" % c.asTable())
        else:

            # set up the modifying iterators
            it = pysam_in.fetch(until_eof=True)

            # function to check if processing should start
            pre_check_f = lambda x: None

            if "unset-unmapped-mapq" in options.methods:
                def unset_unmapped_mapq(i):
                    for read in i:
                        if read.is_unmapped:
                            read.mapq = 0
                        yield read
                it = unset_unmapped_mapq(it)

            if "set-sequence" in options.methods:
                def set_sequence(i):
                    for read in i:
                        # can't get at length of unmapped reads
                        if read.is_unmapped:
                            read.seq = "A"
                            read.qual = "F"
                        else:
                            read.seq = "A" * read.inferred_length
                            read.qual = "F" * read.inferred_length

                        yield read
                it = set_sequence(it)

            if "strip-sequence" in options.methods or "strip-quality" in \
               options.methods:
                def strip_sequence(i):
                    for read in i:
                        read.seq = None
                        yield read

                def check_sequence(reads):
                    if reads[0].seq is None:
                        return 'no sequence present'
                    return None

                def strip_quality(i):
                    for read in i:
                        read.qual = None
                        yield read

                def check_quality(reads):
                    if reads[0].qual is None:
                        return 'no quality information present'
                    return None

                def strip_match(i):
                    for read in i:
                        try:
                            nm = read.opt('NM')
                        except KeyError:
                            nm = 1
                        if nm == 0:
                            read.seq = None
                        yield read

                if options.strip_method == "all":
                    if "strip-sequence" in options.methods:
                        it = strip_sequence(it)
                        pre_check_f = check_sequence
                    elif "strip-quality" in options.methods:
                        it = strip_quality(it)
                        pre_check_f = check_quality
                elif options.strip_method == "match":
                    it = strip_match(it)

            if "unstrip" in options.methods:
                def buildReadDictionary(filename):
                    if not os.path.exists(filename):
                        raise OSError("file not found: %s" % filename)
                    fastqfile = pysam.FastxFile(filename)
                    fastq2sequence = {}
                    for x in fastqfile:
                        if x.name in fastq2sequence:
                            raise ValueError(
                                "read %s duplicate - can not unstrip" % x.name)

                        fastq2sequence[x.name] = (x.sequence, x.quality)
                    return fastq2sequence

                if not options.fastq_pair1:
                    raise ValueError(
                        "please supply fastq file(s) for unstripping")
                fastq2sequence1 = buildReadDictionary(options.fastq_pair1)
                if options.fastq_pair2:
                    fastq2sequence2 = buildReadDictionary(options.fastq_pair2)

                def unstrip_unpaired(i):
                    for read in i:
                        read.seq, read.qual = fastq2sequence1[read.qname]
                        yield read

                def unstrip_pair(i):
                    for read in i:
                        if read.is_read1:
                            read.seq, read.qual = fastq2sequence1[read.qname]
                        else:
                            read.seq, read.qual = fastq2sequence2[read.qname]
                        yield read

                if options.fastq_pair2:
                    it = unstrip_pair(it)
                else:
                    it = unstrip_unpaired(it)

            if "set-nh" in options.methods:
                it = SetNH(it)

            # keep first base of reads by changing the cigarstring to
            # '1M' and, in reads mapping to the reverse strand,
            # changes the pos to aend - 1
            # Needs to be refactored to make it more general
            # (last base, midpoint, ..)
            if "keep_first_base" in options.methods:
                def keep_first_base(i):
                    for read in i:
                        if read.is_reverse:
                            read.pos = read.aend - 1
                            read.cigarstring = '1M'
                        elif not read.is_unmapped:
                            read.cigarstring = '1M'
                        yield read
                it = keep_first_base(it)

            # read first read and check if processing should continue
            # only possible when not working from stdin
            # Refactoring: use cache to also do a pre-check for
            # stdin input.
            if bamfile != "-":
                # get first read for checking pre-conditions
                first_reads = list(pysam_in.head(1))

                msg = pre_check_f(first_reads)
                if msg is not None:
                    if options.force:
                        E.warn('proccessing continues, though: %s' % msg)
                    else:
                        E.warn('processing not started: %s' % msg)
                        pysam_in.close()
                        pysam_out.close()
                        continue

            if "downsample-single" in options.methods:

                if not options.downsample:
                    raise ValueError("Please provide downsample size")

                else:
                    down = SubsetBam(pysam_in=it,
                                     downsample=options.downsample,
                                     paired_end=None,
                                     single_end=True,
                                     random_seed=options.random_seed)
                    it = down.downsample_single()

            if "downsample-paired" in options.methods:

                if not options.downsample:
                    raise ValueError("Please provide downsample size")

                else:
                    down = SubsetBam(pysam_in=it,
                                     downsample=options.downsample,
                                     paired_end=True,
                                     single_end=None,
                                     random_seed=options.random_seed)
                    it = down.downsample_paired()

            # continue processing till end
            for read in it:
                pysam_out.write(read)

        pysam_in.close()
        pysam_out.close()

        if options.inplace:
            # set date and file permissions according to original
            # Note: currently it will not update user and group.
            original = os.stat(bamfile)
            os.utime(tmpfile.name, (original.st_atime, original.st_mtime))
            os.chmod(tmpfile.name, original.st_mode)
            # move new file over original copy
            shutil.move(tmpfile.name, bamfile)
            # re-index
            pysam.index(bamfile)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

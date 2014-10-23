'''
PipelinePreprocess.py - Utility functions for processing short reads
==============================================================

:Author: Tom smith
:Release: $Id$
:Date: |today|
:Tags: Python

UPDATE: Mapping reads is a common task in pipelines. Different pipelines
combine different sources of input (:term:`fastq` files, :term:`sra` files)
of different data (single end, paired end) with different mapping
algorithms (bowtie, tophat, stampy). This module provides utility
functions to abstract some of these variations.

The pipeline does not know what kind of data it gets (a :term:`sra` archive
might contain single end or paired end data or both).

A pipeline might get several input data (:term:`fastq` and :term:`sra`
formatted files at the same time).

The module currently is able to deal with:

   * tophat mapping against genome
   * bowtie mapping against transcriptome, genome and junctions
   * bwa against genome
   * stampy against genome

It implements:
   * .sra: paired-end and single-end
   * .fastq: paired-end and single-end
   * .csfasta: colour-space, single-end

Code
----

'''

import os
import shutil
import glob
import collections
import re
import gzip
import itertools
import CGAT.Pipeline as P
import logging as L
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import CGAT.Fastq as Fastq
import CGAT.IndexedFasta as IndexedFasta
import CGATPipelines.PipelineGeneset as PipelineGeneset
import pysam

SequenceInformation = collections.namedtuple("SequenceInformation",
                                             """paired_end
                                                 filename_first
                                                 filename_second
                                                 readlength_first
                                                 readlength_second
                                                 is_colour""")


def getReadLengthFromFastq(filename):
    '''return readlength from a fasta/fastq file.

    Only the first read is inspected. If there are
    different read lengths in the file, though luck.

    '''

    with IOTools.openFile(filename) as infile:
        record = iterate(infile).next()
        readlength = len(record.seq)
        return readlength


def getReadLengthFromBamfile(filename):
    '''return readlength from a bam file.

    Only the first read is inspected. If there are
    different read lengths in the file, though luck.
    '''

    samfile = pysam.Samfile(filename, "rb")
    record = samfile.fetch().next()
    readlength = record.rlen
    samfile.close()
    return readlength


def getSequencingInformation(track):
    '''glean sequencing information from *track*.'''

    colour = False
    if os.path.exists("%s.fastq.gz" % track):
        first_pair = "%s.fastq.gz" % track
        second_pair = None
    elif os.path.exists("%s.fastq.1.gz" % track):
        first_pair = "%s.fastq.1.gz" % track
        second_pair = "%s.fastq.2.gz" % track
    elif os.path.exists("%s.csfasta.gz" % track):
        first_pair = "%s.csfasta.gz" % track
        second_pair = None
        colour = True

    second_length = None
    if second_pair:
        if not os.path.exists(second_pair):
            raise IOError("could not find second pair %s for %s" %
                          (second_pair, first_pair))
        second_length = getReadLengthFromFastq(second_pair)

    return SequenceInformation._make((second_pair is not None,
                                      first_pair, second_pair,
                                      getReadLengthFromFastq(first_pair),
                                      second_length,
                                      colour))


class Preprocessor(object):

    '''preprocess reads.

    converts the input data, calls processing tools

    All in a single statement to be send to the cluster.
    '''

    datatype = "fastq"

    # compress temporary fastq files with gzip
    compress = False

    # convert to sanger quality scores
    convert = False

    def __init__(self, executable=None,
                 *args, **kwargs):
        if executable:
            self.executable = executable

    def quoteFile(self, filename):
        '''add uncompression for compressed files.
        and programs that expect uncompressed files.

        .. note::
            This will only work if the downstream programs read the
            file only once.
        '''
        if filename.endswith(".gz") and not self.compress:
            return "<( gunzip < %s )" % filename
        else:
            return filename

    def convert(self, infiles, outfile):
        '''build conversion statement

        Build a command line statement that extracts/converts
        various input formats to fastq formatted files.

        returns the statement and the fastq files.
        '''

        assert len(infiles) > 0, "no input files for preprocessing"

        tmpdir_fastq = P.getTempDir()

        # create temporary directory again for nodes
        statement = ["mkdir -p %s" % tmpdir_fastq]
        fastqfiles = []

        # get track by extension of outfile
        track = os.path.splitext(os.path.basename(outfile))[0]

        if self.compress:
            compress_cmd = "| gzip"
            extension = ".gz"
        else:
            compress_cmd = ""
            extension = ""

        for infile in infiles:

            if infile.endswith(".export.txt.gz"):
                # single end illumina export
                statement.append("""gunzip < %(infile)s
                | awk '$11 != "QC" || $10 ~ /(\d+):(\d+):(\d+)/ \
                {if ($1 != "")
                {readname=sprintf("%%%%s_%%%%s:%%%%s:%%%%s:%%%%s:%%%%s",
                $1, $2, $3, $4, $5, $6);}
                else {readname=sprintf("%%%%s:%%%%s:%%%%s:%%%%s:%%%%s",
                $1, $3, $4, $5, $6);}
                printf("@%%%%s\\n%%%%s\\n+\\n%%%%s\\n",readname,$9,$10);}'
                %(compress_cmd)s
                > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" % locals())
                fastqfiles.append(
                    ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),))
            elif infile.endswith(".fa.gz"):
                statement.append(
                    '''gunzip < %(infile)s %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fa''' % locals())
                fastqfiles.append(("%s/%s.fa" % (tmpdir_fastq, track),))
                self.datatype = "fasta"

            elif infile.endswith(".sra"):
                # sneak preview to determine if paired end or single end
                outdir = P.getTempDir()
                # --split-files is present in fastq-dump 2.1.7
                P.execute(
                    """fastq-dump --split-files --gzip -X 1000
                    --outdir %(outdir)s %(infile)s""" % locals())
                # --split-files will create files called prefix_#.fastq.gz
                # where # is the read number.
                # The following cases are:

                # * file cotains paired end data:
                #      output = prefix_1.fastq.gz, prefix_2.fastq.gz
                #
                #    * special case: unpaired reads in a paired end
                #                    run end up in prefix.fastq.gz
                #    * special case: if paired reads are stored in
                #                    a single read, fastq-dump will split.
                #       There might be a joining sequence. The output
                #                 would thus be:
                #       prefix_1.fastq.gz, prefix_2.fastq.gz, prefix_3.fastq.gz
                #
                #      You want files 1 and 3.
                f = sorted(glob.glob(os.path.join(outdir, "*.fastq.gz")))
                ff = [os.path.basename(x) for x in f]
                if len(f) == 1:
                    # sra file contains one read: output = prefix.fastq.gz
                    pass
                elif len(f) == 2:
                    # sra file contains read pairs: output = prefix_1.fastq.gz,
                    # prefix_2.fastq.gz
                    assert ff[0].endswith(
                        "_1.fastq.gz") and ff[1].endswith("_2.fastq.gz")
                elif len(f) == 3:
                    if ff[2].endswith("_3.fastq.gz"):
                        f = glob.glob(os.path.join(outdir, "*_[13].fastq.gz"))
                    else:
                        f = glob.glob(os.path.join(outdir, "*_[13].fastq.gz"))
                E.info("sra file contains the following files: %s" % f)
                shutil.rmtree(outdir)
                fastqfiles.append(
                    ["%s/%s" % (tmpdir_fastq, os.path.basename(x))
                     for x in sorted(f)])
                statement.append(
                    """fastq-dump --split-files --gzip --outdir
                    %(tmpdir_fastq)s %(infile)s""" % locals())

            elif infile.endswith(".fastq.gz"):
                format = Fastq.guessFormat(
                    IOTools.openFile(infile, "r"), raises=False)
                if 'sanger' not in format and self.convert:
                    statement.append("""gunzip < %(infile)s
                    | python %%(scriptsdir)s/fastq2fastq.py
                    --change-format=sanger
                    --guess-format=phred64
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" %
                                     locals())
                    fastqfiles.append(
                        "%s/%s.fastq%s" % (tmpdir_fastq, track, extension))
                else:
                    E.debug("%s: assuming quality score format %s" %
                            (infile, format))
                    fastqfiles.append(infile)

            elif infile.endswith(".csfasta.gz"):
                # single end SOLiD data
                if self.preserve_colourspace:
                    quality = P.snip(infile, ".csfasta.gz") + ".qual.gz"
                    if not os.path.exists(quality):
                        raise ValueError("no quality file for %s" % infile)
                    statement.append("""gunzip < %(infile)s
                    > %(tmpdir_fastq)s/%(track)s.csfasta%(extension)s""" %
                                     locals())
                    statement.append("""gunzip < %(quality)s
                    > %(tmpdir_fastq)s/%(track)s.qual%(extension)s""" %
                                     locals())
                    fastqfiles.append("%s/%s.csfasta%s" %
                                      (tmpdir_fastq, track, extension),
                                      "%s/%s.qual%s" %
                                      (tmpdir_fastq, track, extension))
                    self.datatype = "solid"
                else:
                    quality = P.snip(infile, ".csfasta.gz") + ".qual.gz"

                    statement.append("""solid2fastq
                    <(gunzip < %(infile)s) <(gunzip < %(quality)s)
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fastq%(extension)""" %
                                     locals())
                    fastqfiles.append(
                        "%s/%s.fastq%s" % (tmpdir_fastq, track, extension))

            elif infile.endswith(".csfasta.F3.gz"):
                # paired end SOLiD data
                if self.preserve_colourspace:
                    bn = P.snip(infile, ".csfasta.F3.gz")
                    # order is important - mirrors tophat reads followed by
                    # quals
                    f = []
                    for suffix in ("csfasta.F3", "csfasta.F5",
                                   "qual.F3", "qual.F5"):
                        fn = "%(bn)s.%(suffix)s" % locals()
                        if not os.path.exists(fn + ".gz"):
                            raise ValueError(
                                "expected file %s.gz missing" % fn)
                        statement.append("""gunzip < %(fn)s.gz
                        %(compress_cmd)s
                        > %(tmpdir_fastq)s/%(track)s.%(suffix)s%(extension)s
                        """ %
                                         locals())
                        f.append(
                            "%(tmpdir_fastq)s/%(track)s.%(suffix)s%(extension)s" %
                            locals())
                    fastqfiles.append(f)
                    self.datatype = "solid"
                else:
                    quality = P.snip(infile, ".csfasta.gz") + ".qual.gz"

                    statement.append("""solid2fastq
                    <(gunzip < %(infile)s) <(gunzip < %(quality)s)
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" %
                                     locals())
                    fastqfiles.append(
                        "%s/%s.fastq%s" % (tmpdir_fastq, track, extension))

            elif infile.endswith(".fastq.1.gz"):

                bn = P.snip(infile, ".fastq.1.gz")
                infile2 = "%s.fastq.2.gz" % bn
                if not os.path.exists(infile2):
                    raise ValueError(
                        "can not find paired ended file "
                        "'%s' for '%s'" % (infile2, infile))

                format = Fastq.guessFormat(
                    IOTools.openFile(infile), raises=False)
                if 'sanger' not in format:
                    statement.append("""gunzip < %(infile)s
                    | python %%(scriptsdir)s/fastq2fastq.py
                    --change-format=sanger
                    --guess-format=phred64
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.1.fastq%(extension)s;
                    gunzip < %(infile2)s
                    | python %%(scriptsdir)s/fastq2fastq.py
                    --change-format=sanger
                    --guess-format=phred64
                    --log=%(outfile)s.log
                    %(compress_cmd)s
                    > %(tmpdir_fastq)s/%(track)s.2.fastq%(extension)s
                    """ % locals())
                    fastqfiles.extend(
                        "%s/%s.1.fastq%s" % (tmpdir_fastq, track, extension),
                         "%s/%s.2.fastq%s" % (tmpdir_fastq, track, extension))

                else:
                    E.debug("%s: assuming quality score format %s" %
                            (infile, format))
                    fastqfiles.extend((infile,
                                       infile2))

            else:
                raise NotImplementedError("unknown file format %s" % infile)

        self.tmpdir_fastq = tmpdir_fastq

        assert len(fastqfiles) > 0, "no fastq files for preprocessing"

        return "; ".join(statement) + ";", fastqfiles

    def process(self, infile, list_of_preprocessers, outfile,
                compress=False, save=True):
        '''build prcessing statement on infiles.
        list of processing tools is used to build a command line statement
        in list order
        '''
        total_cmd = ""
        processed_files = []
        print list_of_preprocessers
        file_compress = compress
        save_intermediates = save
        for idx, tool in enumerate(list_of_preprocessers):
            print tool
            print idx
            if idx == len(list_of_preprocessers)-1:
                end = True
            else:
                end = False

            if tool == "cutadapt":
                tool_cmd = ""
            elif tool == "trimmomatic":
                tool_cmd = ""
            elif tool == "trimgalore":
                processer_object = trimgalore(
                    compress=file_compress,
                    save=save_intermediates,
                    final=end)
                tool_cmd, infile = processer_object.build(infile)
            else:
                tool_cmd = "testing_testing"
                # insert error statement here.
                pass
            
            total_cmd += tool_cmd
            processed_files.append(infile,)
        return total_cmd, processed_files

    def postprocess(self, initial_file, preprocesser_list):
        '''collect final output data and postprocess.
        summarise initial fastq and then all processed fastqs'''

        if len(initial_file) == 1:
            num_files = 1
            initial_file = os.path.basename(initial_file[0])
        elif len(initial_file) == 2:
            num_files = 2
            initial_file1, initial_file2 = initial_file
            
        if num_files == 1:
            postprocess_cmd = '''zcat %(initial_file)s |
                          python %%(scriptsdir)s/fastq2summary.py
                          --guess-format=illumina-1.8
                          > processed.dir/%(initial_file)s.summary;
                          ''' % locals()

            processed_file = initial_file
            for tool in preprocesser_list:
                processed_file = tool + "-" + processed_file
                postprocess_cmd += '''zcat processed.dir/%(processed_file)s |
                              python %%(scriptsdir)s/fastq2summary.py
                              --guess-format=illumina-1.8
                              > processed.dir/%(processed_file)s.summary;
                              ''' % locals()

        # edit from here to allow processing of paired end sequences

        elif num_files == 2:
            print "HAVE THEY SPLIT?............."
            print initial_file
            print initial_file1
            print initial_file2
            initial_file1, initial_file2 = initial_file
            postprocess_cmd = '''zcat %(initial_file1)s |
                          python %%(scriptsdir)s/fastq2summary.py
                          --guess-format=illumina-1.8
                          > processed.dir/%(initial_file1)s.summary;
                          zcat %(initial_file2)s |
                          python %%(scriptsdir)s/fastq2summary.py
                          --guess-format=illumina-1.8
                          > processed.dir/%(initial_file2)s.summary;
                          ''' % locals()

            processed_file1 = initial_file1
            processed_file2 = initial_file2
            for tool in preprocesser_list:
                print "THIS IS THE PROBLEM..................."
                print tool
                print processed_file1
                processed_file1 = tool + "-" + processed_file1
                processed_file2 = tool + "-" + processed_file2
                postprocess_cmd += '''zcat processed.dir/%(processed_file1)s |
                              python %%(scriptsdir)s/fastq2summary.py
                              --guess-format=illumina-1.8
                              > processed.dir/%(processed_file1)s.summary;
                              zcat processed.dir/%(processed_file2)s |
                              python %%(scriptsdir)s/fastq2summary.py
                              --guess-format=illumina-1.8
                              > processed.dir/%(processed_file2)s.summary;
                              ''' % locals()

        return postprocess_cmd

    def cleanup(self, outfile):
        '''clean up.'''
        statement = '''rm -rf %s;''' % (self.tmpdir_fastq)
        return statement

    def build(self, infile, outfile, preprocesser_list):
        '''run mapper.'''
        # how best to pass on parameters?
        cmd_preprocess, raw_file = self.convert(infile, outfile)
        cmd_process, processed_files = self.process(raw_file,
                                                    preprocesser_list, outfile)
        cmd_postprocess = self.postprocess(raw_file, preprocesser_list)
        cmd_clean = self.cleanup(outfile)

        assert cmd_preprocess.strip().endswith(";")
        assert cmd_process.strip().endswith(";")
        if cmd_postprocess:
            assert cmd_postprocess.strip().endswith(";")
        if cmd_clean:
            assert cmd_clean.strip().endswith(";")

        statement = " checkpoint; ".join((cmd_preprocess,
                                          cmd_process,
                                          cmd_postprocess,
                                          cmd_clean))

        return statement


class processer(object):
    '''define class attributes for processers'''
    def __init__(self, final=True, compress=False,
                 save=True, *args, **kwargs):
        self.final = final
        self.compress = compress
        if self.final:
            self.save = True
        else:
            self.save = save

    def getfastqAttr(self, infiles):
        num_files = len(infiles)

        if num_files > 1:
            infile1, infile2 = infiles
        else:
            infile1 = infiles[0]

        format = Fastq.guessFormat(IOTools.openFile(infile1), raises=False)
        E.info("%s: format guess: %s" % (infile1, format))
        offset = Fastq.getOffset(format, raises=False)

        return num_files, offset

    def build(self):
        pass


class trimgalore(processer):

    def process(self, infiles):

        num_files, offset = self.getfastqAttr(infiles)

        if self.save:
            outdir = "processed.dir"
        else:
            outdir = P.getTempDir("/ifs/scratch")

        if num_files == 1:
            infile = infiles[0]
            track = os.path.basename(infile)
            outfile = "%(outdir)s/trimgalore-%(track)s" % locals()
            logfile = outdir + "/" + track + ".trim_galore.log"
            trim_out = "%s/%s_trimmed.fq.gz" % (outdir,
                                                P.snip(track, ".fastq.gz"))
            cmd = '''trim_galore %%(trimgalore_options)s
                     -a %%(trimgalore_adapter)s
                     --phred%(offset)s
                     --output_dir %(outdir)s
                     %(infile)s
                     2>>%(logfile)s;
                     mv %(trim_out)s %(outfile)s;
                     ''' % locals()
            outfiles = outfile
            logfile = "%(track)s_trimming_report.txt"

        elif num_files == 2:
            infile1, infile2 = infiles
            track1 = os.path.basename(infile1)
            track2 = P.snip(track1, ".fastq.1.gz") + ".fastq.2.gz"
            outfile1 = "%(outdir)s/trimgalore-%(track1)s" % locals()
            outfile2 = "%(outdir)s/trimgalore-%(track2)s" % locals()

            logfile = outdir + "/" + track1 + ".trim_galore.log"

            cmd = '''trim_galore %%(trimgalore_options)s
                     --paired -a %%(trimgalore_adapter)s
                     --phred%(offset)s --output_dir %(outdir)s
                     %(infile1)s %(infile2)s
                     2>>%(logfile)s;
                     mv %(outdir)s/%(track1)s_val_1.fq.gz %(outfile1)s;
                     mv %(outdir)s/%(track2)s_val_2.fq.gz %(outfile2)s;
                     ''' % locals()
            outfiles = (outfile1, outfile2)

        else:
            cmd = "this shouldn't happen! - throw an error!"

        return cmd, outfiles

    def postprocess(self, logfile):
        return ""

    def build(self, infile):
        process_cmd, fastq_outfiles = self.process(infile)
        post_cmd = self.postprocess(fastq_outfiles)
        statement = " checkpoint; ".join([process_cmd,
                                          post_cmd])

        return statement, fastq_outfiles


    

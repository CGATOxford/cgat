'''
PipelinePreprocess.py - Utility functions for processing short reads
==============================================================

:Author: Tom smith
:Release: $Id$
:Date: |today|
:Tags: Python

UPDATE: Processing reads is a common task before mapping. Different sequencing
technologies and applications require different read processing.
This module provides utility functions to abstract some of these variations

The pipeline does not know what kind of data it gets (a directory
might contain single end or paired end data or both).

The module currently provides modules to perform:
    * hard-trimming (fastx_trimmer, trimmomatic)
    * adapter trimming (trimgalore, trimmomatic)
    * read end quality trimming (trimgalore, trimmomatic)
    * sliding window quality trimming (sickle, trimmomatic)
    * RRBS-specific trimming (trimgalore)

It has been tested with:
   * .fastq: paired-end and single-end


To do
=====
Check this will work with .sra files
How to deal with flash?
e.g this is the only tool that doesn't do one-to-one processing


Code
----

'''

import re
import os
import collections
import CGAT.Pipeline as P
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Fastq as Fastq
import CGATPipelines.PipelineMapping as Mapping

SequenceInformation = collections.namedtuple("SequenceInformation",
                                             """paired_end
                                                 filename_first
                                                 filename_second
                                                 readlength_first
                                                 readlength_second
                                                 is_colour""")


def makeSecond(file1):
    return(re.sub(".fastq.1.gz", ".fastq.2.gz", file1))


def getReadLengthFromFastq(filename):
    '''return readlength from a fasta/fastq file.

    Only the first read is inspected. If there are
    different read lengths in the file, though luck.

    '''

    with IOTools.openFile(filename) as infile:
        record = iterate(infile).next()
        readlength = len(record.seq)
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


class MasterProcessor(Mapping.Mapper):

    '''Processes reads with tools specified.
    inherits from PipelineMapping Mapper class as
    converts the input data, calls processing tools

    All in a single statement to be send to the cluster.
    '''

    datatype = "fastq"

    # compress temporary fastq files with gzip
    compress = False

    # convert to sanger quality scores
    convert = False

    def __init__(self, save=True, summarise=False,
                 threads=1, scriptsdir=None,
                 trimgalore_options=None,
                 trimmomatic_options=None,
                 sickle_options=None,
                 flash_options=None,
                 fastx_trimmer_options=None,
                 *args, **kwargs):
        self.save = save
        self.summarise = summarise
        self.threads = threads
        self.scriptsdir = scriptsdir
        self.trimgalore_opt = trimgalore_options
        self.trimmomatic_opt = trimmomatic_options
        self.sickle_opt = sickle_options
        self.flash_opt = flash_options
        self.fastx_trimmer_opt = fastx_trimmer_options

        if self.save:
            self.outdir = "processed.dir"
        else:
            self.outdir = P.getTempDir("/ifs/scratch")

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

    def process(self, infile, list_of_preprocessers, outfile,
                compress=False, save=True):
        '''build processing statement on infiles.  list of processing
        tools is used to build a command line statement in list order'''

        total_tool_cmd = ""
        total_post_cmd = ""
        processed_files = []
        file_compress = compress
        save_intermediates = save
        first = 1
        f_format = ""
        num_files = ""
        # initiate a dummy processer object with no prefix to summarise
        # the initial input files
        initial_processer_object = process_tool(
            compress=file_compress, save=save_intermediates,
            final=0, f_format=f_format, num_files=num_files,
            first=first, outdir=self.outdir, infiles=infile,
            prefix="")
        tool_cmd, post_cmd, null_infile, null_f_format, null_num_files = (
            initial_processer_object.build(infile))
        total_post_cmd += post_cmd

        for idx, tool in enumerate(list_of_preprocessers):
            if idx == len(list_of_preprocessers)-1:
                end = True
                save_intermediates = True  # if final file must save!
            else:
                end = False
            if idx > 0:
                first = 0
            if tool == "fastx_trimmer":
                fastx_trimmer_object = fastx_trimmer(
                    compress=file_compress, save=save_intermediates,
                    final=end, f_format=f_format, num_files=num_files,
                    first=first, outdir=self.outdir, infiles=infile,
                    prefix="fastx_trimmer-", summarise=self.summarise,
                    options=self.fastx_trimmer_opt, threads=self.threads,
                    scriptsdir=self.scriptsdir)
                tool_cmd, post_cmd, infile, f_format, num_files = (
                    fastx_trimmer_object.build(infile))
            elif tool == "trimmomatic":
                trimmomatic_object = trimmomatic(
                    compress=file_compress, save=save_intermediates,
                    final=end, f_format=f_format, num_files=num_files,
                    first=first, outdir=self.outdir, infiles=infile,
                    prefix="trimmomatic-", summarise=self.summarise,
                    options=self.trimmomatic_opt, threads=self.threads,
                    scriptsdir=self.scriptsdir)
                tool_cmd, post_cmd, infile, f_format, num_files = (
                    trimmomatic_object.build(infile))
            elif tool == "sickle":
                sickle_object = sickle(
                    compress=file_compress, save=save_intermediates,
                    final=end, f_format=f_format, num_files=num_files,
                    first=first, outdir=self.outdir, infiles=infile,
                    prefix="sickle-", summarise=self.summarise,
                    options=self.sickle_opt, threads=self.threads,
                    scriptsdir=self.scriptsdir)
                tool_cmd, post_cmd, infile,  f_format, num_files = (
                    sickle_object.build(infile))
            elif tool == "trimgalore":
                trimgalore_object = trimgalore(
                    compress=file_compress, save=save_intermediates,
                    final=end, f_format=f_format, num_files=num_files,
                    first=first, outdir=self.outdir, infiles=infile,
                    prefix="trimgalore-", summarise=self.summarise,
                    options=self.trimgalore_opt, threads=self.threads,
                    scriptsdir=self.scriptsdir)
                tool_cmd, post_cmd, infile, f_format, num_files = (
                    trimgalore_object.build(infile))
            elif tool == "flash":
                flash_object = flash(
                    compress=file_compress, save=save_intermediates,
                    final=end, f_format=f_format, num_files=num_files,
                    first=first, outdir=self.outdir, infiles=infile,
                    prefix="flash-", summarise=self.summarise,
                    options=self.flash_opt, threads=self.threads,
                    scriptsdir=self.scriptsdir)
                tool_cmd, post_cmd, infile, f_format, num_files = (
                    flash_object.build(infile))
            else:
                E.info("%s is not a supported tool" % tool)

            total_tool_cmd += tool_cmd
            total_post_cmd += post_cmd
            processed_files.append(infile,)
        return total_tool_cmd, total_post_cmd, processed_files

    def cleanup(self, outfile):
        '''clean up.'''
        if self.save:
            statement = 'checkpoint;'
        else:
            statement = 'rm -rf %s;' % self.outdir
        return statement

    def build(self, infile, outfile, processer_list):
        '''run mapper.'''
        cmd_preprocess, raw_file = self.preprocess(infile, outfile)
        raw_file = raw_file[0]
        cmd_process, cmd_post, processed_files = self.process(
            raw_file, processer_list, outfile, save=self.save)
        cmd_clean = self.cleanup(outfile)

        assert cmd_preprocess.strip().endswith(";")
        assert cmd_process.strip().endswith(";")
        assert cmd_post.strip().endswith(";")
        assert cmd_clean.strip().endswith(";")

        statement = " checkpoint; ".join((cmd_preprocess,
                                          cmd_process,
                                          cmd_post,
                                          cmd_clean))
        return statement


class process_tool(object):
    '''defines class attributes for a sequennce utility tool'''

    def __init__(self, first=True, final=True, compress=False,
                 save=True, f_format="", num_files="", prefix="",
                 outdir="", infiles=(), summarise=False, options=None,
                 threads=1, scriptsdir=None, *args, **kwargs):
        self.final = final
        self.compress = compress
        self.first = first
        self.f_format = f_format
        self.num_files = num_files
        self.prefix = prefix
        self.summarise = summarise
        self.processing_options = options
        self.threads = threads
        self.scriptsdir = scriptsdir
        if self.final:
            self.save = True
            self.outdir = "processed.dir"
        else:
            self.save = save
            self.outdir = outdir

    def setfastqAttr(self, infiles):

        if self.first:
            self.num_files = len(infiles)

        if self.num_files == 1:
            infile1 = infiles[0]
        elif self.num_files == 2:
            infile1, infile2 = infiles

        if self.first:
            self.f_format = Fastq.guessFormat(IOTools.openFile(infile1),
                                              raises=False)
            E.info("%s: format guess: %s" % (infile1, self.f_format))
        else:
            E.info("%s: format set as previous: %s" % (
                infile1, self.f_format))
        self.offset = Fastq.getOffset(self.f_format, raises=False)

    def process(self, infiles):
        return ("", infiles)

    def postprocess(self, infiles):

        outdir = self.outdir
        prefix = self.prefix
        scriptsdir = self.scriptsdir
        if self.summarise:
            if self.num_files == 1:
                infile = infiles[0]
                infile_base = os.path.basename(infile)
                postprocess_cmd = '''zcat %(infile)s |
                python %(scriptsdir)s/fastq2summary.py
                --guess-format=illumina-1.8 -v0
                > summary.dir/%(infile_base)s.summary;
                ''' % locals()

            elif self.num_files == 2:
                infile1, infile2 = infiles
                infile_base1, infile_base2, = [
                    os.path.basename(x) for x in infiles]
                postprocess_cmd = '''zcat %(infile1)s |
                python %(scriptsdir)s/fastq2summary.py
                --guess-format=illumina-1.8 -v0
                > summary.dir/%(infile_base1)s.summary;
                zcat %(infile2)s |
                python %(scriptsdir)s/fastq2summary.py
                --guess-format=illumina-1.8 -v0
                > summary.dir/%(infile_base2)s.summary
                ;''' % locals()
        else:
            postprocess_cmd = "checkpoint ;"
        return postprocess_cmd

    def build(self, infiles):
        self.setfastqAttr(infiles)

        process_cmd, fastq_outfiles = self.process(infiles)
        post_cmd = self.postprocess(fastq_outfiles)

        return (process_cmd, post_cmd, fastq_outfiles,
                self.f_format, self.num_files)


class trimgalore(process_tool):

    def process(self, infiles):
        prefix = self.prefix
        offset = self.offset
        outdir = self.outdir
        processing_options = self.processing_options
        # the assigment of infiles is repeated in each process_tool
        # refactor!
        if self.num_files == 1:
            infile = infiles[0]
            track = os.path.basename(infile)
            outfile = "%(outdir)s/%(prefix)s%(track)s" % locals()
            logfile = "log.dir/" + track + ".trim_galore.log"
            trim_out = "%s/%s_trimmed.fq.gz" % (outdir,
                                                P.snip(track, ".fastq.gz"))
            cmd = '''trim_galore %(processing_options)s
                     --phred%(offset)s
                     --output_dir %(outdir)s
                     %(infile)s
                     2>>%(logfile)s;
                     mv %(trim_out)s %(outfile)s;
                     ''' % locals()
            outfiles = (outfile,)

        elif self.num_files == 2:
            infile1, infile2 = infiles
            track1 = os.path.basename(infile1)
            track2 = P.snip(track1, ".fastq.1.gz") + ".fastq.2.gz"
            outfile1 = "%(outdir)s/%(prefix)s%(track1)s" % locals()
            outfile2 = "%(outdir)s/%(prefix)s%(track2)s" % locals()
            logfile = "log.dir/" + track1 + ".trim_galore.log"

            cmd = '''trim_galore %(processing_options)s
                     --paired
                     --phred%(offset)s --output_dir %(outdir)s
                     %(infile1)s %(infile2)s
                     2>>%(logfile)s;
                     mv %(outdir)s/%(track1)s_val_1.fq.gz %(outfile1)s;
                     mv %(outdir)s/%(track2)s_val_2.fq.gz %(outfile2)s;
                     ''' % locals()
            outfiles = (outfile1, outfile2)

        return cmd, outfiles


class sickle(process_tool):

    def process(self, infiles):
        prefix = self.prefix
        offset = self.offset
        outdir = self.outdir
        processing_options = self.processing_options
        rRANGES = {33: 'sanger', 64: 'illumina-1.8', 59: 'solexa'}
        quality = rRANGES[offset]

        if self.num_files == 1:
            infile = infiles[0]
            track = os.path.basename(infile)
            outfile = "%(outdir)s/%(prefix)s%(track)s" % locals()
            logfile = "log.dir/" + track + ".sickle.log"

            cmd = '''sickle se -g %(processing_options)s
                     --qual-type %(quality)s
                     --output-file %(outfile)s
                     --fastq-file %(infile)s
                     2>>%(logfile)s
                     ;''' % locals()
            outfiles = (outfile,)

        elif self.num_files == 2:
            infile1, infile2 = infiles
            track1 = os.path.basename(infile1)
            track2 = P.snip(track1, ".fastq.1.gz") + ".fastq.2.gz"
            outfile1 = "%(outdir)s/%(prefix)s%(track1)s" % locals()
            outfile2 = "%(outdir)s/%(prefix)s%(track2)s" % locals()
            logfile = "log.dir/" + track1 + ".sickle.log"

            cmd = '''sickle pe -g -s %(processing_options)s
                     --qual-type %(quality)s
                     -f %(infile1)s -r %(infile2)s
                     -o %(outfile1)s -p %(outfile2)s
                     2>>%(logfile)s
                     ;''' % locals()
            outfiles = (outfile1, outfile2)

        return cmd, outfiles


class trimmomatic(process_tool):

    def process(self, infiles):
        prefix = self.prefix
        offset = self.offset
        outdir = self.outdir
        threads = self.threads
        processing_options = self.processing_options
        if self.num_files == 1:
            infile = infiles[0]
            track = os.path.basename(infile)
            logfile = "log.dir/" + track + ".trimmomatic.log"
            outfile = "%(outdir)s/%(prefix)s%(track)s" % locals()
            trim_out = "%s/%s_trimmed.fq.gz" % (
                outdir, P.snip(track, ".fastq.gz"))
            cmd = '''trimmomatic SE -threads %(threads)s -phred%(offset)s
                     %(infile)s %(outfile)s
                     %(processing_options)s 2>> %(logfile)s
                     ;''' % locals()
            outfiles = (outfile,)

        elif self.num_files == 2:
            infile1, infile2 = infiles
            track1 = os.path.basename(infile1)
            track2 = re.sub(".fastq.1.gz", ".fastq.2.gz", track1)
            logfile = "log.dir/" + track1 + ".trimmomatic.log"
            outfile1 = "%(outdir)s/%(prefix)s%(track1)s" % locals()
            outfile2 = "%(outdir)s/%(prefix)s%(track2)s" % locals()

            cmd = '''trimmomatic PE -threads %(threads)s -phred%(offset)s
                     %(infile1)s %(infile2)s
                     %(outfile1)s %(outfile1)s.unpaired
                     %(outfile2)s %(outfile2)s.unpaired
                     %(processing_options)s 2>> %(logfile)s
                     ;''' % locals()
            outfiles = (outfile1, outfile2)

        return cmd, outfiles


class fastx_trimmer(process_tool):

    def process(self, infiles):
        prefix = self.prefix
        offset = self.offset
        outdir = self.outdir
        processing_options = self.processing_options

        if self.num_files == 1:
            infile = infiles[0]
            track = os.path.basename(infile)
            logfile = "log.dir/" + track + ".trimmomatic.log"
            outfile = "%(outdir)s/%(prefix)s%(track)s" % locals()
            trim_out = "%s/%s_trimmed.fq.gz" % (outdir,
                                                P.snip(track, ".fastq.gz"))
            cmd = '''zcat %(infile)s | fastx_trimmer -Q%(offset)s
            %(processing_options)s 2> %(logfile)s | gzip > %(outfile)s
            ;''' % locals()

            outfiles = (outfile,)

        else:
            E.info("fastx_trimmer does not support paired end reads")

        return cmd, outfiles


class flash(process_tool):

    def process(self, infiles):
        prefix = self.prefix
        offset = self.offset
        outdir = self.outdir
        processing_options = self.processing_options

        if self.num_files == 2:
            infile1, infile2 = infiles
            track1 = os.path.basename(infile1)
            track2 = makeSecond(track1)
            base_track = P.snip(track1, ".fastq.1.gz")
            track1_single = re.sub(".fastq.1.gz", ".fastq.gz", track1)
            outfile_single = ("%(outdir)s/%(prefix)s%(track1_single)s"
                              % locals())
            outfile1 = "%(outdir)s/%(prefix)s%(track1)s" % locals()
            outfile2 = "%(outdir)s/%(prefix)s%(track2)s" % locals()

            logfile = "log.dir/" + track1 + ".sickle.log"

            cmd = '''flash %(infile1)s %(infile2)s
            -p%(offset)s %(processing_options)s
            -o %(base_track)s -d %(outdir)s
            2>>%(logfile)s; checkpoint;
            cat %(outdir)s/%(base_track)s.extendedFrags.fastq | gzip >
            %(outfile_single)s;
            cat %(outdir)s/%(base_track)s.notCombined_1.fastq | gzip >
            %(outfile1)s;
            cat %(outdir)s/%(base_track)s.notCombined_1.fastq | gzip >
            %(outfile2)s;
            rm -rf %(outdir)s/%(base_track)s.extendedFrags.fastq
            %(outdir)s/%(base_track)s.hist
            %(outdir)s/%(base_track)s.histogram
            %(outdir)s/%(base_track)s.notCombined_1.fastq
            %(outdir)s/%(base_track)s.notCombined_2.fastq
            ;''' % locals()

            outfiles = (outfile_single,)

        else:
            E.info("flash requires paired end reads")

        self.num_files = 1
        # need to set num_file to 1 as only one outfile produced

        return cmd, outfiles

    def postprocess(self, infiles):
        # postprocess used to summarise single end fastq twice so that
        # it gets concatenated at the end of the summary tables for
        # both initial paired end files. if a further single end step
        # is included after flash, currently, it will break at the
        # summarise function in pipeline_preprocess
        scriptsdir = self.scriptsdir

        if self.summarise:
            outdir = self.outdir
            prefix = self.prefix
            infile1 = infiles[0]
            infile_base1 = os.path.basename(infile1)
            infile_base2 = re.sub(".1.fastq.gz", ".2.fastq.gz", infile_base1)
            infile = re.sub(".fastq.1.gz", ".fastq.gz", infile1)
            postprocess_cmd = '''zcat %(infile)s |
            python %(scriptsdir)s/fastq2summary.py
            --guess-format=illumina-1.8 -v0
            > summary.dir/%(infile_base1)s.summary;
            zcat %(infile)s |
            python %(scriptsdir)s/fastq2summary.py
            --guess-format=illumina-1.8 -v0
            > summary.dir/%(infile_base2)s.summary
            ;''' % locals()
        else:
            postprocess_cmd = "checkpoint ;"

        return postprocess_cmd

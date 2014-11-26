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

    def __init__(self, executable=None,
                 *args, **kwargs):
        if executable:
            self.executable = executable

    def getfastqAttr(self, infiles):
        num_files = len(infiles)

        if num_files > 1:
            infile1, infile2 = infiles
        else:
            infile1 = infiles[0]

        print "this is the infile:"
        print infile1
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
        '''build prcessing statement on infiles.
        list of processing tools is used to build a command line statement
        in list order
        '''
        total_cmd = ""
        processed_files = []
        print list_of_preprocessers
        file_compress = compress
        save_intermediates = save
        first = 1
        f_format = ""
        num_files = ""
        for idx, tool in enumerate(list_of_preprocessers):
            print tool
            print idx
            if idx == len(list_of_preprocessers)-1:
                end = True
            else:
                end = False
            if idx > 0:
                first = 0
            if tool == "cutadapt":
                tool_cmd = ""
            elif tool == "sickle":
                sickle_object = sickle(
                    compress=file_compress,
                    save=save_intermediates,
                    final=end,
                    f_format=f_format,
                    num_files=num_files,
                    first=first)
                tool_cmd, infile, f_format, num_files = (
                    sickle_object.build(infile))
            elif tool == "trimgalore":
                trimgalore_object = trimgalore(
                    compress=file_compress,
                    save=save_intermediates,
                    final=end,
                    f_format=f_format,
                    num_files=num_files,
                    first=first)
                tool_cmd, infile, f_format, num_files = (
                    trimgalore_object.build(infile))
                print "just finished first tool"
                print "num_files"
                print num_files
            else:
                tool_cmd = "testing_testing"
                # insert error statement here.
                pass

            total_cmd += tool_cmd
            processed_files.append(infile,)
        return total_cmd, processed_files

    def postprocess(self, initial_file, processer_list):
        '''collect final output data and postprocess.
        summarise initial fastq and then all processed fastqs'''

        if len(initial_file) == 1:
            initial_file = os.path.basename(initial_file[0])
            postprocess_cmd = '''zcat %(initial_file)s |
                          python %%(scriptsdir)s/fastq2summary.py
                          --guess-format=illumina-1.8 -v0
                          > processed.dir/%(initial_file)s.summary;
                          ''' % locals()

            processed_file = initial_file
            for tool in processer_list:
                print "tool"
                print tool
                print "processed_file"
                print processed_file
                processed_file = tool + "-" + processed_file
                postprocess_cmd += '''zcat processed.dir/%(processed_file)s |
                              python %%(scriptsdir)s/fastq2summary.py
                              --guess-format=illumina-1.8 -v0
                              > processed.dir/%(processed_file)s.summary;
                              ''' % locals()

        # edit from here to allow processing of paired end sequences

        elif len(initial_file) == 2:
            initial_file1, initial_file2 = [os.path.basename(x)
                                            for x in initial_file]
            postprocess_cmd = '''zcat %(initial_file1)s |
                          python %%(scriptsdir)s/fastq2summary.py
                          --guess-format=illumina-1.8 -v0
                          > processed.dir/%(initial_file1)s.summary;
                          zcat %(initial_file2)s |
                          python %%(scriptsdir)s/fastq2summary.py
                          --guess-format=illumina-1.8 -v0
                          > processed.dir/%(initial_file2)s.summary;
                          ''' % locals()

            processed_file1 = initial_file1
            processed_file2 = initial_file2
            for tool in processer_list:
                processed_file1 = tool + "-" + processed_file1
                processed_file2 = tool + "-" + processed_file2
                postprocess_cmd += '''zcat processed.dir/%(processed_file1)s |
                              python %%(scriptsdir)s/fastq2summary.py
                              --guess-format=illumina-1.8 -v0
                              > processed.dir/%(processed_file1)s.summary;
                              zcat processed.dir/%(processed_file2)s |
                              python %%(scriptsdir)s/fastq2summary.py
                              --guess-format=illumina-1.8 -v0
                              > processed.dir/%(processed_file2)s.summary;
                              ''' % locals()

        return postprocess_cmd

    def cleanup(self, outfile):
        '''clean up.'''
        statement = '''rm -rf %s;''' % (self.tmpdir_fastq)
        return statement

    def build(self, infile, outfile, processer_list):
        '''run mapper.'''
        cmd_preprocess, raw_file = self.preprocess(infile, outfile)
        raw_file = raw_file[0]
        cmd_process, processed_files = self.process(raw_file,
                                                    processer_list, outfile)
        cmd_postprocess = self.postprocess(raw_file, processer_list)
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
    def __init__(self, first=True, final=True, compress=False,
                 save=True, f_format="", num_files="", *args, **kwargs):
        self.final = final
        self.compress = compress
        self.first = first
        self.f_format = f_format
        self.num_files = num_files
        if self.final:
            self.save = True
        else:
            self.save = save
        print "num files at initiate"
        print self.num_files

    def getfastqAttr(self, infiles):

        if self.first:
            self.num_files = len(infiles)

        if self.num_files == 1:
            infile1 = infiles[0]
        elif self.num_files == 2:
            infile1, infile2 = infiles
        print "num files at get attr"
        print self.num_files
        print self.first

        if self.first:
            self.f_format = Fastq.guessFormat(IOTools.openFile(infile1),
                                              raises=False)
            E.info("%s: format guess: %s" % (infile1, self.f_format))
            print("%s: format guess: %s" % (infile1, self.f_format))
        else:
            E.info("%s: format set as previous: %s" % (infile1, self.f_format))
            print("%s: format set as previous: %s" % (infile1, self.f_format))

    def process(self, infiles):
        pass

    def postprocess(self, logfile):
        return ""

    def build(self, infiles):
        print "num files at get build"
        print self.num_files

        self.getfastqAttr(infiles)
        process_cmd, fastq_outfiles = self.process(infiles)
        post_cmd = self.postprocess(fastq_outfiles)
        statement = " checkpoint; ".join([process_cmd,
                                          post_cmd])

        return statement, fastq_outfiles, self.f_format, self.num_files


class trimgalore(processer):

    def process(self, infiles):

        if self.first:
            if self.num_files == 2:
                infile1, infile2 = infiles
            elif self.num_files == 1:
                infile1 = infiles[0]
        else:
            if self.num_files == 2:
                infile1, infile2 = infiles
            elif self.num_files == 1:
                infile1 = infiles

        offset = Fastq.getOffset(self.f_format, raises=False)
        if self.save:
            outdir = "processed.dir"
        else:
            outdir = P.getTempDir("/ifs/scratch")

        if self.num_files == 1:
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

        elif self.num_files == 2:
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


class sickle(processer):

    def process(self, infiles):

        rRANGES = {33: 'sanger', 64: 'illumina-1.8', 59: 'solexa'}
        offset = Fastq.getOffset(self.f_format, raises=False)
        quality = rRANGES[offset]

        if self.first:
            if self.num_files == 2:
                infile1, infile2 = infiles
            elif self.num_files == 1:
                infile1 = infiles[0]
        else:
            if self.num_files == 2:
                infile1, infile2 = infiles
            elif self.num_files == 1:
                infile1 = infiles

        if self.save:
            outdir = "processed.dir"
        else:
            outdir = P.getTempDir("/ifs/scratch")

        if self.num_files == 1:
            track = os.path.basename(infile1)
            outfile = "%(outdir)s/sickle-%(track)s" % locals()
            logfile = outdir + "/" + track + ".sickle.log"

            cmd = '''sickle se -g %%(sickle_options)s
                     --qual-type %(quality)s
                     --output-file %(outfile)s
                     --fastq-file %(infile1)s
                     2>>%(logfile)s
                     ;''' % locals()
            outfiles = outfile

        elif self.num_files == 2:
            track1 = os.path.basename(infile1)
            track2 = P.snip(track1, ".fastq.1.gz") + ".fastq.2.gz"
            outfile1 = "%(outdir)s/sickle-%(track1)s" % locals()
            outfile2 = "%(outdir)s/sickle-%(track2)s" % locals()

            logfile = outdir + "/" + track1 + ".trim_galore.log"

            cmd = '''sickle pe -g -s %%(sickle_options)s
                     --qual-type %(quality)s
                     -f %(infile1)s
                     -r %(infile2)s
                     -o %(outfile1)s
                     -p %(outfile2)s
                     2>>%(logfile)s
                     ;''' % locals()
            outfiles = (outfile1, outfile2)

        else:
            cmd = "this shouldn't happen! - throw an error!"

        return cmd, outfiles

    def postprocess(self, logfile):
        return ""

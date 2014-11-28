##########################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Tildon Grant Belgard
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
##########################################################################

"""====================
Pre-process pipeline
====================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python

The purpose of this pipeline is to pre-process reads (target
``full``).

Implemented tools and their functionalities are listed below:

+-------------+-------------------+-------------------------------------------+
|tool         |single/paired end  |functionality                              |
+-------------+-------------------+-------------------------------------------+
|fastx_trimmer|single             |hard-trimming                              |
+-------------+-------------------+-------------------------------------------+
|sickle       |both               |sliding-window quality trimming            |
+-------------+-------------------+-------------------------------------------+
|trimgalore   |both               |adapter-trimming                           |
|             |                   |end of read quality trimming               |
|             |                   |RRBS-specfic trimming                      |
+-------------+-------------------+-------------------------------------------+
|trimmomatic  |both               |highly optimisable adapter trimming        |
|             |                   |sliding window quality trimming            |
|             |                   |end/start of read quality trimming         |
|             |                   |hard-trimming                              |
|             |                   |user-specified order of processing         |
+-------------+-------------------+-------------------------------------------+
|flash        |paired             |overlapping paired end merging             |
+-------------+-------------------+-------------------------------------------+

Tools are enabled in the configuration file as a string of comma seperated
tool names

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

No general configuration required.

Input
-----

Reads are imported by placing files or linking to files in the
:term:`working directory`.

The default file format assumes the following convention:

   <sample>-<condition>-<replicate>.<suffix>

``sample`` and ``condition`` make up an :term:`experiment`, while
``replicate`` denotes the :term:`replicate` within an
:term:`experiment`. The ``suffix`` determines the file type.  The
following suffixes/file types are possible:

sra
   Short-Read Archive format. Reads will be extracted using the
   :file:`fastq-dump` tool.

fastq.gz
   Single-end reads in fastq format.

fastq.1.gz, fastq2.2.gz
   Paired-end reads in fastq format. The two fastq files must be
   sorted by read-pair.

.. note::

   Quality scores need to be of the same scale for all input
   files. Thus it might be difficult to mix different formats.

Requirements
------------

* sra-tools
* trimmomatic >= 0.32(optional)
* fastx >= 0.0.13 (optional)
* trim_galore >= 0.3.3 (optional)
* flash >= 1.2.6(optional)
* sickle >= 1.33 (optional)

Pipeline output
===============

By default the pipeline will output only the final processed files.
To save all intermediate files, use save=1 in pipeline.ini

The pipeline also outputs a summary table for each infile

Example
=======

Note: not currently available
Example data is available at
http://www.cgat.org/~andreas/sample_data/pipeline_preprocessing.tgz
To run the example, simply unpack and untar::

   wget http://www.cgat.org/~andreas/sample_data/pipeline_preprocessing.tgz
   tar -xvzf pipeline_preprocessing.tgz
   cd pipeline_preprocessing
   python <srcdir>/pipeline_preprocessing.py make full

To do
=====

Example data
Add threading
Better way to specify adapters?
Remove old code


Code
====

"""

###################################################
###################################################
###################################################
# load modules
###################################################

# import ruffus
from ruffus import *

# import useful standard python modules
import sys
import os
import re
import itertools
import glob
import sqlite3
import string
import pandas as pd
import CGAT.IOTools as IOTools
import CGAT.FastaIterator as FastaIterator
import CGATPipelines.PipelineMapping as PipelineMapping
import CGATPipelines.PipelinePreprocess as PipelinePreprocess
import CGAT.Pipeline as P


###################################################
###################################################
###################################################
# Pipeline configuration
###################################################

# load options from the config file
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])
PARAMS = P.PARAMS


###################################################################
# Helper functions mapping tracks to conditions, etc
###################################################################
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except NameError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]  # not recommended practise


def connect():
    '''connect to database.
    Use this method to connect to additional databases.
    Returns a database connection.
    '''
    dbh = sqlite3.connect(PARAMS["database"])

    return dbh

#########################################################################
#########################################################################
# Read processing
#########################################################################
INPUT_FORMATS = ("*.fastq.1.gz", "*.fastq.gz", "*.sra", "*.csfasta.gz")
REGEX_FORMATS = regex(r"(\S+).(fastq.1.gz|fastq.gz|sra|csfasta.gz)")

SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.fa.gz",
                    "*.sra",
                    "*.export.txt.gz",
                    "*.csfasta.gz",
                    "*.csfasta.F3.gz",
                    )

SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                      for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(r".*/(\S+).(?P<suffix>fastq.1.gz|fastq.gz|fa.gz|\
sra|csfasta.gz|csfasta.F3.gz|export.txt.gz)")

PREPROCESSTOOLS = [tool for tool in
                   P.asList(PARAMS["general_preprocessors"])]
preprocess_prefix = ("-".join(PREPROCESSTOOLS[::-1]) + "-")

print "preprocessing tools: %s",
print PREPROCESSTOOLS
print "\n\n\n\n\n"


@follows(mkdir("processed.dir"),
         mkdir("log.dir"),
         mkdir("summary.dir"))
@transform(SEQUENCEFILES,
           SEQUENCEFILES_REGEX,
           r"processed.dir/%s\1.\g<suffix>" % preprocess_prefix)
def processReads(infile, outfile):
    '''process reads from .fastq or .sra files.

    Tasks specified in PREPROCESSTOOLS are run in order

    '''
    trimmomatic_options = PARAMS["trimmomatic_options"]
    if PARAMS["trimmomatic_adapter"]:
        adapter_options = "ILLUMINACLIP:%s:%s:%s:%s " % (
            PARAMS["trimmomatic_adapter"], PARAMS["trimmomatic_mismatches"],
            PARAMS["trimmomatic_p_thresh"], PARAMS["trimmomatic_c_thresh"])
        trimmomatic_options = adapter_options + trimmomatic_options

    job_threads = PARAMS["general_threads"]
    job_options = "-l mem_free=%s" % PARAMS["general_memory"]
    save = PARAMS["general_save"]

    m = PipelinePreprocess.MasterProcessor(save=save)
    statement = m.build((infile,), outfile, PREPROCESSTOOLS)
    print statement
    P.run()


@transform(processReads,
           SEQUENCEFILES_REGEX,
           r"summary.dir/%s\1.\g<suffix>.processing.tsv" % preprocess_prefix)
def mergeSummaries(infile, outfile):
    '''summarise fastq statistics'''
    initial_file = re.sub(preprocess_prefix, "", os.path.basename(infile))
    initial_file = "summary.dir/" + initial_file + ".summary"
    files = [initial_file]
    print files

    def catSummaries(init, pre_tools_list, files, out):
        init_base = os.path.basename(init)
        stage = ["preprocessing"]
        tool_list = []
        df = pd.read_csv(init, sep="\t", )
        print df
        for tool in pre_tools_list:
            tool_list.append(tool)
            intermediate_file = ("summary.dir/%s-%s" % (
                "-".join(tool_list[::-1]), init_base))
            print "intermediate_file: %s" % intermediate_file
            files.append(intermediate_file)
            temp_df = pd.read_csv(intermediate_file, sep="\t")
            df = pd.concat([df, temp_df], axis=0)
        stage.extend(tool_list)
        df["stage"] = stage
        df["files"] = files
        df.to_csv(out, sep="\t", index=False)

    catSummaries(initial_file, PREPROCESSTOOLS, files, outfile)

    if re.match(".*fastq.1.gz", infile):

        outfile2 = PipelinePreprocess.makeSecond(outfile)
        initial_file2 = PipelinePreprocess.makeSecond(initial_file)
        files2 = [initial_file2]

        catSummaries(initial_file2, PREPROCESSTOOLS, files2, outfile2)


@transform(mergeSummaries,
           suffix(".tsv"),
           ".load")
def loadSummaries(infile, outfile):
    dbh = connect()
    preprocess_prefix2 = ("_".join(PREPROCESSTOOLS[::-1]) + "_")
    tablename = re.sub(preprocess_prefix2, "", P.toTable(outfile))
    scriptsdir = PARAMS["general_scriptsdir"]

    statement = '''cat %(infile)s |
                python %(scriptsdir)s/csv2db.py
                --table %(tablename)s --retry --ignore-empty
                 > %(outfile)s'''

    if re.match(".*fastq.1.gz.processing.tsv", infile):
        infile2 = PipelinePreprocess.makeSecond(infile)
        tablename2 = PipelinePreprocess.makeSecond(tablename)
        outfile2 = PipelinePreprocess.makeSecond(outfile)
        statement += ''';cat %(infile2)s |
                        python %(scriptsdir)s/csv2db.py
                        --table %(tablename2)s --retry --ignore-empty
                        > %(outfile2)s'''
    P.run()

if PARAMS["general_fastqc"]:
    @follows(mkdir(PARAMS["exportdir"]),
             mkdir(os.path.join(PARAMS["exportdir"], "fastqc")))
    @transform(INPUT_FORMATS,
               REGEX_FORMATS,
               r"\1.fastqc")
    def runFastqcInitial(infiles, outfile):
        '''convert sra files to fastq and check mapping qualities are in solexa format.
        Perform quality control checks on reads from .fastq files.'''
        m = PipelineMapping.FastQc(nogroup=PARAMS["readqc_no_group"],
                                   outdir=PARAMS["exportdir"]+"/fastqc")
        statement = m.build((infiles,), outfile)
        print "infiles", infiles, outfile
        P.run()

    @transform(processReads,
               REGEX_FORMATS,
               r"\1.fastqc")
    def runFastqcFinal(infiles, outfile):
        '''convert sra files to fastq and check mapping qualities are in solexa format.
        Perform quality control checks on reads from .fastq files.'''
        m = PipelineMapping.FastQc(nogroup=PARAMS["readqc_no_group"],
                                   outdir=PARAMS["exportdir"]+"/fastqc")
        statement = m.build((infiles,), outfile)
        print "infiles", infiles, outfile
        P.run()
else:
    def runFastqcInitial():
        pass

    def runFastqcFinal():
        pass


########################################################################
# all functions from here need to be removed
# code is being kept temporarily whilst the pipeline is being refactored
########################################################################


#########################################################################
# adaptor trimming
#########################################################################
# these are the adaptors and PCR primers for used in various illumina
# libarary preps see
# https://cgatwiki.anat.ox.ac.uk/xwiki/bin/view/CGAT/Illumina+Sequencing#HIlluminaAdaptors.html
# currently included are primers/adaptors from: TruSeq DNA HT and RNA
# HT Sample Prep Kits; TruSeq DNA v1/v2/LT RNA v1/v2/LT and ChIP
# Sample Prep Kits; Oligonucleotide Sequences for TruSeq Small RNA
# Sample Prep Kits; Oligonucleotide Sequences for Genomic DNA;
# Oligonucleotide Sequences for the v1 and v1.5 Small RNA Sample Prep
# Kits; Paired End DNA Oligonucleotide Sequences; Script Seq Adaptors;
# Oligonucleotide Sequences for the Multiplexing Sample Prep Oligo
# Only Kits.
ILLUMINA_ADAPTORS = {
    "Genomic-DNA-Adaptor": "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    "Genomic/Paired-End/Oligo-Only-Adaptor":
    "ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "Genomic/TruSeq-Universal/PE/OO/ScriptSeq-Adaptor/PCR-Primer":
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "Genomic-PCR-Primer": "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
    "Paired-End-Adaptor": "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
    "Paired-End-PCR-Primer":
    "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
    "TruSeq-HT-Adaptor-I3":
    "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGC",
    "TruSeq-Adaptor-I7": "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGGTTCTATCTCGTAT",
    "TruSeq-Adaptor-I4": "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGCCGGCTATCTCGTAT",
    "TruSeq-HT-Adaptor-I5":
    "AATGATACGGCGACCACCGAGATCTACACNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "TruSeq-HT-Adaptor-I7":
    "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    "TruSeq-LT-Adaptor-I6":
    "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
    "TruSeq-Adaptor-I11": "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGC",
    "TruSeq-Small-RNA-Adaptor": "TGGAATTCTCGGGTGCCAAGG",
    "TruSeq-Small-RNA-RT-Primer": "GCCTTGGCACCCGAGAATTCCA",
    "TruSeq-Small-RNA-PCR-Primer":
    "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",
    "TruSeq-Small-RNA-PCR-Primer-I6":
    "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    "Oligo-Only-Adaptor": "GATCGGAAGAGCACACGTCT",
    "Oligo-Only-PCR-Primer": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "Oligo-Only-PCR-Primer-I7": "CAAGCAGAAGACGGCATACGAGATNNNNNNNTGACTGGAGTTC",
    "Small-RNA-v1-RT-Primer": "CAAGCAGAAGACGGCATACGA",
    "Small-RNA-PCR-Primer": "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
    "ScriptSeq-Adaptor-I6":
    "CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "Exo_WTCHG_V.1_IlluminaWTCHGPrimer1":
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "Exo_WTCHG_V.1_WTCHGIndex1":
    "CAAGCAGAAGACGGCATACGAGATAGTTAACAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "Exo_WTCHG_V.1_ChIP_exo_Adapt1": "AGATCGGAAGA",
    "Exo_WTCHG_V.1_ChIP_exo_Adapt1.1": "TCCCTACACGACGCTCTTCCGATCT",
    "Exo_WTCHG_V.1_ExtPr1": "CCTACACGACGCTCTTCCGATCT",
    "Exo_WTCHG_V.1_ChIP_exo_Adapt2": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "Exo_WTCHG_V.1_ChIP_exo_Adapt2.1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTC",
    "Exo_WTCHG_V.2_IlluminaWTCHGPrimer1":
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "Exo_WTCHG_V.2_WTCHGIndex1":
    "CAAGCAGAAGACGGCATACGAGATAGTTAACAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "Exo_WTCHG_V.2_ChIP_exo_Adapt1": "GATCGGAAGAGCGTCGTGTAGGGA",
    "Exo_WTCHG_V.2_ChIP_exo_Adapt1.1": "TCCCTACACGACGCTCTTCCGATCT",
    "Exo_WTCHG_V.2_ExtPr1": "CCTACACGACGCTCTTCCGATCT",
    "Exo_WTCHG_V.2_ChIP_exo_Adapt2": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "Exo_WTCHG_V.2_ChIP_exo_Adapt2.1": "GATCGGAAGAGCACACGTCTGAACTCCAGTC",
    "SmartIIA": "AAGCAGTGGTATCAACGCAGAGTAC",
    "Illumina-Nextera-v2-Primer1": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "Illumina-Nextera-v2-Primer2": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "Illumina-Paired-End_Primer2":
    "CGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGA",
    "Illumina-Single-End-Adpator2":
    "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGGAAAGGAAGAGCACACG",
    "Nextera-Transposon-End-Sequence": "AGATGTGTATAAGAGACAG",
    "Epicentre-Nextera-Primer1": "AATGATACGGCGACCACCGA",
    "Epicentre-Nextera-Primer2": "CAAGCAGAAGACGGCATACGA",
    "Epicentre-Nextera-Read1": "GCCTCCCTCGCGCCATC",
    "Epicentre-Nextera-Read2": "GCCTTGCCAGCCCGCTC",
    "NETseq-linker-1": "CTGTAGGCACCATCAAT",
    "NETseq-RT-primer_5prime": "ATCTCGTATGCCGTCTTCTGCTTG",
    "NETseq-RT-primer_3prime": "TCCGACGATCATTGATGGTGCCTACAG",
    "NETseq-PCR-primer-oLSC008-bc1":
    "AATGATACGGCGACCACCGAGATCTACACGATCGGAAGAGCA\
CACGTCTGAACTCCAGTCACATGCCATCCGACGATCATTGATGG"
}


@merge(None, "contaminants.fasta")
def outputContaminants(infile, outfile):
    '''output file with contaminants.

    If contamination_reverse_complement is selected, then the reverse
    complement of each sequence is also written to the outfile.

    '''
    outf = IOTools.openFile(outfile, "w")
    for key, value in ILLUMINA_ADAPTORS.iteritems():
        outf.write(">%s\n%s\n" % (key, value))
        if PARAMS["contamination_reverse_complement"]:
            key_rc = key + "_rc"
            value_rc = value[::-1]
            value_rc = value_rc.translate(string.maketrans("ACTGN", "TGACN"))
            outf.write(">%s\n%s\n" % (key_rc, value_rc))
        else:
            continue
    outf.close()


def listAdaptors(infile):
    adaptors = []
    for entry in FastaIterator.FastaIterator(IOTools.openFile(infile)):
        adaptors.append(
            "%s %s" % (PARAMS["contamination_trim_type"], entry.sequence))
    adaptors = " ".join(adaptors)

    return adaptors


@transform([x for x in
            glob.glob("*.fastq.gz") + glob.glob("*.fastq.1.gz") +
            glob.glob("*.fastq.2.gz")
            if not x.startswith("nocontaminants")],
           regex(r"(\S+).(fastq.1.gz|fastq.gz|fastq.2.gz|csfasta.gz)"),
           add_inputs(outputContaminants),
           r"nocontaminants.\1.\2")
def removeContaminants(infiles, outfile):
    '''remove adaptor contamination from fastq files.

    This method uses cutadapt.
    '''

    infile, contaminant_file = infiles

    adaptors = []
    for entry in FastaIterator.FastaIterator(
            IOTools.openFile(contaminant_file)):
        adaptors.append("-a %s" % entry.sequence)

    adaptors = " ".join(adaptors)

    statement = '''
    cutadapt
    %(adaptors)s
    --overlap=%(contamination_min_overlap_length)i
    --format=fastq
    %(contamination_options)s
    <( zcat < %(infile)s )
    2> %(outfile)s.log
    | gzip > %(outfile)s
    '''
    P.run()


def checkPairs(infile):
    '''check for paired read files'''
    if infile.endswith(".fastq.1.gz"):
        infile2 = P.snip(infile, ".fastq.1.gz") + ".fastq.2.gz"
        assert os.path.exists(
            infile2), "second part of read pair (%s) missing" % infile2
    else:
        infile2 = None

    return infile2


def parseCutadapt(lines):
    '''parse cutadapt output.

    Multiple cutadapt outputs are joined.
    '''

    def _chunker(inf):
        chunk = []
        for line in inf:
            if line.startswith("==="):
                if chunk:
                    yield chunk
                chunk = []
            chunk.append(line)

    assert lines[0].startswith("cutadapt")
    results = {}

    del lines[0]
    for x, line in enumerate(lines):
        if not line.strip():
            continue
        if ":" in line:
            if line.strip().startswith("Command"):
                continue
            param, value = line[:-1].split(":")
            param = re.sub(" ", "_", param.strip()).lower()
            value = re.sub("[a-zA-Z ].*", "", value.strip())
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except:
                    pass
            results[param] = value
        else:
            break

    del lines[:x]
    results["unchanged_reads"] = int(
        results["processed_reads"]) - int(results["trimmed_reads"])
    headers = results.keys()

    adapters = {}
    for chunk in _chunker(lines):
        adapter = re.search("=== (.*) ===", chunk[0]).groups()[0]
        length, removed = re.search(
            "Adapter '.*', length (\d+), was trimmed (\d+) times",
            chunk[2]).groups()

        adapters[adapter] = length, removed

    return results, adapters

##################################################################
##################################################################
##################################################################


@transform(processReads,
           suffix(""),
           ".tsv")
def summarizeProcessing(infile, outfile):
    '''build processing summary.'''

    def _parseLog(inf, step):

        inputs, outputs = [], []
        if step == "reconcile":
            for line in inf:
                x = re.search(
                    "first pair: (\d+) reads, second pair: (\d+) \
reads,shared: (\d+) reads", line)
                if x:
                    i1, i2, o = map(int, x.groups())
                    inputs = [i1, i2]
                    outputs = [o, o]
                    break
        elif step == "contaminants":
            lines = inf.readlines()
            assert lines[0].startswith("cutadapt")
            lines = "@@@".join(lines)
            for part in lines.split("cutadapt")[1:]:
                results, adapters = parseCutadapt(
                    ("cutadapt" + part).split("@@@"))
                inputs.append(results["processed_reads"])
                outputs.append(results["unchanged_reads"])
        else:
            for line in inf:
                if line.startswith("Input:"):
                    inputs.append(
                        int(re.match("Input: (\d+) reads.", line).groups()[0]))
                elif line.startswith("Output:"):
                    outputs.append(
                        int(re.match("Output: (\d+) reads.",
                                     line).groups()[0]))

        return zip(inputs, outputs)

    infile2 = checkPairs(infile)
    if infile2:
        track = P.snip(infile, ".fastq.1.gz")
    else:
        track = P.snip(infile, ".fastq.gz")

    outf = IOTools.openFile(outfile, "w")
    outf.write("track\tstep\tpair\tinput\toutput\n")

    for step in "contaminants", "artifacts", "trim", "filter", "reconcile":
        fn = infile + "_%s.log" % step
        if not os.path.exists(fn):
            continue
        for x, v in enumerate(_parseLog(IOTools.openFile(fn), step)):
            outf.write("%s\t%s\t%i\t%i\t%i\n" % (track, step, x, v[0], v[1]))

    outf.close()

#########################################################################
#########################################################################
#########################################################################


@jobs_limit(1, "db")
@transform(summarizeProcessing,
           regex(r"processed.(\S+).fastq.*.gz.tsv"),
           r"\1_processed.load")
def loadProcessingSummary(infile, outfile):
    '''load filtering summary.'''
    P.load(infile, outfile)

#########################################################################
#########################################################################
#########################################################################


@merge(summarizeProcessing, "processing_summary.tsv")
def summarizeAllProcessing(infiles, outfile):
    '''summarize processing information.'''

    outf = IOTools.openFile(outfile, "w")
    data = []
    for infile in infiles:
        inf = IOTools.openFile(infile)
        for line in inf:
            track, step, pair, ninput, noutput = line[:-1].split("\t")
            if track == "track":
                continue
            data.append((track, step, pair, ninput, noutput))

    # sort by track, pair, input
    data.sort(key=lambda x: (x[0], x[2], -int(x[3])))
    first = True
    for key, v in itertools.groupby(data, lambda x: (x[0], x[2])):
        vals = list(v)
        track, pair = key
        ninput = int(vals[0][3])
        outputs = [int(x[4]) for x in vals]
        if first:
            outf.write("track\tpair\tninput\t%s\t%s\t%s\t%s\n" % (
                "\t".join([x[1] for x in vals]),
                "noutput",
                "\t".join(
                    ["percent_%s" % x[1] for x in vals]),
                "percent_output"))
            first = False
        outf.write("%s\t%s\t%i\t%s\t%i\t%s\t%s\n" % (track, pair, ninput,
                                                     "\t".join(
                                                         map(str, outputs)),
                                                     outputs[-1],
                                                     "\t".join(
                                                         ["%5.2f" % (100.0 * x
                                                                     / ninput)
                                                          for x in outputs]),
                                                     "%5.2f" % (100.0 *
                                                                outputs[-1] /
                                                                ninput)))
    outf.close()

#########################################################################
#########################################################################
#########################################################################


@jobs_limit(1, "db")
@transform(summarizeAllProcessing, suffix(".tsv"), ".load")
def loadAllProcessingSummary(infile, outfile):
    P.load(infile, outfile)

#########################################################################
#########################################################################
#########################################################################


@merge(removeContaminants, "filtering.summary.tsv.gz")
def summarizeFiltering(infiles, outfile):
    '''collect summary output from filtering stage.'''

    tracks = {}
    adapters = {}

    for f in infiles:
        track = f[len("nocontaminants."):]
        track = re.sub("[.].*", "", track)
        result, adapter = parseCutadapt(IOTools.openFile(f + ".log"))
        tracks[track] = result
        adapters[track] = adapter
        header = result.keys()

    outf = IOTools.openFile(outfile, "w")
    outf.write("track\t%s\n" % "\t".join(headers))

    for track, results in tracks.iteritems():
        outf.write("%s\t%s\n" %
                   (track, "\t".join(str(results[x]) for x in headers)))
    outf.close()

#########################################################################
#########################################################################
#########################################################################


@jobs_limit(1, "db")
@transform(summarizeFiltering,
           suffix(".summary.tsv.gz"),
           "_summary.load")
def loadFilteringSummary(infile, outfile):
    '''load filtering summary.'''
    P.load(infile, outfile)

#########################################################################
#########################################################################
#########################################################################


@transform([x for x in glob.glob("*.fastq.gz") +
            glob.glob("*.fastq.1.gz") +
            glob.glob("*.fastq.2.gz")],
           regex(r"(\S+).(fastq.1.gz|fastq.gz|fastq.2.gz|csfasta.gz)"),
           r"trim.\1.\2")
def trimReads(infile, outfile):
    '''trim reads to desired length using fastx

    '''

    E.warn("deprecated - use processReads instead")

    to_cluster = True
    statement = '''zcat %(infile)s | fastx_trimmer %(trim_options)s 2>
    %(outfile)s.log | gzip > %(outfile)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################


@transform([x for x in glob.glob("*.fastq.gz") +
            glob.glob("*.fastq.1.gz") +
            glob.glob("*.fastq.2.gz")],
           regex(r"(\S+).(fastq.1.gz|fastq.gz|fastq.2.gz|csfasta.gz)"),
           r"replaced.\1.\2")
def replaceBaseWithN(infile, outfile):
    '''replaces the specified base with N'''

    to_cluster = True
    statement = '''python %(scriptsdir)s/fastq2N.py
    -i %(infile)s %(replace_options)s'''
    P.run()

#########################################################################
#########################################################################
#########################################################################
#########################################################################


@follows(processReads)
def trim():
    pass


@follows(mergeSummaries)
def summarise():
    pass


@follows(runFastqcInitial,
         runFastqcFinal,
         loadSummaries)
def full():
    pass

#########################################################################


@follows()
def publish():
    '''publish files.'''
    P.publish_report()


@follows(mkdir("report"))
def build_report():
    '''build report from scratch.'''

    E.info("starting documentation build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.'''

    E.info("updating documentation")
    P.run_report(clean=False)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

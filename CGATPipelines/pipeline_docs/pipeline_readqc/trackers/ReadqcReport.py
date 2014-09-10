import os
import sys
import re
import types
import itertools
import glob

from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
from collections import OrderedDict as odict

from CGATReport.ResultBlock import ResultBlock, ResultBlocks
from CGATReport import Utils

###################################################################
###################################################################
# parameterization

EXPORTDIR = P.get('readqc_exportdir', P.get('exportdir', 'export'))
DATADIR = P.get('readqc_datadir', P.get('datadir', '.'))
DATABASE = P.get('readqc_backend', P.get('sql_backend', 'sqlite:///./csvdb'))

###################################################################
# cf. pipeline_rnaseq.py
# This should be automatically gleaned from pipeline_rnaseq.py
###################################################################
import CGATPipelines.PipelineTracks as PipelineTracks

TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("%s/*.sra" % DATADIR), "(\S+).sra") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.gz" % DATADIR), "(\S+).fastq.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.1.gz" % DATADIR), "(\S+).fastq.1.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("*.csfasta.gz"), "(\S+).csfasta.gz")

###########################################################################


class ReadqcTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

##############################################################
##############################################################
##############################################################


class TrackerFastQC(ReadqcTracker):
    tracks = ["all"]

    def __call__(self, track, slice=None):

        edir = EXPORTDIR

        toc_text = []
        link_text = []

        tracks = sorted([x.asFile() for x in TRACKS])
        for track in tracks:

            for x, fn in enumerate(glob.glob(os.path.join(EXPORTDIR, "fastqc", "%s*_fastqc" % track))):
                y = x + 1
                toc_text.append("* %(track)s-%(y)i_" % locals())
                link_text.append(
                    ".. _%(track)s-%(y)i: %(fn)s/fastqc_report.html" % locals())

        toc_text = "\n".join(toc_text)
        link_text = "\n".join(link_text)
        rst_text = "\n%(toc_text)s\n\n%(link_text)s\n" % locals()
        return odict((("text", rst_text),))

##########################################################
##########################################################
##########################################################


class FilteringSummary(SingleTableTrackerRows):
    table = "filtering_summary"

##############################################################
##############################################################
##############################################################


class FastQCDetails(ReadqcTracker):
    tracks = ["all"]
    slices = ("duplication_levels",
              "kmer_profiles",
              "per_base_gc_content",
              "per_base_n_content",
              "per_base_quality",
              "per_base_sequence_content",
              "per_sequence_gc_content",
              "per_sequence_quality",
              "sequence_length_distribution")

    def __call__(self, track, slice=None):

        # note there are spaces behind the %(image)s directive to accomodate
        # for path substitution
        block = '''
.. figure:: %(image)s                                     
   :height: 300 
'''

        blocks = ResultBlocks()
        tracks = sorted([x.asFile() for x in TRACKS])

        for track in tracks:

            files = glob.glob(
                os.path.join(EXPORTDIR, "fastqc", "%s*_fastqc" % track))
            for x, fn in enumerate(sorted(files)):
                y = x + 1

                image = os.path.abspath(
                    os.path.join(fn, "Images", "%s.png" % slice))
                if not os.path.exists(image):
                    continue

                blocks.append(ResultBlock(text=block % locals(),
                                          title=os.path.basename(fn)))

        return odict((("rst", "\n".join(Utils.layoutBlocks(blocks, layout="columns-2"))),))


class FastqcSummary(ReadqcTracker):
    pattern = "(.*)_Basic_Statistics"
    slices = ("File type", "Filename", "Encoding",
              "Total Sequences", "Sequence Length", "%GC")

    def __call__(self, track, slice):
        return self.getAll("SELECT * FROM %(track)s_Basic_Statistics WHERE measure = '%(slice)s'")


class FastqcSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "basic_statistics_summary"


class ProcessingDetails(ReadqcTracker):

    '''return summary of the read processing steps.'''
    pattern = "(.*)_processed$"

    def __call__(self, track):
        return self.getAll( """SELECT pair,input,output,pair, 100.0 * output / input as percent 
                              FROM %(track)s_processed""" )


class ProcessingSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "processing_summary"


class CorrelationSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "binned_means_correlation"
    fields = ("sample",)


class GradientSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "binned_means_gradients"
    fields = ("sample",)


class GCContentSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_GC_Content"
    fields = ("GC_Content",)


class LengthSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_length"
    fields = ("length",)


class AASummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_AA"
    fields = ("AA",)


class ATSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_AT"
    fields = ("AT",)


class ACSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_AC"
    fields = ("AC",)


class AGSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_AG"
    fields = ("AG",)


class TASummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_TA"
    fields = ("TA",)


class TTSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_TT"
    fields = ("TT",)


class TCSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_TC"
    fields = ("TC",)


class TGSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_TG"
    fields = ("TG",)


class CASummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_CA"
    fields = ("CA",)


class CTSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_CT"
    fields = ("CT",)


class CCSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_CC"
    fields = ("CC",)


class CGSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_CG"
    fields = ("CG",)


class GASummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_GA"
    fields = ("GA",)


class GTSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_GT"
    fields = ("GT",)


class GCSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_GC"
    fields = ("GC",)


class GGSummary(ReadqcTracker, SingleTableTrackerRows):
    table = "means_binned_GG"
    fields = ("GG",)

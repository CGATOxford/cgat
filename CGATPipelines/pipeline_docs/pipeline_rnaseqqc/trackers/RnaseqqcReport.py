import re
import glob
import pandas as pd
from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
import CGATPipelines.PipelineTracks as PipelineTracks

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


TRACKS = PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
    glob.glob("%s/*.sra" % DATADIR), "(\S+).sra") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.gz" % DATADIR), "(\S+).fastq.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("%s/*.fastq.1.gz" % DATADIR), "(\S+).fastq.1.gz") +\
    PipelineTracks.Tracks(PipelineTracks.Sample).loadFromDirectory(
        glob.glob("*.csfasta.gz"), "(\S+).csfasta.gz")

###########################################################################


class RnaseqqcTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

##############################################################
##############################################################
##############################################################


class CorrelationSummaryA(RnaseqqcTracker):
    table = "binned_means_correlation"
    select = ["AA", "AT", "AC", "AG"]
    select = ",".join(select)

    def getTracks(self, subset=None):
        return ("all")

    def __call__(self, track,  slice=None):
        statement = ("SELECT sample,%(select)s FROM %(table)s")
        # fetch data
        df = pd.DataFrame.from_dict(self.getAll(statement))
        df['sample'] = map(lambda x: re.sub("_quant.sf", "", x), df['sample'])
        df = pd.melt(df, id_vars="sample")
        df2 = pd.DataFrame(map(lambda x: x.split("-"), df['sample']))
        df2.columns = ["id_"+str(x) for x in range(1, len(df2.columns)+1)]
        merged = pd.concat([df, df2], axis=1)
        return merged


class GradientSummaryA(CorrelationSummaryA):
    table = "binned_means_gradients"


class CorrelationSummaryT(CorrelationSummaryA):
    table = "binned_means_correlation"
    select = ["TA", "TT", "TC", "TG"]
    select = ",".join(select)


class GradientSummaryT(CorrelationSummaryT):
    table = "binned_means_gradients"


class CorrelationSummaryC(CorrelationSummaryA):
    table = "binned_means_correlation"
    select = ["CA", "CT", "CC", "CG"]
    select = ",".join(select)


class GradientSummaryC(CorrelationSummaryC):
    table = "binned_means_gradients"


class CorrelationSummaryG(CorrelationSummaryA):
    table = "binned_means_correlation"
    select = ["GA", "GT", "GC", "GG"]
    select = ",".join(select)


class GradientSummaryG(CorrelationSummaryG):
    table = "binned_means_gradients"


class CorrelationSummaryGC(CorrelationSummaryA):
    table = "binned_means_correlation"
    select = ["GC_Content", "length"]
    select = ",".join(select)


class GradientSummaryGC(CorrelationSummaryGC):
    table = "binned_means_gradients"


class BiasFactorPlot(RnaseqqcTracker):
    table = ""
    factor = ""

    def getTracks(self, subset=None):
        return ("all")

    def __call__(self, track, slice=None):
        statement = ("SELECT * FROM %(table)s")
        # fetch data
        df = pd.DataFrame.from_dict(self.getAll(statement))
        df = pd.melt(df, id_vars=self.factor)
        df['variable'] = map(lambda x: re.sub("_quant.sf", "", x),
                             df['variable'])
        df['value'] = ((df['value'] - min(df['value'])) /
                       (max(df['value'])-min(df['value'])))
        df2 = pd.DataFrame(map(lambda x: x.split("_"), df['variable']))
        df2.columns = ["id_"+str(x) for x in range(1, len(df2.columns)+1)]
        merged = pd.concat([df, df2], axis=1)
        return merged


class GCContentSummary(BiasFactorPlot):
    table = "means_binned_GC_Content"
    factor = "GC_Content"


class LengthSummary(BiasFactorPlot):
    table = "means_binned_length"
    factor = "length"


class AASummary(BiasFactorPlot):
    table = "means_binned_AA"
    factor = "AA"


class ATSummary(BiasFactorPlot):
    table = "means_binned_AT"
    factor = "AT"


class ACSummary(BiasFactorPlot):
    table = "means_binned_AC"
    factor = "AC"


class AGSummary(BiasFactorPlot):
    table = "means_binned_AG"
    factor = "AG"


class TASummary(BiasFactorPlot):
    table = "means_binned_TA"
    factor = "TA"


class TTSummary(BiasFactorPlot):
    table = "means_binned_TT"
    factor = "TT"


class TCSummary(BiasFactorPlot):
    table = "means_binned_TC"
    factor = "TC"


class TGSummary(BiasFactorPlot):
    table = "means_binned_TG"
    factor = "TG"


class CASummary(BiasFactorPlot):
    table = "means_binned_CA"
    factor = "CA"


class CTSummary(BiasFactorPlot):
    table = "means_binned_CT"
    factor = "CT"


class CCSummary(BiasFactorPlot):
    table = "means_binned_CC"
    factor = "CC"


class CGSummary(BiasFactorPlot):
    table = "means_binned_CG"
    factor = "CG"


class GASummary(BiasFactorPlot):
    table = "means_binned_GA"
    factor = "GA"


class GTSummary(BiasFactorPlot):
    table = "means_binned_GT"
    factor = "GT"


class GCSummary(BiasFactorPlot):
    table = "means_binned_GC"
    factor = "GC"


class GGSummary(BiasFactorPlot):
    table = "means_binned_GG"
    factor = "GG"

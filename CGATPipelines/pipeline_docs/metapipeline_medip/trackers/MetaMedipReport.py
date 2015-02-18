import os
import glob

from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
import Pipeline
import PipelineTracks

###################################################################
###################################################################
###################################################################
###################################################################
# Run configuration script
EXPORTDIR = P['metamedip_exportdir']
DATADIR = P['metamedip_datadir']
DATABASE = P['metamedip_backend']

###################################################################
# cf. pipeline_medip.py
# This should be automatically gleaned from pipeline_chipseq.py
###################################################################
PARAMS_PIPELINE = Pipeline.peekParameters(".",
                                          "metapipeline_medip.py")

Sample = PipelineTracks.Sample
TRACKS = PipelineTracks.Tracks(Sample).loadFromDirectory(
    glob.glob(os.path.join(DATADIR, "medip_*")),
    "medip_(\S+)")

Sample.setDefault("asTable")

###########################################################################
###########################################################################
###########################################################################


class MetaMedipTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)


class SummaryCalledDMRs(MetaMedipTracker):

    '''return summary stats on differentially called DMRs.'''

    table = "called_dmrs"

    def getTracks(self):
        return self.getValues("SELECT DISTINCT metatrack FROM %s" % self.table)

    def getSlices(self):
        return self.getValues("SELECT DISTINCT test FROM %s" % self.table)

    def __call__(self, track, slice):
        data = self.getAll( """SELECT ntested, nok, nsignificant, n2fold FROM %(table)s 
                                  WHERE metatrack = '%(track)s' AND 
                                        test = '%(slice)s' """ )
        return data


class SummaryMapping(MetaMedipTracker):

    '''return summary stats on mapping.'''

    table = "mapping"

    def getTracks(self):
        return self.getValues("SELECT DISTINCT metatrack FROM %s" % self.table)

    def getSlices(self):
        return self.getValues("SELECT DISTINCT track FROM %s" % self.table)

    def __call__(self, track, slice):
        data = self.getAll( """SELECT * FROM %(table)s
                                  WHERE metatrack = '%(track)s' AND 
                                        track = '%(slice)s' """ )
        return data


class SummaryCpGCoverage(MetaMedipTracker):

    '''return summary stats on differentially called DMRs.'''

    table = "cpg_coverage"

    def getTracks(self):
        return self.getValues("SELECT DISTINCT metatrack FROM %s" % self.table)

    def getSlices(self):
        return self.getValues("SELECT DISTINCT track FROM %s" % self.table)

    def __call__(self, track, slice):
        data = self.getAll( """SELECT coverage, ncovered, pcovered FROM %(table)s 
                                  WHERE metatrack = '%(track)s' AND 
                                        track = '%(slice)s'
                                  ORDER BY coverage""" )
        return data


class ReportsList(MetaMedipTracker):

    '''provide links to reports.'''

    tracks = TRACKS

    def __call__(self, track):
        datadir = os.path.abspath(DATADIR)
        return "`medip_%(track)s <%(datadir)s/medip_%(track)s/report/html/contents.html>`_" % locals()

# TODO:

# number of CpG covered
# use all reads?
# simulations?

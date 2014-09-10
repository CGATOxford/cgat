import os
import sys
import re
import types
import itertools
import glob
import pandas
import pandas.rpy.common

from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P
from collections import OrderedDict as odict

##########################################################################
##########################################################################
# parameterization

EXPORTDIR = P['IDR_exportdir']
DATADIR = P['IDR_datadir']
DATABASE = P['IDR_backend']

###############################################################################


class ProjectTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

##########################################################################
# Base Class... custom __call__ methods defined below
##########################################################################


class GetSummaryTables(TrackerSQL):
    """
    Base class... derive tracks from entries in sqlite table
    """
    table = None

    def getTracks(self, subset=None):
        table = self.table
        return self.getValues("SELECT DISTINCT Condition"
                              " FROM %s" % table)

    def getSlices(self, subset=None):
        table = self.table
        return self.getValues("SELECT DISTINCT Tissue"
                              " FROM %s" % table)

##########################################################################
# Retrieve tables summarizing the number of peaks called for each condition
##########################################################################


class GetPeakCallingTables(GetSummaryTables):
    """
    Return summary of the number of peaks called with lax peak calling
    """
    def __call__(self, track, slice=None):
        table = self.table
        statement = ("SELECT Sample_id,n_peaks"
                     " FROM %(table)s"
                     " WHERE Tissue = '%(slice)s'"
                     " AND Condition = '%(track)s'")
        return self.getAll(statement)


class GetIndividualReplicatePeakCallingSummary(GetPeakCallingTables):
    table = "peakcalling_summary_individual_replicates"


class GetPseudoreplicatePeakCallingSummary(GetPeakCallingTables):
    table = "peakcalling_summary_pseudoreplicates"


class GetPooledPseudoreplicatePeakCallingSummary(GetPeakCallingTables):
    table = "peakcalling_summary_pooled_pseudoreplicates"

##########################################################################
# Retrieve tables summarizing the number of peaks called for each condition
# at IDR thresholds specified in config file
##########################################################################


class GetIDRSummaryTables(GetSummaryTables):
    """
    Return summary of number of peaks called by IDR at specified threshold
    """
    def __call__(self, track, slice=None):
        table = self.table
        statement = ("SELECT sample_1 AS 'Sample_1_id',"
                     "n_peaks AS 'n_peaks_above_idr_threshold'"
                     " FROM %(table)s"
                     " WHERE Tissue = '%(slice)s'"
                     " AND Condition = '%(track)s'")
        return self.getAll(statement)


class GetIndividualReplicateIDRSummary(GetIDRSummaryTables):
    table = "idr_summary_individual_replicates"


class GetPseudoreplicateIDRSummary(GetIDRSummaryTables):
    table = "idr_summary_pseudoreplicates"


class GetPooledPseudoreplicateIDRSummary(GetIDRSummaryTables):
    table = "idr_summary_pooled_pseudoreplicates"

##########################################################################
# Retrieve tables giving the number of peaks called for each condition
# at various IDR thresholds
##########################################################################


class GetIDRTables(GetSummaryTables):
    """
    Return summary of number of peaks called by IDR at various thresholds
    for every pairwise combination.
    """
    table = "peakcalling_summary_individual_replicates"
    stub = None

    def __call__(self, track, slice):
        stub = self.stub
        return self.getAll("SELECT * FROM %(slice)s_%(track)s_%(stub)s")


class GetIndividualReplicateIDRTables(GetIDRTables):
    stub = "npeaks_aboveIDR"


class GetPseudoreplicateIDRTables(GetIDRTables):
    stub = "pseudoreplicates_npeaks_aboveIDR"


class GetPooledPseudoreplicateIDRTables(GetIDRTables):
    stub = "pooled_pseudoreplicates_npeaks_aboveIDR"

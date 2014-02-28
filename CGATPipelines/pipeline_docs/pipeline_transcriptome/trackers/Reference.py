import os
import sys
import re
import types
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats

from SphinxReport.Renderer import *
from SphinxReport.Tracker import *

# for trackers_derived_sets and trackers_master
if not os.path.exists("conf.py"):
    raise IOError("could not find conf.py")

execfile("conf.py")

##########################################################################
##########################################################################
##########################################################################
# Trackers that access reference statistics
##########################################################################


class ReferenceData(TrackerSQL):

    """Base class for Trackers accessing reference table."""
    mPattern = "_annotation$"
    mMaster = trackers_master


class CoverageByTranscripts(ReferenceData):

    """Coverage of reference gene models."""
    mXLabel = "overlap / %"

    def __call__(self, track, slice=None):
        if track == self.mMaster:
            return []
        return odict((("pover1",
                       [ x[0] for x in self.execute( """SELECT pover1 FROM %s_vs_%s WHERE nover > 0""" %
                                                     (self.mMaster,
                                                      track))]), ))


class OverlapByTranscripts(ReferenceData):

    """Number of transcript models overlapping a reference gene model."""
    mXLabel = "number of transcripts"

    def __call__(self, track, slice=None):
        if track == self.mMaster:
            return []
        return odict((("nover1",
                       [ x[0] for x in self.execute( """SELECT nover1 FROM %s_vs_%s WHERE nover > 0""" %
                                                     (self.mMaster,
                                                      track))]), ))


class CoverageStats(CoverageByTranscripts):

    """Coverage of reference gene models - statistics.

    100% covered: number of reference genes completely covered by transcript models
    90% covered: number of reference genes at least 90% covered by transcript models
    1to1: number of reference genes matching to a single transcript model
    1toM: number of reference genes matching to multiple transcript models
    """

    def __call__(self, track, slice=None):
        if track == self.mMaster:
            return []
        columns = Stats.Summary().getHeaders()
        x = CoverageByTranscripts.__call__(self, track, slice)
        stats = Stats.Summary(x["pover1"])
        data = stats.items()
        data.append( ("100% covered", self.getValue( """SELECT COUNT(*) FROM %s_vs_%s WHERE pover1>= 100""" %
                                                     (self.mMaster,
                                                      track))))
        data.append( ("90% covered", self.getValue( """SELECT COUNT(*) FROM %s_vs_%s WHERE pover1>= 90""" %
                                                    (self.mMaster,
                                                     track))))

        data.append( ("1to1", len(list(self.execute( """SELECT gene_id1 FROM %s_vs_%s_ovl GROUP BY gene_id1 HAVING COUNT (gene_id2) = 1""" %
                                                     (self.mMaster,
                                                      track))))))

        data.append( ("1toM", len(list(self.execute( """SELECT gene_id1 FROM %s_vs_%s_ovl GROUP BY gene_id1 HAVING COUNT (gene_id2) > 1""" %
                                                     (self.mMaster,
                                                      track))))))

        return odict(data)


class CoverageVsLengthByReadDepth(TrackerSQL):

    """plot the absolute coverage of a known gene versus its length.
    Dots are colored by read depth.
    """

    mPattern = "_coverage$"
    mMaster = "ensembl"
    mXLabel = "log(length)"

    def __call__(self, track, slice=None):
        master = self.mMaster
        statement = """SELECT AVG(d.exons_sum) AS ref_length,
                                MIN(c.pover1) AS coverage, 
                                AVG(a.mean) AS read_depth
                        FROM %(track)s_coverage AS a,
                                %(master)s_vs_%(track)s_ovl as b,
                                %(master)s_vs_%(track)s as c,
                                %(master)s_annotation as d
                        WHERE   d.gene_id = c.gene_id AND
                                a.gene_id = b.gene_id2 AND
                                c.gene_id = b.gene_id1
                        GROUP BY b.gene_id1"""

        data = [(math.log(x[0]), x[1], math.log(x[2]))
                for x in self.getAll(statement % locals())]

        return odict(zip(("log(length)", "log(coverage), log(read_depth)"), zip(*data)))


class LengthVsReadDepthByLength(TrackerSQL):

    """plot the relative coverage of known gene (1-100%)
    versus the read depth.

    Dots will be colored by the length of the known gene
    (in log scale).
    """

    mPattern = "_coverage$"
    mMaster = "ensembl"
    mXLabel = "log(read_depth)"
    mYLabel = "log(relative coverage)"

    def __call__(self, track, slice=None):
        master = self.mMaster
        statement = """SELECT a.mean AS read_depth, 
                                c.exons_sum AS transcript_length, 
                                d.exons_sum AS ref_length
                        FROM %(track)s_coverage AS a, \
                                %(master)s_vs_%(track)s_ovl as b,
                                %(track)s_annotation as c,
                                %(master)s_annotation as d
                        WHERE   d.gene_id = b.gene_id1 AND
                                a.gene_id = b.gene_id2 AND
                                c.gene_id = b.gene_id2 AND
                                a.mean > 0 AND 
                                c.exons_sum > 0 AND
                                d.exons_sum > 0
                                """

        data = [(math.log(x[0]), math.log(float(x[1]) / x[2]), math.log(x[2]))
                for x in self.getAll(statement % locals()) if x[2] > 0 and x[1] > 0 and x[0] > 0]

        return odict(zip((("log(read_depth)", "log(length_predicted/length_reference), log(length_reference)"), zip(*data))))

# =================================================================
# Coverage
# =================================================================


class MeanVsMaxReadDepth(ReferenceData):

    """maxmimum read depth versus mean read depth of :term:`reference` genes. 
    Dots are coloured by the log(length) of a :term:`reference` gene."""

    mXLabel = "mean read depth"
    mYLabel = "maximum read depth"

    def __call__(self, track, slice=None):
        master = self.mMaster
        statement = "SELECT cov_mean, cov_max, sum FROM %(master)s_vs_%(track)s_readcoverage" % locals(
        )
        data = [(x[0], x[1], math.log(x[2]))
                for x in self.getAll(statement) if x[2] > 0]
        return odict(zip(("mean coverage", "max coverage", "length"), zip(*data)))


class MeanVsMedianReadDepth(ReferenceData):

    """maxmimum read depth versus mean read depth of :term:`reference` genes. 
    Dots are coloured by the log(length) of a :term:`reference` gene."""

    mXLabel = "mean read depth"
    mYLabel = "median read depth"

    def __call__(self, track, slice=None):
        master = self.mMaster
        statement = "SELECT cov_mean, cov_median, sum FROM %(master)s_vs_%(track)s_readcoverage" % locals(
        )
        data = [(x[0], x[1], math.log(x[2]))
                for x in self.getAll(statement) if x[2] > 0]
        return odict(zip(("mean coverage", "median coverage", "length"), zip(*data)))

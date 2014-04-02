from IntervalReport import *

##########################################################################
##########################################################################
##########################################################################
# Annotation of bases with SNPs
##########################################################################


class IntervalsSummary(IntervalTracker):

    """Summary stats of intervals called by the peak finder.
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice=None):
        data = self.getFirstRow(
            "SELECT COUNT(*), AVG(length), SUM(nprobes)  FROM %(track)s_intervals" % locals())
        return odict(zip(("nintervals", "<length>", "nprobes"), data))

##########################################################################
##########################################################################
##########################################################################
# Distribution of interval lengths
##########################################################################


class IntervalLengths(IntervalTracker):

    """Distribution of interval sizes.
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues(
            "SELECT length FROM %(track)s_intervals" % locals())
        return {"length": data}

##########################################################################
##########################################################################
##########################################################################
# Distribution of peak values
##########################################################################


class IntervalPeakValues(IntervalTracker):

    """Distribution of peak values (the number of reads at peak).
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues(
            "SELECT peakval FROM %(track)s_intervals" % locals())
        return {"peakval": data}

##########################################################################
##########################################################################
##########################################################################
# Distribution of peak values
##########################################################################


class EnrichmentOverUnstim(IntervalTracker):

    """For each peakval, present the fold enrichment of a track
    compared to the Unstim set.

    Useful for diagnosing cutoffs.
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice=None):

        cellline, condition, replicate = splitTrack(track)
        if condition == "Unstim":
            return []

        if replicate is None:
            replicate = ""
        other_track = "run" + cellline + "Unstim" + replicate

        data_fg = self.getValues(
            "SELECT peakval FROM %(track)s_intervals ORDER BY peakval DESC" % locals())
        data_bg = self.getValues(
            "SELECT peakval FROM %(other_track)s_intervals ORDER BY peakval DESC" % locals())

        mi = min(data_bg + data_fg)
        ma = max(data_bg + data_fg)

        n_fg = len(data_fg)
        n_bg = len(data_bg)

        hist_fg, bin_edges = numpy.histogram(
            data_fg, bins=numpy.arange(mi, ma + 1.0, 1.0))
        hist_bg, bin_edges = numpy.histogram(
            data_bg, bins=numpy.arange(mi, ma + 1.0, 1.0))

        hist_fg = hist_fg[::-1].cumsum()[::-1]
        hist_bg = hist_bg[::-1].cumsum()[::-1]

        fold = []
        for fg, bg in zip(hist_fg, hist_bg):
            fold.append(float(bg) / fg)

        return odict((("peakval", bin_edges[:-1]),
                      ("fold", fold)))

##########################################################################
##########################################################################
##########################################################################
# Distribution of average values
##########################################################################


class IntervalAverageValues(IntervalTracker):

    """Distribution of average values (the average number of reads within the interval)
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice=None):
        data = self.getValues(
            "SELECT avgval FROM %(track)s_intervals" % locals())
        return {"avgval": data}

##########################################################################
##########################################################################
##########################################################################
# Distribution of average values
##########################################################################


class IntervalLengthVsAverageValue(IntervalTracker):

    """Length vs average value.
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT length, avgval FROM %(track)s_intervals" % locals())
        return odict(zip(("length", "avgval"), zip(*data)))

##########################################################################
##########################################################################
##########################################################################
# Distribution of average values
##########################################################################


class IntervalLengthVsPeakValue(IntervalTracker):

    """Length vs peak value
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT length, peakval FROM %(track)s_intervals" % locals())
        return odict(zip(("length", "peakval"), zip(*data)))

##########################################################################
##########################################################################
##########################################################################
# peak location within interval
##########################################################################


class PeakLocation(IntervalTracker):
    mPattern = "_intervals$"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT (PeakCenter - start) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_intervals" % locals())
        data2 = self.getValues(
            "SELECT (end - PeakCenter) / CAST( Length as FLOAT) - 0.5 FROM %(track)s_intervals" % locals())
        return {"distance": data1 + data2}

##########################################################################
##########################################################################
##########################################################################
# distance of peak to end of interval
##########################################################################


class PeakDistance(IntervalTracker):
    mPattern = "_intervals$"

    def __call__(self, track, slice=None):
        data1 = self.getValues(
            "SELECT PeakCenter - start FROM %(track)s_intervals" % locals())
        data2 = self.getValues(
            "SELECT end - PeakCenter FROM %(track)s_intervals" % locals())
        return {"distance": data1 + data2}

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class IntervalList(IntervalTracker):

    '''list of intervals.'''

    nresults = 10
    mColumnsFixed = ("pos", "length")
    mColumnsVariable = ("peakval", "avgval")
    mPattern = "_intervals$"

    def getSQLStatement(self, track, slice=None):

        statement = '''
          SELECT 
               i.interval_id, i.contig, i.start, i.end, i.peakval, i.avgval
          FROM 
               %(track)s_intervals AS i
          ORDER BY i.peakval DESC''' % locals()

        if self.nresults:
            statement += " LIMIT %i" % self.nresults

        return statement

    def __call__(self, track, slice=None):

        statement = self.getSQLStatement(track, slice)
        data = self.get(statement)
        ucsc_genome = UCSC_GENOME
        n = odict()
        for d in data:
            id, contig, start, end = d[:4]
            pos = "`%(contig)s:%(start)i..%(end)i <http://genome.ucsc.edu/cgi-bin/hgTracks?db=%(ucsc_genome)s&position=%(contig)s:%(start)i..%(end)i>`_" \
                % locals()
            n[str(id)] = odict(
                zip(self.mColumnsFixed + self.mColumnsVariable, (pos, end - start,) + d[4:]))

        return n

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class IntervalListFull(IntervalTracker):

    '''list of all intervals.

    Table for export.
    '''

    nresults = None
    mPattern = "_intervals$"

    def __call__(self, track, slice=None):

        statement = '''
          SELECT 
               i.contig, i.start, i.end, i.peakval, i.avgval
          FROM 
               %(track)s_intervals AS i
          ORDER BY i.peakval DESC''' % locals()

        data = self.get(statement)
        return odict(zip(
            ("contig", "start", "end", "peakval", "avgval"),
            zip(*data)))

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class IntervalListPeakval(IntervalList):

    '''list of intervals.'''

    def getSQLStatement(self, track, slice=None):

        nresults = self.nresults

        statement = '''
          SELECT 
               i.interval_id, i.contig, i.start, i.end, i.peakval, i.avgval, i.length
          FROM 
               %(track)s_intervals AS i
          ORDER BY i.peakval DESC
          LIMIT %(nresults)s''' % locals()

        return statement

##########################################################################
##########################################################################
##########################################################################
# correlations
##########################################################################


class Correlations(IntervalTracker):

    """Correlation between all sets.
    """

    pattern = "(.*)_annotations$"

    def __call__(self, track, slice=None):
        table = "%s_correlation" % self.field
        return self.getValues("SELECT %s AS %s FROM %s ORDER BY id" % (track,
                                                                       self.field,
                                                                       table))


class CorrelationsPeakval(Correlations):
    field = "peakval"


class CorrelationsAvgval(Correlations):
    field = "avgval"


class CorrelationsLength(Correlations):
    field = "length"

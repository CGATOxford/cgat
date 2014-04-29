from PeakcallingReport import *

import CGAT.IOTools


class FoldChangeTracker(TrackerSQL):

    '''the fold change tracker ignores unstimulated tracks.'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

    def getTracks(self, subset=None):
        tracks = TrackerSQL.getTracks(self, subset=subset)
        return [x for x in tracks if TAG_UNSTIM not in x]

##########################################################################
##########################################################################
##########################################################################


class PeaksIntervals:
    pattern = "(.*)_(.*)_peaks$"
    suffix = "peaks"


class RegionsIntervals:
    pattern = "(.*)_(.*)_regions$"
    suffix = "regions"


class SummitsIntervals:
    pattern = "(.*)_(.*)_summits$"
    suffix = "summits"

##########################################################################
##########################################################################
##########################################################################
# Summarise intervals as a table
##########################################################################


class IntervalsSummaryTable(DefaultTracker):

    """Summary stats of intervals called by the peak finder.
    """

    def __call__(self, track, slice):
        table = "%s_%s_%s" % (track, slice, self.suffix)
        if not self.hasTable(table):
            return
        # Get the column containing the recorded qvalues
        fdr_col_headings = ("fdr", "qvalue", "FDR")
        fdr_col_head = None
        for name in fdr_col_headings:
            try:
                test = self.getFirstRow(
                    "SELECT MAX(%s) FROM %s" % (name, table))
                fdr_col_head = name
            except:
                pass

        if fdr_col_head:
            data = self.getFirstRow(
                "SELECT COUNT(*), ROUND(AVG(end-start),2), MIN(end-start), MAX(end-start), MAX(%(fdr_col_head)s) as VARCHAR FROM %(table)s")
        else:
            data = self.getFirstRow(
                "SELECT COUNT(*), ROUND(AVG(end-start),2), MIN(end-start), MAX(end-start), 'unknown' FROM %(table)s")

        # Formatting for prettiness
        if isinstance(data[4], float):
            fdr_value = "%.3g" % data[4]
        else:
            fdr_value = "%s" % data[4]
        fdata = map(str, data[0:4]) + [fdr_value]

        return odict(zip(("nintervals", "avg(length)", "min(length)", "max(length)", "max reported %s" % fdr_col_head), fdata))


class PeaksSummaryTable(PeaksIntervals, IntervalsSummaryTable):
    pass


class RegionsSummaryTable(RegionsIntervals, IntervalsSummaryTable):
    pass


class SummitsSummaryTable(SummitsIntervals, IntervalsSummaryTable):
    pass

##########################################################################
##########################################################################
##########################################################################
# Summarise intervals
##########################################################################


class IntervalsSummary(DefaultTracker):

    """Summary stats of intervals called by the peak finder.
    """

    def __call__(self, track, slice):
        table = "%s_%s_%s" % (track, slice, self.suffix)
        if not self.hasTable(table):
            return
        data = self.getFirstRow(
            "SELECT COUNT(*), AVG(end-start), MIN(end-start), MAX(end-start) as VARCHAR FROM %(table)s")
        return odict(zip(("nintervals", "avg(length)", "min(length)", "max(length)"), data))


class PeaksSummary(PeaksIntervals, IntervalsSummary):
    pass


class RegionsSummary(RegionsIntervals, IntervalsSummary):
    pass


class SummitsSummary(SummitsIntervals, IntervalsSummary):
    pass


##########################################################################
##########################################################################
##########################################################################
# Distribution of interval lengths
##########################################################################
class IntervalLengths(DefaultTracker):

    """Distribution of interval sizes.
    """

    def __call__(self, track, slice=None):
        table = "%s_%s_%s" % (track, slice, self.suffix)
        if not self.hasTable(table):
            return
        
        data = self.getValues(
            "SELECT length FROM %(table)s")
        return {"length": data}


class PeaksLengths(PeaksIntervals, IntervalLengths):
    pass


class RegionsLengths(RegionsIntervals, IntervalLengths):
    pass


class SummitsLengths(SummitsIntervals, IntervalLengths):
    pass

##########################################################################
##########################################################################
##########################################################################
# Distribution of peak values
##########################################################################


class IntervalPeakValues(DefaultTracker):

    """Distribution of peak values (the number of reads at peak).
    """

    def __call__(self, track, slice=None):
        table = "%s_%s_%s" % (track, slice, self.suffix)
        if not self.hasTable(table):
            return
        data = self.getValues(
            "SELECT peakval FROM %(table)s")
        return {"peakval": data}


class PeaksPeakValues(PeaksIntervals, IntervalPeakValues):
    pass


class RegionsPeakValues(RegionsIntervals, IntervalPeakValues):
    pass


class SummitsPeakValues(SummitsIntervals, IntervalPeakValues):
    pass

##########################################################################
##########################################################################
##########################################################################
# Distribution of average values
##########################################################################


class IntervalAverageValues(DefaultTracker):

    """Distribution of average values (the average number of reads within the interval)
    """

    def __call__(self, track, slice=None):
        table = "%s_%s_%s" % (track, slice, self.suffix)
        if not self.hasTable(table):
            return

        data = self.getValues(
            "SELECT avgval FROM %(table)s")
        return {"avgval": data}


class PeaksAverageValues(PeaksIntervals, IntervalAverageValues):
    pass


class RegionsAverageValues(RegionsIntervals, IntervalAverageValues):
    pass


class SummitsAverageValues(SummitsIntervals, IntervalAverageValues):
    pass

##########################################################################
##########################################################################
##########################################################################
# peak location within interval
##########################################################################


class PeakLocation(DefaultTracker):

    def __call__(self, track, slice=None):
        table = "%s_%s_%s" % (track, slice, self.suffix)
        if not self.hasTable(table):
            return

        data1 = self.getValues(
            "SELECT (PeakCenter - start) / CAST( Length as FLOAT) - 0.5 FROM %(table)s")
        data2 = self.getValues(
            "SELECT (end - PeakCenter) / CAST( Length as FLOAT) - 0.5 FROM %(table)s")
        return {"distance": data1 + data2}


class PeaksPeakLocation(PeaksIntervals, PeakLocation):
    pass


class RegionsPeakLocation(RegionsIntervals, PeakLocation):
    pass


class SummitsPeakLocation(SummitsIntervals, PeakLocation):
    pass

##########################################################################
##########################################################################
##########################################################################
# distance of peak to end of interval
##########################################################################


class PeakDistance(DefaultTracker):

    def __call__(self, track, slice=None):
        table = "%s_%s_%s" % (track, slice, self.suffix)
        if not self.hasTable(table):
            return
        data1 = self.getValues(
            "SELECT PeakCenter - start FROM %(table)s")
        data2 = self.getValues(
            "SELECT end - PeakCenter FROM %(table)s")
        return {"distance": data1 + data2}


class PeaksPeakDistance(PeaksIntervals, PeakDistance):
    pass


class RegionsPeakDistance(RegionsIntervals, PeakDistance):
    pass


class SummitsPeakDistance(SummitsIntervals, PeakDistance):
    pass

##########################################################################
##########################################################################
##########################################################################
# Distribution of average values
##########################################################################


class LengthVsAverageValue(DefaultTracker):

    """Length vs average value.
    """

    def __call__(self, track, slice=None):
        table = "%s_%s_%s" % (track, slice, self.suffix)
        if not self.hasTable(table):
            return

        data = self.get(
            "SELECT length, avgval FROM %(table)s")
        return odict(zip(("length", "avgval"), zip(*data)))


class PeaksLengthVsAverageValue(PeaksIntervals, LengthVsAverageValue):
    pass


class RegionsLengthVsAverageValue(RegionsIntervals, LengthVsAverageValue):
    pass


class SummitsLengthVsAverageValue(SummitsIntervals, LengthVsAverageValue):
    pass

##########################################################################
##########################################################################
##########################################################################
# Distribution of average values
##########################################################################


class LengthVsPeakValue(DefaultTracker):

    """Length vs peak value
    """

    def __call__(self, track, slice=None):
        table = "%s_%s_%s" % (track, slice, self.suffix)
        if not self.hasTable(table):
            return

        data = self.get(
            "SELECT length, peakval FROM %(table)s")
        return odict(zip(("length", "peakval"), zip(*data)))


class PeaksLengthVsPeakValue(PeaksIntervals, LengthVsPeakValue):
    pass


class RegionsLengthVsPeakValue(RegionsIntervals, LengthVsPeakValue):
    pass


class SummitsLengthVsPeakValue(SummitsIntervals, LengthVsPeakValue):
    pass

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class IntervalList(DefaultTracker):

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


class IntervalListFull(DefaultTracker):

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
##
##########################################################################


class IntervalListFoldChange(FoldChangeTracker, IntervalList):

    '''list of intervals.'''

    mColumnsVariable = ("fg_anysense_mean", "bg_anysense_mean", "fold_mean",
                        "fg_anysense_max", "bg_anysense_max", "fold_max")
    mMinFold = 1.5
    mPattern = "_readcounts$"

    def getSQLStatement(self, track, slice=None):
        nresults = self.nresults
        minfold = self.mMinFold

        statement = '''
          SELECT 
               i.interval_id, i.contig, i.start, i.end, 
               fg_anysense_mean, bg_anysense_mean,
               fg_anysense_mean / bg_anysense_mean AS fold_mean,
               fg_anysense_max, bg_anysense_max,
               cast(fg_anysense_max as float)/ bg_anysense_max AS fold_max
          FROM 
               %(track)s_intervals AS i,
               %(track)s_readcounts AS c
          WHERE
               i.interval_id = c.gene_id AND 
               cast(fg_anysense_max as float)/ bg_anysense_max > %(minfold)f
          ORDER BY fold_max DESC
          LIMIT %(nresults)s''' % locals()

        return statement

##########################################################################
##########################################################################
##########################################################################
# correlations
##########################################################################


class Correlations(CallingTracker):

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

##########################################################################
##########################################################################
##########################################################################
# fold change
##########################################################################


class FoldChange(FoldChangeTracker):

    """return fold changes for all intervals.
    """
    pattern = "(.*)_readcounts$"

    def __call__(self, track):
        return self.getValues("SELECT CAST(fg_anysense_max AS float)/ bg_anysense_max AS fold FROM %(track)s_readcounts" % locals())

##########################################################################
##########################################################################
##########################################################################
# fold change
##########################################################################


class FoldChangeCounts(FoldChangeTracker):

    """Correlation between all sets.
    """
    pattern = "(.*)_readcounts$"
    mMinFoldChange = 2.0

    def __call__(self, track):
        data = []
        upfold = self.mMinFoldChange
        downfold = 1.0 / upfold
        data.append(("> %5.2f fold" % upfold, self.getValue(
            "SELECT COUNT(*) FROM %(track)s_readcounts WHERE CAST(fg_anysense_max AS float)/ bg_anysense_max > %(upfold)f " % locals())))
        data.append(("unchanged", self.getValue(
            "SELECT COUNT(*) FROM %(track)s_readcounts WHERE CAST(fg_anysense_max AS float)/ bg_anysense_max between %(downfold)f and %(upfold)f" % locals())))
        data.append(("< %5.2f fold" % downfold, self.getValue(
            "SELECT COUNT(*) FROM %(track)s_readcounts WHERE CAST(fg_anysense_max AS float)/ bg_anysense_max < %(downfold)f " % locals())))

        return odict(data)

##########################################################################
##########################################################################
##########################################################################
# peak shape
##########################################################################


class PeakShapeTracker(Tracker):

    '''return peakshape data.

    Only 1000 rows are returned.
    '''

    tracks = [os.path.basename(x)[:-len(".peakshape.tsv.gz")]
              for x in glob.glob(os.path.join(DATADIR, "peakshapes.dir", "*.regions.peakshape.tsv.gz"))]
    slices = ["peak_height", "peak_width"]

    def __call__(self, track, slice=None):

        fn = os.path.join(
            DATADIR, "peakshapes.dir", "%(track)s.peakshape.tsv.gz.matrix_%(slice)s.gz" % locals())
        if not os.path.exists(fn):
            return

        matrix, rownames, colnames = IOTools.readMatrix(IOTools.openFile(fn))

        nrows = len(rownames)
        if nrows < 2:
            return

        if nrows > 1000:
            take = numpy.array(
                numpy.floor(numpy.arange(0, nrows, nrows / 1000)), dtype=int)
            rownames = [rownames[x] for x in take]
            matrix = matrix[take]

        return odict((('matrix', matrix),
                      ('rows', rownames),
                      ('columns', colnames)))


class PeakShapeSummary(Tracker):

    '''summary information about peak shapes.'''
    pattern = "(.*)_peakshape"

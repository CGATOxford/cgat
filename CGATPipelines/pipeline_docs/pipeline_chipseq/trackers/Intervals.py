from ChipseqReport import *


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
# Annotation of bases with SNPs
##########################################################################


class IntervalsSummary(DefaultTracker):

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


class IntervalLengths(DefaultTracker):

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


class IntervalPeakValues(DefaultTracker):

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


class EnrichmentOverUnstim(DefaultTracker):

    """For each peakval, present the fold enrichment of a track
    compared to the Unstim set.

    Useful for diagnosing cutoffs.
    """

    mPattern = "_intervals$"

    def __call__(self, track, slice=None):

        cellline, condition, replicate = splitTrack(track)
        if condition == "Unstim":
            return []

        if replicate == None:
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


class IntervalAverageValues(DefaultTracker):

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


class IntervalLengthVsAverageValue(DefaultTracker):

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


class IntervalLengthVsPeakValue(DefaultTracker):

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


class PeakLocation(DefaultTracker):
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


class PeakDistance(DefaultTracker):
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


class Correlations(ChipseqTracker):

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
              for x in glob.glob(os.path.join(DATADIR, "*.peakshape.tsv.gz"))]
    slices = ["peak_height", "peak_width"]

    def __call__(self, track, slice=None):
        fn = os.path.join(
            DATADIR, "%(track)s.peakshape.tsv.gz.matrix_%(slice)s.gz" % locals())
        if not os.path.exists(fn):
            return

        matrix, rownames, colnames = IOTools.readMatrix(IOTools.openFile(fn))
        nrows = len(rownames)
        if nrows == 0:
            return

        if nrows > 1000:
            take = numpy.array(
                numpy.floor(numpy.arange(0, nrows, nrows / 1000)), dtype=int)
            rownames = [rownames[x] for x in take]
            matrix = matrix[take]

        return odict((('matrix', matrix),
                      ('rows', rownames),
                      ('columns', colnames)))

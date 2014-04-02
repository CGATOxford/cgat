from ChipseqReport import *
import Motifs

##########################################################################
##########################################################################
##########################################################################
# Base class for mast analysis
##########################################################################


class Nubiscan(ChipseqReport.DefaultTracker):
    mPattern = "_nubiscan$"

    def getSlices(self, subset=None):
        if subset is None:
            return "rxrvdr", "nr"
        elif "with-control" in subset:
            return "rxrvdr", "nr", "rxrvdr_shuffled", "nr_shuffled", "rxrvdr_shifted", "nr_shifted"
        else:
            return subset

################################################
################################################
################################################
##
################################################


class NubiscanSummary(Nubiscan):

    """return summary of nubiscan results.
    """

    def __call__(self, track, slice=None):
        data = []
        nintervals = self.getValue(
            "SELECT COUNT(*) FROM %(track)s_intervals" % locals())
        data.append(("ninput", nintervals))
        data.append(("nmotifs",
                     self.getValue("SELECT COUNT(*) FROM %(track)s_nubiscan WHERE motif = '%(slice)s'" % locals())))
        n = self.getValue(
            "SELECT COUNT(DISTINCT id) FROM %(track)s_nubiscan WHERE motif = '%(slice)s'" % locals())
        data.append(("nintervals", "%i" % n))
        data.append(("nintervals / %", "%5.2f" % (100.0 * n / nintervals)))

        nmast = self.getValue( """
        SELECT COUNT(distinct m.id) FROM %(track)s_mast AS m WHERE m.nmatches > 0 AND m.motif = '%(slice)s'""" % locals())
        data.append(("nmast", nmast))
        nshared = self.getValue( """SELECT COUNT(DISTINCT m.id) FROM
                                           %(track)s_mast AS m,
                                           %(track)s_nubiscan AS n
                                   WHERE n.id = m.id AND m.motif = '%(slice)s' AND n.motif = '%(slice)s' AND nmatches > 0""" % locals() )

        data.append(("nshared", nshared))

        return odict(data)

################################################
################################################
################################################
##
################################################


class NubiscanArrangements(Nubiscan):

    """return counts of repeat arrangenments.
    """

    def __call__(self, track, slice=None):
        data = []
        data = self.get(
            "SELECT arrangement, COUNT(*) FROM %(track)s_nubiscan WHERE motif = '%(slice)s' GROUP BY arrangement" % locals())
        data.sort(key=lambda x: (x[0][:2], int(x[0][2:])))
        return odict(data)

################################################
################################################
################################################
##
################################################


class NubiscanArrangementsPerInterval(Nubiscan):

    """return counts of repeat arrangenments.
    """

    def __call__(self, track, slice=None):
        data = []
        data = self.get(
            "SELECT arrangement, COUNT(DISTINCT id) FROM %(track)s_nubiscan WHERE motif = '%(slice)s' GROUP BY arrangement" % locals())
        data.sort(key=lambda x: (x[0][:2], int(x[0][2:])))
        return odict(data)

################################################
################################################
################################################
##
################################################


class NubiscanArrangementsPerIntervalWeighted(Nubiscan):

    """return counts of repeat arrangenments.

    Use weighted counts in the case of multiple possible arrangements.
    """

    def __call__(self, track, slice=None):
        result = collections.defaultdict(float)

        data = self.get(
            "SELECT arrangement, alternatives FROM %(track)s_nubiscan WHERE motif = '%(slice)s' ORDER by arrangement" % locals())
        for arr, alt in data:
            all = [arr]
            if alt:
                all.extend(alt.split(","))
            v = 1.0 / len(all)
            for x in all:
                result[x] += v

        return result

################################################
################################################
################################################
##
################################################


class NubiscanPeakValWithMotif(Nubiscan):

    '''return for each peakval the proportion of intervals
    that have a motif.'''

    def __call__(self, track, slice=None):
        statement = '''
        SELECT distinct i.peakval, i.interval_id, m.motif IS NOT NULL
        FROM %(track)s_intervals AS i
            LEFT JOIN %(track)s_nubiscan AS m ON m.id = i.interval_id AND m.motif = '%(slice)s'
            ORDER BY i.peakval DESC''' % locals()

        data = self.get(statement)

        result = Stats.getSensitivityRecall(
            [(int(x[0]), x[2] > 0) for x in data])
        return odict(zip(("peakval", "proportion with motif", "recall"), zip(*result)))

################################################
################################################
################################################
##
################################################


class NubiscanMastPeakValWithMotif(Nubiscan):

    '''return for each peakval the proportion of intervals
    that have a motif.

    return tracks for both nubiscan and mast
    '''

    def __call__(self, track, slice=None):

        result_nubiscan = NubiscanPeakValWithMotif()(track, slice)
        result_mast = Motifs.MastPeakValWithMotif()(track, slice)

        return odict((("mast", result_mast),
                      ("nubiscan", result_nubiscan)))

################################################
################################################
################################################
##
################################################


class NubiscanROC(Nubiscan):

    '''return a ROC curve. The ROC tests various peak parameters
    whether they are good descriptors of a motif.

    True/false positives are identified by the presence/absence
    of a motif. The presence is tested using the E-value of
    a motif.
    '''

    mPattern = "_nubiscan$"
    mFields = ("peakval", "avgval", "length")

    def __call__(self, track, slice=None):
        data = []

        # retrieve values of interest together with e-value
        rocs = []
        for field in self.mFields:

            statement = """
            SELECT DISTINCT i.interval_id, i.%(field)s, n.motif IS NOT NULL
            FROM %(track)s_intervals as i
            LEFT JOIN %(track)s_nubiscan as n
            ON i.interval_id = n.id AND n.motif = '%(slice)s'
            ORDER BY i.%(field)s DESC""" % locals()

            values = self.get(statement)
            try:
                result = Stats.computeROC([(x[1], x[2] > 0) for x in values])
            except ValueError, msg:
                # ignore results where there are no positives among values.
                continue

            rocs.append(result)

        d = Histogram.Combine(rocs)

        bins = [x[0] for x in d]
        values = zip(*[x[1] for x in d])

        result = odict()
        for f, v in zip(self.mFields, values):
            result[f] = odict((("FPR", bins), (f, v)))
        return result


################################################
################################################
################################################
##
################################################
class NubiscanAUC(NubiscanROC):

    '''return AUC for a ROC curve. The ROC tests various peak parameters
    whether they are good descriptors of a motif.

    True/false positives are identified by the presence/absence
    of a motif.
    '''

    mPattern = "_nubiscan$"
    mFields = ("peakval", "avgval", "length")

    def __call__(self, track, slice=None):

        data = NubiscanROC.__call__(self, track, slice)
        for k, d in data.iteritems():
            data[k] = Stats.getAreaUnderCurve(d['FPR'], d[k])
        return data


################################################
################################################
################################################
##
################################################
class NubiscanMotifLocation(Nubiscan):

    '''plot median position of motifs versus the peak location.'''

    def __call__(self, track, slice=None):

        data = self.getValues( """SELECT ((i.peakcenter - i.start) / 2 - (m.start + (m.end - m.start) / 2)) / ((CAST(i.length AS FLOAT) - (m.end - m.start)) / 2)
                                 FROM %(track)s_nubiscan as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'"""
                               % locals())

        return odict((("distance", data),))

################################################
################################################
################################################
##
################################################


class NubiscanMotifLocationMiddle(Nubiscan):

    '''plot median position of motifs versus the center of the interval.'''

    def __call__(self, track, slice=None):

        data = self.getValues( """SELECT ((200/2) - (m.start + (m.end - m.start) / 2)) / ((CAST(i.length AS FLOAT) - (m.end - m.start))/2)
                                 FROM %(track)s_nubiscan as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'"""
                               % locals())
        return odict((("distance", data),))

################################################
################################################
################################################
##
################################################


class NubiscanControlLocationMiddle(Nubiscan):

    '''plot median position of controls versus the center of the interval.'''

    def __call__(self, track, slice=None):

        data1 = self.getValues( """SELECT ( (m.r_length / 2) - (m.r_start + (m.r_end - m.r_start) / 2) ) / ((CAST( m.r_length as float) - (m.r_end - m.r_start))/2)
                                 FROM %(track)s_nubiscan as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'"""
                                % locals())

        data2 = self.getValues( """SELECT ( (m.l_length / 2) - (m.l_start + (m.l_end - m.l_start) / 2) ) / ((CAST( m.l_length as float) - (m.l_end - m.l_start))/2)
                                 FROM %(track)s_nubiscan as m, %(track)s_intervals as i 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s'"""
                                % locals())

        return odict((("distance", data1 + data2),))

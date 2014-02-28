from ChipseqReport import *

import Motifs
import TSS

###########################################################################
###########################################################################
###########################################################################


class MotifsAndTSS(Motifs.Mast):

    '''number of motifs matching within intervals.'''
    mPattern = "_tss$"

    def __call__(self, track, slice=None):

        statement =  """SELECT COUNT(i.interval_id)
                                 FROM %(track)s_mast as m, %(track)s_intervals as i, %(track)s_tss as d
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s' AND
                                       i.interval_id = d.gene_id AND %%s""" % locals()

        data = []
        data.append(("overlapping with motif", self.getValue(
            statement % "d.is_overlap > 0 AND m.nmatches > 0")))
        data.append(("overlapping without motif", self.getValue(
            statement % "d.is_overlap > 0 AND m.nmatches = 0")))
        data.append(("nonoverlapping with motif", self.getValue(
            statement % "d.is_overlap = 0 AND m.nmatches > 0")))
        data.append(("nonoverlapping without motif", self.getValue(
            statement % "d.is_overlap = 0 AND m.nmatches = 0")))
        return odict(data)

###########################################################################
###########################################################################
###########################################################################


class MotifsOverlappingTSS(Motifs.Mast):

    '''number of motifs matching within intervals.'''
    mPattern = "_tss$"

    def __call__(self, track, slice=None):

        statement =  """SELECT COUNT(i.interval_id)
                                 FROM %(track)s_mast as m, %(track)s_intervals as i, %(track)s_tss as d
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s' AND
                                       i.interval_id = d.gene_id AND %%s""" % locals()

        data = []
        data.append(
            ("with motif", self.getValue(statement % "d.is_overlap > 0 AND m.nmatches > 0")))
        data.append(
            ("without motif", self.getValue(statement % "d.is_overlap > 0 AND m.nmatches = 0")))
        return odict(data)

###########################################################################
###########################################################################
###########################################################################


class MotifsNonOverlappingTSS(Motifs.Mast):

    '''number of motifs matching within intervals.'''
    mPattern = "_tss$"

    def __call__(self, track, slice=None):

        statement =  """SELECT COUNT(i.interval_id)
                                 FROM %(track)s_mast as m, %(track)s_intervals as i, %(track)s_tss as d
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s' AND
                                       i.interval_id = d.gene_id AND %%s""" % locals()

        data = []
        data.append(
            ("with motif", self.getValue(statement % "d.is_overlap = 0 AND m.nmatches > 0")))
        data.append(
            ("without motif", self.getValue(statement % "d.is_overlap = 0 AND m.nmatches = 0")))
        return odict(data)

###########################################################################
###########################################################################
###########################################################################


class MastEValueVersusPeakValueAndDistance(Motifs.Mast):

    '''three way correlation

    between evalue, peak value and distance to TSS
    .'''

    def __call__(self, track, slice=None):

        field = "peakval"
        statement =  """SELECT i.%(field)s, d.closest_dist, m.evalue
                                 FROM %(track)s_mast as m, %(track)s_intervals as i, %(track)s_tss as d
                                 WHERE motif = '%(slice)s' AND 
                                       m.id = i.interval_id AND 
                                       m.id = d.gene_id AND 
                                       d.is_overlap = 0 
                                 ORDER BY i.%(field)s DESC""" % locals()

        data = [(math.log(x[1]), math.log(x[2]), math.log(x[0]))
                for x in self.get(statement % locals()) if x[0] > 00 and x[1] > 0 and x[2] > 0]

        return odict(zip(("log(distance)", "log(evalue)", "log(%s)" % field), zip(*data)))


###########################################################################
###########################################################################
###########################################################################
class MastROC2(Motifs.Mast):

    '''return a ROC curve. The ROC tests various peak parameters
    whether they are good descriptors of a motif.

    True/false positives are identified by the presence/absence
    of a motif.
    '''

    mPattern = "_tss$"
    mFields = ()

    def __call__(self, track, slice=None):
        data = []

        # obtain evalue distribution
        evalues = self.getValues(
            "SELECT evalue FROM %(track)s_mast WHERE motif = '%(slice)s'" % locals())
        bin_edges, with_motifs, explained = Motifs.computeMastCurve(evalues)

        # determine the e-value cutoff as the maximum of "explained"
        cutoff = bin_edges[numpy.argmax(explained)]

        # retrieve values of interest together with e-value
        rocs = []
        values = self.get( """SELECT d.closest_dist, m.evalue 
                                 FROM %(track)s_mast as m, %(track)s_intervals as i, %(track)s_tss AS d 
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s' AND d.gene_id = m.id
                                 AND d.closest_dist > 0
                                 ORDER BY d.closest_dist"""
                           % locals())

        rocs.append(Stats.computeROC([(x[0], x[1] <= cutoff) for x in values]))

        d = Histogram.Combine(rocs)

        bins = [x[0] for x in d]
        values = zip(*[x[1] for x in d])

        result = odict()
        for f, v in zip(self.mFields + ("dist",), values):
            result[f] = odict((("FPR", bins), (f, v)))
        return result

import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *

##################################################################################
class replicatedIntervalsPerContig( cpgTracker ):
    """Summary stats of intervals called by the peak finder. """

    mPattern = "_replicated_intervals$"

    def __call__(self, track, slice = None):
        data = self.getAll( """SELECT i.Contig, g.length as Contig_length, m.mappable_bases, B.repeat_length, COUNT(i.interval_id) as Intervals, A.Predicted_CGIs,  
                              round(COUNT(i.interval_id)/(m.mappable_bases/1000000.0),2) as CAPseq_density,
                              round(AVG(i.length),0) as Mean_length, round(AVG(i.nprobes),0) as Mean_reads 
                              FROM %(track)s_replicated_intervals i, annotations.genome g, annotations.mappable_bases_per_contig m,
                              (select contig, COUNT(id) as Predicted_CGIs from cgi_intervals group by contig) as  A,
                              (select contig, sum(stop-start) as repeat_length from annotations.repeats group by contig) B
                              WHERE i.contig=g.id AND A.contig=i.contig
                              AND B.contig=i.contig AND m.contig=i.contig
                              GROUP BY i.contig ORDER BY g.length desc LIMIT 100;"""  )

        headers = ("Contig length", "CAPseq Intervals", "Predicted CGIs", "CAPseq Density", "Mean Interval Length", "Mean Interval Reads")
        n = odict()
        for d in data:
            contig  = d[:1]
            n[str(contig)] = odict( zip(headers,  d[2:]))

        #result =  zip(headers, zip(*data)) 
        return data



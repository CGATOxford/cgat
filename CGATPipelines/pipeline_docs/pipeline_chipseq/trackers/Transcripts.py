import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy, scipy.stats
import numpy.ma
import Stats
import Histogram
import ChipseqReport


from SphinxReport.Tracker import *

# for trackers_derived_sets and trackers_master
if not os.path.exists("conf.py"):
    raise IOError( "could not find conf.py" )

execfile( "conf.py" )


class CountsTranscripts (TrackerSQL ):
    mPattern = "transcripts"
    mAsTables = True

    def getSlices(self, subset = None ):
        return ("source", "feature", "contig")

    def __call__(self,track, slice = None ):

        
        data = self.get( '''SELECT %(slice)s, count(DISTINCT gene_id), count(DISTINCT transcript_id) 
                                   FROM %(track)s 
                                   GROUP BY %(slice)s''' % locals() )
            
        result = odict()
        for x, genes, transcripts in data:
            result[x] = odict( (('genes', genes),('transcripts',transcripts)))
        return result

class CountsPromotors( CountsTranscripts ):
    mPattern = "promotors"

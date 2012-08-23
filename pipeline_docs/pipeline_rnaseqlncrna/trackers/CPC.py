from SphinxReport.Tracker import *

import IOTools
import GTF
import numpy as np
import scipy.stats
import collections

#################################################
#################################################
#################################################

class CPC(TrackerSQL):

    pattern = "(.*)_result"
    slices = ["total","coding","noncoding"]

    def __call__(self, track, slice = None):
        
        if slice == "total":
            statement = "SELECT COUNT(*) FROM %(track)s_result"
            return self.getValue(statement)
        elif slice == "coding":
            statement = "SELECT COUNT(*) FROM %(track)s_result WHERE C_NC='%(slice)s'"
            return self.getValue(statement)
        elif slice == "noncoding":
            statement = "SELECT COUNT(*) FROM %(track)s_result WHERE C_NC='%(slice)s'"
            return self.getValue(statement)

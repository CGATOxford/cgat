
import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
import cpgReport
import Annotations
import Motifs

from SphinxReport.Tracker import *

##################################################################################
class RepeatOverlap(Annotations.AnnotationsAssociated):
    """Overlap with repeats."""
    mPattern = "_repeats$"
    mColumns = "SUM(CASE WHEN nover>0 THEN 1 ELSE 0 END) as with, SUM(CASE WHEN nover=0 THEN 1 ELSE 0 END) AS without" 
    mTable = "repeats"
    
    def __call__(self, track, slice = None ):
        statement = self.getStatement( track, slice )
        if not statement: return []
        return odict( zip( ("with","without"), self.getFirstRow( statement) ))


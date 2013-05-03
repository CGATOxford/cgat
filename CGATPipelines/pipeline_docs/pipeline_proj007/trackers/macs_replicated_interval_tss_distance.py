import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram
from cpgReport import *
from SphinxReport.Tracker import *
from SphinxReport.odict import OrderedDict as odict

##################################################################################
##################################################################################
##################################################################################
class transcriptTSSOverlap(cpgTracker):
    '''number of transcript TSSs that an interval overlaps.'''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_"+ANNOTATIONS_NAME+"_transcript_tss_distance$"
    mAnnotations = "_annotations"
    mTable = "replicated_"+ANNOTATIONS_NAME+"_transcript_tss_distance"
    mColumn = "d.is_overlap"
    mWhere = "d.is_overlap < 5 "
    def __call__(self, track, slice = None ):

        annotations = self.mAnnotations
        table = self.mTable
        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            data = self.getValues( """SELECT %(column)s FROM %(track)s_%(table)s AS d WHERE %(where)s""" % locals() )
        else:
            data = self.getValues( """SELECT %(column)s FROM %(track)s_%(table)s AS d, %(track)s_%(annotations)s as a 
                                      WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )

        hist, bins = numpy.histogram( data, bins=numpy.arange(0, max(data) + 1, 1) )
        return odict( zip( map(str, bins[:-1]), hist) )

##################################################################################
class TSSClosest(cpgTracker):
    """for each interval, return the distance to the closest TSS."""
    ANNOTATIONS_NAME = P['annotations_name']
    mXLabel = "distance / bases"
    mPattern = "_replicated_"+ANNOTATIONS_NAME+"_transcript_tss_distance$"
    mColumn = "d.closest_dist"
    mWhere = "1"
    mAnnotations = "_annotations"
    mTable = "replicated_"+ANNOTATIONS_NAME+"_transcript_tss_distance"
    def __call__(self, track, slice = None ):

        annotations = self.mAnnotations
        table = self.mTable
        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            statement = """SELECT %(column)s FROM %(track)s_%(table)s AS d WHERE %(where)s""" % locals() 
        else:
            statement = """SELECT %(column)s FROM %(track)s_%(table)s AS d, %(track)s_%(annotations)s as a 
                           WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals()
        #print statement
        data = self.get( statement )
        return data

##################################################################################
class transcriptTSSClosestUpstream(TSSClosest):
    """for each interval, return peakval and the distance to the closest upstream TSS."""
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist5 is null THEN 1000000 ELSE dist5 END as dist5 "

##################################################################################
class transcriptTSSClosestDownstream(TSSClosest):
    """for each interval, return peakval and the distance to the closest downstream TSS."""
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist3 is null THEN 1000000 ELSE dist3 END as dist3 "

##################################################################################
##################################################################################
##################################################################################
class geneTSSOverlap(cpgTracker):
    '''number of gene TSSs that an interval overlaps.'''
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_"+ANNOTATIONS_NAME+"_gene_tss_distance$"
    mAnnotations = "_annotations"
    mTable = "replicated_"+ANNOTATIONS_NAME+"_gene_tss_distance"
    mColumn = "d.is_overlap"
    mWhere = "d.is_overlap < 5 "
    def __call__(self, track, slice = None ):

        annotations = self.mAnnotations
        table = self.mTable
        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            data = self.getValues( """SELECT %(column)s FROM %(track)s_%(table)s AS d WHERE %(where)s""" % locals() )
        else:
            data = self.getValues( """SELECT %(column)s FROM %(track)s_%(table)s AS d, %(track)s_%(annotations)s as a 
                                      WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )

        hist, bins = numpy.histogram( data, bins=numpy.arange(0, max(data) + 1, 1) )
        return odict( zip( map(str, bins[:-1]), hist) )

##################################################################################
class geneTSSClosest(TSSClosest):
    """for each interval, return the distance to the closest TSS."""
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_"+ANNOTATIONS_NAME+"_gene_tss_distance$"
    mColumn = "d.closest_dist"
    mTable = "replicated_"+ANNOTATIONS_NAME+"_gene_tss_distance"

##################################################################################
class geneTSSClosestUpstream(TSSClosest):
    """for each interval, return peakval and the distance to the closest upstream TSS."""
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_"+ANNOTATIONS_NAME+"_gene_tss_distance$"
    mTable = "replicated_"+ANNOTATIONS_NAME+"_gene_tss_distance"
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist5 is null THEN 1000000 ELSE dist5 END as dist5 "

##################################################################################
class geneTSSClosestDownstream(TSSClosest):
    """for each interval, return peakval and the distance to the closest downstream TSS."""
    ANNOTATIONS_NAME = P['annotations_name']
    mPattern = "_replicated_"+ANNOTATIONS_NAME+"_gene_tss_distance$"
    mTable = "replicated_"+ANNOTATIONS_NAME+"_gene_tss_distance"
    mColumn = "CASE WHEN is_overlap>0 THEN 0 WHEN dist3 is null THEN 1000000 ELSE dist3 END as dist3 "

##################################################################################
class transcriptTSSProfile(TrackerImages):
    """TSS and TTS profile per transcript """

##################################################################################
class transcriptTSSProfileCapseq(TrackerImages):
    """TSS and TTS profile per transcript """

##################################################################################
class transcriptTSSProfileNoCapseq(TrackerImages):
    """TSS and TTS profile per transcript """

##################################################################################
class geneTSSProfile(TrackerImages):
    """TSS and TTS profile per gene """
    
##################################################################################
class geneTSSProfileCapseq(TrackerImages):
    """TSS and TTS profile per gene """
    
##################################################################################
class geneTSSProfileNoCapseq(TrackerImages):
    """TSS and TTS profile per gene """
    
    

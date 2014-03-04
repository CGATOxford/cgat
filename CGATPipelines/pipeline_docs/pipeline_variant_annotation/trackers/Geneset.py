import os
import sys
import re
import types
from VariantsReport import *


class GenesetCounts(VariantsTracker):

    '''Count genes, transcripts and proteins in gene annotation set.'''
    tracks = ("annotations.transcript_info",)

    def getSlices(self, subset=None):
        return ("gene_id", "gene_name", "protein_id", "transcript_id", "transcript_name")

    def __call__(self, track, slice=None):
        return odict( (("counts", self.getValue( """SELECT COUNT(DISTINCT %(slice)s) FROM %(track)s""" % locals())),) )

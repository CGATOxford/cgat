import os
import sys
import re
import types
from VariantsReport import *

#####################################################
#####################################################
#####################################################


class SplicingCounts(VariantsTracker):

    '''number of transcripts per gene.'''
    mPattern = "_effects_splicing$"

    def __call__(self, track, slice=None):
        columns = ("nintrons",
                   "ncanonical",
                   "nunchanged",
                   "ndisrupted",
                   "nnonsynonymous",
                   "nnovel",
                   "nnunknown",
                   "nsynonymous",
                   "nframeshifts",
                   "nunchanged_frames",
                   "ncorrected_frames",
                   "nuncorrected_frames")

        result = odict()
        for column in columns:
            result[column] = self.getValue(
                '''SELECT SUM(%(column)s) FROM %(track)s_effects_splicing''' % locals() )
        return result

#####################################################
#####################################################
#####################################################


class FrameShiftCorrection(VariantsTracker):

    '''number of frameshift introns and how many have been corrected.'''
    mPattern = "_effects_splicing$"

    def __call__(self, track, slice=None):
        result = odict(zip(
            ("nframeshifts", "nunchanged", "ncorrected", "nuncorrected"),
            self.getFirstRow(
                '''SELECT SUM(nframeshifts), SUM(nunchanged_frames), SUM(ncorrected_frames), SUM(nuncorrected_frames)
            FROM %(track)s_effects_splicing''' % locals() )))
        return result

#####################################################
#####################################################
#####################################################


class FrameShiftCorrectedTranscripts(VariantsTracker):

    '''return the transcripts that have corrected frameshifts.'''
    mPattern = "_effects_splicing$"

    def __call__(self, track, slice=None):
        return odict(self.get(
            '''SELECT transcript_id, ncorrected_frames
            FROM %(track)s_effects_splicing WHERE ncorrected_frames > 0''' % locals()) )


#####################################################
#####################################################
#####################################################
class VariantSplicingTranscripts(VariantsTracker):

    '''return the transcripts that have abberent splicing due to sequence variants.'''
    mPattern = "^effects$"

    def __call__(self, track, slice=None):

        headers = ("gene_id", "gene_name", "transcript_id",
                   "track", "Genotype", "nvariants")

        statement = '''SELECT
            i.gene_id,
            i.gene_name,
            i.transcript_id,
            e.track,        
            e.splice_genotype,
            e.splice_nvariant_sites
        FROM
            effects e,
            annotations.transcript_info AS i
        WHERE i.transcript_id = e.transcript_id
        AND e.splice_genotype is not null
        ORDER BY i.gene_id''' % locals()
        # return odict( self.get(statement) )
        return odict(zip(headers, zip(*self.get(statement))))

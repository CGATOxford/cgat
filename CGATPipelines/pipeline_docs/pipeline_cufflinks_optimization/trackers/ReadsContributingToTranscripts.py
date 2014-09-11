from CGATReport.Tracker import *


class ReadsContributingToTranscriptsSummary(TrackerSQL):

    '''
    tracker for collecting data for the number of alignments
    and spliced alignments that contribute to transcripts
    '''

    pattern = "(.*)_summary"
    slices = ["alignments", "percent alignments"]

    def __call__(self, track, slice=None):

        if slice is None:
            return self.getRow("SELECT * FROM %(track)s_summary")
        elif slice == "alignments":
            return self.getRow("SELECT aligments_in_transcripts, spliced_alignments_in_transcripts, total_alignments, total_spliced_alignments FROM %(track)s_summary")
        elif slice == "percent alignments":
            return self.getRow("SELECT percent_alignments_in_transcripts, percent_spliced_alignments_in_transcripts FROM %(track)s_summary")

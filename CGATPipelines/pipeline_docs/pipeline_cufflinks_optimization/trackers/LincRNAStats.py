from CGATReport.Tracker import *


class SingleAndMultiExonCounts(TrackerSQL):

    '''
    tracker for retrieving the counts of single and multi exonic lincRNA
    '''
    pattern = "(.*)_count"

    def __call__(self, track, slice=None):

        return self.getRow("SELECT no_multi_exon_transcripts,no_single_exon_transcripts FROM %(track)s_count")


class ProportionSingleExon(TrackerSQL):

    '''
    tracker for retrieving the counts of single and multi exonic lincRNA
    '''
    pattern = "(.*)_count"

    def __call__(self, track, slice=None):

        return self.getRow("SELECT proportion_single FROM %(track)s_count")


class LincRNAExonSummary(TrackerSQL):

    '''
    tracker for collecting data for lincRNA summary
    '''

    pattern = "(.*)_stats"

    def __call__(self, track, slice=None):

        return self.getValue("SELECT avg(no_exons) FROM %(track)s_stats")


class LincRNALengthSummary(TrackerSQL):

    '''
    tracker for collecting data for lincRNA summary
    '''

    pattern = "(.*)_stats"

    def __call__(self, track, slice=None):

        return self.getValue("SELECT avg(transcriptlength) FROM %(track)s_stats")


class LincRNALengthDistribution(TrackerSQL):

    '''
    tracker for collecting data for lincRNA summary
    '''

    pattern = "(.*)_stats"

    def __call__(self, track, slice=None):

        return self.getValues("SELECT transcriptlength FROM %(track)s_stats")

from CGATReport.Tracker import *


class ReadSummary(TrackerSQL):

    '''
    Summarise read counts for each track
    '''

    def __call__(self, track, slice=None):

        return self.getAll("SELECT * FROM reads_summary")

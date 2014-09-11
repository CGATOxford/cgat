from CGATReport.Tracker import *
import sqlite3


class RankedStats(TrackerSQL):

    '''
    tracker fot pulling in the data for the ranked summaries
    '''
    dbh = sqlite3.connect("csvdb")
    cc = dbh.cursor()
    tracks = []
    for track in cc.execute("SELECT track FROM All_stats_combined").fetchall():
        tracks.append(track[0])

    def __call__(self, track, slice=None):

        return self.getAll("""SELECT multi_exon_count
                              single_exon_count, rank_single_exon_count
                              , rank_multi_exon_count
                              , no_exons,rank_no_exons
                              , proportion_reads_to_transcripts
                              , rank_proportion_reads_to_transcripts
                              , proportion_single
                              , transcript_length
                              , rank_transcript_length
                              , rank_sum
                              , rank_all FROM All_stats_combined WHERE track='%(track)s'""")


class PlotRankedStats(TrackerSQL):

    '''
    tracker fot pulling in the data for the ranked summaries
    '''
    dbh = sqlite3.connect("csvdb")
    cc = dbh.cursor()
    tracks = []
    for track in cc.execute("SELECT track FROM All_stats_combined").fetchall():
        tracks.append(track[0])

    def __call__(self, track, slice=None):

        return self.getRow("""SELECT multi_exon_count
                              , no_exons
                              , proportion_reads_to_transcripts
                              , proportion_single
                              , transcript_length
                              FROM All_stats_combined WHERE track='%(track)s'""")

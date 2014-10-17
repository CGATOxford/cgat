

from CGATReport.Tracker import *
import sqlite3


class ReadCountSummary (TrackerSQL):
    pattern = "reads_under_peaks"

    def getTracks(self, subset=None):
        return self.getValues("SELECT track FROM reads_under_peaks ")

    def getSlices(self, subset=None):
        fields = self.getColumns("reads_under_peaks")
        return [x for x in fields if x != "track"]

    def __call__(self, row, col):
        data = self.getValue(
            '''SELECT %(col)s FROM reads_under_peaks WHERE track = '%(row)s' ''' )
        return data


#-------------------------------------------------------------------------

def getNormalisationFactors():
    '''function for retrieving the normalising factors for reads_under_peaks. These are the 
    total number of mapped reads for each sample haveing reads counted i.e. bam files'''

    con = sqlite3.connect('csvdb')
    cursor = con.cursor()
    cursor.execute('''select reads_mapped, track from view_mapping''')
    data2 = cursor.fetchall()
    norms = {}
    for each in data2:
        norms[each[1]] = each[0]

    return norms


class NormalisedTable (ReadCountSummary):

    def __call__(self, row, col):

        norms = getNormalisationFactors()
        data = self.getValue(
            '''SELECT %(col)s FROM reads_under_peaks WHERE track = '%(row)s' ''') / float(norms[row])
        return data

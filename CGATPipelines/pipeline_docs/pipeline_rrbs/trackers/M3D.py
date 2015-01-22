
from SphinxReport.Tracker import *
from rrbsReport import *


class M3DSummaryTable(RrbsTracker):

    # pattern = "^summary_table$"

    def __call__(self, track, slice=None):

        statement = '''
        SELECT group1 as Group_1,  group2 as Group_2,
        total as Total, significant as Significant_Clusters
        FROM summary_table
        ''' % locals()
        return self.getAll(statement)

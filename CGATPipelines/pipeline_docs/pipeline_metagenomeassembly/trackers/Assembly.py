from CGATReport.Tracker import *
import sqlite3
import collections
import numpy as np


class Summary(TrackerSQL):

    '''
    class to provide a summary of the contig
    assembly
    '''
    pattern = "(.*_.*)_summary$"

    def __call__(self, track, slice=None):
        '''
        return sqlite table with summary statistics
        '''
        return self.getAll("""SELECT * FROM %(track)s_summary""")


class LengthDistribution(TrackerSQL):

    '''
    return length distirbution of assembled
    contigs
    '''
    pattern = "(.*)_lengths"

    def __call__(self, track, slice=None):
        x = self.getValues("""SELECT length FROM %(track)s_lengths""")
        return np.log10(x)

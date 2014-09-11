from CGATReport.Tracker import *
import sqlite3
import collections
import numpy as np
import glob
import os


class NumberOfGenes(TrackerSQL):

    '''
    class for collecting the number of 
    genes that were found using metagenemark
    '''

    # use the proteins
    pattern = "(.*)_aa"

    def __call__(self, track, slice=None):

        return self.getValue("""SELECT COUNT(*) FROM %(track)s_aa""")

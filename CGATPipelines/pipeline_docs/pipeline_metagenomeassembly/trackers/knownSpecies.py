from CGATReport.Tracker import *
import sqlite3
import collections
import numpy as np


class SpeciesCount(TrackerSQL):

    def __call__(self, track, slice=None):
        '''
        return the number of reference genomes that 
        are aligned to 
        '''
        genomes = self.execute(
            """SELECT count(*) FROM species_present_fa""").fetchone()[0]
        return {"total_reference_genomes": genomes}


class Species(TrackerSQL):

    def __call__(self, track, slice=None):
        '''
        return the number of reference genomes that 
        are aligned to 
        '''

        return self.getAll("""SELECT * FROM species_present_fa""")


class KnownAlignments(TrackerSQL):

    def __call__(self, track, slice=None):
        '''
        return picard stats results
        '''
        result = {}
        for data in self.execute("""SELECT track, PCT_PF_READS_ALIGNED FROM known_genomes_picard_stats_alignment_summary_metrics"""):
            result[data[0]] = data[1]
        return result

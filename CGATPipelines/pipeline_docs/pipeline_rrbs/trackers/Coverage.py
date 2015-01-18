from SphinxReport.Tracker import *
from rrbsReport import *


class seqStart(RrbsTracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM read_start_summary'''
        return self.getAll(statement)


class coverage(RrbsTracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM coverage'''
        return self.getAll(statement)


class readsRemaining(RrbsTracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM reads_remaining_by_threshold'''
        return self.getAll(statement)


class CpGOverlap(RrbsTracker):
    def __call__(self, track, slice=None):
        statement = '''SELECT * FROM coverage_overlap'''
        return self.getAll(statement)

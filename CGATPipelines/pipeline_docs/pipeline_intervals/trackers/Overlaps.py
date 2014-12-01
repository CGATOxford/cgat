from IntervalReport import *


class OverlapsPercentageExons(SingleTableTrackerEdgeListToMatrix,
                              IntervalTracker):
    table = "overlap"
    row = "set1"
    column = "set2"
    value = "pexons_ovl1"
    value2 = "pexons_ovl2"
    is_square = True


class OverlapsPercentageNucleotides(OverlapsPercentageExons):
    value = "pbases_ovl1"
    value2 = "pbases_ovl2"
    is_square = True

"""remove from a sorted list of links
those which are redundant.
"""

import sys, string, os

class Map:

    def __init__(self ):
        (self.mQueryToken, self.mSbjctToken, self.score,
         self.mQueryFrom, self.mQueryTo, self.mQueryAli,
         self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli ) =\
         ("", "", 0, 0, 0, "", 0, 0, "" )
        self.mIsExpanded = False
        self.mMapQuery2Sbjct = None

    def Read( self, line ):
        
        (self.mQueryToken, self.mSbjctToken, self.score,
         self.mQueryFrom, self.mQueryTo, self.mQueryAli,
         self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli) = line[:-1].split("\t")[:9]

        (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo) = map(
            int, (self.mQueryFrom, self.mQueryTo, self.mSbjctFrom, self.mSbjctTo))

        self.score = float(self.score)

    def __str__( self ):

        return string.join( map(str, (
            self.mQueryToken, self.mSbjctToken, self.score,
            self.mQueryFrom, self.mQueryTo, self.mQueryAli,
            self.mSbjctFrom, self.mSbjctTo, self.mSbjctAli)), "\t")

if __name__ == "__main__":

    last_m = None

    for line in sys.stdin:
        
        if line[0] == "#": continue

        m = Map()
        m.Read( line )

        if not last_m:
            last_m = m
            continue
        
        if m.mQueryToken != last_m.mQueryToken or \
           m.mSbjctToken != last_m.mSbjctToken or \
           m.mQueryFrom > last_m.mQueryTo:
            print str(last_m)
            last_m = m
        else:
            if m.mQueryTo - m.mQueryFrom > last_m.mQueryTo - last_m.mQueryFrom:
                last_m = m

    print str(last_m)
        

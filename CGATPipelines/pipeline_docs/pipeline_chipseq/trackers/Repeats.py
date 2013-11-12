from ChipseqReport import *
import Annotations
import Motifs

##################################################################################
##################################################################################
##################################################################################
## Looking at overlap with repeats
##################################################################################
class RepeatOverlap(Annotations.AnnotationsAssociated):
    """Overlap with repeats."""
    mPattern = "_repeats$"
    mColumns = "SUM(CASE WHEN nover>0 THEN 1 ELSE 0 END) as with, SUM(CASE WHEN nover=0 THEN 1 ELSE 0 END) AS without" 
    mTable = "repeats"
    
    def __call__(self, track, slice = None ):
        statement = self.getStatement( track, slice )
        if not statement: return []
        return odict( zip( ("with","without"), self.getFirstRow( statement) ))

###########################################################################
###########################################################################
###########################################################################
class RepeatsMastEValueVersusPeakValueAndDistance( Motifs.Mast ):
    '''three way correlation.'''
    
    mPattern = "_repeats$"
    def __call__(self, track, slice = None ):

        field = "peakval"
        statement =  """SELECT m.evalue, i.%(field)s, r.pover1 
                                 FROM %(track)s_mast as m, %(track)s_intervals as i, %(track)s_repeats as r
                                 WHERE i.interval_id = m.id AND motif = '%(slice)s' AND
                                       i.interval_id = r.gene_id 
                                 ORDER BY i.%(field)s DESC""" % locals()

        data = [ (x[2], x[1], math.log(x[0])  ) for x in self.get( statement % locals() ) ]

        return odict( zip( ("evalue", field, "log(distance)"), zip(*data) ) )

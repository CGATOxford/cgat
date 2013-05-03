## tracks taken from geneset comparison
## need to integrated

from MappingReport import *

class AnnotationsAssociated:
    pass

##################################################################################
##################################################################################
##################################################################################
## Coverage of transcript models
##################################################################################
class ContaminationCoverage( AnnotationsAssociated ):
    """Check for contamination by listing transcript models with reads ."""
    pattern = "(.*)_coverage$"
    mColumns = "count(*) as total, SUM(CASE WHEN nmatches > 1 THEN 1 ELSE 0 END) AS hico" 
    mTable = "coverage"

    def __call__(self, track, slice = None ):
        statement = self.getStatement( track, slice )
        if not statement: return []
        data = self.getFirstRow( statement) 
        return odict( zip( ("nmatches > 1", "nmatches = 1" ), (data[1], data[0]-data[1]) ))

##################################################################################
##################################################################################
##################################################################################
## Coverage of transcript models
##################################################################################
class PolyATailCounts( AnnotationsAssociated ):
    """Check for contamination by listing transcript models with reads ."""
    pattern = "(.*)_polyA$"
    mColumns = "COUNT(*) AS total, SUM(nmotifs) AS motifs, SUM(tails) AS tails" 
    mTable = "polyA"

    def __call__(self, track, slice = None ):
        statement = self.getStatement( track, slice )
        if not statement: return []
        data = self.getFirstRow( statement) 
        return odict(zip( ("no tail", "with motif", "without motif" ), (data[0] - data[2], data[1], data[2] - data[1]) ) )


##################################################################################
##################################################################################
##################################################################################
## Contamination and repeats
##################################################################################
class ContaminationRepeats( TrackerSQL ):
    """Estimate contamination based on the overlap with repeats.

    repeats
       number of bases in repeats
    genome
       number of base is genome
    prepeats
       proportion of bases in repeats
    repeat_overlap
       number of bases in unknown transcript models overlapping repeats
    length
       number of bases in unknown transcript models
    poverlap
       proportion of bases in unknown transcript models overlapping repeats
    nspliced_ovl
       number of unknown transcript models with introns that overlap repeats
    nspliced_ovl
       number of unknown transcript models with introns
    pspliced
       proportion of unknown transcript models with introns that overlap repeats
    """

    pattern = "(.*)_repeats$"

    def getTracks( self, subset = None ):
        return [ x for x in TrackerSQL.getTracks( self, subset ) if "_vs" not in x]

    def __call__(self, track, slice = None ):
        genome_size = self.getValue( "SELECT SUM(length) FROM repeats_table" )
        repeats_size = self.getValue( "SELECT SUM(nover_bases) FROM repeats_table" )
        
        novl_repeats = self.getValue( "SELECT SUM(nover) FROM %(track)s_repeats as r, %(track)s_annotation as a where a.gene_id = r.gene_id and is_unknown" % locals())
        nlength = self.getValue( "SELECT SUM(exons_sum) FROM %(track)s_repeats as r, %(track)s_annotation as a where a.gene_id = r.gene_id and is_unknown" % locals())
        nspliced = self.getValue( "SELECT COUNT(*) FROM %(track)s_repeats as r, %(track)s_annotation as a where a.gene_id = r.gene_id and is_unknown AND exons_nval > 1" % locals())
        nspliced_ovl_repeats = self.getValue( "SELECT COUNT(*) FROM %(track)s_repeats as r, %(track)s_annotation as a where a.gene_id = r.gene_id and is_unknown AND exons_nval > 1 AND nover > 0" % locals())

        if novl_repeats == None: 
            novl_repeats, nlength, nspliced, nspliced_ovl_repeats = 0, 0, 0, 0

        return odict( ( ("repeats", repeats_size),
                        ("genome", genome_size),
                        ("prepeats", prettyPercent( repeats_size, genome_size)),
                        ("repeat_overlap", novl_repeats),
                        ("length", nlength),
                        ("poverlap", prettyPercent( novl_repeats, nlength ) ),
                        ("pcont", prettyPercent( novl_repeats * genome_size, nlength * repeats_size)),
                        ("nspliced_ovl", nspliced_ovl_repeats ),
                        ("nspliced", nspliced ),
                        ("pspliced_ovl", prettyPercent( nspliced_ovl_repeats, nspliced) ),
                      ) )


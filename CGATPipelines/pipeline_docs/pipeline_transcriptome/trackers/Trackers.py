import os
import sys
import re
import types
import collections
import numpy
import numpy.ma

import sqlalchemy
from SphinxReport.Tracker import *

from AnnotationReport import *

# =======================================================
# =======================================================
# =======================================================
# administrative trackers
# =======================================================


class TrackerSQLCheckTablesGeneId(TrackerSQLCheckTables):

    """Tracker that examines the presence/absence of genes in tables.
    """

    mFields = ["gene_id", ]
    mExcludePattern = "_vs_"
    mPattern = "_annotation$"


class TrackerSQLCheckTableEvol(TrackerSQLCheckTable):

    """Tracker that counts existing entries in the *evol* table.
    """
    mExcludePattern = None
    mPattern = "_evol$"

# =======================================================
# =======================================================
# =======================================================
# data trackers
# =======================================================


class DefaultSlicer:

    """Returns the default slices.

    The option *subset* returns a group of slices. The available
    groups are
       * default: only default slices ("all", "known", "unknown")
       * complete: all slices 
       * derived: only derived slices %s, see conf.py
       * "track": only the derived slices for a "track", see conf.py
       * mixture: combination of default slices and derived slices
       * set.track: combine one of the default slices with all slices of track "track"
       """

    def getSlices(self, subset="default"):

        if type(subset) in (types.ListType, types.TupleType):
            if len(subset) > 1:
                return subset
            else:
                subset = subset[0]

        all_slices = trackers_default_slices + trackers_derived_slices.keys()
        if subset is None or subset == "default":
            return trackers_default_slices
        elif subset == "complete":
            return all_slices
        elif subset == "derived":
            return trackers_derived_slices.keys()
        elif subset == "mixture":
            slices = []
            for x in trackers_default_slices:
                for y in trackers_derived_slices.keys():
                    slices.append("%s.%s" % (x, y))
            return slices
        elif "." in subset:
            slices = []
            set1, track = subset.split(".")
            for y in [x[0] for x in trackers_derived_slices.items() if x[1] == track]:
                slices.append("%s.%s" % (set1, y))
            return slices
        elif subset in all_slices:
            # a single slice
            return [subset, ]
        else:
            # slices by track
            return [x[0] for x in trackers_derived_slices.items() if x[1] == subset]


class AnnotationsAssociated(DefaultSlicer, TrackerSQL):

    """simple join between a data table and table defining slices.

    :attr:`mTable`
       table to join with
    :attr:`mColums`
       columns to output
    """
    pattern = "(.*)_annotation"
    mTable = None
    mColumns = None
    mWhere = "1"
    mSelectAll = "SELECT %(columns)s FROM %(track)s_%(table)s AS t WHERE %(where)s"
    mSelectSubset = "SELECT %(columns)s FROM %(track)s_%(table)s AS t, %(track)s_annotation AS a WHERE a.gene_id = t.gene_id AND a.is_%(slice)s AND %(where)s"
    mSelectSlice = "SELECT %(columns)s FROM %(track)s_%(table)s AS t, %(track)s_%(slice)s AS s WHERE s.gene_id = t.gene_id AND %(where)s"
    mSelectMixture = "SELECT %(columns)s FROM %(track)s_%(table)s AS t, %(subset)s AS s, %(track)s_annotation AS a WHERE a.gene_id = t.gene_id AND a.is_%(slice)s AND s.gene_id = t.gene_id AND %(where)s"

    def getStatement(self, track, slice=None):
        columns = self.mColumns
        table = self.mTable
        where = self.mWhere
        if not table or not columns:
            raise NotImplementedError
        if slice is not None and "." in slice:
            slice, subset = slice.split(".")
            if subset in trackers_derived_slices and track != trackers_derived_slices[subset]:
                return None
            else:
                return self.mSelectMixture % locals()
        elif slice in trackers_derived_slices:
            if track == trackers_derived_slices[slice]:
                return self.mSelectSlice % locals()
            else:
                return None
        elif slice == "all" or slice is None:
            return self.mSelectAll % locals()
        else:
            return self.mSelectSubset % locals()

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################
# note: the order is important!


class TranscriptLengths(DefaultSlicer, TrackerSQL):

    """Lengths of transcript models."""
    pattern = "(.*)_annotation$"

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def __call__(self, track, slice=None):

        table = sqlalchemy.Table(
            "%s_annotation" % track, self.metadata, autoload=True)
        if slice is None or slice == "all":
            data = [x[0]
                    for x in self.execute(sqlalchemy.select(columns=[table.c.exons_sum, ]))]
        elif slice == "linc":
            data = [x[0] for x in self.execute(
                "SELECT exons_sum FROM %s_annotation AS a, %s_%s AS s WHERE s.gene_id = a.gene_id" % (track, track, slice))]
        else:
            data = [x[0] for x in self.execute(
                "SELECT exons_sum FROM %s_annotation WHERE is_%s" % (track, slice))]
        return odict((("length", data),))

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class Annotations(DefaultSlicer, TrackerSQL):

    """Base class for trackers getting info from the annotations tables.

    Derived Trackers should define the two attributes :attr:`mSelect` and 
    :attr:`mColumns`.
    """
    pattern = "(.*)_annotation$"
    mSelect = None
    mColumns = None
    mWhere = "1"

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def __call__(self, track, slice=None):

        where = self.mWhere
        select = self.mSelect

        if slice in trackers_derived_slices:
            if track == trackers_derived_slices[slice]:
                data = self.getFirstRow( """%(select)s
                        FROM %(track)s_annotation AS a, %(track)s_%(slice)s AS s WHERE %(where)s AND a.gene_id = s.gene_id""" % locals() )
            else:
                return []
        elif slice == "all" or slice is None:
            data = self.getFirstRow(
                "%(select)s FROM %(track)s_annotation WHERE %(where)s" % locals())
        else:
            data = self.getFirstRow(
                "%(select)s FROM %(track)s_annotation WHERE %(where)s AND is_%(slice)s" % locals())

        return odict(zip(self.mColumns, data))


class AllAnnotations(Annotations):

    """Annotations of all transcript models."""
    mColumns = ["total", "known", "unknown", "ambiguous",
                "pc", "pseudo", "npc", "utr",
                "intronic", "associated", "intergenic",
                "unclassified",
                "unclassified_known",
                "unclassified_unknown"]

    mSelect = """SELECT count(*) AS total, 
			sum(is_known) AS known, 
                        sum(is_unknown) AS unknown, 
                        sum(is_ambiguous) AS ambiguous, 
			sum( is_pc) AS pc,
                        sum( is_pseudo) AS pseudo, 
                        sum( is_npc) AS npc, 
                        sum(is_utr) AS utr, 
			sum(is_intronic) AS intronic, 
                        sum(is_assoc) AS associated, 
                        sum(is_intergenic) AS intergenic, 
			count(*) - (sum(is_known)+sum(is_unknown)+sum(is_ambiguous)) as unclassified, 
			sum(is_known) - (sum(is_pc)+sum(is_pseudo)+sum(is_npc)+sum(is_utr)) AS unclassified_known, 
			sum(is_unknown) - (sum(is_intronic)+sum(is_intergenic)+sum(is_assoc)) AS unclassified_unknown"""


class KnownAnnotations(Annotations):

    """Annotations of known transcript models."""
    mColumns = ["pc", "pseudo", "npc", "utr", "unclassified"]
    mSelect = """SELECT 
		 sum( is_pc) AS pc,
                 sum( is_pseudo) AS pseudo, 
                 sum( is_npc) AS npc, 
                 sum( is_utr) as utr,
		 sum(is_known) - (sum(is_pc)+sum(is_pseudo)+sum(is_npc)+sum(is_utr)) AS unclassified """


class UnknownAnnotations(Annotations):

    """Annotations of unknown transcript models."""
    mColumns = ["intronic", "associated", "intergenic", "unclassified"]
    mSelect = """SELECT 
			sum(is_intronic) AS intronic, 
                        sum(is_assoc) AS associated, 
                        sum(is_intergenic) AS intergenic,
			sum(is_unknown) - (sum(is_intronic)+sum(is_intergenic)+sum(is_assoc)) AS unclassified"""


class BasicAnnotations(Annotations):

    """Basic annotations of transcript models.

    known: transcript models overlapping reference gene models by more than 95% of its resides.
    novel: transcript models overlapping reference gene models by at most 5% of its resides.
    ambiguous: transcript models partially overlapping reference gene models.
    unclassified: transcript models in none of the categories above.
    """
    mColumns = ["known", "ambiguous", "novel", "unclassified"]
    mSelect = """SELECT sum(is_known) AS known, 
                        sum(is_ambiguous) AS ambiguous, 
                        sum(is_unknown) AS novel, 
			count(*) - (sum(is_known)+sum(is_unknown)+sum(is_ambiguous)) as unclassified"""


class AnnotationsBases(Annotations):

    """Annotations of known transcript models as bases."""
    mColumns = ["total", "CDS", "UTR", "flank", "intronic", "intergenic"]
    mSelect = """SELECT 
                 sum( exons_sum) AS total,
		 sum( nover_CDS ) AS cds,
                 sum( nover_UTR + nover_UTR3 + nover_UTR5 ) AS utr, 
                 sum( nover_flank + nover_3flank + nover_5flank ) AS flank, 
                 sum( nover_intronic) AS intronic,
                 sum( nover_intergenic) AS intergenic
                 """

##########################################################################
##########################################################################
##########################################################################
# Location of UTR transcripts
##########################################################################


class UTRTranscripts(TrackerSQL):

    """return location of UTR transcripts."""

    pattern = "(.*)_annotation$"

    def __call__(self, track, slice=None):

        data = []
        data.append(("utr5", self.getValue(
            "SELECT COUNT(*) FROM %(track)s_annotation WHERE nover1_utr5 > 0" % locals())))
        data.append(("utr3", self.getValue(
            "SELECT COUNT(*) FROM %(track)s_annotation WHERE nover1_utr3 > 0" % locals())))
        data.append(("utr", self.getValue(
            "SELECT COUNT(*) FROM %(track)s_annotation WHERE nover1_utr > 0" % locals())))
        return odict(data)

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class SpliceSites(Annotations):

    """Splice site analysis.

    The columns are:

       total: total number of introns
       with_motif: introns with canonical splice motifs.
    """
    pattern = "(.*)_annotation$"
    mColumns = ("total", "with_motif")
    mSelect = """SELECT 
                 SUM(introns_nval) AS total_transcripts, 
                 SUM(U12_AT_AC+U2_GT_AG+U2_nc_GC_AG) AS with_motif"""


class SpliceMotifs(Annotations):

    """Splice motif analysis.

    The columns are:
       U2_GT/AG: GT/AG introns
       U12_AT/AC: AT/AC introns
       U2_nc_GC/AG: GT/AG introns
       unknown: unknown splice motif
    """
    mPattern = "_annotation$"
    mColumns = ("total", "U2_GT/AG", "U12_AT/AC", "U2_nc_GC/AG,", "unknown")
    mSelect = """SELECT SUM(introns_nval) AS total,                         
                        SUM(U2_GT_AG) AS U2_GT_AG, 
                        SUM(U12_AT_AC) AS U12_AT_AC, 
                        SUM(U2_nc_GC_AG) AS U2_nc_GC_AG,
                        SUM(introns_nval) - SUM(U12_AT_AC+U2_GT_AG+U2_nc_GC_AG) AS unknown"""


class IntronBoundaries(Annotations):

    """Intron boundary analysis.

    The columns are:

       total: introns in transcript models
       found: introns in transcript models present also in reference gene set 
       missed: introns in transcript models not present in reference gene set 
       both: both intron boundaries match
       one:  only one intron boundary matches
       none: no intron boundary matches
       exon_skipping: intron contains an exon in the reference gene set
    """
    pattern = "(.*)_annotation$"
    mColumns = (
        "total", "found", "missed", "both", "one", "none", "exon_skipping")
    mSelect = """SELECT SUM(splice_total) AS total, 
                        SUM(splice_found) AS match, 
                        SUM(splice_missed) AS missed, 
                        SUM(splice_perfect) as both, 
                        SUM(splice_partial) as one,
                        SUM(splice_incomplete) as none,
                        SUM(splice_exon_skipping) as exon_skipping"""

    def getSlices(self, subset=None):
        return ["all", ]


class IntronOverrunLengths(TrackerSQL):

    """Bases of intron overrun within transcript models with intron overrun.
    """
    pattern = "(.*)_overrun$"

    def getSlices(self, subset=None):
        return []

    def getXLabel(self):
        return "bases"

    def __call__(self, track, slice=None):
        data = self.getValues(
            """SELECT nover_intronic from %s_overrun WHERE nover_exonic > 0""" % (track,) )
        return odict((("overrun", data),))


class IntronOverrunCounts(TrackerSQL):

    """Transcript models with intron overrun.
    """
    pattern = "(.*)_overrun$"

    def getSlices(self, subset=None):
        return []

    def __call__(self, track, slice=None):
        data = []
        data.append( ("total", self.getValue( """SELECT COUNT(*)
         FROM %s_overrun WHERE nover_exonic > 0""" % (track,)) ) )
        data.append( ("overrun > 0", self.getValue( """SELECT COUNT(*)
         FROM %s_overrun WHERE nover_intronic > 0 AND nover_exonic > 0""" % (track,)) ) )
        data.append( ("overrun > 10", self.getValue( """SELECT COUNT(*)
         FROM %s_overrun WHERE nover_intronic >= 10 AND nover_exonic > 0""" % (track,)) ) )
        return odict(data)


class IntronOverrunBases(TrackerSQL):

    """Proportion of bases in non CDS sequence of transcripts overlapping CDS.
    """
    pattern = "(.*)_overrun$"

    def getSlices(self, subset=None):
        return []

    def __call__(self, track, slice=None):
        data = self.getValues( '''
             SELECT cast (nover_intronic as float) / length  
             FROM %s_overrun where nover_exonic > 0''' % (track,))
        return odict((("overrun", data),))

##########################################################################
##########################################################################
##########################################################################
# Rates
##########################################################################


class RatesKs(AnnotationsAssociated):

    """Rates of transcript models."""
    pattern = "(.*)_evol$"
    mColumns = "ks"
    mTable = "evol"

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getXLabel(self):
        return self.mColumns

    def __call__(self, track, slice=None):
        data = self.getValues(self.getStatement(track, slice))
        return odict(((self.mColumns, data),))


class RatesKa(RatesKs):

    """Rates of ancestral repeats (ka) around transcript models."""
    mColumns = "ka"


class RatesGC(RatesKs):

    """Percent G+C in transcript models."""
    mColumns = "pgc"


class RatesKi(RatesKs):

    """Rates of introns (ki) within transcript models."""
    mColumns = "ki"


class RatesKsKa(RatesKs):

    """Evolutionary constraint in transcript models. The constraint is computed as the ratio of the substitution rate in the transcript model divided by the median substitution model in ancestral repeats."""
    mColumns = "kska"


class RatesKsKi(RatesKs):

    """Evolutionary constraint in transcript models. The constraint is computed as the ratio of the substitution rate in the transcript model divided by the substitution rate in introns."""
    mColumn = "kski"


class RatesKsVsKa(DefaultSlicer, TrackerSQL):

    """Rate of transcript model versus rate of ancestral repeats."""
    pattern = "(.*)_evol"
    mRate1 = "ks"
    mRate2 = "kska"
    mWhere = "1"

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def __call__(self, track, slice=None):
        if not slice or slice == "all":
            data = list( self.execute( """SELECT %s,%s FROM %s_evol WHERE %s IS NOT NULL AND %s IS NOT NULL AND %s""" %
                                       (self.mRate1, self.mRate2, track, self.mRate1, self.mRate2, self.mWhere)))
        elif slice == "linc":
            data = list( self.execute( """SELECT %s,%s FROM %s_evol AS e, %s_%s AS s WHERE s.gene_id = e.gene_id AND %s IS NOT NULL AND %s IS NOT NULL AND %s""" %
                                       (self.mRate1, self.mRate2, track, track, slice, self.mRate1, self.mRate2, self.mWhere)))
        else:
            data = list( self.execute( """SELECT %s,%s FROM %s_evol AS e, %s_annotation AS a WHERE a.gene_id = e.gene_id AND is_%s AND %s IS NOT NULL AND %s IS NOT NULL AND %s""" %
                                       (self.mRate1, self.mRate2, track, track, slice, self.mRate1, self.mRate2, self.mWhere)))

        return odict(zip((self.mRate1, self.mRate2), zip(*data)))


class RatesKsVsGC(RatesKsVsKa):

    """Substitution rate of transcript model versus G+C content."""
    mRate1 = "ks"
    mRate2 = "pgc"


class RatesKaVsRepeatsGC(RatesKsVsKa):

    """Substitution rate of ancestral repeat rate versus repeat content in repeats."""
    mRate1 = "ka"
    mRate2 = "repeats_gc"


class RatesKsKaVsGC(RatesKsVsKa):

    """Transcript constraint versus G+C content."""
    mRate1 = "kska"
    mRate2 = "pgc"


class RatesKsKaVsRepeatsGC(RatesKsVsKa):

    """Transcript constraint versus G+C content in ancestral repeats."""
    mRate1 = "kska"
    mRate2 = "repeats_gc"


class RatesCoverageVsKs(RatesKsVsKa):

    """Transcript coverage versus substitution rate.

    Transcipt models with coverage of 1 are excluded.
    """
    mRate1 = "meancoverage"
    mRate2 = "ks"
    mWhere = "meancoverage > 1"


class RatesAll(DefaultSlicer, TrackerSQL):

    """All rates."""

    pattern = "(.*)_evol$"
    mColumns = set(("aligned", "ar_aligned", "ir_aligned", "ka", "ki", "ks",
                   "kska", "kski", "length", "meancoverage", "nreads", "pgc", "repeats_gc"))

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def __call__(self, track, slice=None):

        table = sqlalchemy.Table(
            "%s_evol" % track, self.metadata, autoload=True)
        columns = [x.name for x in table.columns if x.name in self.mColumns]

        if not slice or slice == "all":
            data = list(
                self.execute( """SELECT %s FROM %s_evol AS e""" % (",".join( columns), track,)) )
        elif slice == "linc":
            data = list( self.execute( """SELECT e.%s FROM %s_evol AS e, %s_%s AS s WHERE s.gene_id = e.gene_id""" %
                                       (",e.".join(columns), track, track, slice)))
        else:
            data = list( self.execute( """SELECT e.%s FROM %s_evol AS e, %s_annotation AS a WHERE a.gene_id = e.gene_id AND is_%s """ %
                                       (",e.".join(columns), track, track, slice)))

        return odict(zip(columns, zip(*data)))

##########################################################################
##########################################################################
##########################################################################
# Looking at distance
##########################################################################


class ClosestDistance(DefaultSlicer, TrackerSQL):

    """Closest distance of transcript models to gene models in the reference set.
    Transcripts with overlap (closest_dist == 0) are ignored.
    """
    pattern = "(.*)_distances$"
    mXLabel = "distance / bases"
    mColumn = "d.closest_dist"
    mWhere = "d.closest_dist > 0"

    def getTracks(self, subset=None):
        tracks = TrackerSQL.getTracks(self, subset)
        # remove the segment tracks
        tracks = [x for x in tracks if "segments" not in x]
        return tracks

    def __call__(self, track, slice=None):

        column, where = self.mColumn, self.mWhere
        if not slice or slice == "all":
            data = self.getValues(
                """SELECT %(column)s FROM %(track)s_distances AS d WHERE %(where)s""" % locals() )
        elif slice == "linc":
            data = self.getValues( """SELECT %(column)s FROM %(track)s_distances AS d, %(track)s_%(slice)s as s 
                                      WHERE s.gene_id = d.gene_id AND %(where)s""" % locals() )
        else:
            data = self.getValues( """SELECT %(column)s FROM %(track)s_distances AS d, %(track)s_annotation as a 
                                      WHERE d.gene_id = a.gene_id AND a.is_%(slice)s AND %(where)s""" % locals() )
        return odict(((column, data),))


class ClosestDistanceTo5(ClosestDistance):

    """Closest distance of transcript models at 5' end of a gene."""
    mColumn = "d.amin5"
    mWhere = "d.amin5 > 0"


class ClosestDistanceTo3(ClosestDistance):

    """Closest distance of transcript models at 3' end of a gene."""
    mColumn = "d.amin3"
    mWhere = "d.amin3 > 0"


class TranscriptionDensity(DefaultSlicer, TrackerSQL):

    """Closest distance of transcript models to gene models in the reference set."""
    mXLabel = "distance / bases"
    pattern = "(.*)_annotation$"
    mColumn = "d.closest_dist"
    mWhere = "1"
    mMaxSize = 100000
    mReference = trackers_master

    def __call__(self, track, slice=None):

        column, where = self.mColumn, self.mWhere
        reference = self.mReference

        if not slice or slice == "all":
            statement = """SELECT g.start, g.end, e.exons_start, e.exons_end 
                                      FROM %(track)s_distances AS d, 
                                      %(track)s_gtf AS g,
                                      %(reference)s_annotation AS e
                                      WHERE %(where)s AND g.gene_id = d.gene_id AND e.gene_id = d.closest_id""" % locals()

        elif slice == "linc":
            statement = """SELECT g.start, g.end, e.exons_start, e.exons_end 
                                   FROM %(track)s_distances AS d, 
                                   %(track)s_gtf AS g,
                                   %(reference)s_annotation AS e,
                                   %(track)s_%(slice)s as s 
                                   WHERE %(where)s AND g.gene_id = d.gene_id AND e.gene_id = d.closest_id
                                   AND s.gene_id = d.gene_id""" % locals()
        else:
            statement = """SELECT g.start, g.end, e.exons_start, e.exons_end 
                                   FROM %(track)s_distances AS d, 
                                   %(track)s_gtf AS g,
                                   %(track)s_annotation as a, 
                                   %(reference)s_annotation AS e
                                   WHERE %(where)s AND g.gene_id = d.gene_id AND e.gene_id = d.closest_id
                                   AND d.gene_id = a.gene_id AND a.is_%(slice)s""" % locals()

        m = self.mMaxSize
        counts = numpy.zeros(m)

        for seg_start, seg_end, gene_start, gene_end in self.getIter(statement):
            if seg_start > gene_end:
                d = seg_start - gene_end
            elif seg_end < gene_start:
                d = gene_start - seg_end
            l = seg_end - seg_start
            counts[d:min(d + l, m)] += 1

        return odict(zip(("distance", "counts"),
                         (numpy.array(xrange(0, self.mMaxSize))[counts > 0],
                          counts[counts > 0])))

##########################################################################
##########################################################################
##########################################################################
# Looking at the direction of transcription with respect to the closest gene
##########################################################################


class TranscriptionDirectionIntergenic(TrackerSQL):

    """Direction of transcription in novel intergenic transcripts.

    Count by transcripts.
    """
    pattern = "(.*)_annotation$"
    mColumns = "d.closest_strand as strand1, a.exons_strand as strand2, case d.amin3 > 0 when 1 then '3' else '5' END AS pos"
    mWhere = "d.closest_dist < 10000"
    mTable = "distances"

    def getSlices(self, subset=None):
        return []

    def __call__(self, track, slice=None):

        table = self.mTable
        column = self.mColumns
        where = self.mWhere

        statement = """SELECT %(column)s 
                           FROM 
                             %(track)s_%(table)s AS d, 
                             %(track)s_annotation as a 
                           WHERE a.gene_id = d.gene_id AND 
                                 a.is_intergenic AND %(where)s"""

        data = self.get(statement)
        counts = collections.defaultdict(int)
        for x in data:
            counts["".join(x[:3])] += 1

        counts2 = {"same-5": 0, "different-5": 0,
                   "same-3": 0, "different-3": 0,
                   "unknown": 0}

        for key, val in counts.iteritems():
            if "." in key:
                k = "unknown"
            elif key[0] == key[1]:
                k = "same-%s" % key[2]
            else:
                k = "different-%s" % key[2]
            print key, val
            counts2[k] += val

        return odict(list(sorted(counts2.items())))

##########################################################################
##########################################################################
##########################################################################
# Looking at the direction of transcription with respect to the closest gene
##########################################################################


class TranscriptionDirectionIntronic(TrackerSQL):

    """Direction of transcription in novel intronic transcripts.

    Count by transcripts.
    """
    pattern = "%s_vs_(.*)_geneovl$" % trackers_master
    mColumns = "r.exons_strand, a.exons_strand"
    mWhere = "1"
    mReference = trackers_master

    def getSlices(self, subset=None):
        return []

    def __call__(self, track, slice=None):

        column = self.mColumns
        where = self.mWhere
        reference = self.mReference

        if slice == "linc":
            statement = """SELECT %(column)s 
                            FROM %(track)s_annotation AS a, 
                                 %(reference)s_annotation AS r, 
                                 %(reference)s_vs_%(track)s_geneovl AS ovl, 
                                 %(track)s_%(slice)s as s 
                            WHERE a.gene_id = ovl.gene_id2 AND 
                                  r.gene_id = ovl.gene_id1 AND 
                                  a.is_intronic AND 
                                  s.gene_id = a.gene_id AND %(where)s"""
        else:
            statement = """SELECT %(column)s 
                            FROM %(track)s_annotation as a, 
                                    %(reference)s_annotation AS r, 
                                    %(reference)s_vs_%(track)s_geneovl AS ovl
                            WHERE a.gene_id = ovl.gene_id2 AND 
                                  r.gene_id = ovl.gene_id1 AND 
                                  a.is_intronic AND %(where)s"""

        data = self.getAll(statement % locals())

        counts = collections.defaultdict(int)
        for x in data:
            counts["".join(x[:2])] += 1

        counts2 = {"same": 0,
                   "different": 0,
                   "unknown": 0}

        for key, val in counts.iteritems():
            if "." in key:
                k = "unknown"
            elif key[0] == key[1]:
                k = "same"
            else:
                k = "different"
            counts2[k] += val

        return odict(list(sorted(counts2.items())))

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class Segments(TrackerSQL):

    """Return distances between features."""
    mXLabel = "distance"
    pattern = "(.*)_segments_distances$"
    mTable = "_segments_distances"
    mColumn = "distance"

    def getSlices(self, subset=None):
        return []

    def __call__(self, track, slice=None):

        table = self.mTable
        column = self.mColumn
        if not slice or slice == "all":
            data = self.getValues(
                """SELECT d.%(column)s FROM %(track)s%(table)s AS d""" % (locals()) )
        elif slice == "linc":
            data = self.getValues(
                """SELECT d.%(column)s FROM %(track)s%(table)s AS d, %(track)s_%(slice)s as s WHERE s.gene_id = d.gene_id""" % (locals()) )
        else:
            data = self.getValues(
                """SELECT d.%(column)s FROM %(track)s%(table)s AS d, %(track)s_annotation as a WHERE d.gene_id = a.gene_id and a.is_%(slice)s""" % (locals()) )
        return odict(((column, data),))


class SegmentsDistances(Segments):

    """Return distances between features."""
    mXLabel = "distance"
    mTable = "_segments_distances"
    mColumn = "distance"


class SegmentsSizes(Segments):

    """Return sizes between features."""
    mXLabel = "size"
    pattern = "(.*)_segments_sizes$"
    mTable = "_segments_sizes"
    mColumn = "size"


class SegmentsOverlaps(Segments):

    """Return overlaps between features."""
    pattern = "(.*)_segments_overlaps$"
    mXLabel = "overlap"
    mTable = "_segments_overlaps"
    mColumn = "overlap"


class GenesDistances(SegmentsDistances):

    """Return distances between genes."""
    pattern = "(.*)_segments_genes_distances$"
    mTable = "_segments_genes_distances"


class GenesSizes(SegmentsSizes):

    """Return sizes between genes."""
    pattern = "(.*)_segments_genes_sizes$"
    mTable = "_segments_genes_sizes"


class GenesOverlaps(SegmentsOverlaps):

    """Return overlaps between genes."""
    pattern = "(.*)_segments_genes_overlaps$"
    mTable = "_segments_genes_overlaps"

##########################################################################
##########################################################################
##########################################################################
# Introns per transcript
##########################################################################


class IntronsPerTranscript(DefaultSlicer, TrackerSQL):

    """Number of introns per transcript model.
    """
    pattern = "(.*)_annotation$"

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def __call__(self, track, slice=None):
        if not slice or slice == "all":
            data = self.getValues(
                """SELECT introns_nval FROM %s_annotation AS a""" % (track,) )
        elif slice == "linc":
            data = self.getValues(
                """SELECT introns_nval FROM %s_annotation AS a, %s_%s as s WHERE s.gene_id = a.gene_id""" % (track, track, slice) )
        else:
            data = self.getValues(
                """SELECT introns_nval FROM %s_annotation AS a WHERE a.is_%s""" % (track, slice) )

        return odict((("introns_nval", data),))

##########################################################################
##########################################################################
##########################################################################
# Coding potential
##########################################################################


class CodingPotential(AnnotationsAssociated):

    """Coding potential."""
    pattern = "^([^_]+)_coding"
    mColumns = "SUM(t.is_coding) as coding, COUNT(*) - SUM(t.is_coding) AS noncoding"
    mTable = "coding"

    def __call__(self, track, slice=None):
        statement = self.getStatement(track, slice)
        if not statement:
            return []
        return odict(zip(("coding", "non-coding"), self.getFirstRow(statement)))

##########################################################################
##########################################################################
##########################################################################
# Looking at overlap with repeats
##########################################################################


class RepeatOverlap(AnnotationsAssociated):

    """Overlap with repeats."""
    pattern = "(.*)_repeats$"
    mColumns = "SUM(CASE WHEN nover>0 THEN 1 ELSE 0 END) as with, SUM(CASE WHEN nover=0 THEN 1 ELSE 0 END) AS without"
    mTable = "repeats"

    def __call__(self, track, slice=None):
        statement = self.getStatement(track, slice)
        if not statement:
            return []
        return odict(zip(("with", "without"), self.getFirstRow(statement)))

##########################################################################
##########################################################################
##########################################################################
##
##########################################################################


class Overlaps(TrackerSQL):

    """Overlap between sets.

    This tracker returns the overlap between a track and 
    all other tracks. Only one attribute is returned 
    given by :attr:`mColumn`. As the table is not symmetrized,
    the mColumn attribute should be specified without the 
    suffix, i.e. ('nbases_unique' instead of 'nbases_unique1').

    Possible values for mColumn are combinations of 'A_B'.

    A is one of ('nbases', 'nexons', 'ngenes', 'pbases', 'pexons', 'pgenes') 
    where the prefix ''n'' or ''p'' denote the counts or the percent, respectively.

    B is one of ('total', 'uniq', 'ovl')
    """
    mTableName = "overlap_%s"
    mColumn = None
    pattern = "(.*)_annotation$"

    def getSlices(self, subset=None):
        return trackers_default_slices

    def getData(self, track, slice=None):
        if slice == "None":
            slice == "all"
        tablename = self.mTableName % slice
        if self.mColumn is None:
            raise NotImplementedError("mColumn not set in derived class.")
        column = self.mColumn
        data = self.get( "SELECT set2, %(column)s1 FROM %(tablename)s WHERE set1 = '%(track)s'" ) +\
            self.get(
                "SELECT set1, %(column)s2 FROM %(tablename)s WHERE set2 = '%(track)s'")
        return odict(data)


class OverlapsGenesCounts(Overlaps):
    mColumn = "ngenes_ovl"

    def __call__(self, track, slice=None):
        return Overlaps.getData(self, track, slice)


class OverlapsGenesPercent(Overlaps):
    mColumn = "pgenes_ovl"

    def __call__(self, track, slice=None):
        # add the diagonal element of 100%
        return Overlaps.getData(self, track, slice)


class OverlapsBasesCounts(Overlaps):
    mColumn = "nbases_ovl"

    def __call__(self, track, slice=None):
        return Overlaps.getData(self, track, slice)


class OverlapsBasesPercent(Overlaps):
    mColumn = "pbases_ovl"

    def __call__(self, track, slice=None):
        # add the diagonal element of 100%
        return Overlaps.getData(self, track, slice)

##########################################################################
##########################################################################
##########################################################################
# compute transcripts shared between all sets
##########################################################################


class SharedGenes(TrackerSQL):

    """Count the number of shared genes between all sets.

    This counter depends on computing a set merging all transcript
    model (usually called ''merged'').
    ."""
    """Rates of transcript models."""
    pattern = "(.*)_annotation$"
    pattern = "(.*)_annotation"
    mMaster = "merged"

    def getTracks(self, subset=None):
        return ["all", ]

    def getSlices(self, subset=None):
        return trackers_default_slices

    def __call__(self, track, slice=None):
        alltables = TrackerSQL.getTracks(self)
        master = self.mMaster
        tables, data = [], {}

        if slice == "all":
            statement = "SELECT DISTINCT(m.gene_id) FROM %(master)s_vs_%(table)s AS m WHERE nover > 0"
        else:
            statement = "SELECT DISTINCT(m.gene_id) FROM %(master)s_vs_%(table)s AS m, %(master)s_annotation AS a WHERE a.gene_id = m.gene_id AND nover > 0 AND is_%(slice)s"

        # check if set is whithin merged set
        for table in alltables:
            a = self.getValue(
                "SELECT COUNT(DISTINCT gene_id) FROM %(table)s_annotation" % locals())
            try:
                b = self.getValue(
                    "SELECT COUNT(DISTINCT gene_id2) FROM %(master)s_vs_%(table)s_ovl" % locals())
            except SQLError:
                continue

            if b != a:
                continue

            tables.append(table)
            data[table] = set(self.getValues(statement % locals()))

        if len(tables) == 0:
            return

        all = set(data[tables[0]])
        for t in tables[1:]:
            all = all.intersection(set(data[t]))

        return odict([("counts", len(all)), ])

# ########################################################################
# ########################################################################
# ########################################################################
# generic implementations of Annotator results
# ########################################################################
# class Annotator( TrackerSQL ):
#     pattern = "(.*)_annotation"
#     mTableName = None
#     mSelect = "SELECT CASE over WHEN '+' THEN 'over' ELSE 'under' END, category, fold, pvalue, observed, expected FROM %s WHERE fdr "
#     mColumns = ("over", "category", "fold", "pvalue", "observed", "expected")
#     mOrder = "pvalue"
#     def __init__(self, *args, **kwargs ):
#         TrackerSQL.__init__(self, *args, **kwargs)

#     def getTracks( self, subset = None ):
#         if not self.mTableName: raise NotImplementedError("table not specified")
# return self.getValues( "SELECT DISTINCT track FROM %s" % self.mTableName
# )

#     def getSlices( self, subset = None ):
#         if not self.mTableName: raise NotImplementedError("table not specified")
# return self.getValues( "SELECT DISTINCT slice || ':' || subset || ':' ||
# workspace FROM %s" % self.mTableName )

#     def __call__(self, track, slice = None ):

#         if not self.mTableName: raise NotImplementedError("table not specified")

#         select = self.mSelect % self.mTableName
#         order = self.mOrder
#         if slice == "all" or slice is None:
#             data = list( self.execute( """%s AND track = '%s' ORDER BY %s""" % (select, track, order)).fetchone() )
#         else:
#             slice, subset, workspace = slice.split(":")
#             data = self.getAll( """%(select)s AND track = '%(track)s' AND slice = '%(slice)s'
# AND subset='%(subset)s' AND workspace = '%(workspace)s' ORDER BY
# %(order)s""" % locals())

#         return odict( zip(self.mColumns, zip(*data) ) )

# class AnnotatorSummary( Annotator ):
#     mSelect = "SELECT COUNT(*), SUM( CASE over WHEN '+' THEN 1 ELSE 0 END) FROM %s WHERE fdr "
#     mColumns = ("total", "over", "under" )

#     def __call__(self, track, slice = None ):

#         if not self.mTableName: raise NotImplementedError("table not specified")

#         select = self.mSelect % self.mTableName

#         if slice == "all" or slice is None:
#             data = self.getFirstRow( """%s AND track = '%s'""" % (select, track, order))
#         else:
#             slice, subset, workspace = slice.split(":")
#             data = self.getFirstRow( """%(select)s AND track = '%(track)s' AND slice = '%(slice)s'
# AND subset='%(subset)s' AND workspace = '%(workspace)s'""" % locals())

#         data.append( data[0] - data[1] )
#         return odict( zip(self.mColumns, data) )

# class AnnotatorEnrichment( Annotator):
#     mSelect = "SELECT category, 100.0 * (fold - 1), pvalue FROM %s WHERE fdr AND over = '+' "
#     mColumns = ("category", "enrichment/%", "pvalue" )
#     mOrder = "fold DESC"

# class AnnotatorDepletion( Annotator ):
#     mSelect = "SELECT category, 100.0 * (1 - fold), pvalue FROM %s WHERE fdr AND over = '-' "
#     mColumns = ("category", "depletion/%", "pvalue" )
#     mOrder = "fold"

# class AnnotatorPower( TrackerSQL ):
#     pattern = "(.*)_annotation"
#     mTableName = None

#     def __init__(self, *args, **kwargs ):
#         TrackerSQL.__init__(self, *args, **kwargs)

#     def getTracks( self, subset = None ):
#         if not self.mTableName: raise NotImplementedError("table not specified")
# return self.getValues( "SELECT DISTINCT track || ':' || slice || ':' ||
# subset || ':' || workspace FROM %s ORDER BY
# track,slice,subset,workspace" % self.mTableName )

#     def getSlices( self, subset = None ):
#         if not self.mTableName: raise NotImplementedError("table not specified")
#         if not subset: return []
#         return [",".join( subset ) ]

#     def __call__(self, track, slice = None ):

#         if not self.mTableName: raise NotImplementedError("table not specified")
#         if not self.mSelect: raise NotImplementedError("invalid use of base class - only use derived classes")

#         select = self.mSelect % self.mTableName
#         xtrack, xslice, subset, workspace = track.split(":")

#         stmt = """%(select)s AND track = '%(xtrack)s' AND slice = '%(xslice)s'
#                         AND subset='%(subset)s' AND workspace = '%(workspace)s'
#                         ORDER BY category""" % locals()

#         data = self.getValues( stmt )

#         return odict( (("power", data),))

# class AnnotatorPowerEnrichment( AnnotatorPower ):
#     """display power of go analyses.

#     Show are the fold enrichment (in %) that could be detected
#     given the setup of the test (size of sample and size of workspace).
#     """
#     mSelect = "SELECT 100.0 * ( (ci95high / expected) - 1) FROM %s WHERE 1"
#     mXLabel = "Power : fold enrichment / %"

# class AnnotatorPowerDepletion( AnnotatorPower ):
#     """display power of go analyses.

#     Show are the fold depletion (in %) that could be detected
#     given the setup of the test (size of sample and size of workspace).
#     """
#     mSelect = "SELECT 100.0 * ( 1 - (ci95low / expected)) FROM %s WHERE 1"
#     mXLabel = "Power : fold depletion / %"

# class AnnotatorMatrix( TrackerSQL ):
#     """Display GO Annotator results in a matrix.

#     Slicing is possible using a keyword:value syntax. Keyword
#     can be one of 'workspace', 'slice', 'subset' or 'track'. For example,
#     the slices workspace:interegenic,slice:unknown will only display
#     results where workspace is intergenic and slice is unknown.

#     Only results passing the FDR rate test are shown.
#     """

#     mTableName = None
#     mSelect = "SELECT category, CASE over WHEN '+' THEN (100.0 * (fold - 1)) ELSE - (100.0 * (1-fold) ) END FROM %s WHERE fdr"
#     mColumns = ("enrichment", "pvalue")

#     def __init__(self, *args, **kwargs ):
#         TrackerSQL.__init__(self, *args, **kwargs)

#     def getTracks( self, subset = None ):
#         if not self.mTableName: raise NotImplementedError("table not specified")
# return self.getValues( "SELECT DISTINCT track || ':' || slice || ':' ||
# subset || ':' || workspace FROM %s ORDER BY
# track,slice,subset,workspace" % self.mTableName )

#     def getSlices( self, subset = None ):
#         if not self.mTableName: raise NotImplementedError("table not specified")
#         if not subset: return []
#         return [",".join( subset ) ]

#     def __call__(self, track, slice = None ):

#         if not self.mTableName: raise NotImplementedError("table not specified")

#         select = self.mSelect % self.mTableName
#         xtrack, xslice, subset, workspace = track.split(":")

#         if slice:
#             o = collections.defaultdict( list )
#             pairs = slice.split(",")
#             for pair in pairs:
#                 key,value = pair.split(":")
#                 o[key] = value

#             if "workspace" in o and workspace not in o["workspace"]: return []
#             if "subset" in o and subset not in o["subset"]: return []
#             if "slice" in o and xslice not in o["slice"]: return []
#             if "track" in o and xtrack not in o["track"]: return []


#         stmt = """%(select)s AND track = '%(xtrack)s' AND slice = '%(xslice)s'
#                         AND subset='%(subset)s' AND workspace = '%(workspace)s'
#                         ORDER BY category""" % locals()

#         data = self.getAll( stmt )
#         return odict(data)

# =================================================================
# mixin classes to pair annotator results with a table.
# =================================================================
# class AnnotatorGOTerritories(object):
#     mTableName = "goterritories_annotators"

# class AnnotatorGOSlimTerritories(object):
#     mTableName = "goslimterritories_annotators"

# class AnnotatorIntronicGOTerritories(object):
#     mTableName = "intronicgoterritories_annotators"

# class AnnotatorIntronicGOSlimTerritories(object):
#     mTableName = "intronicgoslimterritories_annotators"

# class AnnotatorIntergenicGOTerritories(object):
#     mTableName = "intergenicgoterritories_annotators"

# class AnnotatorIntergenicGOSlimTerritories(object):
#     mTableName = "intergenicgoslimterritories_annotators"

# class AnnotatorUnknownSets(object):
#     mTableName = "unknownsets_annotators"

# class AnnotatorArchitecture(object):
#     mTableName = "architecture_annotators"

# class AnnotatorIntronicSets(object):
#     mTableName = "intronicsets_annotators"

# class AnnotatorIntergenicSets(object):
#     mTableName = "intergenicsets_annotators"

# =================================================================
# specific implementations of Annotator results
# =================================================================
# _annotator_analysis = { "Enrichment" : AnnotatorEnrichment,
#                        "Depletion" : AnnotatorDepletion,
#                        "Summary" : AnnotatorSummary,
#                        "Matrix" : AnnotatorMatrix,
#                        "PowerEnrichment" : AnnotatorPowerEnrichment,
#                        "PowerDepletion" : AnnotatorPowerDepletion}
# _annotator_territories = { "GOTerritories" : AnnotatorGOTerritories,
#                            "GOSlimTerritories" : AnnotatorGOSlimTerritories,
#                            "IntergenicGOTerritories" : AnnotatorIntergenicGOTerritories,
#                            "IntergenicGOSlimTerritories" : AnnotatorIntergenicGOSlimTerritories,
#                            "IntronicGOTerritories" : AnnotatorIntronicGOTerritories,
#                            "IntronicGOSlimTerritories" : AnnotatorIntronicGOSlimTerritories,
#                            "UnknownSets" : AnnotatorUnknownSets,
#                            "Architecture" : AnnotatorArchitecture,
#                            "IntronicSets" : AnnotatorIntronicSets,
#                            "IntergenicSets" : AnnotatorIntergenicSets }


# the order of the base classes is important
# also: make sure that these are new-style classes
# for a, aa in _annotator_analysis.items():
#     for b, bb in _annotator_territories.items():
#         n = "Annotator%s%s" % (a,b)
#         globals()[n] = type( n, (bb,aa), {})

##########################################################################
##########################################################################
##########################################################################
# generic implementations of GO results
##########################################################################
class GO(TrackerSQL):
    pattern = "(.*)_annotation"
    mTableName = None
    mSelect = "SELECT CASE code WHEN '+' THEN 'over' ELSE 'under' END, description, ratio, CASE code WHEN '+' THEN pover ELSE punder END AS pvalue' FROM %s WHERE passed AND code != '?'"
    mColumns = ("over", "category", "fold", "pvalue")
    mOrder = "pvalue"

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getTracks(self, subset=None):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        return self.getValues("SELECT DISTINCT track FROM %s" % self.mTableName)

    def getSlices(self, subset=None):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        return self.getValues("SELECT DISTINCT slice || ':' || subset || ':' || background FROM %s" % self.mTableName)

    def __call__(self, track, slice=None):

        if not self.mTableName:
            raise NotImplementedError("table not specified")

        select = self.mSelect % self.mTableName
        order = self.mOrder
        if slice == "all" or slice is None:
            data = list(self.execute(
                """%s AND track = '%s' ORDER BY %s""" % (select, track, order)).fetchone() )
        else:
            slice, subset, background = slice.split(":")
            data = self.getAll( """%(select)s AND track = '%(track)s' AND slice = '%(slice)s' 
                        AND subset='%(subset)s' AND background = '%(background)s' ORDER BY %(order)s""" % locals())

        return odict(zip(self.mColumns, zip(*data)))


class GOSummary(GO):
    mSelect = "SELECT COUNT(*), SUM( CASE code WHEN '+' THEN 1 ELSE 0 END) FROM %s WHERE passed AND code != '?'"
    mColumns = ("total", "over", "under")

    def __call__(self, track, slice=None):

        if not self.mTableName:
            raise NotImplementedError("table not specified")

        select = self.mSelect % self.mTableName

        if slice == "all" or slice is None:
            data = self.getFirstRow(
                """%s AND track = '%s'""" % (select, track, order))
        else:
            slice, subset, background = slice.split(":")
            data = self.getFirstRow( """%(select)s AND track = '%(track)s' AND slice = '%(slice)s' 
                        AND subset='%(subset)s' AND background = '%(background)s'""" % locals())

        data.append(data[0] - data[1])
        return odict(zip(self.mColumns, data))


class GOEnrichment(GO):
    mSelect = "SELECT description, 100.0 * (ratio - 1), pover FROM %s WHERE passed AND code = '+' "
    mColumns = ("category", "enrichment/%", "pvalue")
    mOrder = "ratio DESC"


class GODepletion(GO):
    mSelect = "SELECT description, 100.0 * (1 - ratio), punder FROM %s WHERE passed AND code = '-' "
    mColumns = ("category", "depletion/%", "pvalue")
    mOrder = "ratio"


class GOPower(TrackerSQL):
    pattern = "(.*)_annotation"
    mTableName = None

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getTracks(self, subset=None):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        return self.getValues("SELECT DISTINCT track || ':' || slice || ':' || subset || ':' || background FROM %s ORDER BY track,slice,subset,background" % self.mTableName)

    def getSlices(self, subset=None):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        if not subset:
            return []
        return [",".join(subset)]

    def __call__(self, track, slice=None):

        if not self.mTableName:
            raise NotImplementedError("table not specified")
        if not self.mSelect:
            raise NotImplementedError(
                "invalid use of base class - only use derived classes")

        select = self.mSelect % self.mTableName
        xtrack, xslice, subset, background = track.split(":")

        stmt = """%(select)s AND track = '%(xtrack)s' AND slice = '%(xslice)s' 
                        AND subset='%(subset)s' AND background = '%(background)s'
                        AND code != '?'
                        ORDER BY category""" % locals()

        data = self.getValues(stmt)

        return odict((("power", data),))


class GOPowerEnrichment(GOPower):

    """display power of go analyses.

    Show are the fold enrichment (in %) that could be detected
    given the setup of the test (size of sample and size of background).
    """
    mSelect = "SELECT 100.0 * ( (ci95upper / mean) - 1) FROM %s WHERE 1"
    mXLabel = "Power : fold enrichment / %"


class GOPowerDepletion(GOPower):

    """display power of go analyses.

    Show are the fold depletion (in %) that could be detected
    given the setup of the test (size of sample and size of background).
    """
    mSelect = "SELECT 100.0 * ( 1 - (ci95lower / mean)) FROM %s WHERE 1"
    mXLabel = "Power : fold depletion / %"


class GOMatrix(TrackerSQL):

    """Display GO results in a matrix.

    Slicing is possible using a keyword:value syntax. Keyword
    can be one of 'background', 'slice', 'subset' or 'track'. For example,
    the slices background:ensembl,slice:unknown will only display
    results where background is ensembl and slice is unknown.
    """
    mTableName = "goslim_goanalysisresults"
    mSelect = "SELECT description, CASE code WHEN '+' THEN (100.0 * (ratio -1)) ELSE - (100.0 * (1-ratio) ) END FROM %s WHERE 1 AND code != '?' "

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, **kwargs)

    def getTracks(self, subset=None):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        return self.getValues("SELECT DISTINCT track || ':' || slice || ':' || subset || ':' || background FROM %s ORDER BY slice,track,subset,background" % self.mTableName)

    def getSlices(self, subset=None):
        if not self.mTableName:
            raise NotImplementedError("table not specified")
        if not subset:
            return []
        return [",".join(subset)]

    def __call__(self, track, slice=None):

        if not self.mTableName:
            raise NotImplementedError("table not specified")
        select = self.mSelect % self.mTableName
        xtrack, xslice, subset, background = track.split(":")

        if slice:
            o = collections.defaultdict(list)
            pairs = slice.split(",")
            for pair in pairs:
                key, value = pair.split(":")
                o[key] = value
            if "background" in o and background not in o["background"]:
                return []
            if "subset" in o and subset not in o["subset"]:
                return []
            if "slice" in o and xslice not in o["slice"]:
                return []
            if "track" in o and xtrack not in o["track"]:
                return []

        stmt = """%(select)s 
                        AND track = '%(xtrack)s' 
                        AND slice = '%(xslice)s' 
                        AND subset='%(subset)s' 
                        AND background='%(background)s'
                        AND passed 
                        AND code != '?'
                        ORDER BY category""" % locals()

        data = self.getAll(stmt)

        return odict(data)

# =================================================================
# mixin classes to pair annotator results with a table.
# =================================================================


class GOGenes(object):
    mTableName = "go_goanalysisresults"


class GOSlimGenes(object):
    mTableName = "goslim_goanalysisresults"


class GOTerritories(object):
    mTableName = "go_territorygoanalysisresults"


class GOSlimTerritories(object):
    mTableName = "goslim_territorygoanalysisresults"

# =================================================================
# specific implementations of GO results
# =================================================================
_go_analysis = {"Enrichment": GOEnrichment,
                "Depletion": GODepletion,
                "Summary": GOSummary,
                "Matrix": GOMatrix,
                "PowerEnrichment": GOPowerEnrichment,
                "PowerDepletion": GOPowerDepletion}

_go_territories = {"GOGenes": GOGenes,
                   "GOSlimGenes": GOSlimGenes,
                   "GOTerritories": GOTerritories,
                   "GOSlimTerritories": GOSlimTerritories}

# the order of the base classes is important
# also: make sure that these are new-style classes
for a, aa in _go_analysis.items():
    for b, bb in _go_territories.items():
        n = "GO%s%s" % (a, b)
        globals()[n] = type(n, (bb, aa), {})

##########################################################################
##########################################################################
##########################################################################


class ReadAssociated(AnnotationsAssociated):

    """simple join between a data table and table defining slices.

    The join works from transcripts to reads.

    :attr:`mTable`
       table to join with
    :attr:`mColums`
       columns to output
    """
    pattern = "(.*)_readmap$"
    mTable = None
    mColumns = None
    mSelectAll = "SELECT %(columns)s FROM %(track)s_%(table)s AS t"
    mSelectSubset = "SELECT %(columns)s FROM %(track)s_%(table)s AS t, %(track)s_annotation AS a, %(track)s_readmap AS m WHERE t.read_id = m.read_id AND m.gene_id = a.gene_id AND a.is_%(slice)s"
    mSelectSlice = "SELECT %(columns)s FROM %(track)s_%(table)s AS t, %(track)s_%(slice)s AS s, %(track)s_readmap AS m WHERE t.read_id = m.read_id AND m.gene_id = s.gene_id"
    mSelectMixture = "SELECT %(columns)s FROM %(track)s_%(table)s AS t, %(subset)s AS s, %(track)s_annotation AS a, %(track)s_readmap AS m  WHERE t.read_id = m.read_id AND m.gene_id = a.gene_id AND a.is_%(slice)s AND s.gene_id = m.gene_id"

# =================================================================
# Read statistics
# =================================================================


class ReadStats(ReadAssociated):

    """Read stats information."""
    pattern = "(.*)_readstats$"
    mColumns = ""
    mTable = "readstats"

    def __init__(self, *args, **kwargs):
        AnnotationsAssociated.__init__(self, *args, **kwargs)

    def getXLabel(self):
        return self.mColumns

    def __call__(self, track, slice=None):
        data = self.getValues(self.getStatement(track, slice))
        return odict(((self.mColumns, data),))


class ReadStatsCoverage(ReadStats):

    """read statistics - coverage of reads."""
    mColumns = "qCov"


class ReadStatsPercentIdentity(ReadStats):

    """read statistics - percent identity of reads."""
    mColumns = "pid"


class ReadStatsCounts(AnnotationsAssociated):

    """simple join between a data table and table defining slices.

    The join works from transcripts to reads.

    :attr:`mTable`
       table to join with
    :attr:`mColums`
       columns to output

    Note: the default slices have been disabled, only known, ambiguous and unknown are returned.
    """
    pattern = "(.*)_readmap$"
    mTable = "readmap"
    mColumns = "COUNT(DISTINCT read_id)"

    def __call__(self, track, slice=None):
        data = []
        data.append(
            ("known", self.getValue(self.getStatement(track, slice="known"))))
        data.append(
            ("ambiguous", self.getValue(self.getStatement(track, slice="ambiguous"))))
        data.append(
            ("novel", self.getValue(self.getStatement(track, slice="unknown"))))
        return odict(data)

# =================================================================
# Read statistics
# =================================================================


class ReadInfo(ReadAssociated):

    """Read stats information."""
    pattern = "(.*)_readinfo$"
    mColumns = ""
    mTable = "readinfo"

    def __init__(self, *args, **kwargs):
        AnnotationsAssociated.__init__(self, *args, **kwargs)

    def getXLabel(self):
        return self.mColumns

    def __call__(self, track, slice=None):
        data = self.getValues(self.getStatement(track, slice))
        return odict(((self.mColumns, data),))


class ReadInfoPercentGC(ReadInfo):

    """read statistics - G+C content of read."""
    mColumns = "pGC"

import os, sys, re, types, itertools
import matplotlib.pyplot as plt
import numpy
import numpy.ma
import Stats
import Histogram

from SphinxReport.Tracker import *
from cpgReport import *
from SphinxReport.odict import OrderedDict as odict

##################################################################################
class genomeCompostionSummary(cpgTracker):
    """Average GC content and CpG density in the genome """

    mPattern = "bam_stats"

    def __call__(self, track, slice = None ):

        data = self.getFirstRow( """SELECT round(sum(length*pGC)/sum(length),3) as mean_GC,
                                    round(sum(length*CpG_ObsExp)/sum(length),3) as mean_CpG_ObsExp,
                                    round(sum(length*pCpG)/sum(length),3) as mean_CpG
                                    FROM annotations.genome
                                    WHERE id <> 'total'
                                    AND id not like 'chrX%%'
                                    AND id not like 'chrY%%'
                                    AND id not like 'chrW%%'
                                    AND id not like 'chrZ%%'
                                    AND id not like 'chrM%%'
                                    AND id not like '%%random%%' """ )
        mColumns = [ "GC content", "CpG Obs/Exp" , "CpG density" ]

        return odict( zip(mColumns, data) )

##################################################################################
class genomeCompostionPerContig(cpgTracker):
    """Average GC content and CpG density per Contig """

    mPattern = "bam_stats"

    def __call__(self, track, slice = None ):

        data = self.getAll( """SELECT g.id as Contig, g.length, m.mappable_bases, a.repeat_length, 
                               g.pGC as GC_content, g.pCpG as CpG_density, g.CpG_ObsExp
                               FROM annotations.genome g, annotations.mappable_bases_per_contig m,
                               (select contig, sum(stop-start) as repeat_length from repeats group by contig) a
                               WHERE g.id=m.contig AND a.contig=g.id AND g.id <> "total"
                               ORDER BY g.length desc LIMIT 100""" )
        return data

##################################################################################
class allTranscriptsByBiotype(cpgTracker):
    """Number of Ensembl transcripts by biotype"""

    mPattern = "bam_stats"

    def __call__(self, track, slice = None ):

        data = self.getAll( """SELECT transcript_biotype, count(transcript_id) as transcripts
                               FROM annotations.transcript_info
                               GROUP BY transcript_biotype
                               ORDER BY transcripts desc """ )
        return data

##################################################################################
class allGenesByBiotype(cpgTracker):
    """Number of Ensembl genes by biotype"""

    mPattern = "bam_stats"

    def __call__(self, track, slice = None ):

        data = self.getAll( """SELECT gene_biotype, count(distinct gene_id) as genes
                               FROM transcript_info
                               GROUP BY gene_biotype
                               ORDER BY genes desc """ )
        return data



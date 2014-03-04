import os, sys, re, types, itertools

from SphinxReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *

class recessives( ExomeTracker ):

    mPattern = "_recessive_table$"

    def __call__(self, track, slice = None):
        data = self.get("SELECT CHROM, POS, CASE WHEN (LENGTH(REF)>6) THEN SUBSTR(REF,1,6)||'...' ELSE REF END AS REF, CASE WHEN (LENGTH(ALT)>6) THEN SUBSTR(ALT,1,6)||'...' ELSE ALT END AS ALT, ID, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID FROM %(track)s_recessive_table WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.01) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.01) AND COMMON is null AND FILTER='PASS' " % locals() )
        return odict( zip( ("CHROM", "POS", "REF", "ALT", "ID", "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_GENE_NAME", "SNPEFF_TRANSCRIPT_ID" ), zip(*data)) )

class dominants( ExomeTracker ):

    mPattern = "_dominant_table$"

    def __call__(self, track, slice = None):
        data = self.get("SELECT CHROM, POS, CASE WHEN (LENGTH(REF)>6) THEN SUBSTR(REF,1,6)||'...' ELSE REF END AS REF, CASE WHEN (LENGTH(ALT)>6) THEN SUBSTR(ALT,1,6)||'...' ELSE ALT END AS ALT, ID, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID, dbNSFP_ESP6500_AA_AF FROM %(track)s_dominant_table WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.001) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.001) AND COMMON is null AND FILTER='PASS' " % locals() )
        return odict( zip( ("CHROM", "POS", "REF", "ALT", "ID", "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_GENE_NAME", "SNPEFF_TRANSCRIPT_ID", "dbNSFP_ESP6500_AA_AF" ), zip(*data)) )

class deNovos( ExomeTracker ):

    mPattern = "_filtered_table$"
    
    def __call__(self, track, slice = None):
        data = self.get("SELECT CHROM, POS, REF, ALT, ID, SNPEFF_EFFECT, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID FROM %(track)s_filtered_table WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.001) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.001) AND COMMON is null AND FILTER='PASS' " % locals() )
        return odict( zip( ("CHROM", "POS", "REF", "ALT", "ID", "SNPEFF_EFFECT", "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_GENE_NAME", "SNPEFF_TRANSCRIPT_ID" ), zip(*data)) )

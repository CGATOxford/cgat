import os
import sys
import re
import types
import itertools

from SphinxReport.Tracker import *
from collections import OrderedDict as odict
from exomeReport import *


class recessives(ExomeTracker):
    
    pattern = "(.*)_recessive_table$"

    def __call__(self, track, slice=None):
        return self.getAll(
            "SELECT CHROM, POS, REF, ALT, ID, SNPEFF_CODON_CHANGE, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_EXON_ID, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID, SNPEFF_EFFECT, SNPEFF_FUNCTIONAL_CLASS, SNPEFF_IMPACT, SNPEFF_GENE_BIOTYPE, dbNSFP_1000Gp1_AF, dbNSFP_ESP6500_AA_AF, dbNSFP_ESP6500_EA_AF, dbNSFP_29way_logOdds, dbNSFP_GERP___NR, dbNSFP_GERP___RS, dbNSFP_Interpro_domain, dbNSFP_Polyphen2_HVAR_pred, dbNSFP_SIFT_score, FILTER, BaseQRankSum, FS, HaplotypeScore, MQ, MQ0, MQRankSum, QD, ReadPosRankSum FROM %(track)s_recessive_table WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.01) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.01) AND COMMON is null AND FILTER='PASS' " % locals())

class dominants(ExomeTracker):

    pattern = "(.*)_dominant_table$"

    def __call__(self, track, slice=None):
        return self.getAll(
            "SELECT CHROM, POS, REF, ALT, ID, SNPEFF_CODON_CHANGE, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_EXON_ID, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID, SNPEFF_EFFECT, SNPEFF_FUNCTIONAL_CLASS, SNPEFF_IMPACT, SNPEFF_GENE_BIOTYPE, dbNSFP_1000Gp1_AF, dbNSFP_ESP6500_AA_AF, dbNSFP_ESP6500_EA_AF, dbNSFP_29way_logOdds, dbNSFP_GERP___NR, dbNSFP_GERP___RS, dbNSFP_Interpro_domain, dbNSFP_Polyphen2_HVAR_pred, dbNSFP_SIFT_score, FILTER, BaseQRankSum, FS, HaplotypeScore, MQ, MQ0, MQRankSum, QD, ReadPosRankSum FROM %(track)s_dominant_table WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.01) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.01) AND COMMON is null AND FILTER='PASS' " % locals())

class deNovos(ExomeTracker):

    pattern = "(.*)_filtered_table$"

    def __call__(self, track, slice=None):
        return self.getAll(
            "SELECT CHROM, POS, REF, ALT, ID, SNPEFF_CODON_CHANGE, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_EXON_ID, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID, SNPEFF_EFFECT, SNPEFF_FUNCTIONAL_CLASS, SNPEFF_IMPACT, SNPEFF_GENE_BIOTYPE, dbNSFP_1000Gp1_AF, dbNSFP_ESP6500_AA_AF, dbNSFP_ESP6500_EA_AF, dbNSFP_29way_logOdds, dbNSFP_GERP___NR, dbNSFP_GERP___RS, dbNSFP_Interpro_domain, dbNSFP_Polyphen2_HVAR_pred, dbNSFP_SIFT_score, FILTER, BaseQRankSum, FS, HaplotypeScore, MQ, MQ0, MQRankSum, QD, ReadPosRankSum FROM %(track)s_filtered_table WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.01) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.01) AND COMMON is null AND FILTER='PASS' " % locals())

#class recessives(ExomeTracker):
#
#    mPattern = "_recessive_table$"
#
#    def __call__(self, track, slice=None):
#        data = self.get(
#            "SELECT CHROM, POS, CASE WHEN (LENGTH(REF)>6) THEN SUBSTR(REF,1,6)||'...' ELSE REF END AS REF, CASE WHEN (LENGTH(ALT)>6) THEN SUBSTR(ALT,1,6)||'...' ELSE ALT END AS ALT, ID, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID FROM %(track)s_recessive_table WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.01) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.01) AND COMMON is null AND FILTER='PASS' " % locals())
#        return odict(zip(("CHROM", "POS", "REF", "ALT", "ID", "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_GENE_NAME", "SNPEFF_TRANSCRIPT_ID"), zip(*data)))
#
#
#class dominants(ExomeTracker):
#
#    mPattern = "_dominant_table$"
#
#    def __call__(self, track, slice=None):
#        data = self.get(
#            "SELECT CHROM, POS, CASE WHEN (LENGTH(REF)>6) THEN SUBSTR(REF,1,6)||'...' ELSE REF END AS REF, CASE WHEN (LENGTH(ALT)>6) THEN SUBSTR(ALT,1,6)||'...' ELSE ALT END AS ALT, ID, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID, dbNSFP_ESP6500_AA_AF FROM %(track)s_dominant_table WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.001) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.001) AND COMMON is null AND FILTER='PASS' " % locals())
#        return odict(zip(("CHROM", "POS", "REF", "ALT", "ID", "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_GENE_NAME", "SNPEFF_TRANSCRIPT_ID", "dbNSFP_ESP6500_AA_AF"), zip(*data)))
#
#
#class deNovos(ExomeTracker):
#
#    mPattern = "_filtered_table$"
#
#    def __call__(self, track, slice=None):
#        data = self.get(
#            "SELECT CHROM, POS, REF, ALT, ID, SNPEFF_EFFECT, SNPEFF_AMINO_ACID_CHANGE, SNPEFF_GENE_NAME, SNPEFF_TRANSCRIPT_ID FROM %(track)s_filtered_table WHERE (dbNSFP_1000Gp1_AF is null OR dbNSFP_1000Gp1_AF<0.001) AND (dbNSFP_ESP6500_AA_AF is null OR dbNSFP_ESP6500_AA_AF<0.001) AND COMMON is null AND FILTER='PASS' " % locals())
#        return odict(zip(("CHROM", "POS", "REF", "ALT", "ID", "SNPEFF_EFFECT", "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_GENE_NAME", "SNPEFF_TRANSCRIPT_ID"), zip(*data)))

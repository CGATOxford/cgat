import os
import sys
import re
import types
import collections
from VariantsReport import *

SNP_CODES = ("E", "O")


class StrainTracker(VariantsTracker):

    '''tracker returning results per strain.'''
    mPattern = "_effects$"

    def getTracks(self, subset=None):

        tracks = VariantsTracker.getTracks(self)

        if subset == None:
            return ["all"] + tracks
        else:
            return tracks

##########################################################################
##########################################################################
##########################################################################
# Annotation of bases with SNPs
##########################################################################


class SNPBaseAnnotation(VariantsTracker):

    """types of snps.
    """

    mPattern = "_annotations$"

    def __call__(self, track, slice=None):
        '''getting data'''
        data = self.get( """SELECT UPPER(code), COUNT(*) FROM %(track)s_annotations 
                              WHERE variant_type IN ('E', 'O')
                              GROUP BY code""" % locals() )
        return odict(data)

##########################################################################
##########################################################################
##########################################################################
# Annotation of bases with SNPs
##########################################################################


class SNPBaseAnnotationWithoutNonCoding(VariantsTracker):

    """types of snps excluding intergenic and intronic
    """

    mPattern = "_annotations$"

    def __call__(self, track, slice=None):

        data = self.get( """SELECT UPPER(code), COUNT(*) FROM %(track)s_annotations 
                            WHERE variant_type IN ('E', 'O')
                            AND code not in ('I','i','G','g') GROUP BY code""" % locals() )
        return odict(data)

##########################################################################
##########################################################################
##########################################################################
# Annotation of bases with SNPs
##########################################################################


class SNPBaseAnnotationsCodingType(VariantsTracker):

    """types of snps excluding intergenic and intronic
    """

    mPattern = "_annotations$"

    def __call__(self, track, slice=None):

        data = self.get( """SELECT UPPER(code), COUNT(*) FROM %(track)s_annotations 
                            WHERE variant_type IN ('E', 'O')
                            AND code not in ('I','i','G','g') GROUP BY code""" % locals() )
        return odict(data)

##########################################################################
##########################################################################
##########################################################################
# Variant types
##########################################################################


class VariantTypeCounts(VariantsTracker):

    """return counts for each variant type.

    The counts are symmetrized, i.e., W,D = D,W
    """

    mPattern = "_annotations$"

    def __call__(self, track, slice=None):

        data = self.get(
            "SELECT variant_type, COUNT(*) FROM %(track)s_annotations GROUP BY variant_type" % locals())
        combined = collections.defaultdict(int)
        for key, counts in data:
            try:
                if "," in key:
                    key = ",".join(sorted(key.split(",")))
            except TypeError:
                pass
            combined[key] += counts
        return combined

##########################################################################
##########################################################################
##########################################################################
# Quality score distribution
##########################################################################


class SNPQuality(VariantsTracker):

    '''Distribution of SNP consensus quality scores.

    '''
    mPattern = "_annotations$"
    mQuality = "snp_quality"

    def __call__(self, track, slice=None):
        quality = self.mQuality
        snpcodes = "','".join(SNP_CODES)
        data = self.getValues(
            "SELECT %(quality)s FROM %(track)s_annotations WHERE variant_type in ('%(snpcodes)s')" % locals())
        return odict((("quality", data),))


class SNPConsensusQuality(SNPQuality):
    mQuality = "consensus_quality"


class SNPRMSMappingQuality(SNPQuality):
    mQuality = "rms_mapping_quality"

##########################################################################
##########################################################################
##########################################################################
# Substitutions
##########################################################################


class SubstitutionMatrixNucleotidesOld(VariantsTracker):

    '''Matrix of single nucleotide substitution polymorphisms.

    Counts both homozygous ([ACGT]) and heterozygous substitutions ([^ACGT])
    '''
    mPattern = "_annotations$"

    def getSlices(self, subset=None):
        return ("A", "C", "G", "T")

    def __call__(self, track, slice=None):
        snpcodes = "','".join(SNP_CODES)
        data = self.get( '''SELECT genotype, COUNT(*) FROM %(track)s_annotations 
                               WHERE reference_base='%(slice)s' AND variant_type IN ('%(snpcodes)s') 
                               GROUP BY genotype''' % locals() )
        d = odict([x for x in data if x[0] in "ACGT"] + [(x[0].lower(), x[1])
                  for x in data if x[0] not in "ACGT"])
        return d

##########################################################################
##########################################################################
##########################################################################
# Substitutions
##########################################################################


class SubstitutionMatrixNucleotides(VariantsTracker):

    '''Matrix of single nucleotide substitution polymorphisms.

    Counts both homozygous ([ACGT]) and heterozygous substitutions ([^ACGT])
    '''
    mPattern = "_annotations_summary$"

    slices = ("all", "coding", "noncoding")

    # def getSlices( self, subset = None):
    # if subset == "all": return ["all"]
    # else: return ["A","B","C","D", "G","I","M","N","P","R","S","U","V","X"]
    #     if subset == None:
    #         return []
    #     else:
    #         return subset

    def __call__(self, track, slice=None):

        if slice == "all":
            where = "1"
        elif slice == "coding":
            where = "upper(code) IN ('A','B','C','D') " % locals()
        elif slice == "noncoding":
            where = "upper(code) IN ('G','I')" % locals()
        else:
            where = "upper(code) = '%(slice)s' " % locals()

        data = self.get( '''SELECT reference_base, genotype, sum(counts)
                               FROM %(track)s_annotations_summary
                               WHERE %(where)s
                               GROUP BY reference_base, genotype''' % locals() )

        d = odict()
        for reference_base, genotype, val in data:
            if reference_base not in d:
                d[reference_base] = odict()
            d[reference_base][genotype] = val
        return d

##########################################################################
##########################################################################
##########################################################################
# Substitutions
##########################################################################


class TransitionTransversionRatios(VariantsTracker):

    '''Matrix of single nucleotide substitution polymorphisms.

    Counts both homozygous ([ACGT]) and heterozygous substitutions ([^ACGT])
    '''
    mPattern = "_annotations_summary$"

    slices = ("all", "coding", "noncoding")

    def __call__(self, track, slice=None):

        if slice == "all":
            where = "1"
        elif slice == "coding":
            where = "upper(code) IN ('A','B','C','D') " % locals()
        elif slice == "noncoding":
            where = "upper(code) IN ('G','I')" % locals()
        else:
            where = "upper(code) = '%(slice)s' " % locals()

        data = self.get( '''SELECT reference_base, genotype, sum(counts)
                               FROM %(track)s_annotations_summary
                               WHERE %(where)s
                               GROUP BY reference_base, genotype''' % locals() )

        bases = "ACGT"

        d = odict([(x, odict([(y, 0) for y in bases])) for x in bases])
        for reference_base, genotype, val in data:
            d[reference_base][genotype] = val

        transitions = d["A"]["G"] + d["C"]["T"] + d["G"]["A"] + d["T"]["C"]
        transversions = d["A"]["C"] + d["A"]["T"] +\
            d["C"]["A"] + d["C"]["G"] +\
            d["G"]["C"] + d["G"]["T"] +\
            d["T"]["A"] + d["T"]["G"]

        if transversions == 0:
            ratio = "na"
        else:
            ratio = float(transitions) / transversions

        return odict((("transitions", transitions),
                      ("transversions", transversions),
                      ("ratio", ratio)))

##########################################################################
##########################################################################
##########################################################################
# Substitutions
##########################################################################


class SubstitutionMatrixAminoAcids(VariantsTracker):

    '''Matrix of amino acid substitution polymorphisms for
    coding sequence.

    Counts only homozygous ([ACGT]) differences. 
    '''
    pattern = "(.*)_annotations$"
    letters = "ACDEFGHIKLMNPQRSTVWXY"
    slices = tuple("ACDEFGHIKLMNPQRSTVWXY")

    def __call__(self, track, slice=None):
        data = self.get(
            "SELECT variant_aa, COUNT(*) FROM %(track)s_annotations WHERE reference_aa='%(slice)s' AND (variant_type='E' OR variant_type='O') and variant_aa != NULL GROUP BY variant_aa" % locals())
        v = odict((x, 0) for x in self.letters)
        v.update(odict([x for x in data if "," not in x[0]]))
        return v


class SubstitutionMatrixCounts(Tracker):

    '''Matrix of the expected effects of changes.'''

    mGeneticCode = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "X",  # stop
        "TAG": "X",  # stop
        "TGT": "C",
        "TGC": "C",
        "TGA": "X",  # stop
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }

    def getTracks(self):
        return ["all", ]

    def __call__(self, track, slice=None):

        nonsyn, stop, syn = 0, 0, 0

        for codon, aa in self.mGeneticCode.items():
            if aa == "X":
                continue
            for pos in range(0, 3):
                other_codon = list(codon)
                for x in "ACGT":
                    if x == codon[pos]:
                        continue
                    other_codon[pos] = x
                    other_aa = self.mGeneticCode["".join(other_codon)]
                    if aa == other_aa:
                        syn += 1
                    elif other_aa == "X":
                        stop += 1
                    else:
                        nonsyn += 1

        return odict((("synonymous", syn),
                      ("non-synonymous", nonsyn),
                      ("stop", stop), )
                     )


class TranscriptCounts(VariantsTracker):

    tracks = ("annotations.transcript_info",)

    def __call__(self, track, slice=None):
        r = odict()
        r["genes"] = self.getValue(
            "SELECT COUNT(DISTINCT gene_id) FROM %(track)s" % locals())
        r["transcripts"] = self.getValue(
            "SELECT COUNT(DISTINCT transcript_id) FROM %(track)s" % locals())
        return r

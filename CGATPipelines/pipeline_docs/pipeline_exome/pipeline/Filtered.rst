=================
Filtered variants
=================

Candidate homozygous recessive variants
=======================================

The following tables present candidate pathogenic variants for
recessive disease.  Please note that this data may be missing if
the project in question contains only dominant families, or
single cases.

.. report:: Filtered.recessives
   :render: xls-table
   :groupby: track
   :force:

For definition of column headings see bottom of page. 

Candidate dominant variants
===========================

The following tables present candidate pathogenic variants for
dominant disease.  Please note that this data may be missing if
the project in question contains only recessive families, or
single cases.

.. report:: Filtered.dominants
   :render: xls-table
   :groupby: track
   :force:

For definition of column headings see bottom of page.

De Novos
========

The following tables present de novo variants detected in trios.
Please note that this data may be missing if the project in question
contains only families with multiple affected individuals or single
cases.

.. report:: Filtered.deNovos
   :render: xls-table
   :groupby: track

Definition of column headings:
CHROM = chromosome
POS = genomic coordinate (hg19)
REF = reference allele
ALT = alternate allele
ID = rs identifier
SNPEFF_CODON_CHANGE = codon change
SNPEFF_AMINO_ACID_CHANGE = amino acid change
SNPEFF_EXON_ID = exon of transcript SNPEFF_TRANSCRIPT_ID in which 
	       the variant resides
SNPEFF_GENE_NAME = gene in which the variant resides
SNPEFF_TRANSCRIPT_ID = transcript upon which annotation is based
SNPEFF_EFFECT = variant effect eg. non-synonymous
SNPEFF_FUNCTIONAL_CLASS = variant functional class eg. none, silent,
			missense, nonsense
SNPEFF_IMPACT = predicted variant impact eg. high, moderate, low,
	      modifier
SNPEFF_GENE_BIOTYPE = if available, eg. protein-coding, pseudogene
dbNSFP_1000Gp1_AF = frequency in 1000 genomes
dbNSFP_ESP6500_AA_AF = frequency in ESP6500 African Americans
dbNSFP_ESP6500_EA_AF = frequency in ESP6500 European Americans
dbNSFP_29way_logOdds = SiPhy score based on 29 mammals genomes. The
		     larger the score, the more conserved the site
dbNSFP_GERP___NR = GERP++ neutral rate
dbNSFP_GERP___RS = GERP++ RS score, the larger the score, the more
		 conserved the site
dbNSFP_Interpro_domain = domain or conserved site on which the variant
		       locates. Domain annotations come from Interpro 
		       database. The number in the brackets following
		       a specific domain is the count of times
		       Interpro assigns the variant position to that 
		       domain, typically coming from different
		       predicting databases
dbNSFP_Polyphen2_HVAR_pred = Polyphen2 prediction based on HumVar, 'D'
			   ('probably damaging'), 'P' ('possibly 
			   damaging') and 'B'('benign'). Multiple
			   entries separated by ';'
dbNSFP_SIFT_score = SIFT score - if a score is smaller than 0.05 the
		  corresponding NS is predicted as 'D(amaging)';
		  otherwise it is predicted as 'T(olerated)'
FILTER = 'PASS' if the variant passes GATK's variant quality score
       recalibration filter
BaseQRankSum = u-based z-approximation from the Mann-Whitney Rank Sum 
	     Test for base qualities(ref bases vs. bases of the
	     alternate allele)
FS = Phred-scaled p-value using Fisher's Exact Test to detect strand
   bias (the variation being seen on only the forward or only the reverse
   strand) in the reads. More bias is indicative of false positive
   calls
HaplotypeScore = Higher scores are indicative of regions with bad
	       alignments, often leading to artifactual SNP and indel
	       calls
MQ = Root Mean Square of the mapping quality of the reads across all
   samples in the family
MQ0 = Total count across all samples of mapping quality zero reads
MQRankSum = u-based z-approximation from the Mann-Whitney Rank Sum
	  Test for mapping qualities (reads with ref bases vs. those
	  with the alternate allele)
QD = Variant confidence (from the QUAL field) / unfiltered depth of
   non-reference samples.  Note that the QD is also normalized by event
   length. Low scores are indicative of false positive calls and
   artifacts
ReadPosRankSum = u-based z-approximation from the Mann-Whitney Rank
	       Sum Test for the distance from the end of the read for
	       reads with the alternate allele. If the alternate
	       allele is only seen near the ends of reads, this is 
	       indicative of error

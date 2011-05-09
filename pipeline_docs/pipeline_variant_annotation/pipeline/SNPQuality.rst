==================
SNP quality scores
==================

SNP quality
=================

The SNP quality is the Phred-scaled probability of the consensus being identical to the reference.

 .. report:: Trackers.SNPQuality
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-bins: arange(0,500,5)

   Distribution of SNP quality scores

Consensus quality
=================

Each consensus genotype is associated with a phred quality that measures the probability that 
the consensus genotype is incorrect. See the maq paper for more details. The quality score takes
into account the expected heterozygosity r and the base qualities. According to the Maq paper:

   Before consensus calling, MAQ first combines mapping quality and base quality. If a read is incorrectly mapped, 
   any sequence differences inferred from the read cannot be reliable. Therefore, the base quality used in SNP calling 
   cannot exceed the mapping quality of the read. MAQ reassigns the quality of each base as the smaller value between 
   the read mapping quality and the raw sequencing base quality.

.. report:: Trackers.SNPConsensusQuality
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-bins: arange(0,500,5)

   Distribution of SNP consensus quality scores

RMS mapping quality
===================

To evaluate the reliability of alignments, MAQ assigns each individual alignment a phred-scaled quality score (capped at 99), 
which measures the probability that the true alignment is not the one found by MAQ. 

Qs = -10 log_10 Pr( Read is wrongly mapped)

The RMS mapping quality is the root mean square 
mapping quality of all reads aligned at the SNP.

.. report:: Trackers.SNPRMSMappingQuality
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-bins: arange(0,100,1)

   Distribution of SNP mapping quality scores

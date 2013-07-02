===============
Panther results
===============

Overview
========

The following data shows global summary statistics over all
predictions. Because there can be several transcripts per gene,
these numbers do include some double counting. For example,
the same SNP in two distinct transcripts will be counted twice.
On the other hand, the same SNP present in several strains will only
be counted once.

The following table provides an overview on how many predictions
were run and how many were successful or unsuccessful.

+------------------------------+---------------------------------------------------------------------------------+
|*column*                      |*content*                                                                        |
+------------------------------+---------------------------------------------------------------------------------+
|total                         |number of SNPs supplied to Panther                                               |
+------------------------------+---------------------------------------------------------------------------------+
|No Panther HMMhit             |Panther could not map and HMM to the SNP's position (no prediction)              |
+------------------------------+---------------------------------------------------------------------------------+
|SNP position in protein does  |SNP maps to an insert state of the HMM.                                          |
|not align to HMM              |                                                                                 |
+------------------------------+---------------------------------------------------------------------------------+
|weak Panther score            |PANTHER refused to predict effect to the insufficient data                       |
+------------------------------+---------------------------------------------------------------------------------+
|predicted                     |number of SNPs with PANTHER prediction                                           |
+------------------------------+---------------------------------------------------------------------------------+

.. report:: Panther.PantherCounts
   :render: table
   :groupby: slice

   Number of snps/peptides for which there are predictions.

Distribution of parameters
==========================

The following plots show the overall distribution of the
SNP effect.

+------------------------------+---------------------------------------------------------------------------------+
|*column*                      |*content*                                                                        |
+------------------------------+---------------------------------------------------------------------------------+
|Pdeleterious                  |probability of a SNP to be called deleterious                                    |
+------------------------------+---------------------------------------------------------------------------------+
|Nic                           |number of independent counts                                                     |
+------------------------------+---------------------------------------------------------------------------------+
|subPSEC                       |subPSEC score (Substitution Position-Specific Evolutionary Conservation)         |
+------------------------------+---------------------------------------------------------------------------------+

.. report:: Panther.PantherDistribution
   :render: line-plot                                  
   :transform: histogram 
   :as-lines:

   Distribution of various columns in the PANTHER table.

In PANTHER, a :term:`subPSEC` score of less than -3 corresponds 
roughly to a probability of 0.5 that this particular SNP is deleterious
(:term:`Pdeleterious`). Such scores correspond to SNPS causing mendelian
diseases (see :pmid:`15492219`).

Effects per strain
==================

A SNP is taken to be deleterious, if the estimated probability of the SNP being 
deleterious greater than :param:`Panther.PantherResults.pdeleterious`. SNPs
are counted by :term:`locus_id`.

+------------------------------+---------------------------------------------------------------------------------+
|*column*                      |*content*                                                                        |
+------------------------------+---------------------------------------------------------------------------------+
|deleterious                   |number of deleterious SNPs                                                       |
+------------------------------+---------------------------------------------------------------------------------+
|neutral                       |number of neutral SNPs                                                           |
+------------------------------+---------------------------------------------------------------------------------+

.. report:: Panther.PantherResults
   :render: table

   Number of SNPS that are predicted to be deleterious per strain. 

.. report:: Panther.PantherResults
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total
	
   Proportion of SNPS that are predicted to be deleterious per strain.

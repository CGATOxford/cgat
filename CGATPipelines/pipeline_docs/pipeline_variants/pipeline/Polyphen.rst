===================================
Polyphen results
===================================

Overview
========

The following data shows global summary statistics over all
predictions. Because there can be several transcripts per gene,
these numbers do include some double counting. For example,
the same SNP in two distinct transcripts will be counted twice.

The following table provides an overview on how many predictions
were run and how many were successful or unsuccessful.

+------------------------------+-----------------------------------------------------+
|*column*                      |*content*                                            |
+------------------------------+-----------------------------------------------------+
|total                         |total number of SNPS/peptides                        |
+------------------------------+-----------------------------------------------------+
|benign                        |SNP is most likely lacking any phenotypic effect     |
|                              |                                                     |
+------------------------------+-----------------------------------------------------+
|possiblydamaging              |SNP is supposed to affect protein function or        |
|                              |structure                                            |
+------------------------------+-----------------------------------------------------+
|probablydamaging              |SNP is with high confidence supposed to affect       |
|                              |protein function or structure                        |
|                              |                                                     |
+------------------------------+-----------------------------------------------------+
|unknown                       |could not predict effect                             |
+------------------------------+-----------------------------------------------------+
|deleterious                   |SNP is in deleterious position                       |
+------------------------------+-----------------------------------------------------+
|neutral                       |SNP is in neutral position                           |
+------------------------------+-----------------------------------------------------+
|none                          |SNP position is unknown                              |
+------------------------------+-----------------------------------------------------+
|with structure                |SNP overlaps a structure                             |
+------------------------------+-----------------------------------------------------+
|with PFAM                     |SNP has a PFAM hit                                   |
+------------------------------+-----------------------------------------------------+

See `polyphen help <http://tux.embl-heidelberg.de/ramensky/doc/pph_help.html>`_
for a more detailed description.

.. report:: Polyphen.PolyphenCounts
   :render: table
   :groupby: track

   Number of snps/peptides/genes for which there are predictions.

+------------------------------+--------------------------------------------------+
|*column*                      |*content*                                         |
+------------------------------+--------------------------------------------------+
|pph2_prob                     |naive Bayes posterior probability that this       |
|                              |mutation is damaging                              |
+------------------------------+--------------------------------------------------+
|dscore                        |difference in PSCI scores                         |
+------------------------------+--------------------------------------------------+
|Nobs                          |number of proteins in multiple alignment          |
+------------------------------+--------------------------------------------------+
|NStruct                       |number of structures mapping to SNP position      |
+------------------------------+--------------------------------------------------+
|dVol                          |difference in volume                              |
+------------------------------+--------------------------------------------------+
|pph2_FDR                      |Polyphen FDR                                      |
+------------------------------+--------------------------------------------------+
|dProp                         |change in accessible surface propensity resulting |
|                              |from the substitution                             |
+------------------------------+--------------------------------------------------+

.. report:: Polyphen.PolyphenDistribution
   :render: line-plot 
   :transform: histogram
   :as-lines:
   :layout: column-2

   Distribution of various columns in the Polyphen table.

The following plot shows the number of snps in a gene versus the gene length. 
Each dot is colored by the proportion of deleterious SNPs found. This plot
shows if there is trend of short or long peptides having an abnormally large
proportion of deleterious SNPs.

.. report:: Polyphen.PolyphenEffectAndLength
   :render: scatter-rainbow-plot
   :mpl-rc: lines.markersize=4;lines.markeredgewith=0
   :palette: Blues

   Logarithm of snps in a gene versus the logarithm of gene length. Each dot
   is coloured by the proportion of deleterious SNPs in the gene.

Effects per strain
==================

The following table lists predicted effects per strain. SNPs are counted by :term:`locus_id`.

+------------------------------+-----------------------------------------------------+
|*column*                      |*content*                                            |
+------------------------------+-----------------------------------------------------+
|benign                        |SNP is most likely lacking any phenotypic effect     |
+------------------------------+-----------------------------------------------------+
|possiblydamaging              |SNP is supposed to affect protein function or        |
|                              |structure                                            |
+------------------------------+-----------------------------------------------------+
|probablydamaging              |SNP is with high confidence supposed to affect       |
|                              |protein function or structure                        |
+------------------------------+-----------------------------------------------------+
|unknown                       |could not predict effect                             |
+------------------------------+-----------------------------------------------------+

.. report:: Polyphen.PolyphenResults
   :render: table

   Number of SNPS that are predicted to be deleterious per strain. 

.. report:: Polyphen.PolyphenResults
   :render: matrix
   :format: %5.2f
   :transform-matrix: normalized-row-total

   Proportion of SNPS that are predicted to be deleterious per strain. 

.. report:: Polyphen.PolyphenResults
   :render: stacked-bar-plot
   :transform-matrix: normalized-row-total
   :layout: columns-2

   Proportion of SNPS that are predicted to be deleterious per strain.

Distribution of SNPs per gene
=============================

The following plot shows the distribution of SNPs with certain
predicted effects per gene. Plotted are the number of genes
with a given number of SNPs of a certain predicted effect. Shown
are the cumulative distributions.

.. report:: Polyphen.PolyphenSnpsPerGeneAndCategory                                                                                                                                                                                          
   :render: line-plot                                                                                                                                                                                                                        
   :transform: histogram                                                                                                                                                                                                                     
   :groupby: track                                                                                                                                                                                                                           
   :as-lines:                                                                                                                                                                                                                                
   :tf-aggregate: cumulative
   :yrange: 0,
   :layout: columns-2

   Number of genes with predicted effects across all strains

The following plot shows the number of genes per strain which carry
at least on SNP with certain predicted effects

.. report:: Polyphen.PolyphenDeleteriousGenesPerStrain
   :render: interleaved-bar-plot
   :layout: columns-2

   Number of genes per strain with predicted effects

Enrichment of SNPs within genes
=================================

The following table and plots show the number of genes
that have a significant enrichment of a :term:`SNP` or a
:term:`dSNP`. The enrichment is computed using a bionomial
distribution and using an FDR threshold of 0.05.

The result is encoded by a three letter code:

+--------------------+------------------------------+
|**position**        |**test**                      |
+--------------------+------------------------------+
|1                   |dSNPS in SNPs per gene        |
+--------------------+------------------------------+
|2                   |SNPs per gene                 |
+--------------------+------------------------------+
|3                   |dSNPs per gene                |
+--------------------+------------------------------+
 
A ``1`` indicates that a test passed, while a ``0``
indicates that a test failed. The code ``---`` signifies 
that no SNP was present in that protein.

.. report:: Polyphen.PolyphenEnrichment
   :render: table

   Proportion of genes significantly enriched in
   various tests.

.. report:: Polyphen.PolyphenEnrichment
   :render: pie-plot
   :layout: columns-2

   Proportion of genes significantly enriched in
   various tests.

The following plots examines if there is a length dependence in the
statistical enrichment tests.

.. report:: Polyphen.PolyphenEnrichmentLength
   :render: box-plot
   :groupby: track
   :yrange: 0,5000
   :layout: columns-2

   Length of genes within the various categories of 
   enriched genes.


Species distribution of deleterious SNPs
========================================

The following plot shows the number of species that a SNP
is found in stratified by the predicted effect of the SNP.

.. report:: Polyphen.PolyphenSpeciesDistribution      
   :render: line-plot
   :transform: histogram  
   :tf-aggregate: normalized-total,cumulative
   :tf-range: ,,1 
   :groupby: track
   :as-lines: 
   :layout: columns-2
   :yrange: 0,1

   Number of strains that a SNP of a certain effect is found in


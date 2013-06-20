*****************************
Analysis of intron overrun
*****************************

Intron overrun can be used as a proxy towards 
genomic contamination. In this analysis the
pipeline counts the number of bases in a
transcript model that fall into regions
of the genome denoted as intronic. Only 
transcript models that have some overlap with
reference gene models are considered.

Possible confounding factors are:

* The mapping pipeline removes reads extending 
   into introns to same extent.

* Missing exons in the reference gene models will be labeled as
   genomic contamination.

Due to nested genes and truncated exons in the ENSEMBL gene set
it can be ambiguous what to call an intron. Hence there is also
usually some intron overrun within the reference gene set.

Intron Overrun
==============

.. report:: Trackers.IntronOverrunCounts
   :render: interleaved-bar-plot
   :transform-matrix: normalized-row-max

   Plot of transcript models with intron overrun.

.. report:: Trackers.IntronOverrunLengths
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total,reverse-cumulative

   Reverse cumulative plot of bases in introns of transcript models that overlap a reference gene model.

.. report:: Trackers.IntronOverrunLengths
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total,reverse-cumulative

   Reverse cumulative plot of bases in introns of transcript models that overlap a reference gene model.

The following plot shows the proportion of bases in transcript models that are believed to be exonic,
but that overlap intronic sequence in the reference gene set. Only transcripts that overlap exons in
the reference gene set are included.

.. report:: Trackers.IntronOverrunBases                                                                                                                                                                                                      
   :render: line-plot                                                                                                                                                                                                                        
   :transform: histogram                                                                                                                                                                                                                     
   :tf-aggregate: normalized-total,reverse-cumulative                                                                                                                                                                                        
   :as-lines:        
                                                                                                                                                                                                                                             
   Proportion of bases mapping to non-CDS sequence in transcript models that overlap CDS sequence.
   The plot is normalized and reverse-cumulative.

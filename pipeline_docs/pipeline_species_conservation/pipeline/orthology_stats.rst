===================
7-way Conservation
===================

Genes Associated with CpG Islands
----------------------------------

The following table presents an overview of the number of protein-coding/conserved genes 
in each species that have a CAPseq interval overlapping the TSS of at least one transcript.

.. report:: orthology.trackSummary
   :render: table

   Orthology Summary

.. report:: orthology.trackSummary
   :render: interleaved-bar-plot
   :slices: proportion_conserved, proportion_with_feature, proportion_conserved_with_feature

   Orthology Summary

Level of Conservation Across Species
------------------------------------

.. report:: orthology.orthologyGroupCounts
   :render: table

   Number of conserved genes with CpG islands

.. report:: orthology.orthologyGroupCounts
   :render: line-plot

   Number of conserved genes with CpG islands


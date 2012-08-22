=======================
Pairwise Conservation
=======================

This heatmap indicates the conservation of NMIs in each species pair. 
The score is calculated as the total number of conserved genes that have an NMI in both species
over the sum of the total conserved genes with an NMI in each species.

.. report:: orthology.pairwiseHeatmap
   :render: matrixNP-plot
   :zrange: 0,1

   Pairwise conservation of NMIs

.. report:: orthology.pairwiseTable
   :render: table

   Pairwise conservation of NMIs
   
.. report:: orthology.pairwiseStats
   :render: table

   Probability of NMIs gene lists overlap
   
Pairwise Conservation Alternative Score
========================================

This heatmap indicates the conservation of NMIs in each species pair. 
The score is calculated as the total number of conserved genes that have an NMI in both species
over the minimum number of conserved genes with an NMI in either species.

.. report:: orthology.pairwiseHeatmap2
   :render: matrixNP-plot
   :zrange: 0,1

   Pairwise conservation of NMIs

.. report:: orthology.pairwiseTable2
   :render: table

   Pairwise conservation of NMIs
   
Threeway Conservation
======================

The following table indicates the numbers of genes with NMIs that are shared between human, mouse and zebrafish.

.. report:: orthology.threewayVenn
   :render: table

   Threeway conservation of NMIs (human, mouse, zebrafish)


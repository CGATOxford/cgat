===================================
Overlapped Genes H3K27Me3 Analysis
===================================

This section looks for enrichment of H3K27Me3 intervals in genes overlapped >90% by NMIs.

Intersection
=======================

.. report:: macs_replicated_overlapped_genes.polycombIntersection
   :render: table
   :groupby: track
   :force:

   Intersection of H3K27Me3 intervals and genes >90% covered by NMIs


.. report:: macs_replicated_overlapped_genes.overlappedGenesH3K27Venn
   :render: gallery-plot
   :glob: plots/*overlapped.genes*H3K27Me3*.png
   :layout: column-2

   Intersection of H3K27Me3 intervals and genes >90% covered by NMIs  
   
GAT
=========

.. report:: macs_replicated_overlapped_genes.polycombGAT
   :render: table
   :groupby: track
   :force:

   GAT analysis


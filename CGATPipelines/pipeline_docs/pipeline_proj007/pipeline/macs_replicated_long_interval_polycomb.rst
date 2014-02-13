======================================
Long Interval Genes H3K27Me3 Analysis
======================================

This section looks for enrichment of H3K27Me3 intervals in genes overlapped by NMIs >3kb in length.

Intersection
=======================

.. report:: macs_replicated_long_interval_genes.longGenesH3K27Venn
   :render: gallery-plot
   :glob: plots/*long*H3K27Me3*.png
   :layout: column-2

   Intersection of H3K27Me3 intervals and genes overlapped by NMIs >3kb in length

.. report:: macs_replicated_long_interval_genes.shortGenesH3K27Venn
   :render: gallery-plot
   :glob: plots/*short*H3K27Me3*.png
   :layout: column-2

   Intersection of H3K27Me3 intervals and genes overlapped by NMIs <2kb in length
      
GAT
=========

.. report:: macs_replicated_long_interval_genes.longPolycombGAT
   :render: table
   :groupby: track
   :force:

   GAT analysis


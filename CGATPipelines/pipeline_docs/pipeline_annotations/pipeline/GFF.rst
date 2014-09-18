GFF formatted intervals
=======================

This page summarizes some basic counts in the GFF formatted
files in the annotation pipeline. GFF files typically contain
genomic features.

.. report:: AnnotationReport.GFFStats
   :render: interleaved-bar-plot
   :slices: features,sources

   Number of ``features`` and ``sources`` in GFF files.


.. report:: AnnotationReport.GFFStats
   :render: interleaved-bar-plot
   :slices: contigs

   Number of chromosomes with annotations in GFF files.


Bases/intervals per source
--------------------------

:term:`GTF` and :term:`GFF` files annotate features
with a ``source``. The following plots show the number of
intervals/bases for each source in the various files:

.. report:: AnnotationReport.GTFSummaryPerSource
   :render: matrix-plot
   :logscale: z
   :mpl-figure: figsize=(20,20)
   :transform: pivot
   :pivot-index: track
   :pivot-column: source
   :pivot-value: intervals

   Number of intervals


.. report:: AnnotationReport.GTFSummaryPerSource
   :render: matrix-plot
   :logscale: z
   :mpl-figure: figsize=(20,20)
   :transform: pivot
   :pivot-index: track
   :pivot-column: source
   :pivot-value: bases

   Number of bases




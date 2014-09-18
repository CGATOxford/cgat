GTF formatted intervals
=======================

This page summarizes some basic counts in the :term:`GTF` formatted
files in the annotation pipeline. GTF files typically contain
gene sets (exons, transcripts, genes).

.. report:: AnnotationReport.GTFStats
   :render: interleaved-bar-plot
   :slices: features,sources

   Number of ``features`` and ``sources`` in GTF files.


.. report:: AnnotationReport.GTFStats
   :render: interleaved-bar-plot
   :slices: contigs

   Number of chromosomes with annotations in GTF files.


Numbers of genes/transcripts/exons
-----------------------------------

.. report:: AnnotationReport.GTFStats
   :render: interleaved-bar-plot
   :slices: genes,transcripts

   Number of genes and transcripts in GTF files.


.. report:: AnnotationReport.GTFStats
   :render: interleaved-bar-plot
   :slices: exon_count_nval

   Number of exons/intervals in GTF files


Size distributions
------------------

.. report:: AnnotationReport.GTFStats
   :render: interleaved-bar-plot
   :slices: exon_count_min, exon_count_median, exon_count_max

   Median, minimum/maximum number of exons in a transcript


.. report:: AnnotationReport.GTFStats
   :render: interleaved-bar-plot
   :slices: transcript_size_min, transcript_size_median, transcript_size_max

   Median, minimum/maximum transcript size


.. report:: AnnotationReport.GTFStats
   :render: interleaved-bar-plot
   :slices: exon_size_min, exon_size_median, exon_size_max

   Median, minimum/maximum exon size


.. report:: AnnotationReport.GTFStats
   :render: interleaved-bar-plot
   :slices: intron_size_min, intron_size_median, intron_size_max

   Median, minimum/maximum of intron sizes


All data
---------

.. report:: AnnotationReport.GTFStats
   :render: table
   :force:

   Summary statistics of GTF files, all files


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






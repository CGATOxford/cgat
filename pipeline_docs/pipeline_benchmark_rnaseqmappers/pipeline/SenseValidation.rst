=====================================
Direction of transcription validation
=====================================

Exons
=====

.. report:: BenchmarkReport.CoverageTotalsExons
   :render: stacked-bar-plot
   :transform: filter
   :tf-fields: sense,antisense

   Number of reads mapped to genes.

.. report:: BenchmarkReport.CoverageProportionsExons
   :render: table
   :transform: stats
   :groupby: all
   
   Proportion of reads in exons.

.. report:: BenchmarkReport.CoverageProportionsExons
   :render: interleaved-bar-plot
   :transform: stats,select
   :tf-fields: mean
   :groupby: all

   Average proportion of antisense reads in exons.

.. report:: BenchmarkReport.CoverageTotalsExons
   :render: table
   :groupby: slice

   Number of reads within exons

.. report:: BenchmarkReport.CoverageTotalsExons
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: anysense_percent,sense_percent,ratio

   Proportion of reads in exons

.. report:: BenchmarkReport.CoverageTotalsExons
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: ratio

   Proportion of reads in exons

Genes
=====

Genes include both exons and introns.

.. report:: BenchmarkReport.CoverageProportionsRegions
   :render: table
   :transform: stats
   :groupby: all
   
   Proportion of reads in genes.

.. report:: BenchmarkReport.CoverageProportionsRegions
   :render: interleaved-bar-plot
   :transform: stats,select
   :tf-fields: mean
   :groupby: all

   Average proportion of antisense reads in genes

.. report:: BenchmarkReport.CoverageTotalsRegions
   :render: table
   :groupby: slice

   Number of reads within genes

.. report:: BenchmarkReport.CoverageTotalsRegions
   :render: interleaved-bar-plot
   :transform: filter
   :tf-fields: anysense_percent,sense_percent,ratio

   Proportion of reads in genes


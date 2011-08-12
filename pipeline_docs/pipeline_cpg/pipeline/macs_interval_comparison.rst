==========================
MACS Interval Comparisons
==========================

Unique Intervals
================

The following table presents the number of intervals that were unique to each dataset:

.. report:: macs_interval_comparison.UniqueIntervals
   :render: table

   Intervals unique to an individual dataset

Shared Intervals
================

The following table presents the number of intervals in each dataset, 
which overlap with one or more intervals in all other datasets:

.. report:: macs_interval_comparison.SharedIntervals
   :render: table

   Intervals shared between all datasets

Pairwise Comparisons
====================

The following table presents the number of intervals that overlap between each dataset pair:

.. report:: macs_interval_comparison.OverlapIntervals
   :render: table
   :force:

   Number of Overlapping intervals for each pair of datasets

Overlap with Predicted CGIs
===========================

The following table presents the number of intervals that overlap with predicted CpG islands for each dataset:

.. report:: macs_interval_comparison.OverlapCpG
   :render: table
   :force:

   Number of intervals overlapping with predicted CpG islands

Overlap with ChIP-seq Datasets
==============================

The following table presents the number of intervals that overlap with revevant ChIP-seq datasets for each dataset:

.. report:: macs_interval_comparison.OverlapChipseq
   :render: table
   :force:

   Number of intervals overlapping with chipseq intervals

Overlap with CAP-seq Datasets
=============================

The following table presents the number of intervals that overlap with revevant CAP-seq datasets for each dataset:

.. report:: macs_interval_comparison.OverlapCAPseq
   :render: table
   :force:

   Number of intervals overlapping with CAPseq intervals

Overlap with Chromatin Marks
============================

The following table presents the number of intervals that overlap with revevant chromatin mark datasets for each dataset:

.. report:: macs_interval_comparison.OverlapChromatinMarks
   :render: table
   :force:

   Number of intervals overlapping with chromatin modification intervals

Genomic Annotation Tester
=========================

The following table presents the correlation of genomic intervals from different datasets.

.. report:: macs_interval_comparison.gatResults
   :render: table

   Genomic Annotation Tester results


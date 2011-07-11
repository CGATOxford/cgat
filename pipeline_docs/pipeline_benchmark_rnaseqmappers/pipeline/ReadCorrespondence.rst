===================
Read correspondence
===================

This section compares the actual locations that 
reads have mapped to.

Mapping versus quality control
==============================

The following plot shows the read quality (as number of 
colors with color quality score less than 20):

.. report:: BenchmarkReport.QCFailedVersusMatched
   :render: matrix-plot
   :palette: Blues
   :max-rows: 100
   :max-cols: 100

   Read quality versus number of numbers being
   able to map a read.

Consistency versus number of hits
=================================

The following plots show the consistency of mapping
versus the number of hits. Consistency here is the
number of mappers that manage to align a particular read.

.. report:: BenchmarkReport.MappedLocationsVersusHits
   :render: matrix-plot
   :palette: Blues
   :layout: column-3
   :width: 300

   Number uniquely matched hits versus discovery in all sets.

Read quality versus number of hits
==================================

The following plots show the read quality (as number of bases with
a quality of less than 20) versus the number of hits a read has.

.. report:: BenchmarkReport.QCFailedVersusHits
   :render: matrix-plot
   :palette: Blues
   :layout: column-3
   :width: 300
   :max-rows: 100
   :max-cols: 100

   Read quality versus number of hits

As box-plots:

.. report:: BenchmarkReport.ReadQualitiesVersusHits
   :render: box-plot
   :layout: column-4
   :width: 200

   Read quality for number of hits

Reads found in each pairwise comparison
=======================================

The following plot shows the number of reads in each pairwise
comparison that are both found.

.. report:: BenchmarkReport.MatrixFound
   :render: matrix-plot
   :transform-matrix: correspondence-analysis
   :palette: Blues

   add caption here


==================
 Feature overview
==================

This part of the pipeline examines the some basic
statistics about the segments and genes input to the analysis.

Segments
========

Distances between segments
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _SegmentsDistances:

.. report:: Trackers.SegmentsDistances
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-bins: dict
   
   Distribution of the distance between segments.

Overlap between segments
~~~~~~~~~~~~~~~~~~~~~~~~

.. _SegmentsOverlaps:

.. report:: Trackers.SegmentsOverlaps
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-bins: dict


   Distribution of the overlaps between neighbouring segments.
   Tracks without overlap appear as solid lines.

Size of segments
~~~~~~~~~~~~~~~~

.. _SegmentsSizes:

.. report:: Trackers.SegmentsSizes
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-bins: dict

   Distribution of the overlaps between neighbouring segments.
   Tracks without overlap appear as solid lines.
 
Segments
========

Distances between genes
~~~~~~~~~~~~~~~~~~~~~~~

.. _GenesDistances:

.. report:: Trackers.GenesDistances
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-bins: dict
   
   Distribution of the distance between genes.

Overlap between genes
~~~~~~~~~~~~~~~~~~~~~~~~

.. _GenesOverlaps:

.. report:: Trackers.GenesOverlaps
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-bins: dict

   Distribution of the overlaps between neighbouring genes.
   Tracks without overlap appear as solid lines.

Size of genes
~~~~~~~~~~~~~~~~

.. _GenesSizes:

.. report:: Trackers.GenesSizes
   :render: line-plot
   :transform: histogram
   :logscale: x
   :tf-aggregate: normalized-total
   :as-lines:
   :tf-bins: dict

   Distribution of the overlaps between neighbouring genes.
   Tracks without overlap appear as solid lines.
 

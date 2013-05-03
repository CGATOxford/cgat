================
Transcript rates
================

Summary Statistics
==================

.. report:: Trackers.RatesKs
   :render: table
   :transform: stats

   Distribution of ks

Plots 
==========

The following plots show the distribution of substitution rates
in transcript models.

.. report:: Trackers.RatesKs
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: @ks_range@
   :as-lines:

   Distribution of ks

Cumulative plots of substitution rates:

.. report:: Trackers.RatesKs
   :render: line-plot
   :transform: histogram
   :tf-bins: @ks_range@
   :tf-aggregate: normalized-total,cumulative
   :as-lines:

   Distribution of ks

Plots per set
=============

The following plots show the distribution of substitution rates
in transcript models grouped by track.

.. report:: Trackers.RatesKs
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: @ks_range@
   :groupby: track
   :as-lines:
   :layout: column-3

   Distribution of ks

Cumulative plots of substitution rates:

.. report:: Trackers.RatesKs
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: @ks_range@
   :groupby: track
   :as-lines:
   :layout: column-3

   Distribution of ks


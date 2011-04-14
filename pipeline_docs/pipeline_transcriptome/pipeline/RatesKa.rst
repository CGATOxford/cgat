*****************************
Rates in ancestral repeats
*****************************

Summary Statistics
==================

.. report:: Trackers.RatesKa
   :render: table
   :transform: stats

   Distribution of ka

Plots 
==========

The following plots show the distribution of substitution rates
in transcript models.

.. report:: Trackers.RatesKa
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: @ka_range@

   Distribution of ka

Cumulative plots of substitution rates:

.. report:: Trackers.RatesKa
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: @ka_range@

   Distribution of ka

Plots per set
=============

The following plots show the distribution of substitution rates
in transcript models.

.. report:: Trackers.RatesKa
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: @ka_range@
   :groupby: track

   Distribution of ka

Cumulative plots of substitution rates:

.. report:: Trackers.RatesKa
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: @ka_range@
   :groupby: track

   Distribution of ka

*****************************
Transcript constraint
*****************************

Methods
=======

Constraint is computed as the ratio of the
substitution rate in the transcript (ks) and the
median substitution rate in ancestral repeats 
surrounding the transcripts.

Summary statistics
==================

.. report:: Trackers.RatesKsKa
   :render: table
   :transform: stats

   Distribution of constraint.

Plots per slice
===============

The following plots show the distribution of constraint
in transcript models grouped by slices.

.. report:: Trackers.RatesKsKa
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total
   :tf-bins: @kska_range@
   :as-lines:
   :layout: column-2   

   Distribution of constraint.

The following are cumulative plots of constraint in transcript models grouped
by slices.

.. report:: Trackers.RatesKsKa
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: @kska_range@
   :as-lines:
   :layout: column-2   

   Distribution of constraint.

Plots per set
=============

The following plots show the distribution of constraint
in transcript models grouped by data set.

.. report:: Trackers.RatesKsKa
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: @kska_range@
   :groupby: track
   :as-lines:
   :layout: column-2   

   Distribution of constraint.

The following plots show the cumulative distribution of constraint
in transcript models grouped by data set.

.. report:: Trackers.RatesKsKa
   :render: line-plot
   :transform: histogram
   :tf-aggregate: normalized-total,cumulative
   :tf-bins: @kska_range@
   :groupby: track
   :as-lines:
   :layout: column-2   
 
   Distribution of constraint.

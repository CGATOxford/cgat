.. _medipsqcplots:

========
QC Plots
========

The Medips package provides a variety of QC plots.

Saturation analysis 
====================

According to Medips, the saturation analysis divides the total set of
available regions into two distinct random sets (A and B) of equal
size. Both sets A and B are again divided into random subsets of equal
size where the number of subsets is deter- mined by the parameter nit
(default=10). For each set, A and B , the saturation analysis
iteratively selects an increasing number of subsets and calculates
short read coverage at genome wide windows where the window sizes are
defined by the window_size parameter. In each iteration step, the
resulting genome wide coverages for the current subsets of A and B are
compared using pearson corre- lation. As the number of considered
reads increases during each iteration step, it is assumed that the
resulting genome wide coverages become more similar, a dependency that
is expressed by an increased correlation.

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: medips.dir/*_saturation.png
   :layout: column-3

   Saturation curves for all data sets.

.. _mepidscpgcoverage:

CpG coverage
============

Medips counts how many CpG dinucleotides are covered by reads.

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: medips.dir/*_pie.png
   :layout: column-3

   Number of CpG covered by certain amount of reads

The histograms below show the number of CpG that have been covered
at a given level.

.. report:: Tracker.TrackerImages
   :render: gallery-plot
   :glob: medips.dir/*_hist.png
   :layout: column-3

   Number of CpG covered by certain amount of reads



.. _deseqspikein:

=========
Spike-Ins
=========

In order to examine the power of the differential methylation
analysis, we added simulated data to the table of observed counts. The
simulated counts were derived from combining different genomic windows
randomly between different experimental conditions. The procedure
worked by permuting rows in the table of observed counts for each
experimental condition independently, keeping replicates within a
condition intact. The counts from the two experimental conditions are
then merged into a simulated counts table. Rows are selected from the
simulated counts table across a range of fold changes and average
read counts such that each bin of fold-change and average read-count
is approximately equally represented. The simulated counts are
added to the observed counts and the statistical tests for
differential tags are performed as above.

Based on the simulated counts, an empirical power can be calculated
given a fixed FDR.

Power analysis
==============

The following plots show the proportion of windows in the observed
data that fall within a certain range of FDR and power.

.. report:: DifferentialMethylation.TrackerDESeqPower
   :render: matrix-plot
   :groupby: track
   :palette: Blues
   :xtitle: power
   :ytitle: FDR
   :layout: column-3

   Power analysis. Each value is the proportion of genomic windows
   that can potentially be detected as differentially methylated
   give a desired level of FDR and desired power.

Spike statistics
================

Number of data points that have been spiked in for each data set. The
actual number is twice the one shown as data points are added in a
symmetrized fashion.

.. report:: DifferentialMethylation.TrackerDESeqSpikeIn
   :render: matrixNP-plot
   :glob: deseq.dir/*.tsv.spiked.gz
   :palette: Blues
   :xtitle: l2fold
   :ytitle: l10counts
   :groupby: track
   :layout: column-3

   Number of data points that have been spiked in.

Spike results
=============

The following plots show the percentage of spiked-in data points that
have been detected at various levels of fold-change and read count
level. The plots shown here are for an FDR level of 10%.

.. report:: DifferentialMethylation.TrackerDESeqSpikeInPercent
   :render: matrix-plot
   :groupby: track
   :palette: BrBG
   :xtitle: l2fold
   :ytitle: l10counts
   :slices: 0.1
   :layout: column-3

   Percentage of spike-in values recovered at 
   certain fold-changes and expression values for FDR=10%

Window distribution
===================

The plots below show the number of windows that fall into the same
categories as the spike-ins in terms of fold-change and read-count
level.

.. report:: DifferentialMethylation.TrackerDESeqSpikeIn
   :render: matrixNP-plot
   :glob: deseq.dir/*.tsv.unspiked.gz
   :palette: Blues
   :xtitle: l2fold
   :ytitle: l10counts
   :groupby: track
   :layout: column-3
   :logscale: z

   Number of data points that have been spiked in.


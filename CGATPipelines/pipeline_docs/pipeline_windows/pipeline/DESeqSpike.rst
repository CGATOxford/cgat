=========
Spike-Ins
=========

Spike statistics
================

Number of data poinst that have been spiked in. The actual number is
twice the one shown as data points are added in a symmetrized fashion.

.. report:: DifferentialMethylation.TrackerDESeqSpikeIn
   :render: matrixNP-plot
   :glob: deseq.dir/*.tsv.spike.gz
   :palette: Blues
   :xtitle: l2fold
   :ytitle: l10counts

   Number of data points that have been spiked in.

Spike results
=============

The following plots show the percentage of spike-ins that have been
detected at various levels of fold-change and read count level.

.. report:: DifferentialMethylation.TrackerDESeqSpikeInPercent
   :render: matrix-plot
   :groupby: track
   :palette: RdBu
   :xtitle: l2fold
   :ytitle: l10counts

   Percentage of spike-in values recovered at 
   certain fold-changes and expression values.


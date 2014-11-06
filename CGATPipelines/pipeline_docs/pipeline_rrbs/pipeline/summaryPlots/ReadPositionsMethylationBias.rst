
Methylation frequencies across the read length
==============================================

The following plots show the average methylation value for each read
position, split by cytosine context, with one plot per sample. The
purpose of these plots is to identify systematic bias in bisulphite
conversion, for example, 3' bias. Note: read positions are taken from
the methylation caller which may be instructed to ignore the first n
bases from each read, in which case the read positions will be shifted
along by n.

.. report:: SummaryPlots.readPositionMethylationBias
   :render: gallery-plot

   Plots can be downloaded using the links below each plot

Difference in methylation frequency across the read length
==========================================================

The following plot(s) presents the difference in methylation between a
pair of treatments groups. It is assumed that most CpGs will remain
unchanged between treatment pairs and the median difference will
therefore be 0 with bioligical and technical variation resulting in a
normal distribution around zero. The purpose of this plot is to
identify instances of unexpectdly high variance which may be due to
systematic technical variation, for example, sequencing errors.

.. report:: SummaryPlots.readPositionMethylationDiff
   :render: gallery-plot

   Plots can be downloaded using the links below each plot

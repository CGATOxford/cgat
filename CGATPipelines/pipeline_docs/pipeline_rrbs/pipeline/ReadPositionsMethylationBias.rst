============================================
Summary Plots - Methylation by read position
============================================

Methylation frequencies across the read length
==============================================

The following plots show the average methylation value for each read
position, split by cytosine context, with one plot per sample. The
purpose of these plots is to identify systematic bias in bisulphite
conversion, for example, 3' bias. Note: read positions are taken from
the methylation caller which may be instructed to ignore the first n
bases from each read, in which case the read positions will be shifted
along by n.

.. report:: rrbsReport.imagesTracker
   :render: gallery-plot
   :glob: plots.dir/*read_position_methylation_bias.png

   High resolution plots can be downloaded using the links below each plot

Difference in methylation frequency across the read length
==========================================================

The following plot(s) presents the difference in methylation between a
pair of treatments groups. It is assumed that most CpGs will remain
unchanged between treatment pairs and the median difference will
therefore be 0 with bioligical and technical variation resulting in a
normal distribution around zero. The purpose of this plot is to
identify instances of unexpectedly high variance which may be due to
systematic technical variation, for example, sequencing errors.

The median for each read position is plotted as a black line, with the
interquartile range (25th-75th percentile) represented by the box and
whiskers containing 99% of the data. Outliers are not shown. The read
position values are obtained by an in-silico MspI digestion of the
reference genome and will therefore not line up with the read
positions from the plots above when the methylation caller is
instructed to ignore the first n positions. Where reads do not align
with MspI fragments, the read position is donoted as "NA"

.. report:: rrbsReport.imagesTracker
   :render: gallery-plot
   :glob: plots.dir/*read_position_meth_diff_*.png 

   High resolution plots can be downloaded using the links below each plot

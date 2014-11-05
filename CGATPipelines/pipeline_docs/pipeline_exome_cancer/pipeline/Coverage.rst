========
Coverage
========

Coverage summary
==========================
The following table and plots present coverage metrics for the
samples.
The number of sequence reads per sample is presented, followed by the
percentage of uniquely mapped reads. The plots after this focus on the
percentage of bases which map to the "target" regions and the coverage
across target regions


The following table presents an overview of the coverage metrics
for each :term:`sample`. pct_aligned = Percentage reads
aligned. pct_on_target = Percentage bases in target
regions. fold_enrich = Fold enrichment of target over
background. pct_2X = Percentage bases in target covered 2X.


.. report:: Mapping.PicardTargetStats
   :render: table
   


The following plot shows the number of sequence reads per
sample. 


.. report:: Mapping.PicardTargetStats
   :render: r-ggplot
   :statement: aes(x=sample,y=reads) +
	       geom_bar(stat='identity',fill='cadetblue3') +
	       xlab('') +
	       ylab('Reads') +
	       theme(
	       axis.text.x=element_text(size=15,angle=90),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.text=element_text(size=20))


The following plot shows the percentage of uniquely aligned reads per sample

.. report:: Mapping.PicardTargetStats
   :render: r-ggplot
   :statement: aes(x=sample,y=100*pct_aligned) +
	       geom_bar(stat='identity',fill='tomato4') +
	       xlab('') +
	       ylab('Uniquely Aligned Reads (%)') +
	       theme(
	       axis.text.x=element_text(size=15,angle=90),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.text=element_text(size=20))


The following plot shows the percentage of bases within the target regions

.. report:: Mapping.PicardTargetStats
   :render: r-ggplot
   :statement: aes(x=sample,y=100*pct_on_target) +
	       geom_bar(stat='identity',fill='olivedrab4') + 
	       xlab('') +
	       ylab('Bases On-Target(%)') +
	       theme(
	       axis.text.x=element_text(size=15,angle=90),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.text=element_text(size=20))

Coverage of Target Regions
==========================
The following plot shows the mean coverage at target regions per sample

.. report:: Mapping.PicardTargetStats
   :render: r-ggplot
   :statement: aes(x=sample,y=mean_coverage) +
	       geom_bar(stat='identity',fill='salmon4') +
	       xlab('') +
	       ylab('Mean Target coverage') +
	       theme(
	       axis.text.x=element_text(size=15,angle=90),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.text=element_text(size=20))


The following plot shows the fold enrichment of target regions over background

.. report:: Mapping.PicardTargetStats
   :render: r-ggplot
   :statement: aes(x=sample,y=fold_enrich) +
	       geom_bar(stat='identity',fill='chocolate3') + 
	       xlab('') +
	       ylab('Fold Enrichment') +
	       theme(
	       axis.text.x=element_text(size=15,angle=90),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.text=element_text(size=20))



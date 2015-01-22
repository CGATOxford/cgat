========
Coverage
========

Coverage across CpGs
==========================

The following plot shows the frequency of CpGs for coverage values
between 1-100.

.. report:: Coverage.coverage
   :render: r-ggplot
   :statement: aes(x=cov,y=freq) +
	       geom_line(size=1,aes(linetype=as.factor(rep), 
	       colour=paste0(sample,'-',condition))) +
	       geom_vline(linetype = "longdash", xintercept = 10) +
	       coord_cartesian(ylim = c(0,80000), xlim = c(0,100)) +
	       ylab('Frequency') +
	       xlab('Coverage') +
	       scale_linetype_discrete(name='Replicate') +
	       scale_colour_discrete(name='Tissue - Treatment') +
	       theme(axis.text.x=element_text(size=20),
	       axis.title.x=element_text(size=20),
	       axis.text.y=element_text(size=20),
	       axis.title.y=element_text(size=20),
	       title=element_text(size=20),
	       legend.text=element_text(size=20))


   The curves show the frequency of CpGs at the the coverage value
   indicated on the x-axis. The y-axis has been trimmed such that the
   frequency of CpGs with very low coverage is not plotted. 
   The dashed black line represents a coverage threshold of 10X which is
   frequently used for making methylation calls. Colours and line
   types are used to indicate the samples and replicate number.



The following plot shows the fraction of sequence reads which would
remain across a range of coverage thresholds

.. report:: Coverage.readsRemaining
   :render: r-ggplot
   :statement: aes(x=threshold,y=percentage,group=file) +
	       geom_line(size=1,aes(linetype=as.factor(rep), 
	       colour=paste0(sample,'-',condition))) + 
	       geom_vline(linetype = "longdash", xintercept = 10) +
	       coord_cartesian(xlim = c(0,100))+
	       ylab('Fraction of sequencing remaining') + 
	       xlab('Coverage threshold') + 
	       scale_linetype_discrete(name='Replicate') +   
	       scale_colour_discrete(name='Tissue - Treatment') + 
	       theme(axis.text.x=element_text(size=20),
	       axis.title.x=element_text(size=20),
	       axis.text.y=element_text(size=20),
	       axis.title.y=element_text(size=20), 
	       title=element_text(size=20),
	       legend.text=element_text(size=20))


   The curves show the fraction of mapped reads which remain using the
   coverage value indicated on the x-axis. 
   The dashed black line represents a coverage threshold of 10X which
   is frequently used for making methylation calls. Colours and line
   types are used to indicate the samples and replicate number.




The following plot shows the number of CpGs which overlap at 10X coverage

.. report:: Coverage.CpGOverlap
   :render: r-ggplot
   :statement: aes(x=overlaps,y=CpGs) +
	       geom_bar(stat="identity")+ 
	       ylab('Number of CpGs') +
	       xlab('Number of samples') +
	       theme(axis.text.x=element_text(size=20),
	       axis.title.x=element_text(size=20),
	       axis.text.y=element_text(size=20),
	       axis.title.y=element_text(size=20), 
	       title=element_text(size=20),
	       legend.text=element_text(size=20))+
	       scale_x_continuous(breaks=seq(1,1000,1))

	     
   The bars represent the number of CpGs which overlap at 10X coverage
   for the number of samples indicated on the x-axis.

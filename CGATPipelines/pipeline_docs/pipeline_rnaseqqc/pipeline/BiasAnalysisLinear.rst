===================================================
Bias analysis results - linear regression - Summary
===================================================

This page presents an overview of the analysis of potential biasing
factors using linear regression. The plot aesthetics are split by the
first identier, e.g tissue.

Genes/transcripts are binned according to their value for each
potential biasing factor (e.g GC content), with each bin containing an
equal number of genes/transcripts.  The mean expression for the
genes/transcripts for each sample is then calculated for each bin. The
relationship between the potential biasing factor and expression level
is explored by computing the Spearman rank correlation and linear
regression. The gradient of the linear regression and rho value of the
correlation are shown below



Summary plots
=========================

.. report:: RnaseqqcReport.CorrelationSummaryGC
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.

.. report:: RnaseqqcReport.CorrelationSummaryA
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.

.. report:: RnaseqqcReport.CorrelationSummaryT
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.

.. report:: RnaseqqcReport.CorrelationSummaryC
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.

.. report:: RnaseqqcReport.CorrelationSummaryG
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=as.factor(sample),
	       colour=as.factor(variable), group=as.factor(variable))+geom_line()+
	       scale_colour_discrete(name=guide_legend(title='biasfactor'))+
	       xlab('')+ylab('Correlation')+
	       theme(axis.text.x=element_text(size=15,hjust=1,angle=90),
	       axis.text.y=element_text(size=15),title=element_text(size=15),
	       legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.
    
.. report:: RnaseqqcReport.GradientSummaryGC
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.

.. report:: RnaseqqcReport.GradientSummaryA
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.


.. report:: RnaseqqcReport.GradientSummaryT
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.

.. report:: RnaseqqcReport.GradientSummaryC
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.

.. report:: RnaseqqcReport.GradientSummaryG
   :render: r-ggplot
   :statement:  aes(y=as.numeric(value), x=as.factor(sample),
		colour=as.factor(variable), group=as.factor(variable))+
		geom_line()+
		scale_colour_discrete(name = guide_legend(title='biasfactor'))+
		xlab('')+ 
		ylab('Gradient')+
		theme(axis.text.x=element_text(size=15,angle=90,hjust=1),
		axis.text.y=element_text(size=15),title=element_text(size=15),
		legend.text=element_text(size=15))

   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.

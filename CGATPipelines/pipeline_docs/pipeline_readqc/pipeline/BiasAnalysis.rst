==============
Bias analysis results
==============

This page presents the results of analysis of potential biasing
factors. 


Summary plots
=========================

.. report:: ReadqcReport.CorrelationSummary
   :render: r-ggplot
   :transform: melt,toframe
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice),group=as.factor(Slice))+geom_line()+scale_colour_discrete(name = guide_legend(title = 'bias factor'))+xlab('')+ylab('Correlation')+theme(axis.text.x=element_text(size=15,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))


   Correlation between gene expression and potential biasing factors
   across all samples.
    
.. report:: ReadqcReport.GradientSummary
   :render: r-ggplot
   :transform: melt,toframe
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice),group=as.factor(Slice))+geom_line()+scale_colour_discrete(name = guide_legend(title = 'bias factor'))+xlab('')+ylab('Gradient')+theme(axis.text.x=element_text(size=15,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))


   Gradient of linear regression between gene expression and potential 
   biasing factors across all samples.

Individual bias factor plots
============================


.. report:: Status.GCContentSummary
   :render: r-ggplot
   :transform: melt
   :statement:
      aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name
      = guide_legend(title = 'Sample'))+xlab('GC Content')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned GC content for each sample

.. report:: Status.LengthSummary
   :render: r-ggplot
   :transform: melt
   :statement:
      aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name
      = guide_legend(title = 'Sample'))+xlab('Log2 length (bp)')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned length for each sample


.. report:: Status.AASummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('AA frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned AA dinucleotide frequency for each sample

.. report:: Status.ATSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('AT frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned AC dinucleotide frequency for each sample

.. report:: Status.ACSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('AC frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned AC dinucleotide frequency for each sample

.. report:: Status.AGSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('AG frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned AG dinucleotide frequency for each sample


.. report:: Status.TASummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('TA frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned TA dinucleotide frequency for each sample

.. report:: Status.TTSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('TT frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned TT dinucleotide frequency for each sample

.. report:: Status.TCSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('TC frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned TC dinucleotide frequency for each sample

.. report:: Status.TGSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('TG frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned TG dinucleotide frequency for each sample


.. report:: Status.CASummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('CA frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned CA dinucleotide frequency for each sample

.. report:: Status.CTSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('CT frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned CT dinucleotide frequency for each sample

.. report:: Status.CCSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('CC frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned CC dinucleotide frequency for each sample

.. report:: Status.CGSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('CG frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned CG dinucleotide frequency for each sample


.. report:: Status.GASummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('GA frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned GA dinucleotide frequency for each sample

.. report:: Status.GTSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('GT frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned GT dinucleotide frequency for each sample

.. report:: Status.GCSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('GC frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned GC dinucleotide frequency for each sample

.. report:: Status.GGSummary
   :render: r-ggplot
   :transform: melt
   :statement: aes(y=as.numeric(Data),x=as.factor(Track),colour=as.factor(Slice))+geom_point()+stat_smooth(aes(group=Slice,colour=Slice),method=lm,se=F)+scale_colour_discrete(name = guide_legend(title = 'Sample'))+xlab('GG frequency')+ylab('Log2 Mean Expression')+theme(axis.text.x=element_text(size=10,angle=90),axis.text.y=element_text(size=15),title=element_text(size=15),legend.text=element_text(size=15))

   Mean expression across binned GG dinucleotide frequency for each sample

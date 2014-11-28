=====================================================================
Bias analysis results - Split by first identifier
=====================================================================

This page presents the analysis of potential biasing factors using
linear regression. The plot aesthetics are split by the first
identier, e.g tissue.

Genes/transcripts are binned according to their value for each
potential biasing factor (e.g GC content), with each bin containing an
equal number of genes/transcripts.  The mean expression for the
genes/transcripts for each sample is then calculated for each
bin. This mean expression is plotted below, along with a linear
regression for each sample.


GC content plots
================

.. report:: RnaseqqcReport.GCContentSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GC_Content, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GC Content (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned GC content for each sample. Linear regression.

.. report:: RnaseqqcReport.GCContentSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GC_Content, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GC Content (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned GC content for each sample. Local
   regression.


Length plots
============

.. report:: RnaseqqcReport.LengthSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=length, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('length (Log 2 bp)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned length for each sample. Linear regression.

.. report:: RnaseqqcReport.LengthSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=length, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('length (Log 2 bp)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned length for each sample. Local
   regression.


AA dinucleotide plots
=====================

.. report:: RnaseqqcReport.AASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AA, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AA dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.AASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AA, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AA dinucleotides for each
   sample. Local regression.


AT dinucleotide plots
=====================

.. report:: RnaseqqcReport.ATSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AT, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AT dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.ATSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AT, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AT dinucleotides for each
   sample. Local regression.


AC dinucleotide plots
=====================

.. report:: RnaseqqcReport.ACSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AC, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AC dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.ACSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AC, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AC dinucleotides for each
   sample. Local regression.

AG dinucleotide plots
=====================

.. report:: RnaseqqcReport.AGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AG, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AG dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.AGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AG, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AG dinucleotides for each
   sample. Local regression.

TA dinucleotide plots
=====================

.. report:: RnaseqqcReport.TASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TA, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TA dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.TASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TA, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TA dinucleotides for each
   sample. Local regression.

TT dinucleotide plots
=====================

.. report:: RnaseqqcReport.TTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TT, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TT dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.TTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TT, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TT dinucleotides for each
   sample. Local regression.

TC dinucleotide plots
=====================

.. report:: RnaseqqcReport.TCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TC, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TC dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.TCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TC, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TC dinucleotides for each
   sample. Local regression.

TG dinucleotide plots
=====================

.. report:: RnaseqqcReport.TGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TG, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TG dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.TGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TG, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TG dinucleotides for each
   sample. Local regression.

CA dinucleotide plots
=====================

.. report:: RnaseqqcReport.CASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CA, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CA dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.CASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CA, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CA dinucleotides for each
   sample. Local regression.

CT dinucleotide plots
=====================

.. report:: RnaseqqcReport.CTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CT, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CT dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.CTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CT, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CT dinucleotides for each
   sample. Local regression.

CC dinucleotide plots
=====================

.. report:: RnaseqqcReport.CCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CC, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CC dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.CCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CC, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CC dinucleotides for each
   sample. Local regression.

CG dinucleotide plots
=====================

.. report:: RnaseqqcReport.CGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CG, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CG dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.CGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CG, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CG dinucleotides for each
   sample. Local regression.

GA dinucleotide plots
=====================

.. report:: RnaseqqcReport.GASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GA, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GA dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.GASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GA, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GA dinucleotides for each
   sample. Local regression.

GT dinucleotide plots
=====================

.. report:: RnaseqqcReport.GTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GT, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GT dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.GTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GT, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GT dinucleotides for each
   sample. Local regression.

GC dinucleotide plots
=====================

.. report:: RnaseqqcReport.GCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GC, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GC dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.GCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GC, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GC dinucleotides for each
   sample. Local regression.

GG dinucleotide plots
=====================

.. report:: RnaseqqcReport.GGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GG, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GG dinucleotides for each
   sample. Linear regression.

.. report:: RnaseqqcReport.GGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GG, colour=id_1)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_1),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='First Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GG dinucleotides for each
   sample. Local regression.

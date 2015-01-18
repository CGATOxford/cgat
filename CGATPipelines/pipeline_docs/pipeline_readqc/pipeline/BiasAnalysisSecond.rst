=====================================================================
Bias analysis results - Split by second identifier
=====================================================================

This page presents the analysis of potential biasing factors using
linear regression. The plot aesthetics are split by the second
identifer, e.g treatment

Genes/transcripts are binned according to their value for each
potential biasing factor (e.g GC content), with each bin containing an
equal number of genes/transcripts.  The mean expression for the
genes/transcripts for each sample is then calculated for each
bin. This mean expression is plotted below, along with a linear
regression for each sample.


GC content plots
================

.. report:: Status.GCContentSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GC_Content, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GC Content (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned GC content for each sample. Linear regression.

.. report:: Status.GCContentSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GC_Content, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.LengthSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=length, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('length (Log 2 bp)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned length for each sample. Linear regression.

.. report:: Status.LengthSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=length, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.AASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AA, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AA dinucleotides for each
   sample. Linear regression.

.. report:: Status.AASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AA, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.ATSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AT, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AT dinucleotides for each
   sample. Linear regression.

.. report:: Status.ATSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AT, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.ACSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AC, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AC dinucleotides for each
   sample. Linear regression.

.. report:: Status.ACSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AC, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.AGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AG, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('AG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage AG dinucleotides for each
   sample. Linear regression.

.. report:: Status.AGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=AG, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.TASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TA, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TA dinucleotides for each
   sample. Linear regression.

.. report:: Status.TASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TA, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.TTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TT, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TT dinucleotides for each
   sample. Linear regression.

.. report:: Status.TTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TT, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.TCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TC, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TC dinucleotides for each
   sample. Linear regression.

.. report:: Status.TCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TC, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.TGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TG, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('TG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage TG dinucleotides for each
   sample. Linear regression.

.. report:: Status.TGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=TG, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.CASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CA, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CA dinucleotides for each
   sample. Linear regression.

.. report:: Status.CASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CA, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.CTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CT, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CT dinucleotides for each
   sample. Linear regression.

.. report:: Status.CTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CT, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.CCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CC, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CC dinucleotides for each
   sample. Linear regression.

.. report:: Status.CCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CC, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.CGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CG, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('CG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage CG dinucleotides for each
   sample. Linear regression.

.. report:: Status.CGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=CG, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.GASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GA, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GA (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GA dinucleotides for each
   sample. Linear regression.

.. report:: Status.GASummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GA, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.GTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GT, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GT (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GT dinucleotides for each
   sample. Linear regression.

.. report:: Status.GTSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GT, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.GCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GC, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GC (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GC dinucleotides for each
   sample. Linear regression.

.. report:: Status.GCSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GC, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
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

.. report:: Status.GGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GG, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=lm,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GG dinucleotides for each
   sample. Linear regression.

.. report:: Status.GGSummary
   :render: r-ggplot
   :statement: aes(y=as.numeric(value), x=GG, colour=id_2)+
	       geom_point()+
	       stat_smooth(aes(group=variable,colour=id_2),method=loess,se=F)+
	       scale_colour_discrete(name=guide_legend(title='Second Identifier'))+
	       scale_y_continuous(limits=c(0,1))+
	       xlab('GG (Fraction)')+
	       ylab('Normalised Expression(Nominal scale)')+
	       theme(axis.text.x=element_text(size=10,angle=90),
	       axis.text.y=element_text(size=15),
	       title=element_text(size=15),
	       legend.text=element_text(size=15))

   Mean expression across binned percentage GG dinucleotides for each
   sample. Local regression.


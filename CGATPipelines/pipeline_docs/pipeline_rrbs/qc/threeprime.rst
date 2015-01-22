=================
Sequence 3' bases
=================

The following plot shows the per sample frequency of three 3' bases.
Note: For paired end samples, only the first read pair is shown

.. report:: Coverage.seqStart
   :render: r-ggplot
   :statement: aes(x=file,fill=sequence,group=1) +
	       geom_bar(position="fill",stat="identity",
	       colour="white",aes(y=reads))+
	       ylab('Fraction') +
	       scale_colour_discrete(label="3' sequence") +
	       theme(axis.text.x=element_text(size=15,angle=90),
	       axis.title.x=element_text(size=20),
	       axis.text.y=element_text(size=20),
	       axis.title.y=element_text(size=20),
	       title=element_text(size=20),
	       legend.text=element_text(size=20))


   3' bases frequency plot. For single end directional RRBS using MspI
   digestion, reads are expected to start with CGG or TGG. For paired
   end non-directional RRBS, reads are expected to start with CGG,
   TGG, CGA or TGA. All remaining sequences are shown as :term:`other`.




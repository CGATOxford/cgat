========
Mutational Signature
========

Mutational Signature
==========================
The following plot presents the mutation signature as represented by
the relative frequency of the 12 possible SNPs with regards to
nucleotide change.


The following plot presents the muatation signature for each :term:`sample`


.. report:: Signature.SnpSummary
   :render: r-ggplot
   :statement: aes(x=patient_id,y=frequency) +
	       geom_bar(position='fill', stat='identity',aes(fill =
	       paste0(ref,">",alt), order=paste0(ref,">",alt))) +
	       xlab('Patient ID') +
	       ylab('Nucleotide change frequency') +
	       scale_fill_discrete(name = '') +
	       theme(
	       axis.text.x=element_text(size=15,angle=90),
	       axis.text.y=element_text(size=15),
	       axis.title.y=element_text(size=15),
	       legend.text=element_text(size=20))


The following table presents the numbers of snp variants called per
base change

.. report:: Signature.SnpSummaryTable
   :render: table
   :large: xls



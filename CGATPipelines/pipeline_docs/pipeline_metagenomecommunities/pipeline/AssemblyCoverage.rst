.. _coverage:



===================
Assembly coverage
===================

To further assess the quality of the assembly we look at the raw read coverage over the
contig assembly. This provides an indication of the coverage of the contigs over the sequenced
metagenome. It also provides potential information on whether it is possible to cluster contigs
based on contig coverage.


All picard statistics
======================

.. report:: AssemblyCoverage.AlignmentSummary
   :render: table


   Picard alignment summary statistics


.. report:: AssemblyCoverage.InsertSizeSummary
   :render: table


   Picard insert size summary statistics


Percent reads aligned
======================

If the assembly was succesfull, the majority of reads in the original sample
should map to the resulting contigs.

.. report:: AssemblyCoverage.ReadsAligned
   :render: table

   percent reads aligned to contigs

.. report:: AssemblyCoverage.ReadsAligned
   :render: interleaved-bar-plot

   percent reads aligned to congtigs

.. report:: AssemblyCoverage.ReadsAlignedInPairs
   :render: table

   percent reads aligned in pairs to contigs

.. report:: AssemblyCoverage.ReadsAlignedInPairs
   :render: interleaved-bar-plot

   percent reads aligned in pairs to contigs

.. report:: AssemblyCoverage.MismatchRate
   :render: table

   mismatch rate

.. report:: AssemblyCoverage.MismatchRate
   :render: interleaved-bar-plot

   mismatch rate


Coverage over contigs
======================

We use Bowtie to assess the coverage over contigs. The average (mean) coverage distributions
are plotted below. This gives an indication as to the number of species that are present if there
is sufficient variablity in abundance between species present.


.. report:: AssemblyCoverage.CoverageMean
   :render: line-plot
   :transform: histogram

   Histogram mean coverage over contigs



Variability in coverage
========================

The presence of highly related species within a metagenomic sample leads to the sequencing of 
highly conserved genomic segments. This can often result in the assembly of chimeric contigs 
whereby reads derived from multiple species are assembled together. On e way in which a chimera
can be assessed without knowing the origin of each read is to assess the uniformity of coverage
over each contig. It is not feasible to plot the coverage over each base for ever individual contig
so we attempt to proxy this by looking at the distribution of the standard deviation of base
coverage over the contigs. This may reveal populations of contigs that show a higher than expected
standard deviation.


.. report:: AssemblyCoverage.CoverageSd
   :render: line-plot
   :transform: histogram

   Histogram of standard deviations of coverage over contigs








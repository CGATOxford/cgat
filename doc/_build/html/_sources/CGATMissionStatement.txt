======================================
Background
======================================

Mission statement
=================

The CGAT code collection has been written over several years in
the context of comparative genomics and more recently next-generation
sequencing analysis.

The aim of the toolkit is to solve practical problems in the analysis
of genomic data. The focus of the toolkit is to facilitate the
interpretation of genomic data in a biological context. Furthermore,
as a training institution our aim is to write code that is well
structured and can serve as an introduction to advanced bioinformatic
scripting for biologists.

Other toolkits with similar functionality
=========================================

The CGAT code collection extends, complements but also overlaps
various other toolkits. As all toolkits, and ours, continue to evolve,
this is a very dynamic relationship. For example, our workflows frequently
use other toolkits, in particular bedtools_ and the `UCSC tools`_, for
high-performance computations. Usage of common genomic file formats
and a command line interface ensures compatibility. 

Below is a list of toolkits with similar or complementarity
functionality to the CGAT code collection and quotes from their
respective web-sites:

* bedtools_
  The BEDTools utilities allow one to address common genomics tasks such
  as finding feature overlaps and computing coverage. The utilities are
  largely based on four widely-used file formats: BED, GFF/GTF, VCF, and
  SAM/BAM. Using BEDTools, one can develop sophisticated pipelines that
  answer complicated research questions by "streaming" several BEDTools
  together.

* samtools_
  SAM Tools provide various utilities for manipulating alignments in
  the SAM format, including sorting, merging, indexing and generating
  alignments in a per-position format.

* `UCSC tools`_
  `Jim Kent's <http://users.soe.ucsc.edu/~kent/>`_ genomic utilities
  for working with genomic features and alignments.

* EMBOSS_
  EMBOSS is "The European Molecular Biology Open Software Suite". EMBOSS
  is a free Open Source software analysis package specially developed
  for the needs of the molecular biology (e.g. EMBnet) user
  community. The software automatically copes with data in a variety of
  formats and even allows transparent retrieval of sequence data from
  the web. Also, as extensive libraries are provided with the package,
  it is a platform to allow other scientists to develop and release
  software in true open source spirit. EMBOSS also integrates a range of
  currently available packages and tools for sequence analysis into a
  seamless whole. EMBOSS breaks the historical trend towards commercial
  software packages. 

* GROK_
  GROK (Genomic Region Operation Toolkit) is "Swiss Army knife" library
  for processing genomic interval data. GROK operates on genomic
  regions, annotated chromosomal intervals that represent sequencing
  short reads, gene locations, ChIP-seq peaks or other genomic
  features. Applications of GROK include file format conversions, set
  operations, overlap queries, and filtering and transformation
  operations. Supported file formats include BAM/SAM, BED, BedGraph,
  CSV, FASTQ, GFF/GTF, VCF and Wiggle. 

* biopieces_
  The Biopieces are a collection of bioinformatics tools that can be
  pieced together in a very easy and flexible manner to perform both
  simple and complex tasks. The Biopieces work on a data stream in such
  a way that the data stream can be passed through several different
  Biopieces, each performing one specific task: modifying or adding
  records to the data stream, creating plots, or uploading data to
  databases and web services.

* fastx-toolkit_
  The FASTX-Toolkit is a collection of command line tools for
  Short-Reads FASTA/FASTQ files preprocessing. 

.. _EMBOSS: http://emboss.sourceforge.net/
.. _GROK: http://csbi.ltdk.helsinki.fi/grok/
.. _biopieces: https://code.google.com/p/biopieces/
.. _fastx-toolkit: http://hannonlab.cshl.edu/fastx_toolkit/

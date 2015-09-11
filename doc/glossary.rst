========
Glossary
========

File formats
=============

.. glossary::

    yaml
    	Language to serialize objects. Used in the CGAT testing
    	framework. (`YAML <http://en.wikipedia.org/wiki/YAML>`_).

    bam
        Format to store genomic alignments in a compressed format.
	(`BAM <http://samtools.sourceforge.net/>`_).

    bed
	File containing genomic intervals. 
	(`BED <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_).
	
    vcf
        `Variant call format
        <http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41>`_.

    gtf
	`General transfer format
        <http://www.ensembl.org/info/website/upload/gff.html>`_.
	Format to store genes and transcripts.
    
    gff
	`General feature format <http://www.ensembl.org/info/website/upload/gff.html>`_.

    bigwig
        Compressed format for displaying numerical values across
        genomic ranges (`BIGWIG <https://genome.ucsc.edu/goldenPath/help/bigWig.html>`_).

    fasta
        Sequence format. 

    wiggle
        Format for displaying numerical values across genomic
        ranges (`Wiggle <https://genome.ucsc.edu/goldenPath/help/wiggle.html>`_).

    psl  
    	Genomic alignment format. The format is described in detail
	`(PSL <https://genome.ucsc.edu/FAQ/FAQformat.html#format2>`_.

    sam
        Format to store genomic alignments
	`(SAM <http://samtools.sourceforge.net/>`_).
	
    gdl
        gdl

    tsv
        Tab separated values. In these tables, records are separated by new-line
        characters and fields by tab characters. Lines with comments are started
        by the ``#`` character and are ignored. The first uncommented line
        should contain the column headers. For example::

	  # This is a comment
	  gene_id	length
	  gene1	1000
	  gene2	2000
	  # Another comment

    svg
        pass

    edge list
        pass

    fastq
        Sequence format containing quality scores, more background is
	`here <http://en.wikipedia.org/wiki/FASTQ_format>`_

    sra
        sra

    axt
        axt

    maf
        maf
   
    agp
        `AGP format <https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/>`_

    rdf
        `Resource description framework <http://en.wikipedia.org/wiki/Resource_Description_Framework>`_

Other terms
===========

.. glossary::

    test directory
        Directory that contains the :file:`test.yaml`, input and
	reference files for testing scripts.
	
    experiment
        experiment

    replicate
        replicate

    graph
	graph

    track
        track

    graph
	graph

    submit host
        pass

    execution host 
        pass

    edge list
        pass

    task
        pass

    sphinxreport
        sphinxreport

    query
        pass

    target
        pass

    code directory
       pass

    go
       pass

    goslim
	pass

    fastq
        pass

    tss
        Transcription start site

    production pipeline
        A pipeline that performs common tasks on a certain type of
        data. The idea of a production pipeline is to provide common
       	preprocessing of data and a first look. A :term:`project
        pipeline` might then take data from one or more
        :term:`production pipeline` to glean biological insight.

    project pipeline
        A pipeline that is project specific. Usually code is developed
	first inside a project pipeline. When it becomes generally
        useful, it may be refactored into a production pipeline.
	 
    stdin
        Unix standard input. Most CGAT tools read data from stdin.

    stdout
        Unix standard output. Most CGAT tools output data to stdout.

    stderr
        Unix standard error. This is where errors go.
  
    loglevel
        Verbosity of logging information. The logging level can be
        determined by the ``--verbose`` option. A
	level of ``0`` means no logging output, while ``1`` is information
	messages only, while ``2`` outputs also debugging information.

 

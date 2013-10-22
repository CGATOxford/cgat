==============
Tool reference
==============

This page summarizes prominent tools within the CGAT Code
collection. The tools are grouped losely by functionality.

Genomic intervals/features
==========================

:doc:`scripts/beds2counts`
    Compute overlap statistics of multiple :term:`bed` files.

:doc:`scripts/bed2fasta`
    Transform interval data in a :term:`bed` formatted file into a
    :term:`fasta` formatted file of sequence data.

:doc:`scripts/bed2gff`
    Convert between interval data. Convert a :term:`bed` formatted
    file to a :term:`gff` or :term:`gtf` formatted file.

:doc:`scripts/gff2gff`
    Work on :term:`gff` formatted files with genomic features. This 
    tools sorts/renames feature files, reconciles chromosome names,
    and more.

:doc:`scripts/bed2bed`
    Filter or merge interval data in a :term:`bed` formatted file.

:doc:`scripts/bed2graph`
    Compare two sets of genomic intervals and output a list of
    overlapping features.

:doc:`scripts/bed2stats`
    Compute summary statistics of genomic intervals.

:doc:`scripts/bed2table`
    Annotate genomic intervals (composition, peak location, overlap, ...)

:doc:`scripts/beds2beds`
    Decompose multiple sets of genomic intervals into various
    intersections and unions.

:doc:`scripts/diff_bed`
    Compare multiple sets of interval data sets. The tools computes
    all-vs-all pairwise overlap summaries. Permits incremental updates
    of similarity table.

:doc:`scripts/gff2bed`
    Convert between formats

:doc:`scripts/split_gff`
    Split a file in gff format into smaller files. The script ensures
    that overlapping intervals remain in the same file.

:doc:`scripts/gff2coverage`
   This script computes the genomic coverage of intervals 
   in a :term:`gff` formatted file. The coverage is computed
   per feature.

:doc:`scripts/gff2fasta`
    Output genomic sequences from intervals.

:doc:`scripts/gff2histogram`
    Compute distributions of interval sizes, intersegmental distances
    and interval ovelap from list of intervals.

:doc:`scripts/gff2stats`
    Summarize features within a :term:`gff` formatted file.	

:doc:`scripts/gff2psl`
    Convert between formats.

Gene sets
=========

:doc:`scripts/gtf2gff`
    Translate a gene set into genomic annotations such as introns,
    intergenic regions, regulatory domains, etc.

:doc:`scripts/gtf2table`
    Annotate transcripts in a :term:`gtf` formatted file. Annotations
    can be in reference to a second gene set (fragments, extensions), 
    aligned reads (coverage, intron overrun, ...) or densities.

:doc:`scripts/gtf2fasta`
    Annotate each base in the genome according to its use within
    a transcript. Outputs lists of junctions.

:doc:`scripts/gtf2gff`
    Derive genomic intervals (intergenic regions, introns) from
    a gene set.

:doc:`scripts/gtf2gtf`
    merge exons/transcripts/genes, filter transcripts/genes, rename
    transcripts/genes, ...

:doc:`scripts/gtf2tsv`
    convert gene set in :term:`gtf` format to tabular format.

:doc:`scripts/gtfs2tsv`
    Compare two gene sets - output common and unique lists of genes.

:doc:`scripts/diff_gtf` 
    Compare multiple gene sets. The tools computes all-vs-all pairwise
    overlap of exons, bases and genes. Permits incremental updates of
    similarity table.

Sequence data
=============

:doc:`scripts/fastqs2fasta`
    Interleave paired reads from two fastq files into a single fasta file.

:doc:`scripts/index_fasta`
    Build an index for a fasta file. Pre-requisite for many CGAT tools.

:doc:`scripts/fasta2kmercontent`
    Count kmer content in a set of :term:`fasta` sequences.

:doc:`scripts/fasta2table`
    Compute features of sequences in :term:`fasta` formatted files

:doc:`scripts/diff_fasta`
    Compare two sets of sequences. Outputs missing, identical
    and fragmented sequences.

:doc:`scripts/fasta2bed`
    Segment sequences based on G+C content, gaps, ...

:doc:`scripts/fastas2fasta`
    Concatentate sequences from multiple files.

:doc:`scripts/fasta2variants`
    In-silico creation of variants of protein coding
    sequences.

NGS data
========

:doc:`scripts/bam2geneprofile`
    Compute meta-gene profiles from aligned reads in a :term:`bam`
    formatted file. Also accepts :term:`bed` or :term:`bigwig`
    formatted files.

:doc:`scripts/bam2bam`
    Operate on :term:`bam` formatted files - filtering, stripping, 
    setting flags.

:doc:`scripts/bam2bed`
    Convert :term:`bam` formatted file of genomic alignments
    into genomic intervals. Permits merging of paired read data
    and filtering by insert-size.

:doc:`scripts/bam2fastq`
    Save sequence and quality information from a :term:`bam` 
    formatted file.

:doc:`scripts/bam2peakshape`
    Compute read densities over a collection of intervals. Also 
    accepts :term:`bed` or :term:`bigwig` formatted files.

:doc:`scripts/bam2stats`
    Compute summary statistics of a :term:`bam` formatted file.

:doc:`scripts/bam2wiggle`
    Convert read coverage in a :term:`bam` formatted file into
    a :term:`wiggle` or :term:`bigwig` formatted file.

:doc:`scripts/bam_vs_gtf`
    Compute stats on exon over-/underrun and spliced reads.

:doc:`scripts/bam_vs_bed`
    Compute coverage of reads within multiple interval types.

:doc:`scripts/bam_vs_bam`
    Outputs side-by-side comparison of residue level counts
    between multiple :term:`bam` formatted files.
	 
:doc:`scripts/fastq2fastq`
    Perform quality score conversion between :term:`fastq` 
    formatted files.

:doc:`scripts/fastqs2fasta`
    Interleave paired end data.

:doc:`scripts/fastq2table`
    Output bases below quality threshold, number of N's, quality score distribution.    

:doc:`scripts/fastqs2fastqs`
    Ensure that paired read :term:`fastq` formatted files are consistent
    after filtering on the individual files.

:doc:`scripts/diff_bam`
    Perform read-by-read comparison of two bam-files.

Variants
========

:doc:`scripts/vcf2vcf`
    Sort a vcf file.

Genomics
========

:doc:`scripts/diff_chains`
    How many residues to the same locations, do different locations,
    etc.

:doc:`scripts/chain2stats`
    Output coverage statistics for a UCSC liftover chain file.





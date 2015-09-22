.. _scripts:

Scripts
=======

This document contains all the scripts for/by CGAT.
Scripts are written to be called from the command line.

Genomics
--------

.. toctree::
   :maxdepth: 1
   
   scripts/bam2geneprofile.rst
   scripts/bed2bed.rst
   scripts/bed2gff.rst

   scripts/gff2bed.rst
   scripts/psl2assembly.rst
   scripts/psl2psl.rst
   scripts/psl2gff.rst
   scripts/psl2map.rst
   scripts/psl2wiggle.rst
   scripts/psl2wiggle_stats.rst
   scripts/diff_gtf.rst
   scripts/fasta2gff.rst
   scripts/gff2psl.rst
   scripts/gff2coverage.rst
   scripts/gff2fasta.rst
   scripts/gtf2alleles.rst
   scripts/gff2gff.rst
   scripts/gff2histogram.rst
   scripts/gff2plot.rst
   scripts/gff2stats.rst
   scripts/gff2view.rst
   scripts/gff_decorate.rst
   scripts/gtf2gff.rst
   scripts/gtf2gtf.rst
   scripts/gtf2reads.rst
   scripts/gtf2table.rst
   scripts/gtfs2graph.rst
   scripts/gtf2fasta.rst
   scripts/snp2table.rst
   scripts/psl2fasta.rst
   scripts/psl2stats.rst
   scripts/psl2table.rst
   scripts/maq2assembly.rst
   scripts/maq2psl.rst
   scripts/analyze_readpositions.rst
   scripts/combine_gff.rst
   scripts/quality2fasta.rst
   scripts/bam2wiggle.rst
   scripts/bed2annotator.rst
   scripts/bed2graph.rst
   scripts/bed2psl.rst
   scripts/chain2psl.rst
   scripts/diff_bed.rst
   scripts/fasta2bed.rst
   scripts/gff2table.rst
   scripts/index2gff.rst
   scripts/maf2psl.rst
   scripts/maq2psl.rst
   scripts/snp2counts.rst
   scripts/snp2maf.rst
   scripts/psl2chain.rst

Trees
-----

.. toctree::
   :maxdepth: 1

   scripts/tree2matrix.rst
   scripts/tree2plot.rst
   scripts/tree2stats.rst
   scripts/tree2taxa.rst
   scripts/tree2tree.rst
   scripts/tree2svg.rst
   scripts/tree_collapse_species.rst
   scripts/tree_diff.rst
   scripts/tree_strain2species.rst
   scripts/trees2sets.rst
   scripts/trees2tree.rst
   scripts/trees2trees.rst
   scripts/tree_species2genes.rst
   scripts/matrix2tree.rst
   scripts/fasta2nj.rst
   scripts/tree2patterns.rst
   scripts/tree_species2genes.rst

Alignment
---------

.. toctree::
   :maxdepth: 1

   scripts/align_transcripts.rst
   scripts/peptides2cds.rst

Visualization
-------------

.. toctree::
   :maxdepth: 1

   scripts/plot_data.rst
   scripts/plot_histogram.rst
   scripts/plot_matrix.rst
   scripts/png2svg.rst
   scripts/go2svg.rst

Graphs
------

.. toctree::
   :maxdepth: 1

   scripts/graph2stats.rst
   scripts/graph_combine_links_redundant.rst
   scripts/graph_filter_links_redundant.rst
   scripts/graph_links2gdl.rst
   scripts/compare_clusters.rst

Sequences and rates
-------------------

.. toctree::
   :maxdepth: 1

   scripts/align_pairs.rst
   scripts/align_all_vs_all.rst
   scripts/jalview.rst
   scripts/mali2bootstrap.rst
   scripts/mali2mali.rst
   scripts/mali2malis.rst
   scripts/mali2summary.rst
   scripts/malis2mali.rst
   scripts/malis2malis.rst
   scripts/malis2profiles.rst
   scripts/sequence2alignment.rst
   scripts/sequences2mali.rst
   scripts/mask_fasta.rst
   scripts/index_fasta.rst
   scripts/mali2kaks.rst
   scripts/mali2rates.rst
   scripts/malis2masks.rst
   scripts/diff_fasta.rst
   scripts/rates2rates.rst
   scripts/introns2rates.rst
   scripts/nr2table.rst
   scripts/extractseq.rst
   scripts/quality2masks.rst
   scripts/map_residues.rst

Matrices and Tables
-------------------

.. toctree::
   :maxdepth: 1

   scripts/combine_histograms.rst
   scripts/csvs2csv.rst
   scripts/csv2csv.rst
   scripts/csv2db.rst
   scripts/csv2xls.rst
   scripts/csv_cut.rst
   scripts/csv_intersection.rst
   scripts/csv_rename.rst
   scripts/csv_set.rst
   scripts/csv_uniq.rst
   scripts/cat_tables.rst
   scripts/matrix2matrix.rst
   scripts/matrix2stats.rst
   scripts/sparse2full.rst
   scripts/filter_tokens.rst
   scripts/table2graph.rst
   scripts/table2table.rst
   scripts/join_tables.rst
   scripts/substitute_tokens.rst

Stats
-----

.. toctree::
   :maxdepth: 1

   scripts/compare_histograms.rst
   scripts/data2roc.rst
   scripts/data2stats.rst
   scripts/data2bins.rst
   scripts/data2histogram.rst
   scripts/data2multiple_anova.rst
   scripts/histogram2histogram.rst
   scripts/merge_tables.rst
   scripts/modify_table.rst
   scripts/r_compare_distributions.rst
   scripts/r_mann_whitney_u.rst
   scripts/r_table2scatter.rst
   scripts/r_test.rst
   scripts/calculate_histogram_2D.rst
   scripts/nmf.rst
   scripts/simulate_function.rst
   scripts/histograms2kl.rst
 
Tools
-----

Cluster and jobs
++++++++++++++++

.. toctree::
   :maxdepth: 1

   scripts/clean.rst
   scripts/split_file.rst

Other
++++++

.. toctree::
   :maxdepth: 1

   scripts/rename_links.rst
   scripts/cgat_script_template.rst
   scripts/convert_time2seconds.rst
   scripts/set_diff.rst

.. Pipelines
.. ---------

.. The following sections lists scripts that are specific to certain
.. pipelines.

.. Gene prediction pipeline
.. +++++++++++++++++++++++++

.. Scripts that are specific for the `gpipe` gene prediction pipeline.

.. .. toctree::
..    :maxdepth: 1


..    scripts/select_transcripts.rst
..    scripts/select_transcripts_per_gene.rst
..    scripts/gff_ensembl2gbrowser.rst
..    scripts/psl2predictions.rst
..    scripts/gff2exons.rst
..    scripts/regions2gff.rst
..    scripts/regions2graph.rst
..    scripts/regions2predictions.rst
..    scripts/predict_genes.rst
..    scripts/cdhit2clusters.rst
..    scripts/prediction2pairs.rst
..    scripts/predictions2assembly.rst
..    scripts/predictions2cds.rst
..    scripts/predictions2disruptions.rst
..    scripts/predictions2genes.rst
..    scripts/predictions2introns.rst
..    scripts/predictions2pseudogenes.rst
..    scripts/predictions2transcripts.rst
..    scripts/compare_predictions.rst
..    scripts/compare_predictions2exons.rst
..    scripts/compare_projects.rst
..    scripts/count_orgs.rst
..    scripts/liftover_predictions.rst
..    scripts/mali_evaluate.rst
..    scripts/analyze_predictions.rst
..    scripts/analyze_queries.rst
..    scripts/analyze_duplications.rst
..    scripts/analyze_genes.rst
..    scripts/analyze_genetrees.rst
..    scripts/setup.rst
..    scripts/split_fasta.rst
..    scripts/split_genome.rst
..    scripts/split_genomic_fasta_file.rst
..    scripts/get_predictions.rst
..    scripts/get_genes.rst
..    scripts/grep_predictions.rst
..    scripts/exonerate2regions.rst
..    scripts/exonerate_combine_regions.rst
..    scripts/exons2clusters.rst
..    scripts/exons2exons.rst
..    scripts/exons2genes.rst
..    scripts/exons2map.rst
..    scripts/exons2stats.rst
..    scripts/gene2gene.rst
..    scripts/genes2quality.rst
..    scripts/pairs2gene_structure.rst
..    scripts/list2regions.rst
..    scripts/cds2codons.rst
..    scripts/assignments2pairs.rst
..    scripts/mali2predictions.rst
..    scripts/components2clusters.rst
..    scripts/extract_regions.rst
..    scripts/id2genes.rst
..    scripts/translate_forward2backward.rst

.. OPTIC pipeline
.. +++++++++++++++

.. Scripts that are specific for the `optic` orthology assignment pipeline.

.. .. toctree::
..    :maxdepth: 1

..    scripts/blast2alignments.rst
..    scripts/prune_multiple_alignment.rst
..    scripts/gff2predictions.rst
..    scripts/gff2transcripts.rst
..    scripts/gff_compare.rst
..    scripts/clade_export.rst
..    scripts/clone_run.rst
..    scripts/gtf2exons.rst
..    scripts/export_aaa.rst
..    scripts/export_all.rst
..    scripts/export_code.rst
..    scripts/export_predictions.rst
..    scripts/extract_clade_data.rst
..    scripts/filter_paralogous_links.rst
..    scripts/mysql_copy_tables.rst
..    scripts/plot_duplications.rst
..    scripts/graph_group_links_by_taxonomy.rst
..    scripts/graph_cluster_by_species.rst
..    scripts/mali2cleaned_mali.rst
..    scripts/create_gbrowser_files.rst
..    scripts/analyze_multiple_orthologs.rst
..    scripts/analyze_orthology.rst
..    scripts/analyze_orthology_multiple.rst
..    scripts/analyze_orthology_pairwise.rst
..    scripts/update_orthology.rst
..    scripts/analyze_cluster_expansion.rst
..    scripts/analyze_clustering.rst
..    scripts/analyze_clusters.rst
..    scripts/plot_multiple_synteny.rst
..    scripts/plot_synteny.rst
..    scripts/prune_fasta.rst
..    scripts/gbrowser_clone_devel.rst
..    scripts/gbrowser_delete_features.rst
..    scripts/analyze_sequences.rst
..    scripts/make2help.rst
..    scripts/cluster_genes_by_category.rst
..    scripts/correlate_fasta_identifier.rst
..    scripts/diff_transcript_sets.rst
..    scripts/mask_fasta.rst
..    scripts/orthologs2genes.rst
..    scripts/orthologs2list.rst
..    scripts/orthologs2transcripts.rst
..    scripts/patch_translations.rst
..    scripts/analyze_synteny.rst
..    scripts/evaluate_mali.rst
..    scripts/evaluate_trees.rst
..    scripts/split_links.rst
..    scripts/links2exons.rst
..    scripts/links2fasta.rst
..    scripts/sequences2graph.rst
..    scripts/transcripts2links.rst
..    scripts/update_blast.rst

.. Codon bias analysis
.. +++++++++++++++++++

.. Scripts relevant to codon bias analysis.

.. .. toctree::
..    :maxdepth: 1

..    scripts/analyze_codonbias.rst
..    scripts/analyze_codonbias_orthology.rst
..    scripts/analyze_codonbias_tables.rst
..    scripts/analyze_ribosomes.rst

.. Unsorted
.. --------

.. .. toctree::
..    :maxdepth: 1

.. Other
.. -----

.. .. toctree::
..    :maxdepth: 1

..    scripts/run_nubiscan.rst
..    scripts/adda2coverage.rst
..    scripts/radar.rst
..    scripts/tbl2veo.rst
..    scripts/check_blast_run.rst
..    scripts/graph2besthits.rst
..    scripts/graph_blast2besthits.rst
..    scripts/graph_blast2pairs.rst
..    scripts/graph_check.rst
..    scripts/graph_check_transitivity.rst
..    scripts/graph_map_links.rst
..    scripts/graph_reweight_links.rst

.. Obsolete
.. --------

.. This section contains obsolete or incomplete scripts.

.. .. toctree::
..    :maxdepth: 1

..    scripts/analyze_go.rst
..    scripts/tree_map_leaves.rst
..    scripts/liftover.rst
..    scripts/go2plot.rst
..    scripts/mali_extract.rst
..    scripts/mali_phylip2fasta.rst
..    scripts/mali_plain2aln.rst
..    scripts/mali_remove_gaps.rst

Unsorted
--------

.. toctree::
   :maxdepth: 1

   scripts/add_random_reads_to_bam.rst
   scripts/align_mali_vs_mali.rst
   scripts/analyze_go.rst
   scripts/bam2UniquePairs.rst
   scripts/bam2bam.rst
   scripts/bam2bed.rst
   scripts/bam2fastq.rst
   scripts/bam2peakshape.rst
   scripts/bam2stats.rst
   scripts/bam2transcriptContribution.rst
   scripts/barplotGo.rst
   scripts/beds2counts.rst
   scripts/bed2fasta.rst
   scripts/bed2stats.rst
   scripts/bed2table.rst
   scripts/beds2beds.rst
   scripts/blast2table.rst
   scripts/chain2stats.rst
   scripts/combine_tables.rst
   scripts/diff_chains.rst
   scripts/fasta2table.rst
   scripts/fasta2variants.rst
   scripts/fastq2N.rst
   scripts/fastq2fastq.rst
   scripts/fastq2solid.rst
   scripts/fastq2table.rst
   scripts/filter_reads.rst
   scripts/genelist_analysis.rst
   scripts/genome_bed.rst
   scripts/go2plot.rst
   scripts/gtf2overlap.rst
   scripts/cgat_import_extensions.rst
   scripts/index2bed.rst
   scripts/intervaltable2bed.rst
   scripts/liftover.rst
   scripts/list_overlap.rst
   scripts/mali2table.rst
   scripts/mali_extract.rst
   scripts/mali_phylip2fasta.rst
   scripts/mali_plain2aln.rst
   scripts/mali_remove_gaps.rst
   scripts/medip_merge_intervals.rst
   scripts/probeset2gene.rst
   scripts/cgat_rebuild_extensions.rst
   scripts/revigo.rst
   scripts/run_function.rst
   scripts/snp2snp.rst
   scripts/solexa2stats.rst
   scripts/tree_map_leaves.rst
   scripts/vcf2vcf.rst
   scripts/vcfstats2db.rst
   scripts/wig2wig.rst
   scripts/bam_vs_bam.rst
   scripts/bam_vs_bed.rst
   scripts/bam_vs_gtf.rst
   scripts/cgat2rdf.rst
   scripts/diff_bam.rst
   scripts/fasta2fasta.rst
   scripts/fasta2kmercontent.rst
   scripts/fastas2fasta.rst
   scripts/fastqs2fasta.rst
   scripts/fastqs2fastqs.rst
   scripts/gtf2tsv.rst
   scripts/gtfs2tsv.rst
   scripts/rnaseq_junction_bam2bam.rst
   scripts/split_gff.rst

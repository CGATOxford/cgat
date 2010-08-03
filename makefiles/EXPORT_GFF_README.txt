Gene predictions for various fly genomes.

Andreas Heger and Chris Ponting (2006)

Citation: This work has not been published yet. 

Method: 

Genes are predicted by homology using transcripts from D. melanogaster and
a pipeline based on exonerate [1].

References:

[1] Slater GS, Birney E. "Automated generation of heuristics for biological sequence comparison."
BMC Bioinformatics. 2005 Feb 15;6(1):31.

Output format:

Gene predictions for each fly genome are in separate files. The files
are in GFF3 format (see http://song.sourceforge.net/gff3.shtml).

Files in this release:

export_aaa_[species]_clean.gff3: 

Cleaned dataset of predictions. Only contains predictions which have conserved
gene structures with respect to the query. Genes can be conserved (CG), 
partially conserved (PG) or a single exon gene (SG). 

export_aaa_[species]_filtered.gff3: 

Non-redundant set of predictions. This dataset contains all predictions after
removing those that result from redundant predictions of the same gene (for 
example, due to prediction by both the ortholog and one or more paralogs) and
those that are dubious. For example, genes that are spanning complete genes
are removed.

export_aaa_[species]_full.gff3: 

All predictions. This is the full dataset of all predictions from the pipeline
before cleaning.

Release:


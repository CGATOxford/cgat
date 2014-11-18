=================================
Testing for Functional enrichment
=================================

This tutorial demonstrates the usage of *gat* with a simple example - 
where does a transcription factor bind in the genome? 

This tutorial uses the SRF data set described in `Valouev et
al. (2008)`_. The data sets used in this tutorial are available at:

http://www.cgat.org/~andreas/sample_data/srf.hg19.bed.gz

This :term:`bed` formatted file contains 556 high confidence peaks
from the analysis of `Valouev et al. (2008)`_ mapped to human
chromosome hg19.

We want to find out, where these binding sites are located in the
genome. First let us download the genomic sequence for hg19 and
index it::

   wget -qO- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
   | cgat index_fasta --file-format=tar.gz hg19 -
   > hg19.log
   	
Next, we need to define where intronic and intergenic regions are
located in the genome. We do this by by obtaining the latest geneset
from ENSEMBL_ and pushing it through a sequence of commands::

   wget -qO- ftp://ftp.ensembl.org/pub/release-72/gtf/homo_sapiens/Homo_sapiens.GRCh37.72.gtf.gz
   | gunzip
   | awk '$2 == "protein_coding"' 
   | cgat gff2gff --genome-file=hg19 --method=sanitize=genome --skip-missing
   | cgat gtf2gtf --method=sort --sort-order=gene
   | cgat gtf2gtf --method=merge-exons --with-utr
   | cgat gtf2gtf --method=filter --filter-method=longest-gene
   | cgat gtf2gtf --method=sort --sort-order=position
   | cgat gtf2gff --genome-file=hg19 --flank-size=5000 --method=genome
   | gzip
   > annotations.hg19.gff.gz

The commands do the following.

1. Reconcile the chromosome names in the gene set (ENSEMBL: 1,2,3)
   with the UCSC convention (chr1,chr2,chr3)::

      | cgat gff2ff --genome-file=hg19 --method=sanitize=genome --skip-missing

2. Sort the gene set by gene making sure that all exons within a gene
   appear in a block::

      | cgat gtf2gtf --method=merge-exons --with-utr

3. Merge overlapping exons from alternative transcripts of the same gene::

      | cgat gtf2gtf --method=merge-exons --with-utr

4. Resolve nested genes. In nested genes a genomic region might be
   both defined intronic and intergenic. Here, we select the longer
   one::

      | cgat gtf2gtf --method=filter --filter-method=longest-gene

5. Sort by genomic position::

      | cgat gtf2gtf --method=sort --sort-order=position

6. Define intronic, intergenic and other gene set based annotations::

      | cgat gtf2gff --genome-file=hg19 --flank-size=5000 --method=genome

The tool gff2stats

We can now use the file :file:`annotations.hg19.gff.gz` to classify
individual peaks with the :doc:`../../scripts/bed2table` tool::

   zcat srf.hg19.bed.gz
   | cgat bed2table --genome-file=hg19 --counter=classify-chipseq --gff-file=annotations.hg19.gff.gz
   | gzip 
   > srf.hg19.tsv.gz

The table :file:`srf.hg19.tsv.gz` contains a row for each interval in
the input file :file:`srf.hg19.bed.gz` describing which genomic
features it overlaps and assigns it to a category such as introning,
intergenic, etc.. We can upload this file into a database or view in 
excel to easily filter and summarize the data.

To get a more global view of where the transcription factor binds,
we make use of the gat_ tool. Gat tests if two sets of genomic
features are overlapping more - or less - than expected by chance
through simulation. Note that gat_ needs to be installed separately 
(``pip install gat``).

For gat, we need the file with genomic annotations
(:file:`annotations.hg19.gff.gz`) we created previously and a workspace - a
set of genomic regions that are accessible for simulation. Here, we
will use the full genome for simulation excluding regions that are
gaps as we do not expect to be able to detect transcription factor
binding sites in those in an NGS experiment.  To get these regions, we
use the :doc:`../../scripts/fasta2bed` tool::

   cat hg19.fasta
   | cgat fasta2bed --method=ungapped --min-gap-size=100
   | awk '$1 ~ /^chr/'
   | cut -f 1,2,3
   | gzip 
   > ungapped.hg19.bed.gz

Gat needs :term:`bed` formatted input files, so let us quickly convert
:file:`annogations.hg19.gff.gz`::

   zcat annogations.hg19.gff.gz
   | cgat gff2bed.py
   | gzip 
   > annotations.bed.gz

We are now ready to run gat::

   gat-run.py 
      --ignore-segment-tracks 
      --segments=srf.hg19.bed.gz
      --annotations=annotations.hg19.bed.gz 
      --workspace-bed-file=ungapped.hg19.bed.gz
      --num-samples=1000 
      --log=gat.log 
   | gzip
   > gat.out

The option `--ignore-segment-tracks` tells *gat* to ignore the fourth
column in the :term:`tracks` file and assume that all intervals in
this file belong to the same :term:`track`. If not given, each
interval would be treated separately. 

The above statement finishes in a few seconds. With large interval
collections or many annotations, *gat* might take a while. It is thus
good practice to always save the output in a file. The option `--log`
tells gat to save information or warning messages into a separate log
file.

The first 11 columns of the output file are the most informative:

+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|track |annotation|observed|expected  |CI95low   |CI95high  |stddev  |fold   |l2fold |pvalue    |qvalue    |
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|telomeric |0       |69.7440   |0.0000    |200.0000  |59.6216 |0.0141 |-6.1445|2.5100e-01|3.9443e-01|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|intergenic|6200    |13909.1770|12989.0000|14800.0000|570.3231|0.4458 |-1.1656|1.0000e-03|2.2000e-03|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|intronic  |8415    |11401.6660|10440.0000|12345.0000|577.7517|0.7381 |-0.4382|1.0000e-03|2.2000e-03|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|UTR3      |284     |305.5370  |114.0000  |500.0000  |120.2095|0.9297 |-0.1051|4.3000e-01|5.2556e-01|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|unknown   |0       |0.0140    |0.0000    |0.0000    |0.3603  |0.9862 |-0.0201|9.9800e-01|9.9800e-01|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|frameshift|0       |0.0050    |0.0000    |0.0000    |0.0947  |0.9950 |-0.0072|9.9700e-01|9.9800e-01|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|3flank    |800     |699.4930  |400.0000  |1045.0000 |187.2328|1.1435 |0.1934 |3.0300e-01|4.1662e-01|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|CDS       |758     |392.1510  |192.0000  |611.0000  |131.0955|1.9306 |0.9490 |3.0000e-03|5.5000e-03|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|flank     |1335    |176.1320  |50.0000   |350.0000  |90.7093 |7.5424 |2.9150 |1.0000e-03|2.2000e-03|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|5flank    |6224    |742.0590  |450.0000  |1071.0000 |191.1824|8.3775 |3.0665 |1.0000e-03|2.2000e-03|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+
|merged|UTR5      |3784    |104.0220  |0.0000    |237.0000  |68.5653 |36.0401|5.1715 |1.0000e-03|2.2000e-03|
+------+----------+--------+----------+----------+----------+--------+-------+-------+----------+----------+

The first two columns contain the name of the :term:`track` and
:term:`annotation` that are being compared. The columns
:term:`observed` and :term:`expected` give the observed and expected
nucleotide overlap, respectively, between the :term:`track` and :term:`annotation`.

The following columns CI95low, CI95high, stddev give 95% confidence
intervals and the standard deviation of the sample distribution,
respectively.

The :term:`fold` column is the fold enrichment or depletion and is 
computed as the ratio of :term:`observed` over :term:`expected`. The
column :term:`l2fold` is the log2 of this ratio.

The column :term:`pvalue` gives the empirical :term:`p-value`, i.e. in what
proportion of samples was a higher enrichment or lower depletion
found than the one that was observed.

The column :term:`qvalue` lists a multiple testing corrected :term:`p-value`.
Setting a qvalue threshold and accepting only those comparisons with a
qvalue below that threshold corresponds to controlling the false discovery
rate at that particular level.

.. _Valouev et al. (2008): http://www.ncbi.nlm.nih.gov/pubmed/19160518
.. _GREAT: http://bejerano.stanford.edu/great/public/html/
.. _MacLean et al. (2010): http://www.ncbi.nlm.nih.gov/pubmed/20436461
.. _Ensembl: http:://www.ensembl.org
.. _GO Gene Ontology: http://www.geneontology.org/
.. _gat: http://code.google.com/p/genomic-association-tester/

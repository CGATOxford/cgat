
Assessing CpG content in long non-coding RNA promoters
=======================================================

The description of pervasive transcription across many mammalian genomes has sparked an interest
in the role of long non-coding RNAs in diverse biological systems. Transcripts derived from non-coding 
loci have been shown to be important in a number of different processes including development and cancer. 
However, some features that are normally associated with protein coding genes are not observed in lncRNAs e.g
they are less conserved. Protein coding gene promoters have a characteristically high GC content and CpG
density. But do lncRNAs display the same bias in their promoters? In this example we show you how to use 
CGAT tools to answer this question. We will be using::

    gtf2gtf.py
    gtf2gff.py
    gff2bed.py
    bed2fasta.py
    fasta2table.py

Our initial input file is a :term:`gtf` formatted file containing genomic coordinates and annotations for
a set of lncRNAs - lncRNA.gtf.gz. We can compute the GC and CpG composition using a single command line
statement using multiple CGAT tools. However, as described in :ref:`quickstart`, we require an CGAT indexed
genome as input to both gtf2gff.py and bed2fasta.py. The first step is therefore to create the indexed genome.

In our example our lncRNA transcript models are from an RNA-seq experiment in human cells. We can index the
human hg19 reference genome by downloading the :term:`fasta` formatted genome from the UCSC website 
and running index_fasta.py::


    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz | index_fasta.py hg19 > hg19.log
    

We can then use this indexed genome as additional input when required. The code to generate a table with GC content and CpG
composition looks like::

    zcat lncRNA.gtf.gz 
    | gtf2gtf.py --sort=gene
    | gtf2gtf.py --merge-transcripts 
    | gtf2gff.py --genome-file=hg19 --method=promotors -p 1500 --sort
    | gff2bed.py
    | bed2fasta.py --genome-file=hg19
    | fasta2table.py --section=cpg 
    | gzip
    > lncRNA_cpg.tsv.gz


The above commands in turn (1) sorts the input file by gene identifier, (2) merges transcripts that have the same gene identifier,
(3) produces a set of lncRNA promoters 1.5Kb upstream of the lncRNA transcription start sites 
(using ``--method=promotors`` in combination with -p 1500), (4) converts gff formatted promoters into :term:`bed` format, 
(5) retrieves sequences from the human hg19 reference genome for lncRNA promoter intervals and (5) outputs statistics related 
to nucleotide composition including CpG content (specified with the ``--section=cpg`` option). 
The output file lncRNA_cpg.tsv.gz will be a tab-delimited text file which will look like:

+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|id                                                      |nC |nG |nA |nT |nN|nUnk|nGC |nAT |nCpG|pC      |pG      |pA      |pT      |pN      |pUnk    |pGC     |pAT     |pCpG    |CpG_ObsExp|
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|ENSG00000224969.1 chr1:948573..950073 (+)               |423|518|277|282|0 |0   |941 |559 |71  |0.282000|0.345333|0.184667|0.188000|0.000000|0.000000|0.627333|0.372667|0.094667|0.486048  |
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|NONCO170 chr1:33464145..33465645 (+)                    |418|396|359|327|0 |0   |814 |686 |37  |0.278667|0.264000|0.239333|0.218000|0.000000|0.000000|0.542667|0.457333|0.049333|0.335291  |
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|NONCO195 chr1:87239820..87241320 (+)                    |354|294|425|427|0 |0   |648 |852 |13  |0.236000|0.196000|0.283333|0.284667|0.000000|0.000000|0.432000|0.568000|0.017333|0.187363  |
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|NONCO55 chr1:108591390..108592890 (+)                   |296|323|471|410|0 |0   |619 |881 |9   |0.197333|0.215333|0.314000|0.273333|0.000000|0.000000|0.412667|0.587333|0.012000|0.141202  |
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|NONCO59 chr1:111181220..111182720 (+)                   |270|380|452|398|0 |0   |650 |850 |9   |0.180000|0.253333|0.301333|0.265333|0.000000|0.000000|0.433333|0.566667|0.012000|0.131579  |
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|NONCO215 chr1:120190857..120192357 (+)                  |350|415|384|351|0 |0   |765 |735 |62  |0.233333|0.276667|0.256000|0.234000|0.000000|0.000000|0.510000|0.490000|0.082667|0.640275  |
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|NONCO66 chr1:121117751..121119251 (+)                   |374|313|340|473|0 |0   |687 |813 |16  |0.249333|0.208667|0.226667|0.315333|0.000000|0.000000|0.458000|0.542000|0.021333|0.205020  |
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|NONCO69 chr1:144569176..144570676 (+)                   |233|299|498|470|0 |0   |532 |968 |21  |0.155333|0.199333|0.332000|0.313333|0.000000|0.000000|0.354667|0.645333|0.028000|0.452151  |
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+
|NONCO70 chr1:144592229..144593729 (+)                   |382|382|361|375|0 |0   |764 |736 |53  |0.254667|0.254667|0.240667|0.250000|0.000000|0.000000|0.509333|0.490667|0.070667|0.544804  |
+--------------------------------------------------------+---+---+---+---+--+----+----+----+----+--------+--------+--------+--------+--------+--------+--------+--------+--------+----------+


The ``--section`` option specifies that we want to include statistics on CpG composition in the output. Alternative options
include::

    length
    na
    aa
    cpg 
    degeneracy
    bias
    codons
    codon-usage
    codon-translator
    sequence  


As the output is in tab separated format it is straight-forward to load into statistical/plotting software such as R and perform further 
downstream analysis. 

==============================================================
RNAseq of Prostate cancer patients pre and post AAT treatment
==============================================================

Introduction
-------------

Prostate Cancer (PCa) is driven by androgen hormones acting thorugh the androgen receptor. As such the mainstay of advanced or metastatic prostate cancer treatments has been androgen ablation thearpy, where the production of androgen hormones such as testosterone is inhibited. This treatment is generally effective initally, but after a period the cancer becomes Castration Resistant (CRPCa). The aim of this project is to compare RNA expression profiles generated from samples taken from patients before and after treatment using RNAseq. This data will be used to investigate hypotheses related to the molecular mechanisms behind CRPCa.
 
This report represents the outcome of a preliminary analysis of the the data provided by the project. Several automated, routine analyses have been run as well as a small number of analyses tailored to the project. This report represents a living document that will be updated as the project continues.

Data Description and quality control
-------------------------------------

The data consists of 10 matched samples from 5 prostate cancer patients taken before and after a course of androgen ablation treatment (AAT). Hereafter each of these samples are referred to as Pre and Post samples for pre- and post-treatment respectively. In each case the samples consist of around 50-55M paired end reads from the Illumina GAII platform.

In all samples the quality score of the reads stayed above 30 for the entire 90bp of the reads. At least 60% of reads mapped to the genome in each sample and for 8 out of the 10 samples more than 70% of reads mapped to the genome (see `mapping status`_).

Detection of fusion transcripts
--------------------------------

Fusion transcripts are common in prostate cancer. In order to detect fusion transcripts in the data, we ran the tophat-fusion algorithm which finds fusion transcripts by finding reads that span junctions and also read pairs that map on both sides. The following fusions were found:


+----------------------------------+------------------------+---------------+
| Left Fragment                    | Right Fragment         | Sample        |
+----------------+------+----------+------+------+----------+---------------+
|Gene            | Chr  | Postion  | Gene | Chr  | Postion  | Sample        |
+----------------+------+----------+------+------+----------+---------------+
| TMPRSS2        | chr21| 42870044 | ERG  | chr21| 39817542 | Patient28_Pre |
|                |      |          |      |      |          |               |
+----------------+------+----------+------+------+----------+---------------+
|                |      |          |      |      |          | Patient29_Pre |
|                |      |          |      |      |          |               |
+----------------+------+----------+------+------+----------+---------------+
| TMPRSS2        | chr21| 42880006 | ERG  | chr21| 39817542 | Patient18_Pre |
|                |      |          |      |      |          |               |
+----------------+------+----------+------+------+----------+---------------+
|                |      |          |      |      |          | Patient29_Pre |
+----------------+------+----------+------+------+----------+---------------+
| IMPDH1         | chr7 | 128037161| KCND2| chr7 | 120336337| Patient12_Pre |
|                |      |          |      |      |          |               |
+----------------+------+----------+------+------+----------+---------------+
|                |      |          |      |      |          | Patient12_Post|
+----------------+------+----------+------+------+----------+---------------+
| HELZ           | chr17| 65110443 | NMT1 | chr17| 43159007 | Patient18_Pre |
+----------------+------+----------+------+------+----------+---------------+
| ENSG00000211651| chr22| 22735711 | IGLL5| chr22| 23235959 | Patient18_Pre |
+----------------+------+----------+------+------+----------+---------------+
| ENSG00000226958| chrX | 108297773| HFM1 | chr1 | 91853138 | Patient12_Pre |
+----------------+------+----------+------+------+----------+---------------+
| RPS3           | chr11| 75115899 | CAGE1| chr6 | 7339050  | Patient11_Pre |
+----------------+------+----------+------+------+----------+---------------+
|                |      |          |      |      |          | Patient11_Post|
+----------------+------+----------+------+------+----------+---------------+
|                |      |          |      |      |          | Patient28_Post|
+----------------+------+----------+------+------+----------+---------------+
| PTPRM          | chr18| 7774268  | RBBP8| chr18| 20531943 | Patient18_Pre |
+----------------+------+----------+------+------+----------+---------------+
|                |      |          |      |      |          | Patient18_Post|
+----------------+------+----------+------+------+----------+---------------+
Table 1: Potential fusions found between transcribed RNAs

It is interesting to note that three of the patients carry a fusion between TMPRSS2 and ERG, a known fusion in prostate cancer (patient 29 may carry two seperate fusions). However, none of the matched post-treatment samples were identified as carrying this fusion, potentially indicating that cells carrying this fusion were eliminated by the treatement. To the best of our knowledge none of the other fusions listed here have been identified previously (they are not listed in the Mitelman Database of Chromosome Aberrations and Gene Fusions in Cancer).

In all cases, apart from that of the RPS3-CAGE1 fusion in Patient28_Post, fusions are either present both before and after treatment or only before, but not only after treatment. 

These results should be interpreted with care. The accuracy of TopHat-fusion is currently unassessed and one can imagine many ways in which false positives might occur. For example, the fusion between RPS3 and CAGE1 is found in three samples, patient 11 both before and after treatment and patient 28 after treatment. However, the CAGE1 primary transcript contains a pseudo-gene from RPS3 within one of its introns. While we will continue to investigate criteria by which we may be more or less confident of these predictions, is important that any novel fusions identified here should be confirmed experimentally.

Genes of Interest
------------------

A number of genes have previously been implicated in the transition to CRPCa. Here we examine the expression of several. While these are not proper statistical measurements of differential expression between two samples, they can help to give an idea of whether a gene is expressed at all, and if particular samples express the gene at a higher or lower level than the others. Values are given in fragments per kilobase per million, which is a measure of expression normalised by the length of the transcript in question. 

One of the primary targets of interest in prostate cancer and its transition to a castration resistant form is the Androgen Receptor (AR) and its isoforms. The cuffdiff package was used to assess the levels of expression for the known isoforms of AR. The levels of each of the isoforms that had non-zero expression are shown below.

.. image:: AR_isoforms.png

As well as expression from the major isoform, ENST00000374690, there is also a much reduced level of expression from three other isoforms in many of the samples. The ENST00000504326 isoform is missing the final 5 exons compared to the major isoform. ENST00000514029 is similar but is degraded by nonsense mediated decay. The final isoform ENST00000544984 has a different start codon to the major isoform, but is otherwise similar in its structure. 

As well as AR, we briefly highlight a few other examples of expression levels.

.. image:: GOI_levels.png

Note that here 'high' and 'low' expression simply refer the scale of the axis rather than rigorous classification of expression levels.

The expression of the interlukins IL6 and IL8 has been previously associated with CRPCa. Here we can see that the interlukins are expressed as a very low level or not at all in many of these samples with some expression of IL8 in a couple of samples. NMB is a mammalian bombesin homolog. Expression of neuro-peptides such as bombesin have been implicated previously in CRPCa.

Differential Expression
-------------------------

We used the software edgeR, which allows the identification of differentially expressed genes using experimental designs based on linear models. This allowed us to do a test that is equivalent to a paired t-test on the expression data. This provides far greater power than other analyses. 

The following plot shows the relationship between expression, fold change and p-values. It is equivalent to an MA or 'smear' plot in microarray analysis. Genes which are significant with an FDR < 0.05 are shown in red.

.. image:: EdgeR_MA.png

This analysis identified 1967 genes that were significantly differentially expressed between Pre and Post samples at an FDR <0.05. Of these, 708 were up-regulated at least two fold and 734 down-regulated at least two fold. 

The following tables show the top twenty up- and down-regulated genes.

.. report:: DESeq.UpReg20
   :render: table
   
   Top 20 up-regulated genes by FDR

.. report:: DESeq.DownReg20
   :render: table

   Top 20 down-regulated genes by FDR

The complete list of differentially expressed genes can be found in the export directory.

 

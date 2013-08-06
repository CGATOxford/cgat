.. _genome_coverage:




================
Genome coverage
================

In a perfect metagenome assembly we would obtain a single contig from each species in the 
community. However, due to the variability in abundance of species and the sampling procedure
associated with the sequencing procedure this is not the case. Genomes from lowly abundant
species will never be fully assembled. We can estimate the expected coverage for each genome
in our simulated data using the formula:

pcoverage = L*N / S

Where L = read length
      N = read number
      S = known genome size

Using contigs that have been aligned back to the genome sequences we can observe how our observed
coverage deviates from what we expect.



.. report:: GenomeCoverage.GenomeCoverage
   :render: scatter-plot

expected and observed genome coverage for each species
 

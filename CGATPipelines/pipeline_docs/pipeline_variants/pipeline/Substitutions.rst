=====================
Substitution matrices
=====================

Nucleotides
===========

All SNPs
--------

The following plot shows the substitution matrices for all SNPs:

.. report:: Trackers.SubstitutionMatrixNucleotides
   :render: matrix-plot
   :transform-matrix: normalized-total,sort
   :groupby: track
   :colorbar-format: %6.4f
   :slices: all
   :layout: column-3
   :width: 200

   Substitution matrix. The reference base is in the rows	
   and consensus bases are in the column.

`Ambiguity codes <http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html#500>`_ 
are used for heterozygous positions.

Coding SNPs
-----------

The following plot shows the substitution matrices for coding SNPs only:

.. report:: Trackers.SubstitutionMatrixNucleotides
   :render: matrix-plot
   :transform-matrix: normalized-total,sort
   :groupby: track
   :colorbar-format: %6.4f
   :slices: coding
   :layout: column-3
   :width: 200

   Substitution matrix. The reference base is in the rows	
   and consensus bases are in the column.

Noncoding SNPs
--------------

The following plot shows the substitution matrices for non-coding SNPs,
i.e., those that are in intergenic or intronic regions:

.. report:: Trackers.SubstitutionMatrixNucleotides
   :render: matrix-plot
   :transform-matrix: normalized-total,sort
   :groupby: track
   :colorbar-format: %6.4f
   :slices: noncoding
   :layout: column-3
   :width: 200

   Substitution matrix. The reference base is in the rows	
   and consensus bases are in the column.

Transition/transversion rations
===============================

All SNPs
--------

.. report:: Trackers.TransitionTransversionRatios 
   :render: table
   :slices: all
                                              
   Transition/transversion ratios for all SNPs

Coding SNPs
-----------

.. report:: Trackers.TransitionTransversionRatios 
   :render: table
   :slices: coding
                                              
   Transition/transversion ratios for coding SNPs

Noncoding SNPs
--------------

.. report:: Trackers.TransitionTransversionRatios 
   :render: table
   :slices: noncoding
                                              
   Transition/transversion ratios for noncoding SNPs

Amino acids
===========

For coding SNPs, the following plots show the substitions between amino
acids. In the following plots, synonymous substitutions are along
the diagonal and non-synonymous substitutions are off-diagonal.

.. report:: Trackers.SubstitutionMatrixAminoAcids
   :render: matrix-plot
   :transform-matrix: normalized-total,sort
   :max-rows: 30
   :groupby: track
   :max-cols: 30
   :colorbar-format: %6.4f
   :layout: column-3
   :width: 300

   Matrix of amino acid substitution frequencies. 
   The reference amino acid is in the rows and the
   and variant amino acid is in the column. Note that 
   only homozygous changes are considered here. Note
   that ``X`` refers to the stop codon.

.. .. report:: Trackers.SubstitutionMatrixAminoAcids
..    :render: matrix
..    :transform-matrix: normalized-total,sort
..    :groupby: track

..    Matrix of amino acid substitution frequencies. 
..    The reference amino acid is in the rows and the
..    and variant amino acid is in the column. Note that 
..    only homozygous changes are considered here. Note
..    that ``X`` refers to the stop codon.


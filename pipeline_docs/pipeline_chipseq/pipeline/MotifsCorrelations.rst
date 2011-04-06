Correlations
------------

The following plot lists correlations between various features
of intervals and motifs. The features are length of a motif,
combined evalue of the motifs, number of motifs found, peak value
of the interval and average value of the interval.

.. needs work
.. .. report:: Motifs.MastAllCorrelations
..    :render: scatter-plot
..    :transform: combine
..    :layout: column-3
..    :width: 300
..    :tf-fields: coefficient

..    Matrix of correlation coefficients between various interval
..    features.

.. .. report:: Motifs.MastAllCorrelations
..    :render: matrix-plot
..    :transform: select,correlation,select
..    :tf-fields: peakval,pvalue
..    :transform-matrix: correspondence-analysis

..    Pearson correlations between various interval features.

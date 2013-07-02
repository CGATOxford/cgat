========================
Transcriptome validation
========================

The following pages summarize the validation via mapping to the transrciptome.
Briefly, reads are mapped to the transcriptome using stringent criteria. 
Read matcheing to both the genome and a known transcript are examined if the
genomic locations between the two matches are consistent.

.. report:: BenchmarkReport.TranscriptomeValidationSummary
   :render: matrix
   :transform-matrix: normalized-row-total
   :slices: removed_mismapped,location_ok,unmapped_transcript

   Proportions of mismapped reads in the various runs

.. report:: BenchmarkReport.TranscriptomeValidationSummary
   :render: pie-plot
   :slices: removed_mismapped,location_ok,unmapped_transcript
   :layout: column-4
   :width: 200

   Proportions of mismapped reads in the various runs

.. report:: BenchmarkReport.TranscriptomeValidationSummary
   :render: pie-plot
   :slices: removed_mismapped,location_ok
   :layout: column-4
   :width: 200

   Proportions of mismapped reads in the various runs (without
   unmapped transcripts).

.. report:: BenchmarkReport.TranscriptomeValidationSummary
   :render: matrix
   :slices: removed_mismapped,location_ok
   :transform-matrix: normalized-row-total

   Proportions of mismapped reads in the various runs (without
   unmapped transcripts).

.. report:: BenchmarkReport.TranscriptomeValidationSummary
   :render: matrix
   :slices: removed_mismapped,location_ok,unmapped_transcript

   Number of mismapped reads in the various runs

.. report:: BenchmarkReport.TranscriptomeValidationSummary
   :render: table
   
   Full results from the transcriptome validation.

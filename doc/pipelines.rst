.. _pipelines:

=========================
CGAT Pipelines
=========================

Best practice for CGAT pipelines:

1. All non-trivial code should be extracted to modules or scripts.
2. Modules should not access PARAMS dictionary directly, but parameters should be passed to the function.
3. Important processing steps where different external tools could potentially be employed the design of the module classes should be carefully considered to ensure consistent input and output file formats for different tools. PipelineMapping provides a good example for this. 
4. All production pipelines should include tests for consistency which can be run automatically.
5. Where appropriate pipelines should include a small test dataset with published results for comparison. This dataset can be run on each pipeline run and included in the pipeline report where it can be used as a pipeline control.
6. Periodic code review meetings where interested parties can agree of major changes to production pipelines and associated modules – to be arranged as required.
7. The best way to manage pipeline improvements is by individuals using pipelines talking responsibility for incremental improvement. As best practice fellows should announce plans to modify particular pipelines and modules on the CGAT members list to avoid duplication of effort. Fellows should log the changes that they make in a change log and document both modules and pipelines in detail. 




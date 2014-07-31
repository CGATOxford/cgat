======
Status
======

MD5 Comparison
==============

The MD5 Comparison compares the files output by
the latest pipeline run against the results of
a reference pipeline run. A run was successfull
if exactly the same files are present and all
files are identical.

.. report:: Status.ComparisonStatus
   :render: status

   Pipeline Status

Completion status
====================

The following report lists the completion status of the pipelines, 
pipelines that have run and completed successfully.

.. report:: TestingReport.PipelineStatus 
   :render: status                  
   
   Status report of pipeline running. The info field shows the last
   time the pipeline was started.

Report status
=============

The following report lists the completion status of the reports that
are associated with each pipelines.

.. report:: TestingReport.ReportTable
   :render: table

   Table with detailed report results.

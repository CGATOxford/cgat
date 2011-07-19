===============
Pipeline status
===============

This page combines various status information.

Mapping
=======

.. report:: Status.MappingStatus
   :render: status        

   Status of mapping stage 

Transcript building
===================

Transcript building requires sufficient read-depth and
a good proportion of spliced reads. If one of them are missing,
ab-initio transrcripts might be incomplete.

.. report:: Status.TranscriptStatus
   :render: status        

   Status of mapping stage 

Differential expression
=======================

The success of estimating differential expression depends on the
read depth. 

.. report:: Status.ExpressionStatus
   :render: status        

   Status of expression stage 


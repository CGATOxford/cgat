===============
Pipeline status
===============

This page combines various status information.

Mapping
=======

.. report:: Status.MappingStatus
   :render: status        

   Status of mapping stage 

For more information, see :ref:`Mapping`.

Transcript building
===================

Transcript building requires sufficient read-depth and
a good proportion of spliced reads. If one of them are missing,
ab-initio transrcripts might be incomplete.

.. report:: Status.TranscriptStatus
   :render: status        

   Status of mapping stage 

For more information, see 
    * :ref:`Reference coverage`,
    * :ref:`FPKM confidence intervals`
    * :ref:`Classification of transcript models`

Differential expression
=======================

The success of estimating differential expression depends on the
read depth. The section below details how accurately expression
of genes (not transcripts) could be measured and how often tests
for differential expression succeeded.

.. report:: Status.ExpressionStatus
   :render: status        

   Status of expression stage for known coding transcripts

.. report:: Status.ExpressionStatusNoncoding
   :render: status        

   Status of expression stage for known non-coding transcripts

.. report:: Status.DifferentialExpressionStatus
   :render: status        

   Status of differential expression stages

For more information, see :ref:`Cuffdiff results`.

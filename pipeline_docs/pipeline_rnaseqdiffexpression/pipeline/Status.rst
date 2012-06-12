===============
Pipeline status
===============

This page combines various status information.

Differential expression
=======================

The success of estimating differential expression depends on the
read depth. The section below details how accurately expression
of genes (not transcripts) could be measured and how often tests
for differential expression succeeded.

.. report:: Status.ExpressionStatus
   :render: status        

   Error bars

.. report:: Status.DifferentialExpressionStatus
   :render: status        

   Status of differential expression stages

For more information, see :ref:`Cuffdiff results`.

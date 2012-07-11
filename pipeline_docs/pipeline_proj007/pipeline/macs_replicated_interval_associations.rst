==================================================================
Association of CAPseq intervals with Genes & Enhancers
==================================================================

CAPseq interval were annotated as within 5kb of a protein-coding or non-coding TSS or overlapping with a H3K4Me1 interval. 
Association was performed hierarchically as follows:

Extended Hierarchy

coding TSS > lncRNA > short RNA TSS > pseudogene > processed transcript > enhancer > RNAseq transcript > intergenic

intermediate Hierarchy

coding TSS > lncRNA > short RNA TSS > enhancer > RNAseq transcript > intergenic

Reduced Hierarchy

coding TSS > noncoding TSS > intergenic

.. report:: macs_replicated_interval_associations.replicatedAssociationsHierarchy
   :render: table

   CAPseq intervals assoctaed with different genomic features (extended hierarchy)

.. report:: macs_replicated_interval_associations.replicatedAssociationsHierarchy3
   :render: table

   CAPseq intervals assoctaed with different genomic features (intermediate hierarchy)
      
.. report:: macs_replicated_interval_associations.replicatedAssociationsHierarchy2
   :render: table

   CAPseq intervals assoctaed with different genomic features (reduced hierarchy)
   
.. report:: macs_replicated_interval_associations.replicatedAssociations
   :render: table

   CAPseq intervals assoctaed with different genomic features



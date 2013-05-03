



=========================
Classification of lncRNA
=========================

LncRNA are classified based on their relative position to protein coding genes. They are
classified into


Classes
=======

* antisense - transcript overlapping protein coding exons on opposite strand                                                                                                                                                               
* antisense_upstream - transcript < 2kb from tss on opposite strand                                                                                                                                                                        
* antisense_downstream - transcript < 2kb from gene end on opposite strand                                                                                                                                                                 
* sense_upstream - transcript < 2kb from tss on same strand                                                                                                                                                                                
* sense_downstream - transcript < 2kb from gene end on same strand                                                                                                                                                                         
* intergenic - >2kb from anyt protein coding gene                                                                                                                                                                                          
* sense_intronic - overlaps protein coding gene intron on same strand (NB. these are not neccessarily 'contained' within an intron)                                                                                                                                                                            
* antisense_intronic - overlaps protein coding intron on opposite strand (NB. these are not neccessarily 'contained' within an intron)




.. report:: Classification.TranscriptClassificationCount
   :render: table
   

   number of lncRNA transcripts for a given calssification


.. report:: Classification.TranscriptClassificationProportion
   :render: pie-plot
   :layout: column-2


   proportion of lncRNA transcripts with a particular classification


.. report:: Classification.GeneClassificationProportion
   :render: pie-plot
   :layout: column-2


   proportion of lncRNA genes with a particular classification






   

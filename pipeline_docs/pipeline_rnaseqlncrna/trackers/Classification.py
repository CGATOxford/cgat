from SphinxReport.Tracker import *

import IOTools
import GTF
import numpy as np
import scipy.stats
import collections

#################################################
#################################################
#################################################
class TranscriptClassificationProportion(TrackerSQL):

    pattern = "(.*)_class"
    slices = ["antisense"
              , "antisense_upstream"
              , "antisense_downstream"
              , "sense_upstream"
              , "sense_downstream"
              , "intergenic" 
              , "sense_intronic" 
              , "antisense_intronic"]


    def __call__(self, track, slice = None):
        
        total = float(self.getValue("SELECT COUNT(*) FROM %(track)s_class"))

        statement = "SELECT COUNT(*) FROM %(track)s_class WHERE class = '%(slice)s'"
        val =  float(self.getValue(statement))
        return (val/total)*100

#-------------------------------------------------------------------
class TranscriptClassificationCount(TranscriptClassificationProportion):

    def __call__(self, track, slice = None):
        statement = "SELECT COUNT(*) FROM %(track)s_class WHERE class = '%(slice)s'"
        return self.getValue(statement)

#-------------------------------------------------------------------
class GeneClassificationProportion(TranscriptClassificationProportion): 
    
    # remove slices
    slices = None

    def __call__(self, track, slice = None):

        total = float(self.getValue("SELECT COUNT(DISTINCT gene_id) FROM %s" % track + "_class"))
        
        # done in a hierarchy so that the 'most important' classification is maintained
        hierarchy = [ "antisense"
                     , "antisense_upstream"
                     , "antisense_downstream"
                     , "sense_upstream"
                     , "sense_downstream"
                     , "intergenic" 
                     , "sense_intronic" 
                     , "antisense_intronic" ]

        result = collections.defaultdict(set)
        for data in self.execute("SELECT gene_id, class FROM %s" % track + "_class").fetchall():
            gene = data[0]
            _class = data[1]
            result[gene].add(_class)
        
        # remove classifications where the gene has multiple classifications
        # e.g antisense trumps antisense_intronic
        genes_multi = {}
        for gene, _class in result.iteritems():
            if len(_class) > 1:
               for h in hierarchy:
                   if h in _class:
                       if gene not in genes_multi:
                           genes_multi[gene] = h
                       else:
                           print "removed %s as %s" % (gene, h)

        for gene, _class in genes_multi.iteritems():
            if gene in result:
                # keep as set for consisteny downstream
                result[gene] = set([_class])
                

        # counts
        counts = collections.defaultdict(int)
        for _class in result.values():
            # all sets so convert to string via list
            c = list(_class)[0]
            if  c not in counts:
                counts[c] = 0
            else:
                counts[c] += 1

        # return dictionary
        for _class, count in counts.iteritems():
            counts[_class] = (count/total)*100
        return dict (counts)







            

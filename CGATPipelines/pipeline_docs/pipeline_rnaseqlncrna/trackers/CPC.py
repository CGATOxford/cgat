from SphinxReport.Tracker import *

import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import numpy as np
import scipy.stats
import collections
import sqlite3
import CGATPipelines.PipelineLncRNA as PipelineLncRNA
from LncRNACounts import *

#################################################
#################################################
#################################################

class CPC(TrackerSQL):

    pattern = "(.*)_result"
    slices = ["total","coding","noncoding"]

    def __call__(self, track, slice = None):
        
        if slice == "total":
            statement = "SELECT COUNT(*) FROM %(track)s_result"
            return self.getValue(statement)
        elif slice == "coding":
            statement = "SELECT COUNT(*) FROM %(track)s_result WHERE C_NC='%(slice)s'"
            return self.getValue(statement)
        elif slice == "noncoding":
            statement = "SELECT COUNT(*) FROM %(track)s_result WHERE C_NC='%(slice)s'"
            return self.getValue(statement)


class CPCScoreCorellation(TrackerSQL):

    tracks = "length vs. CPC score"

    def __call__(self, track):
        
        length = {}
        for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile("gtfs/lncrna_filtered.gtf.gz"))):
            length[transcript[0].transcript_id] = sum([gtf.end - gtf.start for gtf in transcript])
        
        score = {}
        dbh = sqlite3.connect("csvdb")
        cc = dbh.cursor()
        for data in cc.execute("SELECT transcript_id, CP_score FROM lncrna_filtered_cpc_result"):
            score[data[0]] = data[1]
    
        result = {"length": [], "score": []}
        for transcript, value in length.iteritems():
            result["length"].append(np.log10(length[transcript]))
            result["score"].append(score[transcript])
        return result


class CPCLength(AveLength):
    
    def getTracks(self):
        return ["cpc_removed", "lncrna_final"]

class CPCLengthDistribution(LengthDistribution):
    
    def getTracks(self):
        return ["cpc_removed", "lncrna_final"]

class CPCLengthCumFreq(LengthCumFreq):

    def getTracks(self):
        return ["cpc_removed", "lncrna_final"]

class CPCExonNumber(CPCLength):
    
    def getSlices(self):
        
        return None
    
    def __call__(self, track, slice = None):
        
        c_transcript = []
        c_gene = []
        for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(self.getFilename(track)))):
            c_transcript.append(len(transcript))
        for gene in GTF.flat_gene_iterator(GTF.iterator(IOTools.openFile(self.getFilename(track)))):
            c_gene.append(len(gene))

        return odict( ( ("transcript", np.mean(c_transcript)), ("gene",np.mean(c_gene) )) )


class CPCClass(TrackerSQL):
    
    pattern = "(.*)_cpc_result"

    def __call__(self, track, slice = None):

        classes = ["antisense"
              , "antisense_upstream"
              , "antisense_downstream"
              , "sense_upstream"
              , "sense_downstream"
              , "intergenic" 
              , "sense_intronic" 
              , "antisense_intronic"]

        coding_set = {}
        for gtf in GTF.iterator(IOTools.openFile("gtfs/lncrna_filtered.class.gtf.gz")):
            coding_set[gtf.transcript_id] = gtf.source

        result = {"noncoding": {}, "coding":collections.defaultdict(int)}
        total_nc = float(self.getValue("SELECT COUNT(*) FROM %(track)s_cpc_result WHERE C_NC = 'noncoding'"))
        for c in classes:
            result["noncoding"][c] = (float(self.getValue("""SELECT COUNT(*) FROM lncrna_final_class as a, %s_cpc_result as b WHERE a.class = '%s' 
                                                              AND b.C_NC = 'noncoding' 
                                                              AND a.transcript_id = b.transcript_id""" % (track,c)))/total_nc)*100

        
        total_c = len(coding_set.keys())
        for c in classes:
            ids = self.getValues("SELECT transcript_id FROM %(track)s_cpc_result WHERE C_NC = 'coding'")
            for i in ids:
                if i in coding_set.keys():
                    if coding_set[i] == c:
                        result["coding"][c] += 1
            
        for x, y in result["coding"].iteritems():
            result["coding"][x] = (float(y)/total_c)*100
            
        return result










        

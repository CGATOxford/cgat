from SphinxReport.Tracker import *

import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import numpy as np
import scipy.stats
import collections
import os

#################################################
#################################################
#################################################

class Counts(TrackerSQL):

    pattern = "(.*)_stats"

    def __call__(self, track, slice = None):
        
        return self.getRow("SELECT * FROM %(track)s_stats")
        
#-------------------------------------------------------------------
class AveLength(Tracker):
    
    def getFilename(self, track):
        filename = os.path.join("gtfs",track) + ".gtf.gz"
        return filename

    def getTracks(self, subset = None):
        tracks = [ os.path.basename(x)[:-len(".gtf.gz")] for x in glob.glob("gtfs/*_coding.gtf.gz")] + ["lncrna_final"] + ["refcoding"]
        return tracks

    def getSlices(self):
        return ["transcript", "gene"]

    def __call__(self, track, slice = None):
        
        if slice == "transcript":
            lengths_transcripts = []
            for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(self.getFilename(track)))):
                length = sum([gtf.end - gtf.start for gtf in transcript])
                lengths_transcripts.append(length)
            return np.mean(lengths_transcripts)
        
        elif slice == "gene":
            lengths_genes = []
            for gene in GTF.flat_gene_iterator(GTF.iterator(IOTools.openFile(self.getFilename(track)))):
                length = sum([gtf.end - gtf.start for gtf in gene])
                lengths_genes.append(length)
            return np.mean(lengths_genes)

#-------------------------------------------------------------------
class LengthDistribution(AveLength):

    def __call__(self, track, slice = None):
        
        if slice == "transcript":
            lengths_transcripts = []
            for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(self.getFilename(track)))):
                length = sum([gtf.end - gtf.start for gtf in transcript])
                lengths_transcripts.append(length)
            return np.log10(lengths_transcripts)
        
        elif slice == "gene":
            lengths_genes = []
            for gene in GTF.flat_gene_iterator(GTF.iterator(IOTools.openFile(self.getFilename(track)))):
                length = sum([gtf.end - gtf.start for gtf in gene])
                lengths_genes.append(length)
            
            return np.log10(lengths_genes)

#-------------------------------------------------------------------        
class LengthCumFreq(LengthDistribution):

    def __call__(self, track, slice = None):
        
        if slice == "transcript":
            lengths_transcripts = []
            for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(self.getFilename(track)))):
                length = sum([gtf.end - gtf.start for gtf in transcript])
                lengths_transcripts.append(length)
            counts, lower, dx, _ = scipy.stats.cumfreq(lengths_transcripts, numbins=40, defaultreallimits=(0,20000))
            x = np.arange(counts.size) * dx + lower
            return odict( (("length", x), ("cumulative frequency", counts/len(lengths_transcripts))) )

        
        elif slice == "gene":
            lengths_genes = []
            for gene in GTF.flat_gene_iterator(GTF.iterator(IOTools.openFile(self.getFilename(track)))):
                length = sum([gtf.end - gtf.start for gtf in gene])
                lengths_genes.append(length)
            counts, lower, dx, _ = scipy.stats.cumfreq(lengths_genes, numbins=40, defaultreallimits=(0,20000))
            x = np.arange(counts.size) * dx + lower
            return odict( (("length", x), ("cumulative frequency", counts/len(lengths_genes))) )

#-------------------------------------------------------------------        
class TranscriptNumbers(Tracker):
    
    def getFilename(self, track):
        filename = os.path.join("gtfs",track) + ".gtf.gz"
        return filename

    def getTracks(self, subset = None):
        
        # includes the refcoding gene set as a comparison
        tracks = [os.path.basename( glob.glob("gtfs/*_coding.gtf.gz")[0]) [:-len(".gtf.gz")], "lncrna_final", "refcoding"]
        return tracks

    def __call__(self,track, slice = None):
        
        transcript_counts = collections.defaultdict( set )
        counts = []
        for gtf in GTF.iterator(IOTools.openFile(self.getFilename(track))):
            transcript_counts[gtf.gene_id].add(gtf.transcript_id)
        for gene, transcripts in transcript_counts.iteritems():
            counts.append(len(transcripts))
        return np.mean(counts)

#-------------------------------------------------------------------        
            
class TranscriptNumberDistribution(TranscriptNumbers):
    
    def __call__(self,track, slice = None):
        
        transcript_counts = collections.defaultdict( set )
        counts = []
        for gtf in GTF.iterator(IOTools.openFile(self.getFilename(track))):
            transcript_counts[gtf.gene_id].add(gtf.transcript_id)
        for gene, transcripts in transcript_counts.iteritems():
            counts.append(len(transcripts))
        return counts

#-------------------------------------------------------------------                
class TranscriptNumberCumFreq(TranscriptNumbers):    

    def __call__(self,track, slice = None):
        
        transcript_counts = collections.defaultdict( set )
        counts = []
        for gtf in GTF.iterator(IOTools.openFile(self.getFilename(track))):
            transcript_counts[gtf.gene_id].add(gtf.transcript_id)
        for gene, transcripts in transcript_counts.iteritems():
            counts.append(len(transcripts))
        count, lower, dx, _ = scipy.stats.cumfreq(counts, numbins=40, defaultreallimits=(1,15))
        x = np.arange(count.size) * dx + lower
        return odict( (("transcript number", x), ("cumulative frequency", count/len(counts))) )




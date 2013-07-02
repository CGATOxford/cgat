"""
Orthologs.py - tools to deal with Leo's orthology pipeline.
===========================================================

:Author: 
:Release: $Id$
:Date: |today|
:Tags: Python

"""

import string, re

## patch to work around wrong gene assignments for BDGP
global_ignore_bdgp_genes = False

def FilterBDGP( data ):
    """remove genes that are

    transcripts from BDGP and genes from ENSEMBL.

    If only BDGP genes are present, keep them.
    
    """
    rx = re.compile( "CG.+-R." )
    n = []
    
    for d in data:
                    
        if d.mSchema == "dmel_vs_dmel":
            if (rx.match(d.mGene) and not rx.match(d.mTranscript)) or \
                   (not rx.match(d.mGene) and not rx.match(d.mTranscript)):
                continue
        n.append( d )
        
    return n

def GetGenes( transcripts ):
    """from a list of transcripts get genes.
    """
    genes = {}
    for t in transcripts:
        if t.mGene not in genes: genes[t.mGene] = []
        genes[t.mGene].append( t )
        
    return genes

class Transcript:

    mSeparator = "|"
    def __init__( self, a = None):
        if a:
            (self.mSchema, \
             self.mTranscript, \
             self.mGene, \
             self.mQuality) = a.split(self.mSeparator)

    def __str__(self):
        return self.mSeparator.join( (self.mSchema, self.mTranscript, self.mGene, self.mQuality) )
    
def ReadInterpretation( infile, separator,
                        genome1 = None, genome2 = None,
                        filter_restrict_genes1 = {},
                        filter_restrict_genes2 = {},
                        filter_remove_transcripts1 = {},
                        filter_remove_transcripts2 = {},
                        filter_restrict_transcripts1 = {},
                        filter_restrict_transcripts2 = {},
                        ):
    """read interpretation file.
    """
    pairs = []
    neliminated = 0
    for line in infile:
        if line[0] == "#": continue
        
        data = line[:-1].split("\t")
        weight = float(data[0])

        transcripts = []
        ## convert
        # Patch for invalid gene (0). These should not be here in the
        # first place.
        for x in range(1,len(data)):
            t = Transcript( data[x] )
            
            if t.mGene != "0":
                transcripts.append( t )

        if len(transcripts) < 2:
            continue

        first_prefix = transcripts[0].mSchema
        last_prefix  = transcripts[-1].mSchema
        for sep in range(1, len(transcripts)):
            if first_prefix != transcripts[sep].mSchema:
                transcripts1 = transcripts[:sep]
                transcripts2 = transcripts[sep:]
                break
        else:
            continue

        ## revert matches, so that the order is genome1, genome2
        if genome1 and genome2 and first_prefix == genome2 and last_prefix == genome1:
            transcripts2, transcripts1 = transcripts1, transcripts2
            first_prefix, last_prefix = last_prefix, first_prefix

        ## filter if genome don't fit.
        if first_prefix and last_prefix and genome1 and genome2 and\
               first_prefix != genome1 and last_prefix != genome2:
            continue

        if global_ignore_bdgp_genes:
            transcripts1 = FilterBDGP( transcripts1 )
            transcripts2 = FilterBDGP( transcripts2 )            

        ## apply various filters
        if filter_restrict_genes1: transcripts1 = filter( lambda x: x.mGene in filter_restrict_genes1, transcripts1)
        if filter_restrict_genes2: transcripts2 = filter( lambda x: x.mGene in filter_restrict_genes2, transcripts2)
        if filter_remove_transcripts1: transcripts1 = filter( lambda x: x.mTranscript not in filter_remove_transcripts1, transcripts1)
        if filter_remove_transcripts2: transcripts2 = filter( lambda x: x.mTranscript not in filter_remove_transcripts2, transcripts2)

        if filter_restrict_transcripts1: transcripts1 = filter( lambda x: x.mTranscript in filter_restrict_transcripts1, transcripts1)
        if filter_restrict_transcripts2: transcripts2 = filter( lambda x: x.mTranscript in filter_restrict_transcripts2, transcripts2)

        genes1 = GetGenes( transcripts1 )
        genes2 = GetGenes( transcripts2 )        

        if len(transcripts1) == 0 and len(transcripts2) == 0:
            neliminated += 1
            continue
        
        pairs.append( ( transcripts1, transcripts2, genes1, genes2, weight) )

    return pairs

def ReadOrphans( infile, separator,
                 genome1 = None, genome2 = None,
                 filter_restrict_genes1 = {},
                 filter_restrict_genes2 = {},
                 filter_remove_transcripts1 = {},
                 filter_remove_transcripts2 = {},
                 filter_restrict_transcripts1 = {},
                 filter_restrict_transcripts2 = {},
                 ):
    """read interpretation file.
    """
    pairs = []
    neliminated = 0
    for line in infile:
        if line[0] == "#": continue
        
        data = line[:-1].split("\t")

        transcripts = []
        ## convert
        # Patch for invalid gene (0). These should not be here in the
        # first place.
        for x in range(0,len(data)):
            t = Transcript( data[x] )
            
            if t.mGene != "0":
                transcripts.append( t )

        schema = transcripts[0].mSchema
        
        ## apply various filters        
        if schema == genome1:
            if filter_restrict_genes1: transcripts = filter( lambda x: x.mGene in filter_restrict_genes1, transcripts)
            if filter_remove_transcripts1: transcripts = filter( lambda x: x.mTranscript not in filter_remove_transcripts1, transcripts)            
        elif schema == genome2:
            if filter_restrict_genes2: transcripts = filter( lambda x: x.mGene in filter_restrict_genes2, transcripts)
            if filter_remove_transcripts2: transcripts = filter( lambda x: x.mTranscript not in filter_remove_transcripts2, transcripts)            
        else:
            continue

        if len(transcripts) == 0:
            neliminated += 1
            continue

        genes = GetGenes( transcripts )

        if schema == genome1:
            pairs.append( ( transcripts, {}, genes, {}, 0) )
        elif schema == genome2:
            pairs.append( ( {}, transcripts, {}, genes, 0) )
            
    return pairs


##########################################################################
def ClusterOrthologsByGenes( orthologs ):
    """cluster orthologs by genes.

    if an orthologous cluster contains the same genes in
    either species, they are merged.
    """

    map_gene2cluster = {}
    clusters = []

    for t1, t2, g1, g2, w in orthologs:

        ## get previous clusters
        pos = {}
        for g in g1.keys():
            k = "a%s" % g 
            if k in map_gene2cluster: 
                pos[map_gene2cluster[k]] = 1
        for g in g2.keys():
            k = "b%s" % g 
            if k in map_gene2cluster: 
                pos[map_gene2cluster[k]] = 1

        tt1 = t1
        tt2 = t2
        ww  = [w]
        ## add previous clusters to this cluster
        ## and clear previous clusters
        for p in pos:
            tt1 += clusters[p][0]
            tt2 += clusters[p][1]
            ww  += clusters[p][2]
            clusters[p] = None

        ## map previous clusters to new cluster
        n = len(clusters)
        for g in GetGenes( tt1 ):
            map_gene2cluster["a%s" % g] = n
        for g in GetGenes( tt2 ):
            map_gene2cluster["b%s" % g] = n

        ## append to clusters
        clusters.append( [tt1, tt2, ww] )            

    orthologs = []
    for c in clusters:
        if c:
            orthologs.append( (c[0], c[1], GetGenes(c[0]), GetGenes(c[1]), c[2]) )

    return orthologs
    

def GetDegeneracy( t1, t2 ):
    """get degeneracy of orthology assignments.

    given are two lists of transcripts.
    returns a tuple with gene and transcript degeneracy code, respectively.
    """
    g1 = GetGenes( t1 )
    g2 = GetGenes( t2 )        

    return "%i:%i" % (len(g1), len(g2)), "%i:%i" % ( len(t1), len(t2) )

    

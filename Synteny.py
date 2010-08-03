################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
'''
Synteny.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Code
----

'''
import string, re, types
### Tools to analyze synteny
import Genomics

###################################################################
        
class Ortholog:
    def __init__(self, ortholog_id, gene_id, contig, strand, first_res, last_res):
        (self.mOrthologId, self.gene_id, self.contig, self.strand,
         self.mFrom, self.mTo) = ortholog_id, gene_id, contig, strand, first_res, last_res
        self.mRank = 0

    def __str__(self):
        return ";".join(map(str, ( self.mOrthologId, self.gene_id, self.contig, self.strand,
                                   self.mFrom, self.mTo, self.mRank)) )
    
###################################################################
def ReadOrthologs( infile ):
    """read orthologs from file."""
    orthologs = {}
    
    for line in infile:
        if line[0] == "#": continue

        ortholog_id, gene_id, contig, strand, first_res, last_res = line[:-1].split("\t")[:6]

        if ortholog_id not in orthologs:
            orthologs[ortholog_id] = []

        orthologs[ortholog_id].append( Ortholog(ortholog_id,
                                                gene_id, contig, strand,
                                                int(first_res), int(last_res)) )
        
    return orthologs

###################################################################
def ReadOrthologsPerSpecies( infile ):
    """read orthologs from file with species information."""
    orthologs = {}

    read = set()
    
    for line in infile:
        
        if type(line) in (types.TupleType, types.ListType):
            ortholog_id, schema, gene_id, contig, strand, first_res, last_res = line
        else:
            if line[0] == "#": continue
            ortholog_id, schema, gene_id, contig, strand, first_res, last_res = line[:-1].split("\t")[:7]

        key = "%s-%s" % (schema,gene_id)
        if key in read: continue
        read.add(key)

        if schema not in orthologs:
            orthologs[schema] = {}

        o = orthologs[schema]
        if ortholog_id not in o:
            o[ortholog_id] = []

            
        o[ortholog_id].append( Ortholog(ortholog_id,
                                        gene_id, contig, strand,
                                        int(first_res), int(last_res)) )
        
    return orthologs

###################################################################
def FilterJunk( orthologs ):
    """remove assignments to junk contigs.
    """

    for id, oo in orthologs.items():
        oo = filter( lambda x: not Genomics.IsJunk( x.contig ), oo )
        if len(oo) > 0:
            orthologs[id] = oo
        else:
            del orthologs[id]

###################################################################
def CleanOrthologs( orthologs1, orthologs2, filter_junk = False ):
    """only keep orhtologs that are present in both sets."""

    if filter_junk:
        FilterJunk( orthologs1 )
        FilterJunk( orthologs2 )        

    for id, oo in orthologs1.items():
        if id not in orthologs2:
            del orthologs1[id]

    for id, oo in orthologs2.items():
        if id not in orthologs1:
            del orthologs2[id]

###################################################################
def GetContigs( orthologs ):
    """get map of contigs to orthologs.

    An ortholog can be part of only one contig, but the same ortholog_id can
    be part of several contigs.
    """
    
    contigs = {}
    for id, oo in orthologs.items():
        for o in oo:
            if o.contig not in contigs: contigs[o.contig] = []
            contigs[o.contig].append( o )
            
    return contigs

###################################################################
def GetMap( values ):
    map_a2b = {}
    map_b2a = []

    for v in values:
        map_a2b[v] = len(map_b2a)
        map_b2a.append(v)

    return map_a2b, map_b2a

###################################################################
def SetRankToPosition( orthologs, map_contig2orthologs, relative = True ):
    """set rank of ortholog to position within contig."""

    rank = 0    
    for contig, oo in map_contig2orthologs.items():
        oo.sort( lambda x,y: cmp( x.mFrom, y.mFrom) )
        if relative: rank = 0
        for o in oo:
            o.mRank = rank
            rank += 1

###################################################################
def GetContigCounts( orthologs1, orthologs2 ):
    """count number of orthologs per contig pair.
    """
    counts = {}
    
    for id, oo1 in orthologs1.items():
        
        if id not in orthologs2: continue
            
        conts1 = map( lambda x: x.contig, oo1 )
        conts2 = map( lambda x: x.contig, orthologs2[id] )

        for c1 in conts1:
            for c2 in conts2:
                key = (c1, c2)
                if key not in counts: counts[key] = 0
                counts[key] += 1
                
    return counts

###################################################################
def GetContigContents( orthologs1, orthologs2 ):
    """count number of orthologs per contig pair.
    """
    counts = {}
    
    for id, oo1 in orthologs1.items():
        
        if id not in orthologs2: continue
            
        conts1 = map( lambda x: x.contig, oo1 )
        conts2 = map( lambda x: x.contig, orthologs2[id] )

        for c1 in conts1:
            for c2 in conts2:
                key = (c1, c2)
                if key not in counts: counts[key] = []
                counts[key].append( id )
                
    return counts

###################################################################
def SortAssignOrthologs( orthologs1, orthologs2 ):
    """based on the number of shared orthologs, sort the contigs in
    both ortholog sets."""

    contigs1 = GetContigs( orthologs1 )
    contigs2 = GetContigs( orthologs2 )    

    ## Full matrix too large, when I have many contigs.
    ## Thus I use a sparse matrix and do some straightforward
    ## sorting
    counts = GetContigCounts( orthologs1, orthologs2 )

    llinks = []
    for key in counts:
        llinks.append( (counts[key], key) )

    llinks.sort()
    llinks.reverse()

    sorted_list1 = []
    sorted_list2 = []

    sorted1 = {}
    sorted2 = {}

    for n, key in llinks:
        
        ## do not take empty contigs
        if n == 0: continue
        
        a,b = key

        if a not in sorted1:
            sorted1[a] = True
            sorted_list1.append(a)
            
        if b not in sorted2:
            sorted2[b] = True
            sorted_list2.append(b)
            
    return sorted_list1, sorted_list2, contigs1, contigs2

###################################################################
def GetSortedOrthologs( orthologs ):
    """sort orthologs by contig and rank."""

    sorted_orthologs = []
    
    for ortholog_id, oo in orthologs.items():
        sorted_orthologs += oo
        
    sorted_orthologs.sort( lambda x,y: cmp(x.mRank, y.mRank) )
    
    return sorted_orthologs

###################################################################
class Block:
    def __init__(self):
        self.mMembers1 = []
        self.mMembers2 = []        
        self.mBlockId = 0
        self.mBreakers1 = []
        self.mBreakers2 = []        
        
    def __str__(self):
        return "%i\t%i\t%s\t%i\t%i\t%i\t%i\t%i\t%s\t%i\t%i\t%i\t%i" % (self.mBlockId,
                                                               len(self.mMembers1),
                                                               self.contig1,
                                                               self.mFrom1,
                                                               self.mTo1,
                                                               self.mTo1 - self.mFrom1,
                                                               len(self.mBreakers1),                                                               
                                                               len(self.mMembers2),
                                                               self.contig2,
                                                               self.mFrom2,
                                                               self.mTo2,
                                                               self.mTo2 - self.mFrom2,
                                                               len(self.mBreakers2))
    def getHeader( self ):
        return "\t".join(("blockid", "nmem1", "contig1", "from1", "to1", "len1", "nbreaks1",
                         "nmem2", "contig2", "from2", "to2", "len2", "nbreaks2"))

    def cleanMembers( self, members ):
        """remove duplicates from members list."""
        taken = {}
        new_members = []
        for o in members:
            if o.mRank in taken: continue
            
            new_members.append( o )
            taken[o.mRank] = True
                
        return list( new_members )
        
    def cleanBreakers( self, members, breakers ):
        """remove small inversion from the list of synteny breakers."""
        n2 = set( map(lambda x: x.mRank, members) )
        taken = {}
        new_members = []
        for o in breakers:
            if o.mRank in taken: continue
            
            if o.mRank not in n2:
                new_members.append( o )
                taken[o.mRank] = True
                
        return list( new_members )
            
    def update( self ):
        if self.mMembers1:
            self.mMembers1 = self.cleanMembers( self.mMembers1 )
            self.contig1 = self.mMembers1[0].contig
            self.mFrom1 = min(map(lambda x: x.mFrom, self.mMembers1) )
            self.mTo1 = max(map(lambda x: x.mTo, self.mMembers1) )
            self.mBreakers1 = self.cleanBreakers( self.mMembers1, self.mBreakers1 )            
        if self.mMembers2:
            self.mMembers2 = self.cleanMembers( self.mMembers2 )            
            self.contig2 = self.mMembers2[0].contig
            self.mFrom2 = min(map(lambda x: x.mFrom, self.mMembers2) )
            self.mTo2 = max(map(lambda x: x.mTo, self.mMembers2) )
            self.mBreakers2 = self.cleanBreakers( self.mMembers2, self.mBreakers2 )
            
###################################################################
"""The syntenic relationships of the D. melanogaster protein to D. pseudoobscura sequence anchor points described by these BLAST
hits were then manually inspected to define and refine the order of syntenic blocks. Synteny blocks were defined as runs of consecutive
D. melanogaster protein sequence?D. pseudoobscura genomic sequence pairs. Within a syntenic block, gaps were permitted (since there are
genes that fall into sequence gaps) and an occasional gene out of order was also permitted (if it fell within five genes of its expected location).
Gene duplications in one species could be perceived as synteny breaks. Gene duplications were not considered in the derivation of synteny blocks.

Figure 1 shows the synteny blocks of the chromosomes are short and extremely mixed, but the great majority of syntenic sequences are
found on the same Muller element in D. pseudoobscura as they are in D. melanogaster. Thus, as expected, the majority of the chromosomal
rearrangements between the D. pseudoobscura and D. melanogaster lineages have been confined to related chromosome arms. The average number
of D. melanogaster genes in a syntenic block is 10.7, corresponding to ~83 kb. The length distribution of syntenic blocks on different
Muller elements is shown in Supplemental Figure S2.
"""

def GetSyntenicOrthologs( orthologs, last_contig, last_rank,
                          max_synteny_distance = 5 ):
    """get orthologs that are syntenic to previous one."""
    oo = []
    for o in orthologs:
        if o.contig == last_contig and \
               abs(o.mRank - last_rank) <= max_synteny_distance:
            oo.append(o)
    return oo

def GetSyntenyBlocks( orthologs1, orthologs2,
                      max_synteny_distance = 5,
                      max_look_ahead = 1,
                      loglevel = 0 ):
    """compute synteny blocks between orthologs.

    max_synteny_distance: number of rearrangements and deletions that are to be
    tolerated.

    This is the distance between adjacent orthologs, minimum should be 1.
    The ones causeing the rearrangement/deletion are saved in
    block.mRearrangements.

    max_look_ahead: if there is a synteny break, look ahead these number of
    orthologs in species 1 to see if the block continues.
    """

    ## sort orthologs1 by contig and rank and group into
    ## micro-synteny blocks

    sorted_orthologs = GetSortedOrthologs( orthologs1 )
    sorted_orthologs2 = GetSortedOrthologs( orthologs2 )    

    blocks = []

    last_contig1 = None
    last_contig2 = None

    last_rank2 = None

    block = None

    no1 = 0
    while no1 < len(sorted_orthologs):

        o1 = sorted_orthologs[no1]
        test_oo2 = orthologs2[o1.mOrthologId]
        
        if loglevel >= 2:
            print "# last=", last_contig1, last_contig2, last_rank2
            print "# o1=", str(o1)
            print "# o2=", map(str, orthologs2[o1.mOrthologId])
        
        # find the syntenic ones within the other orthologs
        # In case of local duplication, this might be several
        found_oo2 = GetSyntenicOrthologs( test_oo2,
                                          last_contig2, last_rank2,
                                          max_synteny_distance )
        ## if no match has been found
        if o1.contig != last_contig1 or \
               len(found_oo2) == 0:

            ## check for lookahead:
            for x in range(no1 + 1, min(no1 + 1 + max_look_ahead, len(sorted_orthologs) )):
                
                new_oo2 = GetSyntenicOrthologs( orthologs2[sorted_orthologs[x].mOrthologId],
                                                last_contig2, last_rank2,
                                                max_synteny_distance )
                ## if something is found
                ## save orthologs causing the breaks
                ## advance index
                if len(new_oo2) > 0:
                    for o2 in test_oo2:
                        if loglevel >= 3:
                            print "## added breaker to 2: %s" % str(o2)
                        block.mBreakers2.append( o2 )
                    no1 = x - 1
                    found_oo2 = new_oo2
                    last_rank2 = found_oo2[0].mRank                    
                    break
                
            else:
                ## if lookahead does not continue break, continue
                if block:
                    if loglevel >= 2:
                        print "############saving cluster %i#####################################" % block.mBlockId
                    
                    block.update()
                    blocks.append( block )

                block = Block()
                block.mBlockId = len(blocks)
                    
                found_oo2 = orthologs2[o1.mOrthologId]
                last_rank2 = found_oo2[0].mRank
        else:
            ## save rank for synteny breakers that are skipped over.
            ## find minimum distance to last_rank2
            min_rank = min(map(lambda x: x.mRank, found_oo2 ))
            max_rank = max(map(lambda x: x.mRank, found_oo2 ))

            if last_rank2 < min_rank:
                for x in range( last_rank2 + 1, min_rank):
                    if loglevel >= 3:
                        print "# added skippers to 2: ", sorted_orthologs2[x]
                    block.mBreakers2.append( sorted_orthologs2[x] )
                last_rank2 = max_rank
            elif last_rank2 > max_rank:
                for x in range( max_rank+1, last_rank2):
                    if loglevel >= 3:
                        print "# added skippers to 2: ", sorted_orthologs2[x]
                    block.mBreakers2.append( sorted_orthologs2[x] )
                last_rank2 = min_rank
            else:
                ## do nothing, if rank within middle (happens in m:m clusters)
                pass

        block.mMembers1.append( o1 )
        for o2 in found_oo2:
            block.mMembers2.append( o2 )        
            
        last_contig1 = o1.contig
        last_contig2 = o2.contig
        no1 += 1

    return blocks

##--------------------------------------------------------------------------------
def Blocks2Orthologs( blocks ):
    """from a list of blocks, return two lists of orthologs.
    """
    orthologs1 = []
    orthologs2 = []
    
    for block in blocks:
        o1 = Ortholog( block.mBlockId, block.mBlockId, block.contig1, "+", block.mFrom1, block.mTo1)
        orthologs1.append(o1)
        o2 = Ortholog( block.mBlockId, block.mBlockId, block.contig2, "+", block.mFrom2, block.mTo2)
        orthologs2.append(o2)

    return orthologs1, orthologs2


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
mali2cleaned_mali.py - remove split genes from a multiple alignment
===================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Remove split genes from a multiple alignment. Split genes confuse
trees as both parts might not cluster together in the tree. This then
causes problems on tree-reconciliation.

Usage
-----

Example::

   python mali2cleaned_mali.py --help

Type::

   python mali2cleaned_mali.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import optparse

USAGE="""python %s [OPTIONS] < graph.in > graph.out

Version: $Id: mali2cleaned_mali.py 2781 2009-09-10 11:33:14Z andreas $



""" % sys.argv[0]

import CGAT.Experiment as E
import CGAT.Mali as Mali
import CGAT.Exons as Exons

import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import networkx

###########################################################
def AddSynonym( rep_transcript, mem_transcript, map_synonyms, map_to_synonyms):

    if rep_transcript not in map_synonyms:
        map_synonyms[rep_transcript] = []

    if mem_transcript in map_synonyms:            

        for s in map_synonyms[mem_transcript]:
            map_to_synonyms[s] = rep_transcript

        map_synonyms[rep_transcript] += map_synonyms[mem_transcript]                    
        del map_synonyms[mem_transcript]

    map_synonyms[rep_transcript].append( mem_transcript )
    map_to_synonyms[mem_transcript] = rep_transcript

###########################################################
def filterJoiningTranscripts( links, options ):
    """cluster links in infile by transcript."""
    
    map_priority = {}
    for x in range( len(options.quality_priority) ):
        map_priority[options.quality_priority[x]] = x
        
    transcripts = {}
    ninput, noutput, nskipped = 0, 0, 0

    ## build matrix
    matrix = {}
    
    for transcript1, transcript2, weight in links:

        weight = float(weight)
        ninput += 1
        
        if options.min_weight != None and weight < options.min_weight:
            nskipped += 1
            continue

        if options.max_weight != None and weight > options.max_weight:
            nskipped += 1
            continue

        species1 = transcript1.split(options.separator)[0]
        species2 = transcript2.split(options.separator)[0]

        if transcript1 > transcript2: transcript1, transcript2 = transcript2, transcript1
        
        if species1 not in transcripts: transcripts[species1] = []
        if species2 not in transcripts: transcripts[species2] = []
        
        transcripts[species1].append( transcript1 )
        transcripts[species2].append( transcript2 )    

        if transcript1 not in matrix: matrix[transcript1] = {}
        if transcript2 not in matrix: matrix[transcript2] = {}        
        
        matrix[transcript1][transcript2] = weight
        matrix[transcript2][transcript1] = weight        

        noutput += 1
        
    ## select each transcripts for each species and go through list
    ## would be more efficient if I had the length

    
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped==%i\n" % (ninput, noutput, nskipped))

##-------------------------------------------------------------------------------
def hasGenomicOverlap( key1, key2, exons ):
    """returns true if two transcripts overlap on the genome."""
    ee = exons[key1]
    contig1, strand1, start1, end1 = ee[0].mSbjctToken, ee[0].mSbjctStrand, ee[0].mGenomeFrom, ee[-1].mGenomeTo
    ee = exons[key2]
    contig2, strand2, start2, end2 = ee[0].mSbjctToken, ee[0].mSbjctStrand, ee[0].mGenomeFrom, ee[-1].mGenomeTo

    return contig1 == contig2 and strand1 == strand2 and min(end1,end2) - max(start1, start2) > 0

##-------------------------------------------------------------------------------
def isSynchronous( key1, key2, mali, exons ):
    """check if two transcripts are on synchronous positions.

    Two transcripts are synchroneous, if their position in the
    alignment and on the genomic sequence are both a < b.

    Both front and end are taken into account. Note that they
    can still be overlapping.
    """

    def getStartEnd( s ):
        start, end = 0, len(s)
        x = re.match("([-.]+)", s )
        if x: start = len(x.groups()[0])
        x = re.search("([-.]+)$", s )
        if x: end -= len(x.groups()[0])        
        
        return start, end

    mstart1, mend1 = getStartEnd( mali[key1] )
    mstart2, mend2 = getStartEnd( mali[key2] )    
    
    gstart1, gend1 = exons[key1][0].mGenomeFrom, exons[key1][-1].mGenomeTo
    gstart2, gend2 = exons[key2][0].mGenomeFrom, exons[key2][-1].mGenomeTo    

    return (mstart1 < mstart2 and mend1 < mend2 and gstart1 < gstart2 and gend1 < gend2) or \
           (mstart2 < mstart1 and mend2 < mend1 and gstart2 < gstart1 and gend2 < gend1)

##-------------------------------------------------------------------------------
def getGenomicDistance( key1, key2, exons ):
    """returns the number of nucleotides between two genes.

    return negative value if they are on different contigs/strands or
    are overlapping.
    """
    ee = exons[key1]
    contig1, strand1, start1, end1 = ee[0].mSbjctToken, ee[0].mSbjctStrand, ee[0].mGenomeFrom, ee[-1].mGenomeTo
    ee = exons[key2]
    contig2, strand2, start2, end2 = ee[0].mSbjctToken, ee[0].mSbjctStrand, ee[0].mGenomeFrom, ee[-1].mGenomeTo

    if contig1 == contig2 and strand1 == strand2:
        return -(min(end1,end2) - max(start1, start2))
    else:
        return -1

##-------------------------------------------------------------------------------
def getPercentOverlap( seq1, seq2, gap_chars = (".", "-") ):
    """calculate overlap distance between two sequences."""
    naligned, nunaligned = 0, 0
    for c1, c2 in zip( seq1.upper(), seq2.upper()):
        if c1 in gap_chars or c2 in gap_chars:
            if not (c1 in gap_chars and c2 in gap_chars):
                nunaligned += 1
        else:
            naligned += 1
        
    return 100.0 * float( naligned) / float(nunaligned + naligned)

##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------    
def removeJoiningTranscripts( mali, exons, map_species2transcripts, options ):
    
    ###############################################################
    ###############################################################
    ###############################################################        
    ## identify all joining transcripts
    ###############################################################
    joining_transcripts = set()
    
    for species, transcripts in map_species2transcripts.items():

        transcripts.sort()
        transcripts.reverse()

        for t1 in range(len(transcripts)-1):
            transcript1 = transcripts[t1][1]

            for t2 in range(t1+1, len(transcripts)):
                transcript2 = transcripts[t2][1]
                
                if not hasGenomicOverlap( transcript1, transcript2, exons ):
                    continue
                
                for t3 in range(t2+1, len(transcripts)):
                    transcript3 = transcripts[t3][1]

                    if not hasGenomicOverlap( transcript1, transcript3, exons ):
                        continue

                    if hasGenomicOverlap( transcript2, transcript3, exons ):
                        continue
                    
                    overlap12 = getPercentOverlap( mali[transcript1], mali[transcript2] )
                    overlap13 = getPercentOverlap( mali[transcript1], mali[transcript3] )                    
                    overlap23 = getPercentOverlap( mali[transcript2], mali[transcript3] )

                    if options.loglevel >= 2:
                        options.stdlog.write("# transcript: %s connects %s and %s - mali overlap: %i %i %i\n" % (transcript1, transcript2, transcript3, overlap12, overlap13, overlap23) )
                        joining_transcripts.add( (transcript1, transcript2, "joining transcript") )
                        
    return joining_transcripts

##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------        
def removeSplitTranscripts( mali, exons, map_species2transcripts, options ):
    """remove split transcripts from mulitple alignment.

    This procedure finds two genes of the same species and if they do not overlap
    on the multiple alignment will remove the shorter one. The algorithm works like this

    1. Sort all transcripts per species by length
    2. In descending length:
       Check overlap with all other transcripts of same species and shorter length
       If there is no overlap, remove the shorter fragment.

    Further filters include:
       * remove the shorter gene, if two genes overlap on the genome
         This is like a sanity check.
       
    """

    ###############################################################
    ###############################################################
    ###############################################################        
    ## identify all split transcripts
    ###############################################################
    removed_transcripts = set()

    if options.loglevel >= 2:
        options.stdlog.write("# removeSplitTranscripts started.\n" )
    
    for species, transcripts in map_species2transcripts.items():

        transcripts.sort()
        transcripts.reverse()

        for t1 in range(len(transcripts)-1):
            transcript1 = transcripts[t1][1]

            for t2 in range(t1+1, len(transcripts)):
                transcript2 = transcripts[t2][1]
                
                overlap12 = getPercentOverlap( mali[transcript1], mali[transcript2] )

                if options.loglevel >= 2:
                    options.stdlog.write("# %s %s: overlap=%5.2f\n" % (transcript1,transcript2,overlap12))
                
                if exons and hasGenomicOverlap( transcript1, transcript2, exons ):
                    if options.loglevel >= 1:
                        options.stdlog.write("# transcript: %s removed because genomic overlap with %s\n" % (transcript2, transcript1) )
                    removed_transcripts.add( (transcript2, transcript1, "genomic overlap") )
                    continue
                
                if exons:
                    
                    genomic_distance = getGenomicDistance( transcript1, transcript2, exons )
                    is_synchronous = isSynchronous( transcript1, transcript2, mali, exons )
                    
                    if genomic_distance > 0 and genomic_distance < options.min_genomic_distance and \
                           overlap12 < options.min_percent_overlap and \
                           is_synchronous:
                        if options.loglevel >= 1:
                            options.stdlog.write("# transcript: %s removed because genomic proximity to %s: %i/%5.2f\n" % (transcript2, transcript1, genomic_distance, overlap12) )
                        removed_transcripts.add( (transcript2, transcript1, "genomic proximity: %i/%5.2f%%" % (genomic_distance, overlap12) ) )
                        continue
                    
                if overlap12 < options.max_percent_overlap:
                    if options.loglevel >= 1:
                        options.stdlog.write("# transcript: %s removed because no overlap with %s\n" % (transcript2, transcript1) )
                    removed_transcripts.add( (transcript2, transcript1, "no mali overlap" ) )
                    continue
                    
    return removed_transcripts

##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------
##-------------------------------------------------------------------------------    
if __name__ == "__main__":

    parser = E.OptionParser( version = "%prog version: $Id: mali2cleaned_mali.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-m", "--genome-master", dest="genome_master", type="string",
                      help="genome to use as master." )

    parser.add_option("-s", "--filename-removed", dest="filename_removed", type="string",
                      help="output filename for deleted entries." )

    parser.add_option("-e", "--filename-exons", dest="filename_exons", type="string",
                      help="filename on where to exon information."  )

    parser.add_option("-u", "--filename-summary", dest="filename_summary", type="string",
                      help="output filename of component summary." )

    parser.add_option("-c", "--filename-components", dest="filename_components", type="string",
                      help="output filename for components." )

    parser.add_option("--min-percent-overlap", dest="min_percent_overlap", type="float",
                      help="minimum percent overlap for splitting multiple alignment into components.")

    parser.add_option("--max-percent-overlap", dest="max_percent_overlap", type="float",
                      help="maximum percent overlap for split genes.")

    parser.add_option("--min-genomic-distance", dest="min_genomic_distance", type="int",
                      help="minimum genomic distance for adjacent genes to be considered dodgy.")

    parser.add_option("-o", "--mode", dest="mode", type="choice", 
                      choices=("joining", "split"),
                      help="""how to filter the alignment.
                      joining: remove joining transcripts (spindly genes)
                      split:  remove split transcripts""")

    parser.add_option("-g", "--gene-mode", dest="gene_mode", action="store_true",
                      help="""the aligned sequences are genes. This forces the exon boundaries to
                      collated by genes.""")
    
    parser.set_defaults( \
        genome_master = None,
        filename_removed = None,
        filename_components = None,        
        filename_summary = None,
        filename_exons = None,
        mode="joining",
        input_format = "fasta",
        output_format = "fasta", 
        max_percent_overlap = 0,
        min_percent_overlap = 0,
        gene_mode = False,
        separator = "|")
        
    (options, args) = E.Start( parser )

    ###############################################################
    ###############################################################
    ###############################################################        
    ## input
    ###############################################################                


    mali = Mali.Mali()
    mali.readFromFile( sys.stdin, format = options.input_format )
    all_identifiers = mali.getIdentifiers()

    if options.filename_exons:
        ## read exon boundaries and keep forward coordinates
        
        if options.gene_mode:
            exons = Exons.ReadExonBoundaries( open(options.filename_exons, "r"),
                                              from_zero = True)
            
            gene_exons = {}
            for id, ee in exons.items():
                data = id.split(options.separator)
                new_id = options.separator.join((data[0], data[2]))
                if new_id not in gene_exons: gene_exons[new_id] = []
                for e in ee:
                    e.mQueryToken = new_id
                gene_exons[new_id] += ee
            for id, ee in gene_exons.items():
                ee.sort( lambda x,y: cmp( x.mGenomeFrom, y.mGenomeFrom ) )
            exons = gene_exons

        else:
            exons = Exons.ReadExonBoundaries( open(options.filename_exons, "r"),
                                              filter=set(all_identifiers),
                                              from_zero = True)
            
    else:
        exons = {}

    ###############################################################
    ###############################################################
    ###############################################################        
    ## collect all transcripts for a species together with their
    ## aligned length
    ###############################################################                
    map_species2transcripts = {}
    
    for id in mali.getIdentifiers():
        data = id.split(options.separator)

        species = data[0]

        if exons:
            l = exons[id][-1].mGenomeTo - exons[id][0].mGenomeFrom
        else:
            l = len(mali.getEntry(id).getSequence())
        
        try:
            map_species2transcripts[species].append( (l, id ) )
        except KeyError:
            map_species2transcripts[species]= [( l, id) ]

    if options.mode == "joining":
        mapped_transcripts = removeJoiningTranscripts( mali, exons, map_species2transcripts, options )
        
    elif options.mode == "split":
        mapped_transcripts = removeSplitTranscripts( mali, exons, map_species2transcripts, options )

    ###############################################################
    ###############################################################
    ###############################################################        
    ## now build overlap graph of remaining sequences split multiple
    ## alignment in components.
    ## Compute reciprocal best match graph
    ###############################################################                
    graph = networkx.Graph()

    removed_transcripts = set(map(lambda x: x[0], mapped_transcripts ))

    for t in all_identifiers:
        if t not in removed_transcripts:
            graph.add_node( t )

    for t1 in range(len(all_identifiers)-1):
        transcript1 = all_identifiers[t1]
        if transcript1 in removed_transcripts: continue
        
        for t2 in range(t1+1, len(all_identifiers)):
            transcript2 = all_identifiers[t2]            
            if transcript2 in removed_transcripts: continue

            overlap = getPercentOverlap( mali[transcript1], mali[transcript2] )
            if overlap > 5:
                graph.add_edge( transcript1, transcript2 )
                
    ## compute components
    components = networkx.connected_components( graph )

    ###############################################################
    ###############################################################
    ###############################################################        
    ## output
    ###############################################################                
    if options.filename_components:
        n = 1
        outfile = open(options.filename_components, "w" )
        
        outfile.write( "id\tcomponent\n" )    
        for component in components:
            for c in component:
                outfile.write( "%s\t%i\n" % (c, n))
            n += 1
        outfile.close()
        
    if options.filename_removed and len(removed_transcripts) > 0:
        outfile = open(options.filename_removed, "w")
        outfile.write("removed\trepresentative\treason\n")
        for removed_transcript, rep_transcript, reason in mapped_transcripts:
            outfile.write("%s\t%s\t%s\n" % (removed_transcript, rep_transcript, reason) )
        outfile.close()

    if options.filename_summary:
        n = 1
        outfile = open(options.filename_summary, "w")
        outfile.write("component\tsize\tnspecies\tnmaster\n" )
        for component in components:
            species = map( lambda x: x.split(options.separator)[0], component)
            outfile.write( "%i\t%i\t%i\t%i\t%i\n" % ( n,
                                                      len(component),
                                                      len(species),
                                                      len(filter( lambda x: x == options.genome_master, species))))

            n += 1

    for transcript in removed_transcripts:
        mali.deleteEntry( transcript )

    new_identifiers = mali.getIdentifiers()
    
    mali.removeGaps( minimum_gaps = len( new_identifiers) )

    mali.writeToFile( options.stdout, format = options.output_format )
        
    if options.loglevel >= 1:
        options.stdlog.write( "# input=%i, output=%i, removed=%i, ncomponents=%i\n" % (len(all_identifiers), len(new_identifiers), len(removed_transcripts), len(components) ))
        options.stdlog.write( "# final component sizes: %s\n" % ",".join(map( lambda x: str(len(x)), components)))

    E.Stop()

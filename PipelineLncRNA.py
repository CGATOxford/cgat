'''
functions and classes for use with the lincRNA
pipeline
'''

import re, sys, os
import GTF
import IOTools
import gzip
import collections
import IndexedGenome
import Pipeline as P
import Experiment as E
import sqlite3

########################################################
# gene set building
########################################################
def buildCodingGeneSet(abinitio_coding, reference, outfile):
    '''
    takes the output from cuffcompare of a transcript
    assembly and filters for annotated protein coding
    genes. 
    
    NB "pruned" refers to nomenclature in the transcript
    building pipeline - transcripts that appear in at least
    two samples.
    
    Because an abinitio assembly will often contain
    fragments of known transcripts and describe them as 
    novel, the default behaviour is to produce a set that
    is composed of 'complete' transcripts
    '''
    inf = IOTools.openFile(abinitio_coding)
    outf = gzip.open(outfile, "w")

    coding = {}
    coding["protein_coding"] =  GTF.readAndIndex( GTF.iterator_filtered( GTF.iterator(IOTools.openFile(reference))
                                                                        , source="protein_coding" )
                                                                        , with_value = False  )

    for gtf in GTF.iterator(inf):
        if coding["protein_coding"].contains(gtf.contig, gtf.start, gtf.end):
            if gtf.class_code == "=":
                outf.write("%s\n" % str(gtf))
    outf.close()

#-------------------------------------------------------------------------------
def buildRefcodingGeneSet(coding_set, refcoding_set, outfile):
    '''
    takes genes from an ab initio assembly and filters a reference coding set
    for these genes. Allows for comparisons of known transcripts for those genes
    that are assembled ab initio. Does this by gene name
    '''
    keep_genes = set()
    for gtf in GTF.iterator(IOTools.openFile(coding_set)):
        keep_genes.add(gtf.gene_name)

    outf = gzip.open(outfile, "w")
    for gtf in GTF.iterator(IOTools.openFile(refcoding_set)):
        if gtf.gene_name in keep_genes:
            outf.write("%s\n" % gtf)
    outf.close()
    
#-------------------------------------------------------------------------------
def buildRefnoncodingGeneSet(reference, outfile):
    '''
    filter the refnoncoding geneset for things that are described in ensembl
    as being:
    Ambiguous_orf
    Retained_intron
    Sense_intronic
    antisense
    Sense_overlapping
    Processed transcript
    '''

    statement = '''zcat %(reference)s 
                   | awk '$2 == "lincRNA" || $2 == "non_coding" || $2 == "3prime_overlapping_ncrna" || $2 == "ncRNA_host"' | gzip > %(outfile)s'''                                                                                                                                             
    P.run()

#-------------------------------------------------------------------------------
def buildLncRNAGeneSet( abinitio_lincrna, reference, refnoncoding, pseudogenes_gtf, numts_gtf, outfile , min_length):
    '''
    build lncRNA gene set. 
    
    In contrast to the pipeline, this lincRNA set does not contain
    the reference noncoding gene set. It is transcripts in the abinitio set that
    do not overlap at any protein coding, processed or pseudogene transcripts 
    (exons+introns) in a reference gene set.
    lincRNA genes are often expressed at low level and thus the resultant transcript
    models are fragmentory. To avoid some double counting in downstream analyses
    transcripts overlapping on the same strand are merged - this does not neccessarily seem to work well if 
    not using the reference set.
    
    Transcripts need to have a length of at least 200 bp.

    '''
    
    infile_abinitio, reference_gtf, refnoncoding_gtf, pseudogenes_gtf, numts_gtf = \
    abinitio_lincrna, reference, refnoncoding, pseudogenes_gtf, numts_gtf

    E.info( "indexing geneset for filtering" )

    input_sections = ("protein_coding", 
                      "processed_pseudogene",
                      "unprocessed_pseudogene",
                      "nonsense_mediated_decay",
                      "retained_intron")
                    
    indices = {}
    for section in input_sections:
        indices[section] = GTF.readAndIndex( 
            GTF.iterator_filtered(  GTF.iterator( IOTools.openFile( reference_gtf ) ),
                                   source = section ),
            with_value = True )
        
    E.info( "built indices for %i features" % len(indices))

    indices["numts"] = GTF.readAndIndex( GTF.iterator( IOTools.openFile( numts_gtf) ), with_value = True )

    E.info( "added index for numts" )

    indices["pseudogenes"] = GTF.readAndIndex( GTF.iterator( IOTools.openFile( pseudogenes_gtf) ), with_value = True )

    E.info( "added index for pseudogenes" )

    noncoding = {}
    noncoding["noncoding"] = GTF.readAndIndex(GTF.iterator(IOTools.openFile(refnoncoding_gtf)), with_value = False)

    E.info("created index for known noncoding exons to avoid filtering")

    sections = indices.keys()

    total_transcripts, remove_transcripts = set(), collections.defaultdict( set )
    transcript_length = collections.defaultdict( int )
    inf = GTF.iterator( IOTools.openFile( infile_abinitio ) )

    E.info( "collecting genes to remove" )

    min_length = min_length

    for gtfs in GTF.transcript_iterator( inf ): 
        l = sum([x.end-x.start for x in gtfs])
        if l < min_length:
            remove_transcripts[gtfs[0].transcript_id].add("length")

        for section in sections:
            for gtf in gtfs:
                transcript_id = gtf.transcript_id   
                total_transcripts.add( transcript_id )
                if not noncoding["noncoding"].contains(gtf.contig, gtf.start, gtf.end):
                    if indices[section].contains( gtf.contig, gtf.start, gtf.end):
                        
                        # retain antisense transcripts
                        for gtf2 in indices[section].get(gtf.contig, gtf.start, gtf.end):
                            if gtf.strand == gtf2[2].strand:
                                remove_transcripts[transcript_id].add( section )
                                                    
                
    E.info( "removing %i out of %i transcripts" % (len(remove_transcripts), len(total_transcripts)) )
    
    outf = open("lincrna_removed.tsv", "w")
    outf.write("transcript_id" + "\t" + "removed" + "\n" )
    for x, y in remove_transcripts.iteritems():
        outf.write("%s\t%s\n" % (x, ",".join(y)))

    # write out transcripts that are not in removed set
    outf = gzip.open(outfile, "w")
    for entry in GTF.iterator(IOTools.openFile(infile_abinitio)):
        if entry.transcript_id in remove_transcripts: continue
        outf.write("%s\n" % str(entry))
    outf.close()

#-------------------------------------------------------------------------------
def buildFilteredLncRNAGeneSet(flagged_gtf, refnoncoding_gtf, outfile, geneset_previous = None):
    '''
    creates a final lincRNA geneset. This geneset will not include any
    single exon lincRNA unless they have been seen previously i.e. it overlaps
    a previously identified lincRNA
    '''
    previous = IndexedGenome.IndexedGenome()
    previous_all = IndexedGenome.IndexedGenome()
    if geneset_previous:
        for transcript in GTF.transcript_iterator( GTF.iterator(IOTools.openFile(geneset_previous)) ):
            # use an indexed genome to assess novelty of lncRNA
            for gtf in transcript:
                previous_all.add(gtf.contig, gtf.start, gtf.end, gtf.strand)
            # only use if single exon
            if len(transcript) == 1:
                previous.add(transcript[0].contig, transcript[0].start, transcript[0].end, [transcript[0].strand, transcript[0].transcript_id])


    # same for refnoncoding geneset        
    refnoncoding = IndexedGenome.IndexedGenome()
    refnoncoding_all = IndexedGenome.IndexedGenome()
    for transcript in GTF.transcript_iterator( GTF.iterator(IOTools.openFile(refnoncoding_gtf)) ):
        # use an indexed genome to assess novelty of lncRNA
        for gtf in transcript:
            refnoncoding_all.add(gtf.contig, gtf.start, gtf.end, gtf.strand)
        # only use if single exon
        if len(transcript) == 1:
            refnoncoding.add(transcript[0].contig, transcript[0].start, transcript[0].end, [transcript[0].strand, transcript[0].transcript_id])

    outf = gzip.open(outfile, "w")
    keep = set()
    for gene in GTF.flat_gene_iterator(GTF.iterator(IOTools.openFile(flagged_gtf))):
        gene_id = gene[0].gene_id
        if gene[0].exon_status == "m":
            # check if there is overlap in the known sets
            for gtf in gene:
                if previous:
                    if refnoncoding_all.contains(gtf.contig, gtf.start, gtf.end):
                        for gtf2 in refnoncoding_all.get(gtf.contig, gtf.start, gtf.end):
                            if gtf.strand == gtf[2]:
                                # add novelty flag
                                gtf.setAttribute("gene_status", "known")
                                outf.write("%s\n" % str(gtf))
                            else:
                                gtf.setAttribute("gene_status", "novel")
                                outf.write("%s\n" % str(gtf))
                    
                    elif previous_all.contains(gtf.contig, gtf.start, gtf.end):
                        for gtf2 in previous_all.get(gtf.contig, gtf.start, gtf.end):
                            if gtf.strand == gtf[2]:
                                # add novelty flag
                                gtf.setAttribute("gene_status", "known")
                                outf.write("%s\n" % str(gtf))
                    else:
                        gtf.setAttribute("gene_status", "novel")
                        outf.write("%s\n" % str(gtf))
                else:
                    if refnoncoding_all.contains(gtf.contig, gtf.start, gtf.end):
                        for gtf2 in refnoncoding_all.get(gtf.contig, gtf.start, gtf.end):
                            if gtf.strand == gtf[2]:
                                # add novelty flag
                                gtf.setAttribute("gene_status", "known")
                                outf.write("%s\n" % str(gtf))
                            else:
                                gtf.setAttribute("gene_status", "novel")
                                outf.write("%s\n" % str(gtf))
                    else:
                        gtf.setAttribute("gene_status", "novel")
                        outf.write("%s\n" % str(gtf))
                    
        # check for single exon status
        elif entry.exon_status == "s":
            if refnoncoding.contains(entry.contig, entry.start, entry.end):                               
                
                # retain strandedness
                for gtf in refnoncoding.get(entry.contig, entry.start, entry.end):                               
                    if entry.strand == gtf[2][0]:
                        keep.add(gtf[2][1])
                        
            # if geneset is not included then this won't be used
            elif previous.contains(entry.contig, entry.start, entry.end):                                     
                for gtf in previous.get(entry.contig, entry.start, entry.end):                                
                    if entry.strand == gtf[2][0]:
                        keep.add(gtf[2][1])

    # write out ones to keep
    if geneset_previous:
        for gtf in GTF.iterator(IOTools.openFile(geneset_previous)):
            if gtf.transcript_id in keep:
                outf.write("%s\n" % gtf)
        for gtf in GTF.iterator(IOTools.openFile(refnoncoding_gtf)):
            if gtf.transcript_id in keep:
                outf.write("%s\n" % gtf)
    else:
        for gtf in GTF.iterator(IOTools.openFile(refnoncoding_gtf)):
            if gtf.transcript_id in keep:
                outf.write("%s\n" % gtf)

    outf.close()
#-------------------------------------------------------------------------------
def buildFinalLncRNAGeneSet(filteredLncRNAGeneSet, cpc_table, outfile):
    '''
    filters lncRNA set based on the coding potential as output from 
    the CPC
    '''
    # get the transcripts that are designated as coding
    coding_set = set()
    dbh = sqlite3.connect("csvdb")
    cc = dbh.cursor()
    for transcript_id in cc.execute("SELECT transcript_id from %s WHERE C_NC='coding'" % cpc_table):
        coding_set.add(transcript_id[0])

    remove = set()
    for gtf in GTF.iterator(IOTools.openFile(filteredLncRNAGeneSet)):
        if gtf.transcript_id in coding_set:
            remove.add(gtf.gene_id)

    outf = gzip.open(outfile, "w")
    for gtf in GTF.iterator(IOTools.openFile(filteredLncRNAGeneSet)):
        if gtf.gene_id in remove: continue
        outf.write("%s\n" % gtf)
    outf.close()

########################################################
# counter classes
########################################################
class CounterExons:
    '''
    various classes for counting features in a gtf file
    '''

    def __init__(self, gtffile):
        
        self.gtffile = GTF.iterator(IOTools.openFile(gtffile))

    def count(self):
        c = 0
        for gtf in self.gtffile:
            c += 1
        return c
            
#-------------------------------------------------------------------------------
class CounterTranscripts(CounterExons):

    def count(self):
        c = 0
        for transcript in GTF.transcript_iterator(self.gtffile):
            c += 1

        return c

#-------------------------------------------------------------------------------
class CounterGenes(CounterExons):

    def count(self):
        c = 0
        for gtf in GTF.flat_gene_iterator(self.gtffile):
            c += 1
        return c

#-------------------------------------------------------------------------------
class CounterExonsPerTranscript(CounterExons):
    '''
    returns the average number of exons per transcript
    '''

    def count(self):
        
        no_exons = []
        for transcript in GTF.transcript_iterator(self.gtffile):
            no_exons.append(len(transcript))
        return float(sum(no_exons))/len(no_exons)

#-------------------------------------------------------------------------------
class CounterExonsPerGene(CounterExons):
    '''
    returns the average number of exons per transcript
    '''

    def count(self):
        
        no_exons = []
        for gene in GTF.flat_gene_iterator(self.gtffile):
            no_exons.append(len(gene))
        return float(sum(no_exons))/len(no_exons)

#-------------------------------------------------------------------------------
class CounterSingleExonTranscripts(CounterExons):
    
    def count(self):
        c = 0
        for transcript in GTF.transcript_iterator(self.gtffile):
            if len(transcript) == 1:
                c += 1
        return c

#-------------------------------------------------------------------------------
class CounterMultiExonTranscripts(CounterExons):
    
    def count(self):
        c = 0
        for transcript in GTF.transcript_iterator(self.gtffile):
            if len(transcript) > 1:
                c += 1
        return c

#-------------------------------------------------------------------------------
class CounterSingleExonGenes(CounterExons):
    
    def count(self):
        
        gene_ids = set()
        for transcript in GTF.transcript_iterator(self.gtffile):
            if len(transcript) == 1:
                gene_ids.add(transcript[0].gene_id)
        return len(gene_ids)

#-------------------------------------------------------------------------------
class CounterMultiExonGenes(CounterExons):
    
    def count(self):
        
        # note that this approach work swhen there are also single
        # exon gtf entries from the same gene (doesn't work to use flat_gene_iterator)
        gene_ids = set()
        for transcript in GTF.transcript_iterator(self.gtffile):
            if len(transcript) > 1:
                gene_ids.add(transcript[0].gene_id)
        return len(gene_ids)
#----------------------------------------------------------------------
def flagExonStatus(gtffile, outfile):
    '''
    adds an attribute to a gtf based on the multi exonic status of a
    transcript - "s" = single, "m" = multi
    '''
    infile = GTF.iterator(IOTools.openFile(gtffile))
    outf = gzip.open(outfile, "w")

    for transcript in GTF.transcript_iterator(infile):
        if len(transcript) > 1:
            for gtf in transcript:
                gtf.setAttribute("exon_status", "m")
                outf.write("%s\n" % gtf)
        elif len(transcript) == 1:
            for gtf in transcript:
                gtf.setAttribute("exon_status", "s")
                outf.write("%s\n" % gtf)

#############################################################
# classifying lincRNA relative to protein coding transcripts
#############################################################

def classifyLncRNA(lincRNA_gtf, reference, outfile, dist = 2):
    '''
    function to classify lincRNA in terms of their relationship to 
    protein coding genes - creates indices for intervals on the 
    fly - mayb should be creating additional annotations:

    antisense - transcript overlapping protein coding exons on opposite strand
    antisense_upstream - transcript < 2kb from tss on opposite strand
    antisense_downstream - transcript < 2kb from gene end on opposite strand
    sense_upstream - transcript < 2kb from tss on same strand
    sense_downstream - transcript < 2kb from gene end on same strand
    antisense_span - start and end of transcript span and extend beyond protein coding gene on opposite strand
    sense_span - start and end of transcript span and extend beyond protein coding gene on same strand
    intergenic - >2kb from anyt protein coding gene
    intronic - overlaps protein coding gene intron on same strand
    antisense_intronic - overlaps protein coding intron on opposite strand
    '''
    
    # index the reference geneset
    ref = {}
    ref["ref"] = GTF.readAndIndex(GTF.iterator(IOTools.openFile(reference)), with_value = True)

    # create index for intronic intervals
    intron = IndexedGenome.IndexedGenome()
    
    # create index for up and downstream intervals
    plus_up = IndexedGenome.IndexedGenome()
    plus_down = IndexedGenome.IndexedGenome()
    minus_up =  IndexedGenome.IndexedGenome()
    minus_down = IndexedGenome.IndexedGenome()


    # iterate over reference transcripts
    for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(reference))):
        start = transcript[0].end
        for i in range(1,len(transcript)):
            intron.add( transcript[i].contig, start, transcript[i].start, transcript[i].strand)

        # create up and downstream intervals on plus strand
        if transcript[0].strand == "+":
            plus_up.add(transcript[0].contig, transcript[0].start - (dist*1000), transcript[0].start, transcript[0].strand)
            plus_down.add(transcript[0].contig, transcript[len(transcript) - 1].end,  transcript[len(transcript) - 1].end + (dist*1000), transcript[0].strand)

        # create up and downstream intervals on minus strand
        elif transcript[0].strand == "-":
            minus_up.add(transcript[0].contig, transcript[len(transcript) - 1].end,  transcript[len(transcript) - 1].end + (dist*1000), transcript[0].strand)
            minus_down.add(transcript[0].contig, transcript[0].start - (dist*1000), transcript[0].start, transcript[0].strand)

        else:
            print "no strand specified for %s" % transcript[0].transcript_id

    # iterate over lincRNA transcripts
    transcript_class = {}
    for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(lincRNA_gtf))):
        transcript_id = transcript[0].transcript_id
        for gtf in transcript:
            
            # antisense to protein coding gene
            if ref["ref"].contains(gtf.contig, gtf.start, gtf.end):
                for gtf2 in ref["ref"].get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand != gtf2[2].strand:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "antisense"
    
            # upstream sense and antisense - if up or downstream of a protein 
            # coding gene then that classification is prioritised over intronic
            # classification. upstream is prioritised over downstream
            # sense upstream is prioritised over antisense upstream
            elif plus_up.contains(gtf.contig, gtf.start, gtf.end):
                for gtf2 in plus_up.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand == gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "sense_upstream"
                    elif gtf.strand != gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "antisense_upstream"
            
            elif minus_up.contains(gtf.contig, gtf.start, gtf.end):
                for gtf2 in minus_up.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand == gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "sense_upstream"
                    elif gtf.strand != gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "antisense_upstream"
            
            # downstream sense and antisense - downstream antisense is prioritised as
            # less likely to be part of the coding transcript
            elif plus_down.contains(gtf.contig, gtf.start, gtf.end):
                for gtf2 in plus_down.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand != gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "antisense_downstream"
                    elif gtf.strand == gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "sense_downstream"

            elif minus_down.contains(gtf.contig, gtf.start, gtf.end):
                for gtf2 in minus_down.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand != gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "antisense_downstream"
                            
                    elif gtf.strand == gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "sense_downstream"

            # intronic sense and antisense - intronic antisense is prioritised
            elif intron.contains(gtf.contig, gtf.start, gtf.end):
                for gtf2 in intron.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand != gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "antisense_intronic"
                    elif gtf.strand == gtf2[2]:
                        if transcript_id not in transcript_class:
                            transcript_class[transcript_id] = "sense_intronic"
            
            elif transcript_id not in transcript_class:
                transcript_class[transcript_id] = "intergenic"

    outf = gzip.open(outfile, "w")
    outf_unclassified = gzip.open("lncrna_unclassified.gtf.gz", "w")
    for gtf in GTF.iterator(IOTools.openFile(lincRNA_gtf)):
        if gtf.transcript_id in transcript_class:
            gtf.source = transcript_class[gtf.transcript_id]
            outf.write("%s\n" % gtf)
        else:
             outf_unclassified.write("%s\n" % gtf)
    outf.close()
        

    








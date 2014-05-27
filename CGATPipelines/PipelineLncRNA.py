'''
PipelineLncRNA.py - functions and classes for use with the lincRNA pipeline
===========================================================================

'''

import re
import sys
import os
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import gzip
import collections
import CGAT.IndexedGenome as IndexedGenome
import CGAT.IndexedFasta as IndexedFasta
import CGAT.Pipeline as P
import CGAT.Experiment as E
import sqlite3
import CGAT.Experiment as E

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
    coding["protein_coding"] = GTF.readAndIndex(GTF.iterator_filtered(
        GTF.iterator(IOTools.openFile(reference)), source="protein_coding"), with_value=False)

    for gtf in GTF.iterator(inf):
        if coding["protein_coding"].contains(gtf.contig, gtf.start, gtf.end):
            if gtf.class_code == "=":
                outf.write("%s\n" % str(gtf))
    outf.close()

#-------------------------------------------------------------------------


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

#-------------------------------------------------------------------------


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

#-------------------------------------------------------------------------


def buildLncRNAGeneSet(abinitio_lincrna, reference, refnoncoding, pseudogenes_gtf, numts_gtf, outfile, min_length):
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

    E.info("indexing geneset for filtering")

    # 17/11/2012 Jethro added Ig genes to list. (NB 31 entries in $2 of reference.gtf.gz)
    # 10/12/2012 Nick added CDS. This will then filter anything that overlaps with a CDS
    # NB this is the annotation for RefSeq and so requires the user to have created a
    # reference annotation that includes this annotation
    input_sections = ("protein_coding",
                      "processed_pseudogene",
                      "unprocessed_pseudogene",
                      "nonsense_mediated_decay",
                      "retained_intron",
                      "IG_V_gene",
                      "IG_J_gene",
                      "IG_C_gene",
                      "IG_D_gene",
                      "IG_LV_gene",
                      "TR_V_gene",
                      "CDS")

    # create a dictionary containing a separate index for each input section
    indices = {}
    for section in input_sections:
        indices[section] = GTF.readAndIndex(
            GTF.iterator_filtered(GTF.iterator(IOTools.openFile(reference_gtf)),
                                  source=section),
            with_value=True)
    E.info("built indices for %i features" % len(indices))

    # add psuedogenes and numts to dictionary of indices
    indices["numts"] = GTF.readAndIndex(
        GTF.iterator(IOTools.openFile(numts_gtf)), with_value=True)
    E.info("added index for numts")

    indices["pseudogenes"] = GTF.readAndIndex(
        GTF.iterator(IOTools.openFile(pseudogenes_gtf)), with_value=True)
    E.info("added index for pseudogenes")

    # iterate through assembled transcripts, identify those that intersect
    # with indexed sections
    total_transcripts = set()
    remove_transcripts = collections.defaultdict(set)
    E.info("collecting genes to remove")
    for gtf in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(infile_abinitio))):
        # remove those transcripts too short to be classified as lncRNAs
        l = sum([x.end - x.start for x in gtf])
        if l < min_length:
            remove_transcripts[gtf[0].transcript_id].add("length")
        # remove transcripts where one or more exon intersects one or more
        # sections
        for section in indices.iterkeys():
            for exon in gtf:
                transcript_id = exon.transcript_id
                total_transcripts.add(transcript_id)
                if indices[section].contains(
                        exon.contig, exon.start, exon.end):
                    # check if any of the intersected intervals are on same
                    # stand as query exon
                    for interval in indices[section].get(
                            exon.contig, exon.start, exon.end):
                        # only remove transcripts where an exon is on
                        # the same strand as the interval it
                        # intersects
                        if exon.strand == interval[2].strand:
                            remove_transcripts[transcript_id].add(section)

    E.info("removing %i out of %i transcripts" %
           (len(remove_transcripts), len(total_transcripts)))

    # Jethro - removed the automatic retention of transcripts that
    # intersect known non-coding intervals regardless of
    # strand. Instead, any removed transcript that also intersects a
    # lncRNA on the same strand is now output to a separate gtf.
    noncoding = GTF.readAndIndex(
        GTF.iterator(IOTools.openFile(refnoncoding_gtf)), with_value=True)
    rej_gtf = os.path.join(os.path.dirname(outfile),
                           "lncrna_removed_nc_intersect.gtf.gz")
    rej_gtf = IOTools.openFile(rej_gtf, "w")
    for gtf in GTF.transcript_iterator(GTF.iterator(
            IOTools.openFile(infile_abinitio))):
        if gtf[0].transcript_id in remove_transcripts.keys():
            for exon in gtf:
                if noncoding.contains(exon.contig, exon.start, exon.end):
                    if exon.strand in [x[2].strand for x in
                                       list(noncoding.get(exon.contig,
                                                          exon.start,
                                                          exon.end))]:
                        for exon2 in IOTools.flatten(gtf):
                            rej_gtf.write(str(exon2) + "\n")
                        break
    rej_gtf.close()

    outf = open("lncrna_removed.tsv", "w")
    outf.write("transcript_id" + "\t" + "removed" + "\n")
    for x, y in remove_transcripts.iteritems():
        outf.write("%s\t%s\n" % (x, ",".join(y)))
    outf.close()

    # write out transcripts that are not in removed set
    temp = P.getTempFile(".")
    for entry in GTF.iterator(IOTools.openFile(infile_abinitio)):
        if entry.transcript_id in remove_transcripts:
            continue
        temp.write("%s\n" % str(entry))
    temp.close()

    filename = temp.name
    statement = '''cat %(filename)s | python %(scriptsdir)s/gtf2gtf.py
                   --sort=gene 
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

    os.unlink(temp.name)


#-------------------------------------------------------------------------
def buildFilteredLncRNAGeneSet(flagged_gtf,
                               outfile,
                               genesets_previous,
                               filter_se="transcripts"):
    '''
    creates a filtered lincRNA geneset. This geneset will not include any
    single exon lincRNA unless they have been seen previously i.e. it overlaps
    a previously identified lincRNA

    At this point we add a flag for whether the gene is a novel lncRNA or not
    NB this is a large function and should be modified in the future

    note genesets_previous provided as a list - priority is
    placed on the first in the list
    '''

    # keep single and multi exonic lncRNA separate
    previous_single = IndexedGenome.IndexedGenome()
    previous_multi = IndexedGenome.IndexedGenome()

    E.info("indexing previously identified lncRNA")
    for prev in genesets_previous:
        inf = IOTools.openFile(prev)
        for transcript in GTF.transcript_iterator(GTF.iterator(inf)):
            # use an indexed genome to assess novelty of lncRNA
            if len(transcript) > 1:
                for gtf in transcript:
                    previous_multi.add(
                        gtf.contig, gtf.start, gtf.end, [transcript[0].strand, transcript[0].gene_id])
            # add single exons
            elif len(transcript) == 1:
                previous_single.add(transcript[0].contig, transcript[0].start, transcript[
                                    0].end, [transcript[0].strand, transcript[0].gene_id])

    # create sets for keeping and discarding genes
    temp = P.getTempFile(dir=".")
    keep = set()
    known = set()
    novel = set()

    # iterate over the flagged GTF - flagged for exon status
    E.info("checking for overlap with previously identified sets")
    for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(flagged_gtf))):
        gene_id = transcript[0].gene_id
        # check if there is overlap in the known sets in order to add
        # gene_status attribute
        if transcript[0].exon_status == "m":
            if previous_multi.contains(transcript[0].contig,
                                       min([gtf.start for gtf in transcript]),
                                       max([gtf.end for gtf in transcript])):
                known.add(transcript[0].gene_id)
                # add previous multi exonic lincRNA gene id
                # to known set - to avoid addin gto gtf later
                for gtf2 in previous_multi.get(transcript[0].contig,
                                               min([
                                                   gtf.start for gtf in transcript]),
                                               max([gtf.end for gtf in transcript])):
                    known.add(gtf2[2][1])
            else:
                novel.add(transcript[0].gene_id)
        elif filter_se == "locus" and transcript[0].exon_status_locus == "m":
            if previous_multi.contains(transcript[0].contig,
                                       min([gtf.start for gtf in transcript]),
                                       max([gtf.end for gtf in transcript])):
                known.add(transcript[0].gene_id)
                # add previous multi exonic lincRNA gene id
                # to known set - to avoid addin gto gtf later
                for gtf2 in previous_multi.get(transcript[0].contig,
                                               min([
                                                   gtf.start for gtf in transcript]),
                                               max([gtf.end for gtf in transcript])):
                    known.add(gtf2[2][1])
            else:
                novel.add(transcript[0].gene_id)
        else:
            continue

    E.info("writing filtered lncRNA geneset")
    E.info("writing %i assembled but known" % len(known))
    E.info("writing %i assembled novel" % len(novel))

    # write out ones to keep from assembled data
    for gtf in GTF.iterator(IOTools.openFile(flagged_gtf)):
        if gtf.gene_id in known and gtf.exon_status == "m":
            gtf.setAttribute("gene_status", "known")
            temp.write("%s\n" % gtf)
        elif gtf.gene_id in novel and gtf.exon_status == "m":
            gtf.setAttribute("gene_status", "novel")
            temp.write("%s\n" % gtf)

    # write out ones to keep - from previous evidence
    # i.e. anything single exon or anything
    # multi-exonic non-overlapping the assembled set

    # hierarchically done so add to 'done' if found as we iterate
    # over the previous sets
    done = IndexedGenome.IndexedGenome()
    known_count = 0
    for prev in genesets_previous:
        inf = IOTools.openFile(prev)
        for gene in GTF.flat_gene_iterator(GTF.iterator(inf)):
            gene_id = gene[0].gene_id

            # dont' write out ones that we build
            if gene_id in known:
                continue
            elif done.contains(gene[0].contig, min([gtf.start for gtf in gene]), max([gtf.end for gtf in gene])):
                continue
            else:
                for gtf in gene:
                    gtf.setAttribute("gene_status", "known")
                    temp.write("%s\n" % gtf)
                known_count += 1
                done.add(gene[0].contig, min([gtf.start for gtf in gene]), max(
                    [gtf.end for gtf in gene]), gene_id)
    E.info("written %i from %s " % (known_count, ",".join(genesets_previous)))
    temp.close()

    filename = temp.name
    statement = '''cat %(filename)s | python %(scriptsdir)s/gtf2gtf.py 
                   --sort=transcript 
                   --log=%(outfile)s.log
                   | gzip > %(outfile)s'''
    P.run()

#-------------------------------------------------------------------------


def buildFinalLncRNAGeneSet(filteredLncRNAGeneSet,
                            cpc_table,
                            outfile,
                            filter_cpc=False,
                            cpc_threshold=1,
                            rename_lncRNA=False):
    '''filters lncRNA set based on the coding potential as output from
    the CPC

    '''

    if filter_cpc:
        # set threshold for filtering transcripts on coding potential
        cpc_thresh = float(cpc_threshold)

        # get the transcripts that are designated as coding
        coding_set = set()
        dbh = sqlite3.connect("csvdb")
        cc = dbh.cursor()
        for transcript_id in cc.execute("SELECT transcript_id"
                                        " FROM %s"
                                        " WHERE CP_score > %f"
                                        % (cpc_table, cpc_thresh)):
            coding_set.add(transcript_id[0])

        remove = set()
        outf_coding = gzip.open("gtfs/cpc_removed.gtf.gz", "w")
        for gtf in GTF.iterator(IOTools.openFile(filteredLncRNAGeneSet)):
            if gtf.transcript_id in coding_set:
                remove.add(gtf.gene_id)
                outf_coding.write("%s\n" % gtf)
        outf_coding.close()
    else:
        # create empty set
        remove = set()

    # get temporary file for built lncrna
    temp = P.getTempFile(".")

    # get temporary file for known lncrna
    temp2 = P.getTempFile(".")

    for gtf in GTF.iterator(IOTools.openFile(filteredLncRNAGeneSet)):
        if gtf.gene_id in remove:
            continue
        if gtf.transcript_id.find("TCONS") != -1:
            # output known and built transcripts separately
            temp.write("%s\n" % gtf)
        else:
            temp2.write("%s\n" % gtf)
    temp.close()
    temp2.close()

    filename = temp.name
    filename2 = temp2.name

    if rename_lncRNA:
        filename3 = P.getTempFilename(".")
        statement = ("cat %(filename)s |"
                     " python %(scriptsdir)s/gtf2gtf.py"
                     "  --sort=gene"
                     "  --log=%(outfile)s.log |"
                     " python %(scriptsdir)s/gtf2gtf.py"
                     "  --renumber-genes=NONCO%%i"
                     "  --log=%(outfile)s.log |"
                     " python %(scriptsdir)s/gtf2gtf.py"
                     "  --sort=gene"
                     "  --log=%(outfile)s.log"
                     " > %(filename3)s;"
                     " cat %(filename2)s %(filename3)s |"
                     " python %(scriptsdir)s/gtf2gtf.py"
                     "  --sort=contig+gene"
                     "  --log=%(outfile)s.log |"
                     " gzip > %(outfile)s")
        P.run()
        os.unlink(filename3)
    else:
        statement = ("cat %(filename)s %(filename2)s |"
                     " python %(scriptsdir)s/gtf2gtf.py"
                     "  --sort=contig+gene"
                     "  --log=%(outfile)s.log |"
                     " gzip > %(outfile)s")
        P.run()

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

#-------------------------------------------------------------------------


class CounterTranscripts(CounterExons):

    def count(self):
        c = 0
        for transcript in GTF.transcript_iterator(self.gtffile):
            c += 1

        return c

#-------------------------------------------------------------------------


class CounterGenes(CounterExons):

    def count(self):
        c = 0
        for gtf in GTF.flat_gene_iterator(self.gtffile):
            c += 1
        return c

#-------------------------------------------------------------------------


class CounterExonsPerTranscript(CounterExons):

    '''
    returns the average number of exons per transcript
    '''

    def count(self):

        no_exons = []
        for transcript in GTF.transcript_iterator(self.gtffile):
            no_exons.append(len(transcript))
        return float(sum(no_exons)) / len(no_exons)

#-------------------------------------------------------------------------


class CounterExonsPerGene(CounterExons):

    '''
    returns the average number of exons per transcript
    '''

    def count(self):

        no_exons = []
        for gene in GTF.flat_gene_iterator(self.gtffile):
            no_exons.append(len(gene))
        return float(sum(no_exons)) / len(no_exons)

#-------------------------------------------------------------------------


class CounterSingleExonTranscripts(CounterExons):

    def count(self):
        c = 0
        for transcript in GTF.transcript_iterator(self.gtffile):
            if len(transcript) == 1:
                c += 1
        return c

#-------------------------------------------------------------------------


class CounterMultiExonTranscripts(CounterExons):

    def count(self):
        c = 0
        for transcript in GTF.transcript_iterator(self.gtffile):
            if len(transcript) > 1:
                c += 1
        return c

#-------------------------------------------------------------------------


class CounterSingleExonGenes(CounterExons):

    def count(self):

        gene_ids = set()
        for gene in GTF.flat_gene_iterator(self.gtffile):
            if gene[0].exon_status_locus == "s":
                gene_ids.add(gene[0].gene_id)
        return len(gene_ids)

#-------------------------------------------------------------------------


class CounterMultiExonGenes(CounterExons):

    def count(self):

        # note that this approach works when there are also single
        # exon gtf entries from the same gene (doesn't work to use
        # flat_gene_iterator)
        gene_ids = set()
        for gene in GTF.flat_gene_iterator(self.gtffile):
            if gene[0].exon_status_locus == "m":
                gene_ids.add(gene[0].gene_id)
        return len(gene_ids)

#----------------------------------------------------------------------


def flagExonStatus(gtf_file, outfile):
    '''
    Adds two attributes to a gtf, the first species transcript exon status,
    the second specifies gene exon status. 
    It is possible for genes to contain single-exon transcripts but not the 
    reciprocal.
    '''

    tmpf1 = P.getTempFilename(".")
    tmpf2 = P.getTempFilename(".")
    outf = IOTools.openFile(tmpf2, "w")

    # sort infile
    statement = ("zcat %(gtf_file)s |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --sort=gene"
                 "  --log=%(outfile)s.log"
                 " > %(tmpf1)s")
    P.run()

    # create dictionary where key is transcript_id,
    # value is a list of gtf_proxy objects
    for gene in GTF.gene_iterator(GTF.iterator(IOTools.openFile(tmpf1))):
        trans_dict = collections.defaultdict(list)
        # load current gene into dictionary
        for transcript in gene:
            for exon in transcript:
                trans_dict[exon.transcript_id].append(exon)

        # set exon status for transcripts
        for transcript in trans_dict.iterkeys():
            if len(trans_dict[transcript]) == 1:
                exon_status = "s"
            else:
                exon_status = "m"

            for exon in trans_dict[transcript]:
                exon.setAttribute("exon_status", exon_status)

        # collate transcript exon status for a gene
        transcript_status = set()
        for exons in trans_dict.itervalues():
            transcript_status.update([exon.exon_status for exon in exons])
        # set gene_exon_status
        if "m" in transcript_status:
            gene_exon_status = "m"
        else:
            gene_exon_status = "s"

        # write gene model to outfile
        for transcript in trans_dict.itervalues():
            for exon in transcript:
                exon.setAttribute("exon_status_locus", gene_exon_status)
                outf.write(str(exon) + "\n")

    outf.close()

    statement = ("cat %(tmpf2)s |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --sort=transcript"
                 "  --log=%(outfile)s.log |"
                 " gzip > %(outfile)s")
    P.run()

    os.unlink(tmpf1)
    os.unlink(tmpf2)


#############################################################
# classifying lincRNA TRANSCRIPT level relative to protein
# coding transcripts
#############################################################


def classifyLncRNA(lincRNA_gtf, reference, outfile, dist=2):
    '''Classify lincRNA in terms of their proximity to protein coding
    genes - creates indices for intervals on the fly - maybe should be
    creating additional annotations:

    antisense - GENE overlapping protein coding exons or introns on opposite strand
    antisense_upstream - GENE < Xkb from tss on opposite strand
    antisense_downstream - GENE < Xkb from gene end on opposite strand
    sense_upstream - GENE < Xkb from tss on same strand
    sense_downstream - GENE < Xkb from gene end on same strand
    intergenic - >Xkb from any protein coding gene
    intronic - overlaps protein coding gene intron on same strand
    antisense_intronic - overlaps protein coding intron on opposite strand

    '''

    # index the reference geneset
    ref = {}
    ref["ref"] = GTF.readAndIndex(
        GTF.iterator(IOTools.openFile(reference)), with_value=True)

    # create index for intronic intervals
    intron = IndexedGenome.IndexedGenome()

    # create index for up and downstream intervals
    plus_up = IndexedGenome.IndexedGenome()
    plus_down = IndexedGenome.IndexedGenome()
    minus_up = IndexedGenome.IndexedGenome()
    minus_down = IndexedGenome.IndexedGenome()

    # iterate over reference transcripts and create intervals in memory
    outf = open("introns.bed", "w")
    for transcript in GTF.transcript_iterator(GTF.iterator(
            IOTools.openFile(reference))):
        start = transcript[0].end
        for i in range(1, len(transcript)):
            intron.add(
                transcript[i].contig, start,
                transcript[i].start, transcript[i].strand)
            start = transcript[i].end

        # create up and downstream intervals on plus strand
        if transcript[0].strand == "+":
            plus_up.add(transcript[0].contig, transcript[
                        0].start - (dist * 1000), transcript[0].start, transcript[0].strand)
            plus_down.add(transcript[0].contig, transcript[
                          len(transcript) - 1].end,  transcript[len(transcript) - 1].end + (dist * 1000), transcript[0].strand)

        # create up and downstream intervals on minus strand
        elif transcript[0].strand == "-":
            minus_up.add(transcript[0].contig, transcript[
                         len(transcript) - 1].end,  transcript[len(transcript) - 1].end + (dist * 1000), transcript[0].strand)
            minus_down.add(transcript[0].contig, transcript[
                           0].start - (dist * 1000), transcript[0].start, transcript[0].strand)
        else:
            print "WARNING: no strand specified for %s" % transcript[0].transcript_id

    # iterate over lincRNA transcripts
    transcript_class = {}
    for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(lincRNA_gtf))):
        transcript_id = transcript[0].transcript_id
        for gtf in transcript:

            # antisense to protein coding gene
            if ref["ref"].contains(gtf.contig, gtf.start, gtf.end):
                for gtf2 in ref["ref"].get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand != gtf2[2].strand:
                        transcript_class[transcript_id] = "antisense"

            # upstream sense and antisense - if up or downstream of a protein
            # coding gene then that classification is prioritised over intronic
            # classification. upstream is prioritised over downstream
            # sense upstream is prioritised over antisense upstream
            elif plus_up.contains(transcript[0].contig, transcript[0].start, transcript[len(transcript) - 1].end) and not intron.contains(transcript[0].contig, transcript[0].start, transcript[len(transcript) - 1].end):
                for gtf2 in plus_up.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand == gtf2[2]:
                        transcript_class[transcript_id] = "sense_upstream"
                    elif gtf.strand != gtf2[2]:
                        transcript_class[transcript_id] = "antisense_upstream"

            # and not intron.contains(transcript[0].contig,
            # transcript[0].start, transcript[len(transcript) -1].end):
            elif minus_up.contains(transcript[0].contig, transcript[0].start, transcript[len(transcript) - 1].end):
                for gtf2 in minus_up.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand == gtf2[2]:
                        transcript_class[transcript_id] = "sense_upstream"
                    elif gtf.strand != gtf2[2]:
                        transcript_class[transcript_id] = "antisense_upstream"

            # downstream sense and antisense - downstream antisense is prioritised as
            # less likely to be part of the coding transcript
            elif plus_down.contains(transcript[0].contig, transcript[0].start, transcript[len(transcript) - 1].end) and not intron.contains(transcript[0].contig, transcript[0].start, transcript[len(transcript) - 1].end):
                for gtf2 in plus_down.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand != gtf2[2]:
                        transcript_class[
                            transcript_id] = "antisense_downstream"
                    elif gtf.strand == gtf2[2]:
                        transcript_class[transcript_id] = "sense_downstream"

            elif minus_down.contains(transcript[0].contig, transcript[0].start, transcript[len(transcript) - 1].end) and not intron.contains(transcript[0].contig, transcript[0].start, transcript[len(transcript) - 1].end):
                for gtf2 in minus_down.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand != gtf2[2]:
                        transcript_class[
                            transcript_id] = "antisense_downstream"
                    elif gtf.strand == gtf2[2]:
                        transcript_class[transcript_id] = "sense_downstream"

            # intronic sense and antisense - intronic antisense is prioritised
            elif intron.contains(transcript[0].contig, transcript[0].start, transcript[len(transcript) - 1].end):
                for gtf2 in intron.get(gtf.contig, gtf.start, gtf.end):
                    if gtf.strand != gtf2[2]:
                        transcript_class[transcript_id] = "antisense_intronic"

                    elif gtf.strand == gtf2[2]:
                        transcript_class[transcript_id] = "sense_intronic"

            # we have known lncRNA that overlap protein coding exons so
            # classify these as well
            elif ref["ref"].contains(gtf.contig, gtf.start, gtf.end):
                if gtf.gene_status == "known":
                    transcript_class[transcript_id] = "sense_protein_coding"

        # catch intergenic transcripts - based on the
        # assumption that anything not falling into the above classes
        # is intergenic
        if transcript_id not in transcript_class:
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

#############################################################
# classifying lincRNA GENE level relative to protein coding
# transcripts
#############################################################


def classifyLncRNAGenes(lincRNA_gtf, reference, outfile, dist=2):

      # index the reference geneset
    ref = {}
    ref["ref"] = GTF.readAndIndex(
        GTF.iterator(IOTools.openFile(reference)), with_value=True)

    # create index for intronic intervals
    intron = IndexedGenome.IndexedGenome()

    # create index for up and downstream intervals
    plus_up = IndexedGenome.IndexedGenome()
    plus_down = IndexedGenome.IndexedGenome()
    minus_up = IndexedGenome.IndexedGenome()
    minus_down = IndexedGenome.IndexedGenome()

    # iterate over reference transcripts and create intervals in memory
    for transcript in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(reference))):
        start = transcript[0].end
        for i in range(1, len(transcript)):
            intron.add(
                transcript[i].contig, start, transcript[i].start, transcript[i].strand)
            start = transcript[i].end

        # create up and downstream intervals on plus strand
        if transcript[0].strand == "+":
            plus_up.add(transcript[0].contig, transcript[
                        0].start - (dist * 1000), transcript[0].start, transcript[0].strand)
            plus_down.add(transcript[0].contig, transcript[
                          len(transcript) - 1].end,  transcript[len(transcript) - 1].end + (dist * 1000), transcript[0].strand)

        # create up and downstream intervals on minus strand
        elif transcript[0].strand == "-":
            minus_up.add(transcript[0].contig, transcript[
                         len(transcript) - 1].end,  transcript[len(transcript) - 1].end + (dist * 1000), transcript[0].strand)
            minus_down.add(transcript[0].contig, transcript[
                           0].start - (dist * 1000), transcript[0].start, transcript[0].strand)
        else:
            print "WARNING: no strand specified for %s" % transcript[0].transcript_id

    # iterate over lincRNA genes
    outf_introns = os.path.join(os.path.dirname(outfile),
                                "sense_intronic_removed.gtf.gz")
    outf_introns = gzip.open(outf_introns, "w")
    gene_class = {}

    for gtf in GTF.merged_gene_iterator(GTF.iterator(
            IOTools.openFile(lincRNA_gtf))):
        gene_id = gtf.gene_id

        # the first classification resolves any gene
        # that overlaps a gene. We don't mind whether it
        # overlaps protein coding gene exons or introns
        if ref["ref"].contains(gtf.contig, gtf.start, gtf.end):
            for gtf2 in ref["ref"].get(gtf.contig, gtf.start, gtf.end):
                if gtf.strand != gtf2[2]:
                    gene_class[gene_id] = "antisense"
                else:
                    gene_class[gene_id] = "sense"

        # remove intronic sense transcripts at this point
        elif intron.contains(gtf.contig, gtf.start, gtf.end):
            for gtf2 in intron.get(gtf.contig, gtf.start, gtf.end):
                if gtf.strand == gtf2[2]:
                    outf_introns.write("%s\n" % gtf)
                else:
                    gene_class[gene_id] = "antisense"

        # the second classification resolves sense and antisense genes up and
        # downstream of protein coding genes - nb having some problems with the
        # merged gene iterator
        elif plus_up.contains(gtf.contig, gtf.start, gtf.end):
            for gtf2 in plus_up.get(gtf.contig, gtf.start, gtf.end):
                if gtf.strand != gtf2[2]:
                    if gene_id in gene_class:
                        continue
                    gene_class[gene_id] = "antisense_upstream"
                else:
                    if gene_id in gene_class:
                        continue
                    gene_class[gene_id] = "sense_upstream"
        elif minus_up.contains(gtf.contig, gtf.start, gtf.end):
            for gtf2 in minus_up.get(gtf.contig, gtf.start, gtf.end):
                if gtf.strand != gtf2[2]:
                    if gene_id in gene_class:
                        continue
                    gene_class[gene_id] = "antisense_upstream"
                else:
                    if gene_id in gene_class:
                        continue
                    gene_class[gene_id] = "sense_upstream"
        elif plus_down.contains(gtf.contig, gtf.start, gtf.end):
            for gtf2 in plus_down.get(gtf.contig, gtf.start, gtf.end):
                if gtf.strand != gtf2[2]:
                    if gene_id in gene_class:
                        continue
                    gene_class[gene_id] = "antisense_downstream"
                else:
                    if gene_id in gene_class:
                        continue
                    gene_class[gene_id] = "sense_downstream"
        elif minus_down.contains(gtf.contig, gtf.start, gtf.end):
            for gtf2 in minus_down.get(gtf.contig, gtf.start, gtf.end):
                if gtf.strand != gtf2[2]:
                    if gene_id in gene_class:
                        continue
                    gene_class[gene_id] = "antisense_downstream"
                else:
                    if gene_id in gene_class:
                        continue
                    gene_class[gene_id] = "sense_downstream"

        # the third classification assumes all genes have been classified leaving
        # intergenic genes
        else:
            gene_class[gene_id] = "intergenic"

    outf = gzip.open(outfile, "w")
    for gtf in GTF.iterator(IOTools.openFile(lincRNA_gtf)):
        if gtf.gene_id in gene_class:
            gtf.source = gene_class[gtf.gene_id]
            outf.write("%s\n" % gtf)
    outf.close()
    outf_introns.close()


##########################################################################
##########################################################################
# An alternative method for classifying LncRNAs relative to supplied geneset
##########################################################################


def write_to_temp(tempfile, interval_list, transcript, check_strand=True):
    if check_strand:
        for interval in interval_list:
            if interval[2][6] == transcript[0].strand:
                tempfile.write(transcript[0].gene_id + "\t"
                               + str(interval[0]) + "\t"
                               + str(interval[1]) + "\t"
                               + str(transcript[0]) + "\t"
                               + "\t".join(interval[2]) + "\n")
    else:
        for interval in interval_list:
            tempfile.write(transcript[0].gene_id + "\t"
                           + str(interval[0]) + "\t"
                           + str(interval[1]) + "\t"
                           + str(transcript[0]) + "\t"
                           + "\t".join(interval[2]) + "\n")


def reClassifyLncRNAGenes(lncRNA_gtf,
                          reference_gtf,
                          outfile,
                          upstr_dist=5,
                          dstr_dist=5,
                          wdir="."):
    """
    This re-write of classifyLncRNAGenes() does not throw out intronic loci, but 
    labels them as either sense-intronic or sense-overlap.
    It also fixes the bug that cause sense gene-models encompassing reference 
    exons to be output as antisense. 
    Because lncRNA boundaries intersect multiple intervals in indexes, rather
    than classifying each lncRNA multiple times, lncRNA strand is instead 
    compared to a list of interval strand values, if any of these are sense, 
    then the lncRNA is classified as sense etc. 

    Sense-intronic: when lncRNA loci start and end are contained within a single
    intron. 
    Sense-overlap: when lncRNA loci start and end are in different introns. Note 
    that different introns do not necessarily come from different gene-models. 
    """

    # index exons in the reference gene-set
    ref_index = IndexedGenome.IndexedGenome()
    for exon in GTF.iterator(IOTools.openFile(reference_gtf)):
        ref_index.add(exon.contig, exon.start, exon.end, str(exon).split())

    # create index for all other intervals to be classified
    intron = IndexedGenome.IndexedGenome()
    plus_up = IndexedGenome.IndexedGenome()
    plus_down = IndexedGenome.IndexedGenome()
    minus_up = IndexedGenome.IndexedGenome()
    minus_down = IndexedGenome.IndexedGenome()

    # iterate over reference transcripts and create intervals in memory
    ref_file = IOTools.openFile(reference_gtf)
    for transcript in GTF.transcript_iterator(GTF.iterator(ref_file)):
        start = transcript[0].end
        for i in range(1, len(transcript)):
            intron.add(transcript[i].contig,
                       start,
                       transcript[i].start,
                       str(transcript[i]).split())
            start = transcript[i].end

        # create up and downstream intervals on plus strand
        if transcript[0].strand == "+":
            plus_up.add(transcript[0].contig,
                        transcript[0].start - (upstr_dist * 1000),
                        transcript[0].start,
                        str(transcript[0]).split())
            plus_down.add(transcript[0].contig,
                          transcript[len(transcript) - 1].end,
                          transcript[
                              len(transcript) - 1].end + (dstr_dist * 1000),
                          str(transcript[len(transcript) - 1]).split())

        # create up and downstream intervals on minus strand
        elif transcript[0].strand == "-":
            minus_up.add(transcript[0].contig,
                         transcript[len(transcript) - 1].end,
                         transcript[len(transcript) - 1].end +
                         (upstr_dist * 1000),
                         str(transcript[len(transcript) - 1]).split())
            minus_down.add(transcript[0].contig,
                           transcript[0].start - (dstr_dist * 1000),
                           transcript[0].start,
                           str(transcript[0]).split())
        else:
            E.warn("WARNING: no strand specified for %s" %
                   transcript[0].transcript_id)

    # create single representative transcript for each lncRNA gene_id
    merged_lncRNA_gtf = P.getTempFilename(wdir)
    to_cluster = False
    statement = ("zcat %(lncRNA_gtf)s |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --sort=gene"
                 "  --log=%(outfile)s.log |"
                 " python %(scriptsdir)s/gtf2gtf.py"
                 "  --merge-exons"
                 "  --log=%(outfile)s.log"
                 " > %(merged_lncRNA_gtf)s")
    P.run()

    # create a temp directory containing the indexed intervals used to classify
    # the lncRNA transcripts created (for debugging purposes)
    # create a temporary count of # of gene_models in each category
    tempdir = P.getTempDir(wdir)
    E.info("intersecting intervals are being written to %s"
           % os.path.abspath(tempdir))
    temp_file_names = ["sense",
                       "sense_intronic",
                       "sense_overlap",
                       "antisense",
                       "sense_downstream",
                       "sense_upstream",
                       "antisense_downstream",
                       "antisense_upstream",
                       "intergenic"]
    temp_files = {}
    temp_count = {}
    for handle in temp_file_names:
        temp_count[handle] = 0
        temp_files[handle] = IOTools.openFile(os.path.join(tempdir,
                                                           handle), "w")

    # iterate through the representative (i.e. merged) lncRNA transcripts
    # each lncRNA transcript is classified only once.
    # In situations where a lncRNA fits > 1 classification, priority is:
    # (sense > antisense)
    # & (overlap_exons > overlap_introns > downstream > upstream > intergenic)
    lnc_file = IOTools.openFile(merged_lncRNA_gtf)
    gene_class = {}  # dictionary of gene_id : classification
    input_transcripts = 0  # keep track of # transcripts in lncRNA_gtf
    for transcript in GTF.transcript_iterator(GTF.iterator(lnc_file)):
        input_transcripts += 1
        gene_id = transcript[0].gene_id
        strand = transcript[0].strand

        # create lists of indexed intervals that intersect transcript exons
        overlap_list = []
        intron_list = []
        plus_down_list = []
        minus_down_list = []
        plus_up_list = []
        minus_up_list = []
        for exon in transcript:
            if exon.contig in ref_index.mIndex.keys():
                overlap_list.extend([x for x in list(ref_index.get(exon.contig,
                                                                   exon.start,
                                                                   exon.end))])
            else:
                E.warn("Contig %s not in reference exon index "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in intron.mIndex.keys():
                intron_list.extend([x for x in list(intron.get(exon.contig,
                                                               exon.start,
                                                               exon.end))])
            else:
                E.warn("Contig %s not in reference intron index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in plus_down.mIndex.keys():
                plus_down_list.extend([x for x in list(plus_down.get(exon.contig,
                                                                     exon.start,
                                                                     exon.end))])
            else:
                E.warn("Contig %s not in plus downstream index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in minus_down.mIndex.keys():
                minus_down_list.extend([x for x in list(minus_down.get(exon.contig,
                                                                       exon.start,
                                                                       exon.end))])
            else:
                E.warn("Contig %s not in minus downstream index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in plus_up.mIndex.keys():
                plus_up_list.extend([x for x in list(plus_up.get(exon.contig,
                                                                 exon.start,
                                                                 exon.end))])
            else:
                E.warn("Contig %s not in plus upstream index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))
            if exon.contig in minus_up.mIndex.keys():
                minus_up_list.extend([x for x in list(minus_up.get(exon.contig,
                                                                   exon.start,
                                                                   exon.end))])
            else:
                E.warn("Contig %s not in minus upstream index, "
                       "failed to retrieve intervals for %s" % (exon.contig,
                                                                exon.gene_id))

        # check if any exon in lncRNA intersects an reference exon
        if overlap_list:
            # if the intersecting exons are on the same strand,
            # classify lncRNA as sense.
            if strand in [x[2][6] for x in overlap_list]:
                gene_class[gene_id] = "sense"
                write_to_temp(temp_files[gene_class[gene_id]],
                              overlap_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # otherwise check if lncRNA has sense overlap with a reference
            # intron
            elif intron_list and strand in [x[2][6] for x in intron_list]:
                last = len(transcript) - 1
                start_list = [(x[0], x[1]) for x in list(intron.get(transcript[
                    0].contig, transcript[0].start, transcript[0].start + 1)) if x[2][6] == strand]
                end_list = [(x[0], x[1]) for x in list(intron.get(transcript[
                    last].contig, transcript[last].end, transcript[last].end + 1)) if x[2][6] == strand]
                # if start and end of transcript are within the same sense
                # introns, then lncRNA is classified as 'sense_intronic'
                if set(start_list) == set(end_list):
                    gene_class[gene_id] = "sense_intronic"
                    write_to_temp(temp_files[gene_class[gene_id]],
                                  intron_list,
                                  transcript)
                    temp_count[gene_class[gene_id]] += 1

                # if start/end are within different sense introns,
                # then lncRNA is classified as 'sense overlap'
                else:
                    gene_class[gene_id] = "sense_overlap"
                    write_to_temp(temp_files[gene_class[gene_id]],
                                  intron_list,
                                  transcript)
                    temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the plus strand...
            elif plus_down_list and strand in [x[2][6] for x in plus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the minus strand...
            elif minus_down_list and strand in [x[2][6] for x in minus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the plus strand...
            elif plus_up_list and strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the minus strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...if none of the above... classify as antisense
            else:
                gene_class[gene_id] = "antisense"
                write_to_temp(temp_files[gene_class[gene_id]],
                              overlap_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # if lncRNA doesn't intersect a reference exon,
        # check if it overlaps a reference intron
        elif intron_list:
            if strand in [x[2][6] for x in intron_list]:
                last = len(transcript) - 1
                start_list = [(x[0], x[1]) for x in list(intron.get(transcript[
                    0].contig, transcript[0].start, transcript[0].start + 1)) if x[2][6] == strand]
                end_list = [(x[0], x[1]) for x in list(intron.get(transcript[
                    last].contig, transcript[last].end, transcript[last].end + 1)) if x[2][6] == strand]
                # if start and end of transcript are within the same sense
                # introns, then lncRNA is classified as 'sense_intronic'
                if set(start_list) == set(end_list):
                    gene_class[gene_id] = "sense_intronic"
                    write_to_temp(temp_files[gene_class[gene_id]],
                                  intron_list,
                                  transcript)
                    temp_count[gene_class[gene_id]] += 1

                # if start/end are within different sense introns,
                # then lncRNA is classified as 'sense overlap'
                else:
                    gene_class[gene_id] = "sense_overlap"
                    write_to_temp(temp_files[gene_class[gene_id]],
                                  intron_list,
                                  transcript)
                    temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the plus strand...
            elif plus_down_list and strand in [x[2][6] for x in plus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the minus strand...
            elif minus_down_list and strand in [x[2][6] for x in minus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the plus strand...
            elif plus_up_list and strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream in the minus strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...if none of the above, lncRNAs intersecting introns on
            # the opposite strand are classified as antisense
            else:
                gene_class[gene_id] = "antisense"
                write_to_temp(temp_files[gene_class[gene_id]],
                              intron_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # if lncRNA doesn't intersect reference introns or exons...
        # check if it's downstream on the plus strand...
        elif plus_down_list:
            # ... check if lncRNA is sense downstream...
            if strand in [x[2][6] for x in plus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense downstream on the minus strand...
            elif minus_down_list and strand in [x[2][6] for x in minus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense usptream on the plus strand...
            elif plus_up_list and strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the pluse strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # if none of the above, lncRNA is classified as
            # antisense_downstream
            else:
                gene_class[gene_id] = "antisense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_down_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # check if lncRNA is downstream on the minus strand...
        elif minus_down_list:
            # check if lncRNA is sense downstream
            if strand in [x[2][6] for x in minus_down_list]:
                gene_class[gene_id] = "sense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the plus strand...
            elif plus_up_list and strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the minus strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # if none of the above, lncRNA is classified as
            # antisense_downstream
            else:
                gene_class[gene_id] = "antisense_downstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_down_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # check if lncRNA is upstream on the plus strand...
        elif plus_up_list:
            # check if lncRNA is sense upstream...
            if strand in [x[2][6] for x in plus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # ...check if lncRNA is sense upstream on the plus strand...
            elif minus_up_list and strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # if none of the above, lncRNA is classified as
            # antisense upstream
            else:
                gene_class[gene_id] = "antisense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              plus_up_list,
                              transcript,
                              check_strand=False)
                temp_count[gene_class[gene_id]] += 1

        # check if lncRNA is upstream on the minus strand...
        elif minus_up_list:
            # check if lncRNA is sense upstream...
            if strand in [x[2][6] for x in minus_up_list]:
                gene_class[gene_id] = "sense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

            # otherwise classify as antisense upstream
            else:
                gene_class[gene_id] = "antisense_upstream"
                write_to_temp(temp_files[gene_class[gene_id]],
                              minus_up_list,
                              transcript)
                temp_count[gene_class[gene_id]] += 1

        # lncRNA that do not fall into any of the above categories
        # are classified as intergenic
        else:
            gene_class[gene_id] = "intergenic"
            temp_files[gene_class[gene_id]].write(str(transcript[0]) + "\n")
            temp_count[gene_class[gene_id]] += 1

    # check that all the numbers add up
    E.info("Number of lncRNA loci falling into each category are as follows:")
    for key, value in temp_count.iteritems():
        print(key + "\t" + str(value))
    total_classified = sum(temp_count.values())
    E.info("Total number of lncRNA loci classified: %i" % total_classified)
    E.info("Total number of lncRNA loci in input gtf: %i" % input_transcripts)

    # sanity check:
    assert total_classified == input_transcripts, "Not all lncRNAs in input gtf were successfully classified"

    # close the tempfiles
    for handle in temp_file_names:
        temp_files[handle].close()

    # write the genes plus their classification to the outfile
    outf = IOTools.openFile(outfile, "w")
    for gtf in GTF.iterator(IOTools.openFile(lncRNA_gtf)):
        if gtf.gene_id in gene_class:
            gtf.source = gene_class[gtf.gene_id]
            outf.write(str(gtf) + "\n")
        else:
            E.info("Warning the gene_id %s is not classified" % gtf.gene_id)
    outf.close()

    os.unlink(merged_lncRNA_gtf)

    return tempdir


##########################################################################
##########################################################################
##########################################################################
# Extract pairwise MAF alignments
##########################################################################
# This section of the pipeline makes use of galaxy's maf_utilties (written by
# Dan Blankenberg - see galaxy-dist/lib/galaxy/tools/util) for indexing maf files
# and for retrieving maf blocks that intersect gene models.
# Because maf_utilities.py is not available outside of galaxy, required functions
# have been copied directly (below). These functions have not been altered in the
# hope that maf_utilities will one day be available as a stand-alone module.

# The following classes/functions have been lifted directly from maf_utilities:
# RegionAlignment()
# GenomicRegionAlignment()
# SplicedAlignment()
# build_maf_index_species_chromosomes()
# build_maf_index()
# component_overlaps_region()
# chop_block_by_region()
# orient_block_by_region()
# iter_blocks_split_by_species()
# reduce_block_by_primary_genome()
# fill_region_alignment()
# get_spliced_region_alignment()
# get_starts_ends_fields_from_gene_bed()
# iter_components_by_src()

import tempfile
import string
from copy import deepcopy
import bx.intervals.io
import bx.align.maf
import bx.intervals
import bx.interval_index_file

GAP_CHARS = ['-']
SRC_SPLIT_CHAR = '.'


def src_split(src):
    fields = src.split(SRC_SPLIT_CHAR, 1)
    spec = fields.pop(0)
    if fields:
        chrom = fields.pop(0)
    else:
        chrom = spec
    return spec, chrom

# an object corresponding to a reference layered alignment


class RegionAlignment(object):

    DNA_COMPLEMENT = string.maketrans("ACGTacgt", "TGCAtgca")
    MAX_SEQUENCE_SIZE = sys.maxint  # Maximum length of sequence allowed

    def __init__(self, size, species=[]):
        assert size <= self.MAX_SEQUENCE_SIZE, "Maximum length allowed for an individual sequence has been exceeded (%i > %i)." % (
            size, self.MAX_SEQUENCE_SIZE)
        self.size = size
        self.sequences = {}
        if not isinstance(species, list):
            species = [species]
        for spec in species:
            self.add_species(spec)

    # add a species to the alignment
    def add_species(self, species):
        # make temporary sequence files
        self.sequences[species] = tempfile.TemporaryFile()
        self.sequences[species].write("-" * self.size)

    # returns the names for species found in alignment, skipping names as
    # requested
    def get_species_names(self, skip=[]):
        if not isinstance(skip, list):
            skip = [skip]
        names = self.sequences.keys()
        for name in skip:
            try:
                names.remove(name)
            except:
                pass
        return names

    # returns the sequence for a species
    def get_sequence(self, species):
        self.sequences[species].seek(0)
        return self.sequences[species].read()

    # returns the reverse complement of the sequence for a species
    def get_sequence_reverse_complement(self, species):
        complement = [base for base in self.get_sequence(
            species).translate(self.DNA_COMPLEMENT)]
        complement.reverse()
        return "".join(complement)

    # sets a position for a species
    def set_position(self, index, species, base):
        if len(base) != 1:
            raise Exception("A genomic position can only have a length of 1.")
        return self.set_range(index, species, base)
    # sets a range for a species

    def set_range(self, index, species, bases):
        if index >= self.size or index < 0:
            raise Exception(
                "Your index (%i) is out of range (0 - %i)." % (index, self.size - 1))
        if len(bases) == 0:
            raise Exception(
                "A set of genomic positions can only have a positive length.")
        if species not in self.sequences.keys():
            self.add_species(species)
        self.sequences[species].seek(index)
        self.sequences[species].write(bases)

    # Flush temp file of specified species, or all species
    def flush(self, species=None):
        if species is None:
            species = self.sequences.keys()
        elif not isinstance(species, list):
            species = [species]
        for spec in species:
            self.sequences[spec].flush()


class GenomicRegionAlignment(RegionAlignment):

    def __init__(self, start, end, species=[]):
        RegionAlignment.__init__(self, end - start, species)
        self.start = start
        self.end = end


class SplicedAlignment(object):

    DNA_COMPLEMENT = string.maketrans("ACGTacgt", "TGCAtgca")

    def __init__(self, exon_starts, exon_ends, species=[]):
        if not isinstance(exon_starts, list):
            exon_starts = [exon_starts]
        if not isinstance(exon_ends, list):
            exon_ends = [exon_ends]
        assert len(exon_starts) == len(
            exon_ends), "The number of starts does not match the number of sizes."
        self.exons = []
        for i in range(len(exon_starts)):
            self.exons.append(
                GenomicRegionAlignment(exon_starts[i], exon_ends[i], species))

    # returns the names for species found in alignment, skipping names as
    # requested
    def get_species_names(self, skip=[]):
        if not isinstance(skip, list):
            skip = [skip]
        names = []
        for exon in self.exons:
            for name in exon.get_species_names(skip=skip):
                if name not in names:
                    names.append(name)
        return names

    # returns the sequence for a species
    def get_sequence(self, species):
        sequence = tempfile.TemporaryFile()
        for exon in self.exons:
            if species in exon.get_species_names():
                sequence.write(exon.get_sequence(species))
            else:
                sequence.write("-" * exon.size)
        sequence.seek(0)
        return sequence.read()

    # returns the reverse complement of the sequence for a species
    def get_sequence_reverse_complement(self, species):
        complement = [base for base in self.get_sequence(
            species).translate(self.DNA_COMPLEMENT)]
        complement.reverse()
        return "".join(complement)

    # Start and end of coding region
    @property
    def start(self):
        return self.exons[0].start

    @property
    def end(self):
        return self.exons[-1].end


def build_maf_index_species_chromosomes(filename, index_species=None):
    species = []
    species_chromosomes = {}
    indexes = bx.interval_index_file.Indexes()
    blocks = 0
    try:
        maf_reader = bx.align.maf.Reader(open(filename))
        while True:
            pos = maf_reader.file.tell()
            block = maf_reader.next()
            if block is None:
                break
            blocks += 1
            for c in block.components:
                spec = c.src
                chrom = None
                if "." in spec:
                    spec, chrom = spec.split(".", 1)
                if spec not in species:
                    species.append(spec)
                    species_chromosomes[spec] = []
                if chrom and chrom not in species_chromosomes[spec]:
                    species_chromosomes[spec].append(chrom)
                if index_species is None or spec in index_species:
                    forward_strand_start = c.forward_strand_start
                    forward_strand_end = c.forward_strand_end
                    try:
                        forward_strand_start = int(forward_strand_start)
                        forward_strand_end = int(forward_strand_end)
                    except ValueError:
                        # start and end are not integers, can't add component
                        # to index, goto next component
                        continue
                        # this likely only occurs when parse_e_rows is True?
                        # could a species exist as only e rows? should the
                    if forward_strand_end > forward_strand_start:
                        # require positive length; i.e. certain lines have
                        # start = end = 0 and cannot be indexed
                        indexes.add(
                            c.src, forward_strand_start, forward_strand_end, pos, max=c.src_size)
    except Exception, e:
        # most likely a bad MAF
        log.debug('Building MAF index on %s failed: %s' % (filename, e))
        return (None, [], {}, 0)
    return (indexes, species, species_chromosomes, blocks)

# builds and returns ( index, index_filename ) for specified maf_file


def build_maf_index(maf_file, species=None):
    indexes, found_species, species_chromosomes, blocks = build_maf_index_species_chromosomes(
        maf_file, species)
    if indexes is not None:
        fd, index_filename = tempfile.mkstemp()
        out = os.fdopen(fd, 'w')
        indexes.write(out)
        out.close()
        return (bx.align.maf.Indexed(maf_file, index_filename=index_filename, keep_open=True, parse_e_rows=False), index_filename)
    return (None, None)


def component_overlaps_region(c, region):
    if c is None:
        return False
    start, end = c.get_forward_strand_start(), c.get_forward_strand_end()
    if region.start >= end or region.end <= start:
        return False
    return True


def chop_block_by_region(block, src, region, species=None, mincols=0):
    # This chopping method was designed to maintain consistency with how start/end padding gaps have been working in Galaxy thus far:
    #   behavior as seen when forcing blocks to be '+' relative to src sequence (ref) and using block.slice_by_component( ref, slice_start, slice_end )
    #   whether-or-not this is the 'correct' behavior is questionable, but this will at least maintain consistency
    # comments welcome
    slice_start = block.text_size  # max for the min()
    slice_end = 0  # min for the max()
    old_score = block.score  # save old score for later use
    # We no longer assume only one occurance of src per block, so we need to
    # check them all
    for c in iter_components_by_src(block, src):
        if component_overlaps_region(c, region):
            if c.text is not None:
                rev_strand = False
                if c.strand == "-":
                    # We want our coord_to_col coordinates to be returned from
                    # positive stranded component
                    rev_strand = True
                    c = c.reverse_complement()
                start = max(region.start, c.start)
                end = min(region.end, c.end)
                start = c.coord_to_col(start)
                end = c.coord_to_col(end)
                if rev_strand:
                    # need to orient slice coordinates to the original block
                    # direction
                    slice_len = end - start
                    end = len(c.text) - start
                    start = end - slice_len
                slice_start = min(start, slice_start)
                slice_end = max(end, slice_end)

    if slice_start < slice_end:
        block = block.slice(slice_start, slice_end)
        if block.text_size > mincols:
            # restore old score, may not be accurate, but it is better than 0
            # for everything?
            block.score = old_score
            if species is not None:
                block = block.limit_to_species(species)
                block.remove_all_gap_columns()
            return block
    return None


def orient_block_by_region(block, src, region, force_strand=None):
    # loop through components matching src,
    # make sure each of these components overlap region
    # cache strand for each of overlaping regions
    # if force_strand / region.strand not in strand cache, reverse complement
    # we could have 2 sequences with same src, overlapping region, on
    # different strands, this would cause no reverse_complementing
    strands = [c.strand for c in iter_components_by_src(
        block, src) if component_overlaps_region(c, region)]
    if strands and (force_strand is None and region.strand not in strands) or (force_strand is not None and force_strand not in strands):
        block = block.reverse_complement()
    return block

# split a block into multiple blocks with all combinations of a species
# appearing only once per block


def iter_blocks_split_by_species(block, species=None):
    def __split_components_by_species(components_by_species, new_block):
        if components_by_species:
            # more species with components to add to this block
            components_by_species = deepcopy(components_by_species)
            spec_comps = components_by_species.pop(0)
            for c in spec_comps:
                newer_block = deepcopy(new_block)
                newer_block.add_component(deepcopy(c))
                for value in __split_components_by_species(components_by_species, newer_block):
                    yield value
        else:
            # no more components to add, yield this block
            yield new_block

    # divide components by species
    spec_dict = {}
    if not species:
        species = []
        for c in block.components:
            spec, chrom = src_split(c.src)
            if spec not in spec_dict:
                spec_dict[spec] = []
                species.append(spec)
            spec_dict[spec].append(c)
    else:
        for spec in species:
            spec_dict[spec] = []
            for c in iter_components_by_src_start(block, spec):
                spec_dict[spec].append(c)

    empty_block = bx.align.Alignment(score=block.score, attributes=deepcopy(
        block.attributes))  # should we copy attributes?
    empty_block.text_size = block.text_size
    # call recursive function to split into each combo of spec/blocks
    for value in __split_components_by_species(spec_dict.values(), empty_block):
        # restore original component order
        sort_block_components_by_block(value, block)
        yield value

# reduces a block to only positions exisiting in the src provided


def reduce_block_by_primary_genome(block, species, chromosome, region_start):
    # returns ( startIndex, {species:texts}
    # where texts' contents are reduced to only positions existing in the
    # primary genome
    src = "%s.%s" % (species, chromosome)
    ref = block.get_component_by_src(src)
    start_offset = ref.start - region_start
    species_texts = {}
    for c in block.components:
        species_texts[c.src.split('.')[0]] = list(c.text)
    # remove locations which are gaps in the primary species, starting from
    # the downstream end
    for i in range(len(species_texts[species]) - 1, -1, -1):
        if species_texts[species][i] == '-':
            for text in species_texts.values():
                text.pop(i)
    for spec, text in species_texts.items():
        species_texts[spec] = ''.join(text)
    return (start_offset, species_texts)


def fill_region_alignment(alignment, index, primary_species, chrom, start, end, strand='+', species=None, mincols=0, overwrite_with_gaps=True):
    region = bx.intervals.Interval(start, end)
    region.chrom = chrom
    region.strand = strand
    primary_src = "%s.%s" % (primary_species, chrom)

    # Order blocks overlaping this position by score, lowest first
    blocks = []
    for block, idx, offset in index.get_as_iterator_with_index_and_offset(primary_src, start, end):
        score = float(block.score)
        for i in range(0, len(blocks)):
            if score < blocks[i][0]:
                blocks.insert(i, (score, idx, offset))
                break
        else:
            blocks.append((score, idx, offset))

    #gap_chars_tuple = tuple( GAP_CHARS )
    gap_chars_str = ''.join(GAP_CHARS)
    # Loop through ordered blocks and layer by increasing score
    for block_dict in blocks:
        # need to handle each occurance of sequence in block seperately
        for block in iter_blocks_split_by_species(block_dict[1].get_at_offset(block_dict[2])):
            if component_overlaps_region(block.get_component_by_src(primary_src), region):
                block = chop_block_by_region(
                    block, primary_src, region, species, mincols)  # chop block
                block = orient_block_by_region(
                    block, primary_src, region)  # orient block
                start_offset, species_texts = reduce_block_by_primary_genome(
                    block, primary_species, chrom, start)
                for spec, text in species_texts.items():
                    # we should trim gaps from both sides, since these are not
                    # positions in this species genome (sequence)
                    text = text.rstrip(gap_chars_str)
                    gap_offset = 0
                    # python2.4 doesn't accept a tuple for .startswith()
                    while True in [text.startswith(gap_char) for gap_char in GAP_CHARS]:
                    # while text.startswith( gap_chars_tuple ):
                        gap_offset += 1
                        text = text[1:]
                        if not text:
                            break
                    if text:
                        if overwrite_with_gaps:
                            alignment.set_range(
                                start_offset + gap_offset, spec, text)
                        else:
                            for i, char in enumerate(text):
                                if char not in GAP_CHARS:
                                    alignment.set_position(
                                        start_offset + gap_offset + i, spec, char)
    return alignment


def get_spliced_region_alignment(index,
                                 primary_species,
                                 chrom,
                                 starts,
                                 ends,
                                 strand='+',
                                 species=None,
                                 mincols=0,
                                 overwrite_with_gaps=True):
    """
    Returns a filled spliced region alignment for specified region with start
    and end lists.
    """
    # create spliced alignment object
    if species is not None:
        alignment = SplicedAlignment(starts, ends, species)
    else:
        alignment = SplicedAlignment(starts, ends, [primary_species])
    for exon in alignment.exons:
        fill_region_alignment(exon,
                              index,
                              primary_species,
                              chrom,
                              exon.start,
                              exon.end,
                              strand,
                              species,
                              mincols,
                              overwrite_with_gaps)
    return alignment

# read a GeneBed file, return list of starts, ends, raw fields


def get_starts_ends_fields_from_gene_bed(line):
    # Starts and ends for exons
    starts = []
    ends = []

    fields = line.split()
    # Requires atleast 12 BED columns
    if len(fields) < 12:
        raise Exception("Not a proper 12 column BED line (%s)." % line)
    chrom = fields[0]
    tx_start = int(fields[1])
    tx_end = int(fields[2])
    name = fields[3]
    strand = fields[5]
    if strand != '-':
        strand = '+'  # Default strand is +
    cds_start = int(fields[6])
    cds_end = int(fields[7])

    # Calculate and store starts and ends of coding exons
    region_start, region_end = cds_start, cds_end
    exon_starts = map(int, fields[11].rstrip(',\n').split(','))
    exon_starts = map((lambda x: x + tx_start), exon_starts)
    exon_ends = map(int, fields[10].rstrip(',').split(','))
    exon_ends = map((lambda x, y: x + y), exon_starts, exon_ends)
    for start, end in zip(exon_starts, exon_ends):
        start = max(start, region_start)
        end = min(end, region_end)
        if start < end:
            starts.append(start)
            ends.append(end)
    return (starts, ends, fields)


def iter_components_by_src(block, src):
    for c in block.components:
        if c.src == src:
            yield c


def sort_block_components_by_block(block1, block2):
    # orders the components in block1 by the index of the component in block2
    # block1 must be a subset of block2
    #occurs in-place
    return block1.components.sort(cmp=lambda x, y: block2.components.index(x) - block2.components.index(y))


##########################################################################
# CGAT functions for extracting MAF alignments
##########################################################################

def gtfToBed12(infile, outfile, model):
    """
    Convert a gtf file to bed12 format.
    """
    model = "gene"
    outfile = IOTools.openFile(outfile, "w")

    for all_exons in GTF.transcript_iterator(GTF.iterator(IOTools.openFile(infile, "r"))):
        chrom = all_exons[0].contig
        # GTF.iterator returns start co-ordinates as zero-based
        start = str(all_exons[0].start)
        end = str(all_exons[len(all_exons) - 1].end)
#        if model == "gene":
#            name = all_exons[0].gene_id
#        elif model == "transcript":
#            name = all_exons[0].transcript_id
        name = all_exons[0].gene_id + "__" + all_exons[0].transcript_id
#        else:
#            raise ValueError( "model must either be gene or transcript" )
        score = "0"
        strand = all_exons[0].strand
        thickStart = start
        thickEnd = end
        colourRGB = "0"
        blockCount = str(len(all_exons))

        sizes = []
        starts = []
        for exon in all_exons:
            blockSize = str(exon.end - (exon.start))
            sizes.append(blockSize)
            # start + blockStart should return a zero-based co-ordinate for the
            # exon start
            blockStart = str((exon.start) - int(start))
            starts.append(str(blockStart))
        blockSizes = ','.join(sizes)
        blockStarts = ','.join(starts)

        outfile.write(chrom + "\t"
                      + start + "\t"
                      + end + "\t"
                      + name + "\t"
                      + score + "\t"
                      + strand + "\t"
                      + thickStart + "\t"
                      + thickEnd + "\t"
                      + colourRGB + "\t"
                      + blockCount + "\t"
                      + blockSizes + "\t"
                      + blockStarts + "\n")

    outfile.close()


def filterMAF(infile, outfile, removed, filter_alignments=False):
    """
    Iterates through the MAF file. If filter_alignments == int, then will remove MAF
    blocks for which length < int.
    """
    inf = bx.align.maf.Reader(open(infile))
    outf = bx.align.maf.Writer(open(outfile, "w"))
    outf_rj = bx.align.maf.Writer(open(removed, "w"))

    removed = 0
    included = 0
    if filter_alignments:
        for block in inf:
            if len([x for x in block.column_iter() if x[0] != "-"]) < int(filter_alignments):
                removed += 1
                outf_rj.write(block)
            else:
                included += 1
                outf.write(block)
    else:
        for block in inf:
            outf.write(block)

    outf.close()
    outf_rj.close()

    return (removed, included)


def extractGeneBlocks(bedfile,
                      maf_file,
                      outfile,
                      primary_species,
                      secondary_species):
    """
    Is based on Dan Blankenburg's interval_to_maf_merged_fasta.py 
    Receives a bed12 containing intervals of interest, and maf containing 
    pairwise genomic alignment. Iterates through bed intervals and outputs 
    intervals in fasta format.
    See comments in maf_utilities.py: maf blocks are always extracted in a 
    positive strand orientation relative to the src alignment, reverse 
    complementing is then done using method AlignedSequence method 
    get_sequence_reverse_complement.
    MAF file is indexed on the fly using bx.align.maf.MultiIndexed
    """
    index, index_filename = build_maf_index(maf_file)
    output = IOTools.openFile(outfile, "w")
    regions_extracted = 0

    # iterate through intervals
    for line_count, line in enumerate(IOTools.openFile(bedfile).readlines()):
        try:
            # retrieve exon starts & ends, plus all gtf fields
            starts, ends, fields = get_starts_ends_fields_from_gene_bed(line)

            # create spliced alignment (N.B. strand always +ve)
            alignment = get_spliced_region_alignment(index,
                                                     primary_species,
                                                     fields[0],
                                                     starts,
                                                     ends,
                                                     strand='+',
                                                     species=[primary_species,
                                                              secondary_species],
                                                     mincols=0,
                                                     overwrite_with_gaps=False)
            primary_name = secondary_name = fields[3]
            alignment_strand = fields[5]
        except Exception, e:
            print "Error loading exon positions from input line %i: %s" % (line_count, e)
            break

        # write the stiched sequence to outfile in the correct orientation
        output.write(">%s.%s\n" % (primary_species, primary_name))
        if alignment_strand == "-":
            output.write(
                alignment.get_sequence_reverse_complement(primary_species))
        else:
            output.write(alignment.get_sequence(primary_species))
        output.write("\n")

        output.write(">%s.%s\n" % (secondary_species, secondary_name))
        if alignment_strand == "-":
            output.write(
                alignment.get_sequence_reverse_complement(secondary_species))
        else:
            output.write(alignment.get_sequence(secondary_species))
        output.write("\n")

        regions_extracted += 1

    output.close()

    return regions_extracted

##########################################################################
##########################################################################


def complement(template):
    """ Generates the reverse complement of a template sequence """
    # Turn the string around using slicing
    backward_template = template[::-1]
    # Create the complementary sequence from the backward template
    reverse_template = backward_template.translate(
        string.maketrans("ACTGactg", "TGACtgac"))

    return reverse_template


def extractMAFGeneBlocks(bedfile,
                         maf_file,
                         genome_file,
                         outfile,
                         primary_species,
                         secondary_species,
                         keep_gaps=True):
    """
    Is based on Dan Blankenburg's interval_to_maf_merged_fasta.py 
    Receives a bed12 containing intervals of interest, and maf containing 
    pairwise genomic alignment. Iterates through bed intervals and outputs 
    intervals in fasta format.
    See comments in maf_utilities.py: maf blocks are always extracted in a 
    positive strand orientation relative to the src alignment, reverse 
    complementing is then done using method AlignedSequence method 
    get_sequence_reverse_complement.
    MAF file is indexed on the fly using bx.align.maf.MultiIndexed
    """
    index, index_filename = build_maf_index(maf_file)
    output = IOTools.openFile(outfile, "w")
    regions_extracted = 0

    # iterate through intervals
    for line_count, line in enumerate(IOTools.openFile(bedfile).readlines()):
        try:
            # retrieve exon starts & ends, plus all gtf fields
            starts, ends, fields = get_starts_ends_fields_from_gene_bed(line)

            # create spliced alignment (N.B. strand always +ve)
            alignment = get_spliced_region_alignment(index,
                                                     primary_species,
                                                     fields[0],
                                                     starts,
                                                     ends,
                                                     strand='+',
                                                     species=[primary_species,
                                                              secondary_species],
                                                     mincols=0,
                                                     overwrite_with_gaps=False)
            primary_name = secondary_name = fields[3]
            alignment_strand = fields[5]
        except Exception, e:
            print "Error loading exon positions from input line %i: %s" % (line_count, e)
            break

        if keep_gaps:
            # write the stiched sequence to outfile in the correct orientation
            output.write(">%s.%s\n" % (primary_species, primary_name))
            if alignment_strand == "-":
                output.write(
                    alignment.get_sequence_reverse_complement(primary_species))
            else:
                output.write(alignment.get_sequence(primary_species))
            output.write("\n")

            output.write(">%s.%s\n" % (secondary_species, secondary_name))
            if alignment_strand == "-":
                output.write(
                    alignment.get_sequence_reverse_complement(secondary_species))
            else:
                output.write(alignment.get_sequence(secondary_species))
            output.write("\n")

        else:
            # create indexed fasta
            fasta = IndexedFasta.IndexedFasta(genome_file)
            # retrieve exon sequence
            exon_list = []
            for exon in alignment.exons:
                exon_list.append(
                    fasta.getSequence(fields[0], "+", exon.start, exon.end))

            # write the stitched sequence to outfile in the correct orientation
            output.write(">%s.%s\n" % (primary_species, primary_name))
            if alignment_strand == "-":
                output.write("".join([complement(x) for x in exon_list]))
            else:
                output.write("".join([x for x in exon_list]))
            output.write("\n")

            output.write(">%s.%s\n" % (secondary_species, secondary_name))
            if alignment_strand == "-":
                output.write(
                    alignment.get_sequence_reverse_complement(secondary_species))
            else:
                output.write(alignment.get_sequence(secondary_species))
            output.write("\n")

        regions_extracted += 1

    output.close()

    return regions_extracted

##########################################################################
##########################################################################


def splitAlignedFasta(infile, out_stub, name_dict):
    """
    Receives a fasta file containing multiple sequence alignments. Splits into 
    multiple outfiles, each containing sequence with the same interval_id.
    Identifiers must be in format >species.interval_id
    name_dict specifies the format for the outfile identifiers
    """
    infile = IOTools.openFile(infile).readlines()

    current_id = ""
    line_number = 0
    for line in infile:
        line = line.rstrip()
        line_number += 1
        if line_number == 1:
            current_id = line.split(".")[1]
            out = open(os.path.join(out_stub, current_id + ".fasta"), "w")
        if line_number % 2 == 1:
            gene_id = line.split(".")[1]
            if gene_id == current_id:
                species = line.split(".")[0]
                out.write(name_dict[species] + "\n")
            else:
                out.close()
                current_id = gene_id
                out = os.path.join(out_stub, current_id + ".fasta")
                if os.path.exists(out):
                    raise IOError("There are two transcript with gene_id %s"
                                  " and transcript_id %s" % current_id.split("__"))
                else:
                    out = open(out, "w")
                    out.write(name_dict[line.split(".")[0]] + "\n")
        else:
            out.write(line + "\n")
    out.close()


def removeGapsFromAlignedFasta(in_dir, out_dir, min_length=0):
    """
    Receives a fasta file containing two or more aligned sequences, removes gaps
    from the reference (first) alignment in file and removes corresponding gaps
    intervals from subsequence sequence.
    Any region of sequence flanked by gaps that is smaller than a specified min
    length may also be removed 
    WARNING: using this function will likely cause frame shifts in resulting 
    sequence.
    """
    if not out_dir.endswith("/"):
        out_dir = out_dir + "/"

    fasta_files = os.listdir(in_dir)
    X = 0
    for fasta in fasta_files:
        X += 1
        if X > 5:
            break
        print fasta
        print os.path.abspath(fasta)
        print os.path.join(in_dir, fasta)
        file_name = P.snip(os.path.basename(fasta),  ".fasta")
        lines = [line.strip() for line in open(fasta).readlines()]
        # reference sequence name
        name = lines[0]
        # the reference sequence
        seq = lines[1]

        # return a list of start, end for ungapped alignment regions in
        # reference
        regex = re.compile("\w+")
        intervals = [(x.start(), x.end()) for x in regex.finditer(seq)]

        x = 1
        for interval in intervals:
            outfile = "".join(file_name, "__", str(x), ".fasta")
            # remove intervals below specified length
            if interval[1] - interval[0] <= int(min_length):
                continue
            else:
                # remove gapped regions from reference and subsequent
                # alignments
                out = open(os.path.join(out_dir, outfile), "w")
                out.write(name + "\n")
                out.write(seq[interval[0]:interval[1]])
                for line in enumerate(lines[2:], 1):
                    if line[0] % 2 == 1:
                        out.write(line + "\n")
                    else:
                        out.write(line[interval[0]:interval[1]])
                out.close()
                x += 1


def runPhyloCSF(in_fasta, tmp_dir, outfile):
    statement = ("zcat %(in_fasta)s |"
                 " %(scriptsdir)s/farm.py"
                 "  --split-at-regex='\n(\s*)\n'"
                 "  --chunksize=1"
                 "  --output-header"
                 "  --tmpdir=%(tmp_dir)s"
                 "  --use-cluster"
                 "  --log=%(outfile)s.log"
                 " PhyloCSF 29mammals"
                 "  --frames=3"
                 "  --removeRefGaps"
                 "  --species=%(species)s"
                 " > %(outfile)s")

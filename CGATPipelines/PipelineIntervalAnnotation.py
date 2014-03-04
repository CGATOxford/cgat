'''
PipelineIntervalAnnotation.py - Tasks associated with annotation of genomic intervals
=====================================================================================
'''

import sys
import tempfile
import optparse
import shutil
import itertools
import csv
import math
import random
import re
import glob
import os
import shutil
import collections
import gzip
import sqlite3
import cStringIO
import fileinput
import CGAT.IOTools as IOTools
import CGAT.IndexedFasta as IndexedFasta
import CGAT.IndexedGenome as IndexedGenome
import CGAT.FastaIterator as FastaIterator
import CGAT.Genomics as Genomics
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import pysam
import numpy
import CGAT.Experiment as E
import CGAT.Pipeline as P

############################################################
############################################################
############################################################
# Pipeline configuration
P.getParameters(["%s.ini" % __file__[:-len(".py")], "pipeline.ini"])
PARAMS = P.PARAMS

############################################################
############################################################
############################################################


def exportIntervalsAsBed(database, query, outfile):
    '''export intervals from SQlite database as bed files. '''

    dbhandle = sqlite3.connect(database)
    cc = dbhandle.cursor()
    cc.execute(query)

    outs = IOTools.openFile(outfile, "w")
    for result in cc:
        contig, start, end, interval_id, score = result
        outs.write("%s\t%i\t%i\t%s\t%s\n" %
                   (contig, start, end, str(interval_id), str(score)))
    cc.close()
    outs.close()

############################################################
############################################################
############################################################


def BedFileVenn(infiles, outfile):
    '''merge :term:`bed` formatted *infiles* by intersection
    and write to *outfile*.

    Only intervals that overlap in all files are retained.
    Interval coordinates are given by the first file in *infiles*.

    Bed files are normalized (overlapping intervals within 
    a file are merged) before intersection. 

    Intervals are renumbered starting from 1.
    '''
    bed1, bed2 = infiles
    liver_name = P.snip(os.path.basename(liver), ".replicated.bed")
    testes_name = P.snip(os.path.basename(testes), ".replicated.bed")
    to_cluster = True

    statement = '''cat %(liver)s %(testes)s | mergeBed -i stdin | awk 'OFS="\\t" {print $1,$2,$3,"CAPseq"NR}' > replicated_intervals/liver.testes.merge.bed;
                   echo "Total merged intervals" > %(outfile)s; cat replicated_intervals/liver.testes.merge.bed | wc -l >> %(outfile)s; 
                   echo "Liver & testes" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.merge.bed -b %(liver)s -u | intersectBed -a stdin -b %(testes)s -u > replicated_intervals/liver.testes.shared.bed; cat replicated_intervals/liver.testes.shared.bed | wc -l >> %(outfile)s; 
                   echo "Testes only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.merge.bed -b %(liver)s -v > replicated_intervals/%(testes_name)s.liver.testes.unique.bed; cat replicated_intervals/%(testes_name)s.liver.testes.unique.bed | wc -l >> %(outfile)s; 
                   echo "Liver only" >> %(outfile)s; intersectBed -a replicated_intervals/liver.testes.merge.bed -b %(testes)s -v > replicated_intervals/%(liver_name)s.liver.testes.unique.bed; cat replicated_intervals/%(liver_name)s.liver.testes.unique.bed | wc -l >> %(outfile)s;                   
                   sed -i '{N;s/\\n/\\t/g}' %(outfile)s; '''

    if len(infiles) == 1:
        shutil.copyfile(infiles[0], outfile)

    elif len(infiles) == 2:

        if P.isEmpty(infiles[0]) or P.isEmpty(infiles[1]):
            P.touch(outfile)
        else:
            statement = '''
        intersectBed -u -a %s -b %s 
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %%(outfile)s 
        ''' % (infiles[0], infiles[1])
            P.run()

    else:

        tmpfile = P.getTempFilename(".")

        # need to merge incrementally
        fn = infiles[0]
        if P.isEmpty(infiles[0]):
            P.touch(outfile)
            return

        statement = '''mergeBed -i %(fn)s > %(tmpfile)s'''
        P.run()

        for fn in infiles[1:]:
            if P.isEmpty(infiles[0]):
                P.touch(outfile)
                os.unlink(tmpfile)
                return

            statement = '''mergeBed -i %(fn)s | intersectBed -u -a %(tmpfile)s -b stdin > %(tmpfile)s.tmp; mv %(tmpfile)s.tmp %(tmpfile)s'''
            P.run()

        statement = '''cat %(tmpfile)s
        | cut -f 1,2,3,4,5 
        | awk 'BEGIN { OFS="\\t"; } {$4=++a; print;}'
        > %(outfile)s '''
        P.run()

        os.unlink(tmpfile)


############################################################
############################################################
############################################################
def makeIntervalCorrelation(infiles, outfile, field, reference):
    '''compute correlation of interval properties between sets
    '''

    dbhandle = sqlite3.connect(PARAMS["database"])

    tracks, idx = [], []
    for infile in infiles:
        track = P.snip(infile, ".bed")
        tablename = "%s_intervals" % P.quote(track)
        cc = dbhandle.cursor()
        statement = "SELECT contig, start, end, %(field)s FROM %(tablename)s" % locals(
        )
        cc.execute(statement)
        ix = IndexedGenome.IndexedGenome()
        for contig, start, end, peakval in cc:
            ix.add(contig, start, end, peakval)
        idx.append(ix)
        tracks.append(track)
    outs = IOTools.openFile(outfile, "w")
    outs.write("contig\tstart\tend\tid\t" + "\t".join(tracks) + "\n")

    for bed in Bed.iterator(infile=open(reference, "r")):

        row = []
        for ix in idx:
            try:
                intervals = list(ix.get(bed.contig, bed.start, bed.end))
            except KeyError:
                row.append("")
                continue

            if len(intervals) == 0:
                peakval = ""
            else:
                peakval = str((max([x[2] for x in intervals])))
            row.append(peakval)

        outs.write(str(bed) + "\t" + "\t".join(row) + "\n")

    outs.close()

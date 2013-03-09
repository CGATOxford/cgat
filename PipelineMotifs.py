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
pipeline_vitaminD_motifs.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python pipeline_vitaminD_motifs.py --help

Type::

   python pipeline_vitaminD_motifs.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys, re, os, tempfile, collections, shutil

import logging as L
import Experiment as E
import Pipeline as P
import sqlite3
import IndexedFasta
import Masker
import Glam2Scan
import MAST
import IOTools
import Bed
import Bioprospector
import FastaIterator
import glob

PARAMS = P.getParameters()

############################################################
############################################################
############################################################
def filterMotifsFromMEME( infile, outfile, selected ):
    '''select motifs from a MEME file and save into outfile
    '''

    outs = open(outfile, "w" )
    if len(selected) == 0:
        outs.close()
        return

    keep = True

    for line in open(infile, "r"):
        if line.startswith( "MOTIF" ):
            motif = re.match( "MOTIF\s+(\d+)", line ).groups()[0]
            if motif in selected:
                keep = True
            else:
                keep = False
        if keep: outs.write( line )
    outs.close()

############################################################
############################################################
############################################################
def maskSequences( sequences, masker = None):
    '''return a list of masked sequence.

    *masker* can be one of
        dust/dustmasker * run dustmasker on sequences
        softmask        * use softmask to hardmask sequences
    '''

    if masker in ("dust", "dustmasker"):
        masker_object = Masker.MaskerDustMasker()
    else:
        masker_object = None

    if masker == "softmask":
        # the genome sequence is repeat soft-masked
        masked_seq = sequences
    elif masker in ("dust", "dustmasker"):
        # run dust
        masked_seq = masker_object.maskSequences( [ x.upper() for x in sequences ] )
    elif masker == None:
        masked_seq = [x.upper() for x in sequences ]
    else:
        raise ValueError("unknown masker %s" % masker )

    # hard mask softmasked characters
    masked_seq = [re.sub( "[a-z]","N", x) for x in masked_seq ]

    return masked_seq

############################################################
############################################################
############################################################
def exportSequencesFromBedFile( infile, outfile, masker = None, mode = "intervals" ):
    '''export sequences for intervals in :term:`bed`-formatted *infile* 
    to :term:`fasta` formatted *outfile*
    '''

    track = P.snip( infile, ".bed.gz" )

    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] ) )
    outs = IOTools.openFile( outfile, "w")

    ids, seqs = [], []
    for bed in Bed.setName(Bed.iterator( IOTools.openFile(infile) )):
        lcontig = fasta.getLength( bed.contig )

        if mode == "intervals":
            seqs.append( fasta.getSequence( bed.contig, "+", bed.start, bed.end) )
            ids.append( "%s_%s %s:%i..%i" % (track, bed.name, bed.contig, bed.start, bed.end) )

        elif mode == "leftright":
            l = bed.end - bed.start

            start, end = max(0,bed.start-l), bed.end-l
            ids.append( "%s_%s_l %s:%i..%i" % (track, bed.name, bed.contig, start, end) )
            seqs.append( fasta.getSequence( bed.contig, "+", start, end) )
            
            start, end = bed.start+l, min(lcontig,bed.end+l)
            ids.append( "%s_%s_r %s:%i..%i" % (track, bed.name, bed.contig, start, end) )
            seqs.append( fasta.getSequence( bed.contig, "+", start, end) )
            
    masked = maskSequences( seqs, masker )
    outs.write("\n".join( [ ">%s\n%s" % (x,y) for x,y in zip(ids, masked) ] ) )

    outs.close()

############################################################
############################################################
############################################################
def writeSequencesForIntervals( track, filename,
                                dbhandle,
                                full = False,
                                halfwidth = None,
                                maxsize = None,
                                proportion = None,
                                masker = [],
                                offset = 0,
                                shuffled = False,
                                min_sequences = None ):
    '''build a sequence set for motif discovery. Intervals are taken 
    from the table <track>_intervals in the database *dbhandle* and 
    save to *filename* in :term:`fasta` format.

    If num_shuffles is set, shuffled copies are created as well with
    the shuffled number appended to the filename.

    The sequences are masked before shuffling (is this appropriate?)

    If *full* is set, the whole intervals will be output, otherwise
    only the region around the peak given by *halfwidth*

    If *maxsize* is set, the output is truncated at *maxsize* characters.

    If proportion is set, only the top *proportion* intervals are output
    (sorted by peakval).

    *masker* can be a combination of 
        * dust, dustmasker: apply dustmasker
        * softmask: mask softmasked genomic regions

    '''

    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] ) )

    cc = dbhandle.cursor()

    if PARAMS["score"] == "peakval":
        orderby = " ORDER BY peakval DESC"
    elif PARAMS["score"] == "max":
        orderby = " ORDER BY score DESC"
    elif PARAMS["score"] == "min":
        orderby = " ORDER BY score ASC"
    else:
        raise ValueError("Unknown value passed as score parameter, check your ini file")
         
    tablename = "%s_intervals" % P.quote( track )
    if full:
        statement = '''SELECT start, end, interval_id, contig 
                       FROM %(tablename)s 
                       ''' % locals() + orderby
    elif halfwidth:
        statement = '''SELECT peakcenter - %(halfwidth)s, peakcenter + %(halfwidth)s,
                       interval_id, contig 
                       FROM %(tablename)s 
                       ''' % locals() + orderby
    else:
        raise ValueError("either specify full or halfwidth" )
    
    cc.execute( statement )
    data = cc.fetchall()
    cc.close()

    if proportion:
        cutoff = int(len(data) * proportion) + 1
        if min_sequences:
            cutoff = max( cutoff, min_sequences )
    else:
        cutoff = len(data)
        L.info( "writeSequencesForIntervals %s: using at most %i sequences for pattern finding" % (track, cutoff) )

    L.info( "writeSequencesForIntervals %s: masker=%s" % (track,str(masker)))

    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"]) )

    sequences = []
    current_size, nseq = 0, 0
    new_data = []
    for start, end, id, contig in data[:cutoff]:
        lcontig = fasta.getLength( contig )
        start, end = max(0, start + offset), min(end + offset, lcontig)
        if start >= end:
            L.info( "writeSequencesForIntervals %s: sequence %s is empty: start=%i, end=%i, offset=%i - ignored" % (track, id, start, end, offset))
            continue
        seq = fasta.getSequence( contig, "+", start, end )
        sequences.append( seq )
        new_data.append( (start, end, id, contig) )
        current_size += len(seq)
        if maxsize and current_size >= maxsize: 
            L.info( "writeSequencesForIntervals %s: maximum size (%i) reached - only %i sequences output (%i ignored)" % \
                        (track, maxsize, nseq, len(data) - nseq ) )
            break
        nseq += 1
        
    data = new_data
            
    if shuffled:
        # note that shuffling is done on the unmasked sequences
        # Otherwise N's would be interspersed with real sequence
        # messing up motif finding unfairly. Instead, masking is
        # done on the shuffled sequence.
        sequences = [ list(x) for x in sequences ]
        for sequence in sequences: random.shuffle(sequence)
        sequences = maskSequences( ["".join(x) for x in sequences ], masker )
        
    c = E.Counter()
    outs = IOTools.openFile(filename, "w" )
    for masker in masker:
        sequences = maskSequences( sequences, masker )

    for sequence, d in zip( sequences, data ):
        c.input += 1
        if len(sequence) == 0: 
            c.empty += 1 
            continue
        start, end, id, contig = d
        id = "%s_%s %s:%i..%i" % (track, str(id), contig, start, end)
        outs.write( ">%s\n%s\n" % (id, sequence ) )
        c.output += 1
    outs.close()
    
    E.info("%s" % c )

    return c.output

############################################################
############################################################
############################################################
def runRegexMotifSearch( infiles, outfile ):
    '''run a regular expression search on sequences.
    compute counts.
    '''

    motif = "[AG]G[GT]T[CG]A"
    reverse_motif = "T[GC]A[CA]C[TC]"

    controlfile, dbfile  = infiles
    if not os.path.exists( controlfile ):
        raise P.PipelineError( "control file %s for %s does not exist" % (controlfile, dbfile))

    motifs = []
    for x in range( 0, 15):
        motifs.append( ("DR%i" % x, re.compile( motif + "." * x + motif, re.IGNORECASE ) ) )
    for x in range( 0, 15):
        motifs.append( ("ER%i" % x, re.compile( motif + "." * x + reverse_motif, re.IGNORECASE ) ) )

    db_positions = Motifs.countMotifs( IOTools.openFile( dbfile, "r" ), motifs )
    control_positions = Motifs.countMotifs( IOTools.openFile( controlfile, "r" ), motifs )

    db_counts, control_counts = Motifs.getCounts( db_positions), Motifs.getCounts( control_positions )
    db_seqcounts, control_seqcounts = Motifs.getOccurances( db_positions), Motifs.getCounts( control_positions )
    
    ndb, ncontrol = len( db_positions), len( control_positions )
    outf = IOTools.openFile( outfile, "w" )
    outf.write( "motif\tmotifs_db\tmotifs_control\tseq_db\tseq_db_percent\tseq_control\tseq_control_percent\tfold\n" )
    for motif, pattern in motifs:
        try:
            fold = float(db_seqcounts[motif]) * ncontrol / (ndb * control_seqcounts[motif])
        except ZeroDivisionError:
            fold = 0
            
        outf.write( "%s\t%i\t%i\t%i\t%s\t%i\t%s\t%5.2f\n" % \
                    (motif,
                     db_counts[motif],
                     control_counts[motif],
                     db_seqcounts[motif],
                     IOTools.prettyPercent( db_seqcounts[motif], ndb),
                     control_seqcounts[motif],
                     IOTools.prettyPercent( control_seqcounts[motif], ncontrol),
                     fold) )


############################################################
############################################################
############################################################
def runGLAM2SCAN( infiles, outfile ):
    '''run glam2scan on all intervals and motifs.
    '''

    to_cluster = True
    # only use new nodes, as /bin/csh is not installed
    # on the old ones.
    # job_options = "-l mem_free=8000M"

    controlfile, dbfile, motiffiles  = infiles
    controlfile = dbfile[:-len(".fasta")] + ".controlfasta"
    if not os.path.exists( controlfile ):
        raise P.PipelineError( "control file %s for %s does not exist" % (controlfile, dbfile))

    if os.path.exists(outfile): os.remove( outfile )
    
    for motiffile in motiffiles:
        of = IOTools.openFile(outfile, "a")
        motif, x = os.path.splitext( motiffile )
        of.write(":: motif = %s ::\n" % motif )
        of.close()
        
        statement = '''
        cat %(dbfile)s %(controlfile)s | %(execglam2scan)s -2 -n %(glam2scan_results)i n %(motiffile)s - >> %(outfile)s
        '''
        P.run()

############################################################
############################################################
############################################################
def loadGLAM2SCAN( infile, outfile ):
    '''parse mast file and load into database.

    Parse several motif runs and add them to the same
    table.
    '''
    tablename = outfile[:-len(".load")]
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    tmpfile.write( "motif\tid\tnmatches\tscore\tscores\tncontrols\tmax_controls\n" )

    lines = IOTools.openFile(infile).readlines()
    chunks = [x for x in range(len(lines)) if lines[x].startswith("::") ]
    chunks.append( len(lines) )

    for chunk in range(len(chunks)-1):

        # use real file, as parser can not deal with a
        # list of lines
        
        try:
            motif = re.match( ":: motif = (\S+) ::", lines[chunks[chunk]]).groups()[0]
        except AttributeError:
            raise P.PipelineError("parsing error in line '%s'" % lines[chunks[chunk]])

        if chunks[chunk]+1 == chunks[chunk+1]:
            L.warn( "no results for motif %s - ignored" % motif )
            continue
        
        tmpfile2 = tempfile.NamedTemporaryFile(delete=False)        
        tmpfile2.write( "".join( lines[chunks[chunk]+1:chunks[chunk+1]]) )
        tmpfile2.close()
        glam = Glam2Scan.parse( IOTools.openFile(tmpfile2.name, "r") )
                    
        os.unlink( tmpfile2.name )        

        # collect control data
        full_matches = collections.defaultdict( list )
        controls = collections.defaultdict( list )
        for match in glam.matches:
            m = match.id.split("_")
            track, id = m[:2]
            if len(m) == 2:
                full_matches[id].append( match )
            else:
                controls[id].append( match.score )

        for id, matches in full_matches.iteritems():
            
            nmatches = len(matches)
            scores = [x.score for x in matches ]
            score = max(scores)
            # move to genomic coordinates
            #contig, start, end = re.match( "(\S+):(\d+)..(\d+)", match.id).groups()
            #start, end = int(start), int(end)
            #match.start += start
            #match.end += start
            contig = ""

            if id not in controls:
                P.warn( "no controls for %s - increase evalue?" % id )
            
            c = controls[id]
            if len(c) == 0: mmax = "" 
            else: mmax = max(c)

            tmpfile.write( "\t".join( map(str,
                    (motif, id, 
                     nmatches,
                     score,
                     ",".join(map(str,scores)),
                     len(c),
                     mmax))) + "\n" )

    tmpfile.close()
    tmpfilename = tmpfile.name

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
              -b sqlite \
              --index=id \
              --index=motif \
              --index=id,motif \
              --table=%(tablename)s \
              --map=base_qualities:text \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run()
    os.unlink( tmpfile.name )


############################################################
############################################################
############################################################
def loadMAST( infile, outfile ):
    '''parse mast file and load into database.

    Parse several motif runs and add them to the same
    table.

    Add columns for the control data as well.
    '''

    tablename = P.toTable( outfile )

    tmpfile = tempfile.NamedTemporaryFile(delete=False)

    tmpfile.write( MAST.Match().header +\
                   "\tmotif\tcontig" \
                   "\tl_evalue\tl_pvalue\tl_nmatches\tl_length\tl_start\tl_end" \
                   "\tr_evalue\tr_pvalue\tr_nmatches\tr_length\tr_start\tr_end" \
                   "\tmin_evalue\tmin_pvalue\tmax_nmatches" + "\n" )

    lines = IOTools.openFile(infile).readlines()
    chunks = [x for x in range(len(lines)) if lines[x].startswith("::") ]
    chunks.append( len(lines) )

    def readChunk( lines, chunk ):
        # use real file, as MAST parser can not deal with a
        # list of lines
        tmpfile2 = tempfile.NamedTemporaryFile(delete=False)        
        try:
            motif, part = re.match( ":: motif = (\S+) - (\S+) ::", lines[chunks[chunk]]).groups()
        except AttributeError:
            raise P.PipelineError("parsing error in line '%s'" % lines[chunks[chunk]])

        E.info( "reading %s - %s" % (motif, part))

        tmpfile2.write( "".join( lines[chunks[chunk]+1:chunks[chunk+1]]) )
        tmpfile2.close()

        mast = MAST.parse( IOTools.openFile(tmpfile2.name, "r") )

        os.unlink( tmpfile2.name )        
    
        return motif, part, mast
    
    def splitId( s, mode ):
        '''split background match id

        has three parts: track _ id _ pos

        track might contain '_'.
        '''
        d = match.id.split("_")
        if mode == "bg":
            return "_".join(d[:-2]), d[-2], d[-1]
        elif mode == "fg":
            return "_".join(d[:-1]), d[-1]

    for chunk in range(0, len(chunks)-1, 2):

        motif_fg, part, mast_fg = readChunk( lines, chunk )
        assert part == "foreground"
        motif_bg, part, mast_bg = readChunk( lines, chunk + 1 )
        assert part == "background"
        assert motif_fg == motif_bg

        # index control data
        controls = collections.defaultdict( dict )
        for match in mast_bg.matches:
            track, id, pos = splitId( match.id, "bg" )
            controls[id][pos] = (match.evalue, match.pvalue, match.nmotifs, match.length, match.start, match.end )

        for match in mast_fg.matches:
            # remove track and pos
            track, match.id = splitId( match.id, "fg" )
            # move to genomic coordinates
            contig, start, end = re.match( "(\S+):(\d+)..(\d+)", match.description).groups()
            if match.nmotifs > 0:
                start, end = int(start), int(end)
                match.start += start
                match.end += start
                match.positions = [ x + start for x in match.positions ]

            id = match.id
            if id not in controls:
                P.warn( "no controls for %s - increase MAST evalue" % id )

            if "l" not in controls[id]: controls[id]["l"] = (float(PARAMS["mast_evalue"]), 1, 0, 0, 0, 0)
            if "r" not in controls[id]: controls[id]["r"] = (float(PARAMS["mast_evalue"]), 1, 0, 0, 0, 0)

            min_evalue = min( controls[id]["l"][0], controls[id]["r"][0])
            min_pvalue = min( controls[id]["l"][1], controls[id]["r"][1])
            max_nmatches = max( controls[id]["l"][2], controls[id]["r"][2])

            tmpfile.write( str(match) + "\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
                               (motif_fg, contig,
                                "\t".join( map(str, controls[id]["l"] )),
                                "\t".join( map(str, controls[id]["r"] )),
                                str(min_evalue),
                                str(min_pvalue),
                                str(max_nmatches),
                                ) + "\n" )
            
    tmpfile.close()
    tmpfilename = tmpfile.name

    statement = '''
    python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
              -b sqlite 
              --index=id 
              --index=motif 
              --index=id,motif 
              --table=%(tablename)s 
              --allow-empty
              --map=base_qualities:text 
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run()
    os.unlink( tmpfile.name )

############################################################
############################################################
############################################################
def runBioProspector( infiles, outfile, dbhandle ):
    '''run bioprospector for motif discovery.

    Bioprospector is run on only the top 10% of peaks.
    '''

    # bioprospector currently not working on the nodes
    to_cluster = False

    # only use new nodes, as /bin/csh is not installed
    # on the old ones.
    # job_options = "-l mem_free=8000M"

    tmpfasta = P.getTempFilename( "." )
    track = outfile[:-len(".bioprospector")]
    nseq = writeSequencesForIntervals( track,
                                       tmpfasta,
                                       dbhandle,
                                       full = True,
                                       masker = "dust",
                                       proportion = PARAMS["bioprospector_proportion"] )

    if nseq == 0:
        E.warn( "%s: no sequences - bioprospector skipped" % track )
        P.touch( outfile )
    else:
        statement = '''
    BioProspector -i %(tmpfasta)s %(bioprospector_options)s -o %(outfile)s > %(outfile)s.log
    '''
        P.run()

    os.unlink( tmpfasta )

############################################################
############################################################
############################################################
##
############################################################
def loadBioProspector( infile, outfile ):
    '''load results from bioprospector.'''

    tablename = outfile[:-len(".load")]
    target_path = os.path.join( os.path.abspath( PARAMS["exportdir"] ), "bioprospector" )

    try:
        os.makedirs( target_path )
    except OSError: 
        pass

    track = infile[:-len(".bioprospector")]

    results = Bioprospector.parse( IOTools.openFile(infile, "r") )

    tmpfile = P.getTempFile()
    tmpfile.write( "id\tmotif\tstart\tend\tstrand\tarrangement\n" )
    
    for x, motifs in enumerate( results ):
        outname = os.path.join( target_path, "%s_%02i.png" % (track, x ) )
        Bioprospector.build_logo( [y.sequence for y in motifs.matches],
                                  outname )

        for match in motifs.matches:

            distance = abs(match.start + match.width1 - (match.end - match.width2))

            if match.strand in ( "+-", "-+"):
                arrangement = "ER" 
            elif match.strand in ("++", "--"):
                arrangement = "DR"
            else:
                arrangement = "SM"
                distance = 0
                
            arrangement += "%i" % distance
            strand = match.strand[0]

            id = re.sub( ".*_", "", match.id )
            tmpfile.write( "%s\t%i\t%i\t%i\t%s\t%s\n" % \
                           (id,
                            x,
                            match.start,
                            match.end,
                            strand,
                            arrangement) )
    tmpfile.close()
    tmpfilename = tmpfile.name

    statement = '''
   python %(scriptsdir)s/csv2db.py %(csv2db_options)s \
    --allow-empty \
    -b sqlite \
    --index=id \
    --index=motif \
    --index=id,motif \
    --table=%(tablename)s \
    < %(tmpfilename)s > %(outfile)s
    '''

    P.run()


############################################################
############################################################
############################################################
def runMAST( infiles, outfile ):
    '''run mast on all intervals and motifs.

    Collect all results for an E-value up to 10000 so that
    all sequences are output and MAST curves can be computed. 

    10000 is a heuristic.
    '''
    to_cluster = True

    # job_options = "-l mem_free=8000M"

    controlfile, dbfile, motiffiles  = infiles

    if IOTools.isEmpty( dbfile ):
        P.touch( outfile )
        return

    if not os.path.exists( controlfile ):
        raise P.PipelineError( "control file %s for %s does not exist" % (controlfile, dbfile))

    # remove previous results
    if os.path.exists(outfile): os.remove( outfile )
    
    tmpdir = P.getTempDir( "." )
    tmpfile = P.getTempFilename( "." )

    for motiffile in motiffiles:
        if IOTools.isEmpty( motiffile ):
            L.info( "skipping empty motif file %s" % motiffile )
            continue
        
        of = IOTools.openFile(tmpfile, "a")
        motif, x = os.path.splitext( motiffile )
        of.write(":: motif = %s - foreground ::\n" % motif )
        of.close()
        
        # mast bails if the number of nucleotides gets larger than
        # 2186800982?
        # To avoid this, run db and control file separately.
        statement = '''
        cat %(dbfile)s 
        | mast %(motiffile)s - -nohtml -oc %(tmpdir)s -ev %(mast_evalue)f %(mast_options)s >> %(outfile)s.log 2>&1;
        cat %(tmpdir)s/mast.txt >> %(tmpfile)s 2>&1
        '''
        P.run()

        of = IOTools.openFile(tmpfile, "a")
        motif, x = os.path.splitext( motiffile )
        of.write(":: motif = %s - background ::\n" % motif )
        of.close()

        statement = '''
        cat %(controlfile)s 
        | mast %(motiffile)s - -nohtml -oc %(tmpdir)s -ev %(mast_evalue)f %(mast_options)s >> %(outfile)s.log 2>&1;
        cat %(tmpdir)s/mast.txt >> %(tmpfile)s 2>&1
        '''
        P.run()

    statement = "gzip < %(tmpfile)s > %(outfile)s" 
    P.run()

    shutil.rmtree( tmpdir )
    os.unlink( tmpfile )


############################################################
############################################################
############################################################
def runGLAM2( infile, outfile, dbhandle ):
    '''run glam2 on all intervals and motifs.

    In order to increase the signal/noise ratio,
    MEME is not run on all intervals but only the 
    top 10% of intervals (peakval) are used. 
    Also, only the segment of 200 bp around the peak
    is used and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker
    '''
    to_cluster = True

    target_path = os.path.join( os.path.abspath( PARAMS["exportdir"] ), "glam2", outfile )
    track = infile[:-len(".fasta")]

    tmpdir = tempfile.mkdtemp()
    tmpfasta =  os.path.join( tmpdir, "in.fa")
    
    nseq = PipelineMotifs.writeSequencesForIntervals( track, tmpfasta,
                                                      dbhandle,
                                                      full = False,
                                                      halfwidth = int(PARAMS["meme_halfwidth"]),
                                                      maxsize = int(PARAMS["meme_max_size"]),
                                                      proportion = PARAMS["meme_proportion"] )

    min_sequences = int(nseq / 10.0)
    statement = '''
    %(execglam2)s -2 -O %(tmpdir)s %(glam2_options)s -z %(min_sequences)i n %(tmpfasta)s > %(outfile)s.log
    '''
    P.run()

    # copy over results
    try:
        os.makedirs( os.path.dirname( target_path ) )
    except OSError: 
        # ignore "file exists" exception
        pass

    if os.path.exists( target_path ): shutil.rmtree( target_path )
    shutil.move( tmpdir, target_path )

    shutil.copyfile( os.path.join(target_path, "glam2.txt"), outfile)


############################################################
############################################################
############################################################
def collectMEMEResults( tmpdir, target_path, outfile ):
    '''collect output from a MEME run in tmpdir
    and copy all over to target_path

    convert images output by MEME (.eps files) to 
    .png files.'''

    # copy over results
    try:
        os.makedirs( os.path.dirname( target_path ) )
    except OSError: 
        # ignore "file exists" exception
        pass

    if os.path.exists( target_path ): shutil.rmtree( target_path )
    shutil.move( tmpdir, target_path )

    shutil.copyfile( os.path.join(target_path, "meme.txt"), outfile)

    # convert images to png
    epsfiles = glob.glob( os.path.join( target_path, "*.eps" ) )

    for epsfile in epsfiles:
        b, ext = os.path.splitext( epsfile )
        pngfile = b + ".png" 
        statement = '''convert %(epsfile)s %(pngfile)s" '''
        P.run()
    
############################################################
############################################################
############################################################
def runMEME( track, outfile, dbhandle ):
    '''run MEME to find motifs.

    In order to increase the signal/noise ratio,
    MEME is not run on all intervals but only the 
    top 10% of intervals (peakval) are used. 
    Also, only the segment of 200 bp around the peak
    is used and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker

    This method is deprecated - use runMEMEOnSequences instead.
    '''
    to_cluster = True
    # job_options = "-l mem_free=8000M"

    target_path = os.path.join( os.path.abspath(PARAMS["exportdir"]), "meme", outfile )

    fasta = IndexedFasta.IndexedFasta( os.path.join( PARAMS["genome_dir"], PARAMS["genome"] ) )

    tmpdir = P.getTempDir( "." )
    tmpfasta =  os.path.join( tmpdir, "in.fa")
    
    nseq = writeSequencesForIntervals( track, tmpfasta,
                                       dbhandle,
                                       full = False,
                                       masker = PARAMS['motifs_masker'],
                                       halfwidth = int(PARAMS["meme_halfwidth"]),
                                       maxsize = int(PARAMS["meme_max_size"]),
                                       proportion = PARAMS["meme_proportion"],
                                       min_sequences = PARAMS["meme_min_sequences"] )

    if nseq == 0:
        E.warn( "%s: no sequences - meme skipped" % outfile)
        P.touch( outfile )
    else:
        statement = '''
        meme %(tmpfasta)s -dna -revcomp -mod %(meme_model)s -nmotifs %(meme_nmotifs)s -oc %(tmpdir)s -maxsize %(meme_max_size)s %(meme_options)s > %(outfile)s.log
        '''
        P.run()

        collectMEMEResults( tmpdir, target_path, outfile )

############################################################
############################################################
############################################################
def runMEMEOnSequences( infile, outfile ):
    '''run MEME to find motifs.

    In order to increase the signal/noise ratio,
    MEME is not run on all intervals but only the 
    top 10% of intervals (peakval) are used. 
    Also, only the segment of 200 bp around the peak
    is used and not the complete interval.

    * Softmasked sequence is converted to hardmasked
      sequence to avoid the detection of spurious motifs.

    * Sequence is run through dustmasker
    '''
    to_cluster = True
    # job_options = "-l mem_free=8000M"

    nseqs = int(FastaIterator.count( infile ))
    if nseqs == 0:
        E.warn( "%s: no sequences - meme skipped" % outfile)
        P.touch( outfile )
        return
    
    target_path = os.path.join( os.path.abspath(PARAMS["exportdir"]), "meme", outfile )
    tmpdir = P.getTempDir( "." )

    statement = '''
        meme %(infile)s -dna -revcomp 
                        -mod %(meme_model)s 
                        -nmotifs %(meme_nmotifs)s 
                        -oc %(tmpdir)s 
                        -maxsize %(motifs_max_size)s 
                        %(meme_options)s 
       > %(outfile)s.log
    '''

    P.run()

    collectMEMEResults( tmpdir, target_path, outfile )

############################################################
############################################################
############################################################
def runTomTom( infile, outfile ):
    '''compare ab-initio motifs against tomtom.'''

    tmpdir = P.getTempDir( "." )
    
    to_cluster = True
    databases = " ".join(P.asList( PARAMS["tomtom_databases"]))

    target_path = os.path.join( os.path.abspath(PARAMS["exportdir"]), "tomtom", outfile )
    
    if IOTools.isEmpty( infile ):
        E.warn( "input is empty - no computation performed" )
        P.touch( outfile )
        return

    statement = '''
           tomtom %(tomtom_options)s -oc %(tmpdir)s %(infile)s %(databases)s > %(outfile)s.log
    '''

    P.run()
    
    # copy over results
    try:
        os.makedirs( os.path.dirname( target_path ) )
    except OSError: 
        # ignore "file exists" exception
        pass

    if os.path.exists( target_path ): shutil.rmtree( target_path )
    shutil.move( tmpdir, target_path )

    shutil.copyfile( os.path.join(target_path, "tomtom.txt"), outfile)
    


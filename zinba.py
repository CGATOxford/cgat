############################################################
############################################################
############################################################
## Find intervals using SICER with input
@follows( normaliseBAMs )
@files( [ (("%s/bam/%s.norm.bam" % (x, x.asFile()), "%s/bam/%s.norm.bam" % (getControl(x), getControl(x).asFile())), 
           "%s/sicer/%s.sicer" % (x, x.asFile()) ) for x in TRACKS ] )
def runSICER( infiles, outfile ):
    '''Run SICER for peak detection.'''
    infile, controlfile = infiles
    to_cluster = False

    track = P.snip( os.path.basename(infile), ".norm.bam" )
    control = P.snip( os.path.basename(controlfile), ".norm.bam" )
    inputdir = os.path.dirname(outfile)
    try: os.mkdir( track )
    except OSError: pass
    try: os.mkdir( '''%(track)s/sicer''' % locals() )
    except OSError: pass

    # convert bam to bed
    statement = '''bamToBed -i %(infile)s > %(track)s/sicer/%(track)s.bed; 
                   bamToBed -i %(controlfile)s > %(track)s/sicer/%(control)s.bed; '''

    # Run SICER
    statement += '''cd %(inputdir)s; SICER.sh . %(track)s.bed %(control)s.bed . %(genome)s %(sicer_params)s >& %(track)s.sicer''' 
    P.run() 
    
############################################################
@transform(runSICER, regex(r"(\S+)/sicer/(\S+).sicer"), r"\1/sicer/\2.sicer.load" )
def loadSICER( infile, outfile ):
    track = P.snip( os.path.basename(infile), ".sicer" )
    sicerdir = os.path.dirname(infile)
    window = PARAMS["sicer_window"]
    gap = PARAMS["sicer_gap"]
    FDR = PARAMS["sicer_fdr"]
    bedfile = sicerdir + "/" + track + "-W" + str(window) + "-G" + str(gap) + "-islands-summary-FDR" + str(FDR)
    tablename = "%s_sicer_intervals" % track
    headers="contig,start,stop,chip_reads,control_reads,pvalue,fold,fdr"

    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --header=%(headers)s
                       --index=contig,start
                       --table=%(tablename)s
                       --allow-empty 
                   < %(bedfile)s > %(outfile)s'''
    P.run()
    
    
############################################################
@merge( runSICER, "sicer.summary" )
def summarizeSICER( infiles, outfile ):
    '''run SICER for peak detection.'''
    def __get( line, stmt ):
        x = line.search(stmt )
        if x: return x.groups() 

    map_targets = [
        ("Window average: (\d+)", "window_mean",()),
        ("Minimum num of tags in a qualified window:  (\d+)", "window_min",()),
        ("The score threshold is:  (\d+)", "score_threshold",()),
        ("Total number of islands:  (\d+)","total_islands", ()),
        ("chip library size   (\d+)", "chip_library_size", ()),
        ("control library size   (\d+)", "control_library_size", ()),
        ("Total number of chip reads on islands is:  (\d+)", "chip_island_reads", ()),
        ("Total number of control reads on islands is:  (\d+)", "control_island_reads", ()),
        ("Given significance 0.01 ,  there are (\d+) significant islands",  "significant_islands", ()) ]

    mapper, mapper_header = {}, {}
    for x,y,z in map_targets: 
        mapper[y] = re.compile( x )
        mapper_header[y] = z

    keys = [ x[1] for x in map_targets ]

    outs = open(outfile,"w")

    headers = []
    for k in keys:
        if mapper_header[k]:
            headers.extend( ["%s_%s" % (k,x) for x in mapper_header[k] ])
        else:
            headers.append( k )
    outs.write("track\t%s" % "\t".join(headers) + "\n" )

    for infile in infiles:
        results = collections.defaultdict(list)
        with open( infile ) as f:
            for line in f:
                if "diag:" in line: break
                for x,y in mapper.items():
                    s = y.search( line )
                    if s: 
                        results[x].append( s.groups()[0] )
                        break
                
        row = [ P.snip( os.path.basename(infile), ".sicer" ) ]
        for key in keys:
            val = results[key]
            if len(val) == 0: v = "na"
            else: 
                c = len(mapper_header[key])
                if c >= 1: assert len(val) == c, "key=%s, expected=%i, got=%i, val=%s, c=%s" %\
                   (key,
                    len(val),
                    c,
                    str(val), mapper_header[key])
                v = "\t".join( val )
            row.append(v)
        outs.write("\t".join(row) + "\n" )

    outs.close()


############################################################
@transform( summarizeSICER, suffix(".summary"), "_summary.load" )
def loadSICERSummary( infile, outfile ):
    '''load sicer summary.'''
    
    table = P.snip( os.path.basename(outfile), ".load" )
    statement = '''python %(scriptsdir)s/csv2db.py %(csv2db_options)s
                      --index=track 
                      --table=%(table)s
                   < %(infile)s > %(outfile)s'''
    P.run()

############################################################
@transform( loadSICER, regex(r"(\S+)/sicer/(\S+).sicer.load"), r"\1/sicer/\2.sicer.bed" )
def exportSicerAsBed( infile, outfile ):
    '''export locations for all intervals.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )
    track = P.snip( os.path.basename(infile), ".sicer.load" ).replace("-","_")

    cc = dbhandle.cursor()
    statement = "SELECT contig, start, stop FROM %s_sicer_intervals ORDER by contig, start" % track
    cc.execute( statement )

    outs = open( outfile, "w")
    for result in cc:
        contig, start, stop = result
        outs.write( "%s\t%i\t%i\n" % (contig, start, stop) )
    cc.close()
    outs.close()

############################################################
############################################################
############################################################
## Run Zinba to call peaks
@follows( dedup )
@files( [ (("%s/bam/%s.norm.bam" % (x, x.asFile()), "%s/bam/%s.norm.bam" % (getControl(x), getControl(x).asFile())), 
           "%s/zinba/%s.peaks" % (x, x.asFile()) ) for x in TRACKS ] )
def runZinba( infiles, outfile ):
    '''Run Zinba for peak detection.'''
    infile, controlfile = infiles
    to_cluster = False

    track = P.snip( os.path.basename(infile), ".norm.bam" )
    control = P.snip( os.path.basename(controlfile), ".norm.bam" )
    inputdir = os.path.dirname(outfile)
    frag_len = PARAMS['zinba_fragment_size']
    mappability = PARAMS['zinba_mappability']
    genome = PARAMS['zinba_genome']

    try: os.mkdir( track )
    except OSError: pass
    try: os.mkdir( '''%(track)s/zinba''' % locals() )
    except OSError: pass
    try: os.mkdir( '''%(track)s/zinba/map_ext%(frag_len)s''' % locals() )
    except OSError: pass

    # convert bam to bed
    statement = '''bamToBed -i %(infile)s > %(track)s/zinba/%(track)s.bed; 
                   bamToBed -i %(controlfile)s > %(track)s/zinba/%(control)s.bed; '''
    P.run()

    # Run Zinba
    R.library( 'zinba' )
    R( '''generateAlignability( mapdir='%(mappability)s', outdir='%(track)s/zinba/map_ext%(frag_len)s', athresh=1, extension=%(frag_len)s, twoBitFile='%(genome)s' )''' % locals() )
    R( '''basealigncount( inputfile='%(track)s/zinba/%(track)s.bed', outputfile='%(track)s/zinba/%(track)s.basecount', extension=%(frag_len)s, filetype='bed', twoBitFile='%(genome)s' )'''  % locals() )
    R( '''zinba( refinepeaks=1, seq='%(track)s/zinba/%(track)s.bed', input='%(track)s/zinba/%(control)s.bed', filetype='bed',  align='%(track)s/zinba/map_ext%(frag_len)s', twoBit='%(genome)s', outfile='%(track)s/zinba/%(track)s', extension=%(frag_len)s, basecountfile='%(track)s/zinba/%(track)s.basecount', numProc=4, threshold=0.01, broad=FALSE, printFullOut=1, interaction=FALSE, mode='peaks', FDR=TRUE) '''  % locals() )

############################################################
@transform(runZinba, regex(r"(\S+)/zinba/(\S+).peaks"), r"\1/zinba/\2.zinba.load" )
def loadZinba( infile, outfile ):
    track = P.snip( os.path.basename(infile), ".peaks" )
    #zinbadir = os.path.dirname(infile)

    tablename = "%s_zinba_intervals" % track
    headers="peakid,contig,start,stop,strand,sig,maxloc,maxval,pstart,pstop,median,qvalue"

    statement = '''cat | sed -v 'PEAKID'
                   | python %(scriptsdir)s/csv2db.py %(csv2db_options)s 
                       --header=%(headers)s
                       --index=contig,start
                       --table=%(tablename)s
                       --allow-empty 
                   > %(outfile)s'''
    P.run()

############################################################
@transform( loadZinba, regex(r"(\S+)/zinba/(\S+).zinba.load"), r"\1/zinba/\2.zinba.bed" )
def exportZinbaAsBed( infile, outfile ):
    '''export locations for all intervals.'''

    dbhandle = sqlite3.connect( PARAMS["database"] )

    track = P.snip( os.path.basename(infile), ".zinba.load" ).replace("-","_")

    cc = dbhandle.cursor()
    statement = "SELECT contig, pstart, pstop FROM %s_zinba_intervals ORDER by contig, start" % track
    cc.execute( statement )

    outs = open( outfile, "w")

    for result in cc:
        contig, start, stop = result
        outs.write( "%s\t%i\t%i\n" % (contig, start, stop) )
    cc.close()
    outs.close()



############################################################
############################################################
############################################################
## Compare intervals from different peak callers
@transform( mergeIntervals, regex(r"(\S+)/intervals/(\S+).merged.bed"), r"\1/intervals/\2.macs.sicer.bed")
def getSicerOverlap(infile, outfile):
    '''identify intervals overlapping SICER intervals for each datasets'''
    track = P.snip( os.path.basename( infile ), ".merged.bed")
    macsdir = os.path.dirname(infile)
    sicer = track + ".sicer.bed"
    sicerdir = macsdir.replace("macs","sicer")
    statement = '''intersectBed -a %(infile)s -b %(sicerdir)s/%(sicer)s -u > %(outfile)s; '''
    P.run()

############################################################
@transform( getSicerOverlap, regex(r"(\S+)/intervals/(\S+).macs.sicer.bed"), r"\1/intervals/\2.macs.sicer.load")
def loadSicerIntervals(infile, outfile):
    '''Load intervals overlapping SICER intervals into database '''
    track = P.snip( os.path.basename( infile ), ".macs.sicer.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_macs_sicer_intervals_shared
                      --header=%(header)s
                      --ignore-empty
                   > %(outfile)s '''
    P.run()

############################################################
@transform( mergeIntervals, regex(r"(\S+)/intervals/(\S+).merged.bed"), r"\1/intervals/\2.macs.zinba.bed")
def getZinbaOverlap(infile, outfile):
    '''identify intervals overlapping ZINBA intervals for each datasets'''
    track = P.snip( os.path.basename( infile ), ".merged.bed")
    macsdir = os.path.dirname(infile)
    zinba = track + ".zinba.bed"
    zinbadir = macsdir.replace("macs","zinba")
    statement = '''intersectBed -a %(infile)s -b %(zinbadir)s/%(zinba)s -u > %(outfile)s; '''
    P.run()

############################################################
@transform( getZinbaOverlap, regex(r"(\S+)/intervals/(\S+).macs.zinba.bed"), r"\1/macs/\2.intervals.zinba.load")
def loadZinbaIntervals(infile, outfile):
    '''Load intervals overlapping ZINBA intervals into database '''
    track = P.snip( os.path.basename( infile ), ".macs.zinba.bed" ).replace(".","_").replace("-","_")
    header = "contig,start,stop"
    statement = '''cat %(infile)s | python %(scriptsdir)s/csv2db.py
                      --table=%(track)s_macs_zinba_intervals_shared
                      --header=%(header)s
                      --ignore-empty
                   > %(outfile)s '''
    P.run()

@follows( runSICER, loadSICER, 
          summarizeSICER, loadSICERSummary, 
          exportSicerAsBed )
def buildIntervalsSicer():
    '''Find peaks using SICER using control sample'''
    pass

@follows( runZinba, loadZinba, exportZinbaAsBed )
def buildIntervalsZinba():
    '''Find peaks using ZINBA using control sample'''
    pass


@follows( getSicerOverlap, loadSicerIntervals )
def compareCallers():
    '''Compare intervals from different peak callers'''
    pass

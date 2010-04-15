import os, sys, re, csv
import sqlite3

import Stats, IOTools
import Experiment as E
import Pipeline as P

PARAMS = P.getParameters()

############################################################
############################################################
############################################################
def splitExpressionTrack( track ):
    '''split track from expression experiments into fields cellline, condition, replicate'''
    assert track.startswith( "exp" )
    cellline, timepoint, stimulus = track[3:].split("_")
    return cellline, timepoint, stimulus

def getExpressionMatch( infile ):
    '''from an expression-diff get track, method and control.'''
    assert infile.endswith( ".expdiff" )
    track, method = infile[:-len(".expdiff")].split(".")
    control = getExpressionControl( track )
    return track, method, control

def getExpressionControl( track ):
    '''return appropriate control for a track.'''
    cellline, timepoint, stimulus = splitExpressionTrack( track )
    # all stimulus have D3 as control
    return "exp" + "_".join( (cellline, "00", "D3") )

def getExpressionReplicates( track ):
    
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
    cc.execute( "PRAGMA table_info(%(track)s_levels)" % locals() )
    replicates_sample = [ x[1] for x in cc ]
    cc.close()
    return sorted([x for x in replicates_sample if re.match( "R\d", x) ])

############################################################
############################################################
############################################################
##
############################################################
def getExpressionMeasurements( track ):
    '''return a tuple (probesets, treatments, controls)
    where probesets is an array of n probsets and
    treatments/controls are tuples of n-length arrays.'''

    control = getExpressionControl( track )
    if track == control:
        raise P.PipelineError( "track (%s) == control (%s)" % (track, control) )

    replicates_sample = getExpressionReplicates( track )
    assert len(replicates_sample) > 0, "no replicates for sample %s" % track
    
    replicates_control = getExpressionReplicates( control )
    assert len(replicates_control) > 0, "no replicates for control %s" % control
    
    dbhandle = sqlite3.connect( PARAMS["database"] )
    cc = dbhandle.cursor()
 
    track_samples = ",".join( ["a.%s" % x for x in replicates_sample ] )
    control_samples = ",".join( ["b.%s" % x for x in replicates_control ] )

    cc = dbhandle.cursor()
    statement = """SELECT a.cluster_id, 
                        %(track_samples)s, 
                        %(control_samples)s
                        FROM %(track)s_levels AS a, 
                             %(control)s_levels AS b 
                        WHERE a.cluster_id = b.cluster_id""" % locals()

    cc.execute( statement )

    r = zip(*cc.fetchall())
    nreplicates_sample = len(replicates_sample)
    return (control, r[0], r[1:nreplicates_sample+1], r[nreplicates_sample+1:])


############################################################
############################################################
############################################################
def executeGO( outfile,
               statement_foreground,
               statement_background,
               go_file  ):
    '''check for GO enrichment within differentially expressed genes.'''

    to_cluster = True
    
    dbhandle = sqlite3.connect( PARAMS["database"] )

    cc = dbhandle.cursor()
    fg = set( [x[0] for x in cc.execute( statement_foreground).fetchall() ] )
    bg = set( [x[0] for x in cc.execute( statement_background).fetchall() ] )

    if len(fg) == 0:
        P.touch( outfile )
        return

    outdir = os.path.abspath( outfile + ".dir" )

    try: os.makedirs( outdir )
    except OSError: pass

    fg_file = os.path.join( outdir, "foreground" )
    bg_file = os.path.join( outdir, "background" )
    outf = open( fg_file, "w")
    outf.write("\n".join( map(str, fg ) ) + "\n" )
    outf.close()
    outf = open( bg_file, "w")
    outf.write("\n".join( map(str, bg ) ) + "\n" )
    outf.close()
    
    statement = '''
    python %(scriptsdir)s/GO.py \
    --filename-input=%(go_file)s \
    --genes=%(fg_file)s \
    --background=%(bg_file)s \
    --sample=1000 \
    --fdr \
    --filename-ontology=%(filename_ontology)s \
    --output-filename-pattern='%(outdir)s/%%(go)s.%%(section)s' \
    > %(outfile)s'''

    P.run( **dict( locals().items() + PARAMS.items() ) )    

############################################################
############################################################
############################################################
##
############################################################
def buildExpressionTracks( infile, outfiles, map_exp2columns, suffix ):
    '''build expression tracks.

    read the analysis from FILENAME_EXPRESSION
    
    ..note::
       The file A589_Data_RMA.csv does NOT always contain the probeset_id 
       in the first column, but instead it might be the transcript_cluster_id.
       A possible explanation is that if several probesets map to the same
       transcript cluster, the transcript cluster is normalized.
       
       The set of cluster_id and probeset ids are completely non-overlapping.

    Hence, the :term:`cluster_id` will be used.
    '''

    E.info( "importing expression data from %s" % infile )
    
    dbhandle = sqlite3.connect( PARAMS["database"] )

    cc = dbhandle.cursor()
    statement = "SELECT DISTINCT probeset, cluster_id, transcript_id FROM probeset2transcript"
    cc.execute( statement )
    map_cluster2transcript, map_probeset2cluster = {}, {}
    for probeset, cluster, transcript_id in cc.fetchall():
        map_probeset2cluster[probeset] = cluster
        map_cluster2transcript[cluster] = transcript_id

    reader = csv.reader( open(infile,"rU") )

    first = True
    # do not delete old files as this function is called several times
    output_files = IOTools.FilePool( output_pattern = "exp%s.data", force = False )

    headers = (
        ("Probe Set ID", "cluster_id"),
        ("Gene Symbol", "genesymbol"),
        ("mRna - Description", "description"),
        ('mRNA Accession',  'mrna_id'),
        ('mRNA  Source', 'source' ),
        ('mRNA - xhyb', 'xhyb'),
        ('GO Biological Process ID', 'go_biol_id'),
        ('GO Biological Process Term', 'go_biol_term'),    
        ('GO Cellular Component ID', 'go_cell_id'),  
        ('GO Cellular Component Term', 'go_cell_term'), 
        ('GO Molecular Function ID', 'go_mol_id'),
        ('GO Molecular Function Term', 'go_mol_term'),
        ('Pathway Source', 'pw_source' ),       
        ('Pathway Name', 'pw_name' ) )

    old_headers = set( [x[0] for x in headers] )
    new_headers = [x[1] for x in headers]
    take = []
    index_soure, index_accession, index_probeset = None, None, None
    counts = E.Counter()
    found = set()

    outf = open( outfiles[0] + suffix, "w")
    outf.write( "# %s\n" % infile )
    outs = open( outfiles[1] + suffix, "w")
    outs.write( "# %s\n" % infile )
    
    writer = csv.writer( outf )

    for row in reader:
        if first:
            first = False
            writer.writerow( row )

            for x, old_header in enumerate(row ):
                if old_header == "mRNA  Source": index_source = len(take)
                if old_header == "mRNA Accession": index_accession = len(take)
                if old_header == "Probe Set ID": index_probeset = len(take)
                if old_header in old_headers: take.append( x )

            # write headers to all files
            outs.write("\t".join(new_headers)+ "\n")

            for exp,columns in map_exp2columns.items():
                output_files.write( exp, 
                                    "\t".join( ("cluster_id",
                                                Stats.Summary().getHeader(),
                                                "\t".join(["R%i" % i for i in range(len(columns))])))+ "\n")
        else:
            new_row = []
            for x in take:
                if row[x].strip() != "---":
                    new_row.append(row[x].strip())
                else:
                    new_row.append("")

            probeset = new_row[index_probeset].strip()
            if probeset in map_probeset2cluster:
                probeset = map_probeset2cluster[probeset]
                counter.mapped_to_cluster += 1
                
            if probeset not in map_cluster2transcript:
                writer.writerow( row )
                counts.skipped += 1
                continue 
            else:
                if probeset in found:
                    counts.duplicates += 1
                counts.output += 1
                found.add(probeset)

            outs.write("\t".join( new_row )+ "\n")

            for exp,cols in map_exp2columns.items():
                data = [row[x] for x in cols ]
                output_files.write( exp, "\t".join( (probeset,
                                                     str(Stats.Summary([float(x) for x in data ])),
                                                     "\t".join( data ) )) + "\n" )

                
    outf.close()
    if counts.duplicates > 0:
        P.warn( "duplicate probeset/clusters" )

    P.info( "probeset source information: %s" % str(counts) )
    output_files.close()

############################################################
############################################################
############################################################
##
############################################################
def runGO( infile, outfile, go_file ):

    assert infile.endswith( "_expdiff.import" )
    track, method, control = getExpressionMatch(
        infile[:-len("_expdiff.import")]+".expdiff")

    if track == control:
        P.touch( outfile )
        return

    tablename = "%s_vs_%s_%s" % (track,control,method)

    if "sam" in infile:
        # use <=, as qvalues has been set artificially
        # to the threshold as siggenes sometimes reports
        # the qvalue to be higher than the threshold.
        qvalue = PARAMS["sam_fdr"]
        # do not use DISTINCT as it is much slower
        statement_fg = '''SELECT m.gene_id FROM
        %(tablename)s AS e,
        probeset2transcript AS m
        WHERE e.cluster_id = m.cluster_id AND
        e.qvalue <= %(qvalue)f''' %\
        (locals())
    elif "ttest" in infile:
        pvalue = PARAMS["ttest_pvalue_cutoff"]
        statement_fg = '''SELECT m.gene_id FROM
        %(tablename)s AS e,
        probeset2transcript AS m
        WHERE e.cluster_id = m.cluster_id AND
        e.pvalue < %(pvalue)f''' %\
        (locals())

    statement_bg = '''SELECT m.gene_id FROM
    %(tablename)s AS e,
    probeset2transcript AS m
    WHERE e.cluster_id = m.cluster_id''' %\
    (locals())
    
    executeGO( outfile,
           statement_fg,
           statement_bg,
           go_file )
    
############################################################
############################################################
############################################################
##
############################################################
def importGO( infile, outfile, suffix ):
    '''import GO results into a table.'''

    x = "_expdiff.%s" % suffix 
    assert infile.endswith( x )
    track, method, control = getExpressionMatch(
        infile[:-len(x)] + ".expdiff" )

    if track == control: return
    
    tablename = "%(track)s_vs_%(control)s_%(method)s_%(suffix)s" % locals()

    indir = infile + ".dir"

    if not os.path.exists( indir ):
        P.touch( outfile )
        return

    statement = '''
    python %(toolsdir)s/cat_tables.py %(indir)s/*.overall |\
    csv2db.py %(csv2db_options)s \
              --allow-empty \
              --index=category \
              --index=goid \
              --table=%(tablename)s \
    > %(outfile)s
    '''
    P.run( **dict( locals().items() + PARAMS.items() ) )

    
             

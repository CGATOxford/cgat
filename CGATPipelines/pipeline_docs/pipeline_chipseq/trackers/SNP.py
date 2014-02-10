from ChipseqReport import *

# for trackers_derived_sets and trackers_master
if not os.path.exists("conf.py"):
    raise IOError( "could not find conf.py" )

execfile( "conf.py" )

class SNPSlicer( ChipseqReport.DefaultTracker ):
    def getSlices( self, subset = None):
        return self.getValues( "SELECT DISTINCT snp FROM snps_of_interest" )

class SNPTable( TrackerSQL ):

    mPattern = "snps_of_interest"
    
    def getTracks( self, subset = None ):
        return ["all"]

    def __call__(self, track, slice = None ):
        columns = ("snp", "contig", "pos", "reference", "alleles", "url")
        data = self.get( """
        SELECT snp, contig, pos, reference, alleles, url
        FROM snps_of_interest
        """)
        
        return odict( zip( (columns),
                           zip(*data ) ) )

########################################################
########################################################
########################################################        
class SNPPileup( SNPSlicer ):
    '''create table with SNP variants from the unnormalized read files.'''
    mPattern = "_(?<!norm_)snpcoverage$"
    suffix = "snpcoverage"

    def __call__(self, track, slice = None ):

        columns = ("snp",
                   "nreads",
                   "contig",
                   "position",
                   "reference_base",
                   "consensus_base",
                   "consensus_quality",
                   "snp_quality",
                   "rms_mapping_quality",
                   "coverage",
                   "read_bases" )

        # can't use - gives rest parsing errors
        #                   "base_qualities" )
        
        fields = ",".join(columns)
        suffix = self.suffix
        statement = '''
        SELECT %(fields)s FROM %(track)s_%(suffix)s WHERE snp = '%(slice)s'
        ''' % locals()
        data = self.get(statement)
        
        data = [ x[:-2] + ('``%s``' % x[-1],) for x in data ]

        return odict( zip( (columns),
                           zip(*data) ) )

########################################################
########################################################
########################################################        
class SNPPileupNormalized( SNPPileup ):
    '''create table with SNP variants from the normalized read files.'''
    mPattern = "_norm_snpcoverage$"
    suffix = "norm_snpcoverage"

########################################################
########################################################
########################################################        
class SNPIntervals( SNPSlicer ):
    '''return intervals that SNPs overlap with.'''

    def __call__(self, track, slice = None ):

        columns = (
            "s.contig",
            "i.start",
            "i.end",
            "s.snp",
            "s.pos",
            "i.interval_id",
            "i.peakval",
            "m.motif",
            "m.nmatches",
            "m.evalue",
            "m.diagram",
            )

        fields = ",".join( columns )
        statement = '''
        SELECT
        %(fields)s
        FROM snps_of_interest AS s,
             %(track)s_intervals AS i,
             %(track)s_mast AS m
        WHERE
             s.snp = '%(slice)s' AND 
             i.contig = s.contig AND
             m.id = i.interval_id AND 
             s.pos BETWEEN i.start AND i.end-1
             '''

        data = self.get( statement % locals() )


        # check for overlap and add link
        data = [ (ChipseqReport.linkToUCSC( x[0], x[1], x[2]),) + x for x in data ]

        return odict( zip( (["link"] + [x[2:] for x in columns]),
                           zip(*data) ) )

########################################################
########################################################
########################################################        
class SNPMotifs( SNPSlicer ):
    '''return intervals that SNPs overlap with that contain a motif that the
    SNP ovelaps with.'''
    motif_width = 15

    def __call__(self, track, slice = None ):

        columns = (
            "s.contig",
            "i.start",
            "i.end",
            "s.snp",
            "s.pos",
            "i.interval_id",
            "i.peakval",
            "m.motif",
            "m.nmatches",
            "m.evalue",
            "m.diagram",
            "m.positions",
            )

        fields = ",".join( columns )
        statement = '''
        SELECT
        %(fields)s
        FROM snps_of_interest AS s,
             %(track)s_intervals AS i,
             %(track)s_mast AS m
        WHERE
             s.snp = '%(slice)s' AND 
             i.contig = s.contig AND
             m.id = i.interval_id AND 
             s.pos BETWEEN i.start AND i.end-1
             '''

        data = self.get( statement % locals() )
        motif_width = self.motif_width
        found = []
        for d in data:
            snp_position = d[4]
            
            try:
                motif_positions = map(int, d[-1].split(","))
            except AttributeError:
                motif_positions = (d[-1],)
                
            for p in motif_positions:
                if p == None: continue
                if p <= snp_position < p + motif_width:
                    found.append( d[:-1] )
        
        data = [ (ChipseqReport.linkToUCSC( x[0], x[1], x[2]),) + x for x in found ]
        
        return odict( zip( (["link"] + [x[2:] for x in columns[:-1]]),
                           zip(*data) ) )
    

             

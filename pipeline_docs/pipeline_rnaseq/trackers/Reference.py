import os, sys, re, types, itertools, math

from RnaseqReport import *
from SphinxReport.ResultBlock import ResultBlock, ResultBlocks

##################################################################################
##################################################################################
##################################################################################
## Trackers that access reference statistics
##################################################################################

class ReferenceData(RnaseqTracker):
    """Base class f or Trackers accessing reference table."""
    pattern = "(.*)_transcript_counts$" 
    reference = "refcoding"

class TranscriptCoverage(ReferenceData):
    """Coverage of reference transcripts."""
    mXLabel = "overlap / %"
    def __call__(self, track, slice = None ):
        data = self.getValues( """SELECT coverage_pcovered FROM %(track)s_transcript_counts WHERE coverage_nval > 0""" )
        return odict( (("covered", data ) ,) )

class GeneCoverage(ReferenceData):
    '''Coverage of reference genes - max transcript coverage per gene.'''
    mXLabel = "number of transcripts"
    def __call__(self, track, slice = None ):
        data = self.getValues( """SELECT max(c.coverage_pcovered) FROM 
                                            %(track)s_transcript_counts as c,
                                            %(reference)s_transcript2gene as i
                                         WHERE c.coverage_nval > 0
                                   AND i.transcript_id = c.transcript_id 
                                   GROUP BY i.gene_id """ )
        return odict( (("covered", data ) ,) )

class CoverageVsLengthByReadDepth(ReferenceData):
    """plot the absolute coverage of a known gene versus its length.
    Dots are colored by read depth.
    """

    mXLabel = "log(length)"

    def __call__(self, track, slice = None):
        reference = self.reference
        statement = """SELECT AVG(exons_sum) AS ref_length,
                                MIN(c.coverage_pcovered) AS coverage, 
                                AVG(c.coverage_mean) AS read_depth
                        FROM %(track)s_transcript_counts AS c,
                                %(reference)s_transcript2gene as i
                        WHERE i.transcript_id = c.transcript_id AND 
                                c.coverage_nval > 0
                        GROUP BY i.gene_id"""

        data = [ (math.log(x[0]), x[1], math.log(x[2]) ) for x in self.get( statement % locals() ) ]
        r = odict( zip( ("log(length)", "log(coverage)", "log(read_depth)"), zip(*data)))
        return odict( zip( ("log(length)", "log(coverage)", "log(read_depth)"), zip(*data)))

##=================================================================
## Coverage
##=================================================================
class MeanVsMaxReadDepth( ReferenceData ):
    """maxmimum read depth versus mean read depth of :term:`reference` genes. 
    Dots are coloured by the log(length) of a :term:`reference` gene."""

    mXLabel = "mean read depth"
    mYLabel = "maximum read depth"

    def __call__(self, track, slice = None ):
        reference = self.reference
        statement = "SELECT coverage_mean, coverage_max, exons_sum FROM %(track)s_transcript_counts" % locals()
        data = [ (x[0], x[1], math.log( x[2]) ) for x in self.get( statement) if x[2] > 0 ]
        return odict( zip( ("mean coverage", "max coverage", "length" ), zip(*data) ) )

class MeanVsMedianReadDepth( ReferenceData ):
    """maxmimum read depth versus mean read depth of :term:`reference` genes. 
    Dots are coloured by the log(length) of a :term:`reference` gene."""

    mXLabel = "mean read depth"
    mYLabel = "median read depth"

    def __call__(self, track, slice = None ):
        reference = self.reference
        statement = "SELECT coverage_mean, coverage_median, exons_sum FROM %(track)s_transcript_counts" % locals()
        data = [ (x[0], x[1], math.log( x[2]) ) for x in self.get( statement) if x[2] > 0 ]
        return odict( zip( ("mean coverage", "median coverage", "length" ), zip(*data) ) )



##=================================================================
## Directionality
##=================================================================

class ReadDirectionality(RnaseqTracker):
    '''return antisense / sense direction of reads in introns/genes.

    +1 is added as pseudo-count.
    '''
    
    pattern = "(.*)_intron_counts$" 

    slices = ("intron", "gene")

    def __call__(self, track, slice = None ):
        data = self.getValues( \
                    """SELECT CAST( (antisense_unique_counts + 1) AS FLOAT) / (sense_unique_counts + 1)  
                      FROM %(track)s_%(slice)s_counts """ )
        return odict( ( ( "direction", data) ,) )


class IntronicExonicReadDepth(RnaseqTracker):
    '''return the maximum read depth in introns
    and exons of a gene.

    +1 is added as pseudo-count.
    '''
    pattern = "(.*)_intron_counts$" 

    slices = ( "anysense", "antisense", "sense")
    min_coverage = 0

    def __call__(self, track, slice = None ):
        data = self.getAll( \
            """SELECT e.coverage_%(slice)s_max + 1 AS exon, i.coverage_%(slice)s_max + 1 as intron
               FROM %(track)s_gene_counts as e, %(track)s_intron_counts as i 
               WHERE e.gene_id = i.gene_id
                     AND e.coverage_%(slice)s_max >= %(min_coverage)i""" )
        return data


##############################################################
##############################################################
##############################################################
class UTRReadDensityInPlace( Tracker ):
    tracks = [ x.asFile() for x in TRACKS ]
    slices = ( "raw", "scaled", "fit" )
    def __call__(self, track, slice = None ):
        edir = EXPORTDIR
        method = "utr_extension"

        blocks = ResultBlocks()

        filepath = "%(edir)s/%(method)s/%(track)s.readextension_%(region)s_%(direction)s.%(slice)s.png"

        block = \
'''
.. figure:: %(filename)s
   :height: 300 
'''
        # append spaces for file extension
        block = "\n".join( [ x + " " * 200 for x in block.split("\n") ] )

        for region, direction in itertools.product( ("downstream", "upstream"),
                                                    ("sense", "antisense", "anysense" )):
            
            filename = filepath % locals()

            if os.path.exists( filename ):
                blocks.append( ResultBlock( text = block % locals(),
                                            title = "%(track)s %(region)s %(direction)s" % locals() ) )
            # else:
            #     blocks.append( ResultBlock( "",
            #                                 title = "%(track)s %(region)s %(direction)s" % locals() ) )

        return odict( (("rst", "\n".join(Utils.layoutBlocks( blocks, layout="columns-3" ))),))


class UTRReadDensityTable( Tracker ):
    tracks = [ x.asFile() for x in TRACKS ]
    slices = ( "raw", "scaled", "fit" )
    def __call__(self, track, slice = None ):
        edir = EXPORTDIR
        method = "utr_extension"

        filepath = "%(edir)s/%(method)s/%(track)s.readextension_%(region)s_%(direction)s.%(slice)s.png"

        toc_text = []
        link_text = []

        for region, direction in itertools.product( ("downstream", "upstream"),
                                                    ("sense", "antisense", "anysense" )):
        
            filename = filepath % locals()
            if not os.path.exists(filename): continue

            linktext = "%(track)s-%(region)s-%(direction)s-%(slice)s" % locals()

            toc_text.append( "* %(linktext)s_" % locals()) 
            link_text.append( ".. _%(linktext)s: %(filename)s" % locals() )

        toc_text = "\n".join(toc_text)
        link_text =  "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

        return odict( (("text", rst_text),) )

class UTRExtension( RnaseqTracker ):
    
    pattern = "(.*)_extension_counts_utr$" 
    slices = ( "5utr", "3utr" )

    def __call__(self, track, slice = None ):
        return self.getAll( '''SELECT old_%(slice)s_length AS known,
                                      new_%(slice)s_length AS predicted
                                      FROM %(track)s_extension_counts_utr
                                      WHERE new_%(slice)s_length IS NOT NULL''' )
    

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
PipelineMapping.py - Utility functions for mapping
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Mapping reads is a common task in pipelines. Different pipelines
combine different sources of input (:term:`fastq` files, :term:`sra` files)
of different data (single end, paired end) with different mapping
algorithms (bowtie, tophat, stampy). This module provides utility 
functions to abstract some of these variations.

The pipeline does not know what kind of data it gets (a :term:`sra` archive
might contain single end or paired end data or both).

A pipeline might get several input data (:term:`fastq` and :term:`sra`
formatted files at the same time).

The module currently is able to deal with:

   * tophat mapping against genome
   * bowtie mapping against transcriptome

It implements:
   * .sra: paired-end and single-end
   * .fastq: paired-end and single-end
   * .csfasta: colour-space, single-end

Code
----

'''

import os, sys, shutil, glob
import Pipeline as P
import IOTools

class Mapper( object ):
    '''map reads.

    preprocesses the input data, calls mapper and post-process the output data.
    
    All in a single statement to be send to the cluster.
    '''
    
    datatype = "fastq"

    def __init__(self):
        pass

    def preprocess( self, infiles, outfile ):
        '''build preprocessing statement

        returns statement and fastq files to map.

        Mapping qualities are changed to solexa format.
        '''

        assert len(infiles) > 0, "no input files for mapping"

        tmpdir_fastq = P.getTempDir()
        #print outfile

        # create temporary directory again for nodes
        statement = [ "mkdir -p %s" % tmpdir_fastq ]
        fastqfiles = []
        for infile in infiles:

            if infile.endswith( ".export.txt.gz"):
                # single end illumina export
                track = P.snip( os.path.basename( infile ), ".export.txt.gz" )
                statement.append( """gunzip < %(infile)s 
                     | awk '$11 != "QC" || $10 ~ /(\d+):(\d+):(\d+)/ \
                        { if ($1 != "") 
                             { readname=sprintf( "%%%%s_%%%%s:%%%%s:%%%%s:%%%%s:%%%%s", $1,$2,$3,$4,$5,$6);}
                        else { readname=sprintf( "%%%%s:%%%%s:%%%%s:%%%%s:%%%%s", $1,$3,$4,$5,$6); }
                       printf("@%%%%s\\n%%%%s\\n+\\n%%%%s\\n",readname,$9,$10);}'
                     > %(tmpdir_fastq)s/%(track)s.fastq""" % locals() )
                fastqfiles.append( ("%s/%s.fastq" % (tmpdir_fastq, track),) )

            elif infile.endswith( ".sra"):
                track = P.snip( infile, ".sra" )

                # sneak preview to determine if paired end or single end
                outdir = P.getTempDir()
                P.execute( "fastq-dump -X 1000 --outdir %(outdir)s %(infile)s" % locals() )
                f = glob.glob( os.path.join( outdir, "*.fastq" ) )
                if len(f) == 3:
                    f = glob.glob( os.path.join( outdir, "*_[12].fastq" ) )
                shutil.rmtree( outdir )
                fastqfiles.append( [ "%s/%s" % (tmpdir_fastq, os.path.basename( x )) for x in sorted(f) ] )
                statement.append( "fastq-dump --outdir %(tmpdir_fastq)s %(infile)s" % locals() )
                
            elif infile.endswith( ".fastq.gz" ):
                track = P.snip( os.path.basename( infile ), ".fastq.gz" )
                statement.append(  """gunzip < %(infile)s 
                                      | python %%(scriptsdir)s/fastq2fastq.py --change-format=sanger --guess-format=phred64 --log=%(outfile)s.log
                                      > %(tmpdir_fastq)s/%(track)s.fastq""" % locals() )
                fastqfiles.append( ("%s/%s.fastq" % (tmpdir_fastq, track),) )

            elif infile.endswith( ".csfasta.gz" ):
                track = P.snip( os.path.basename( infile ), ".csfasta.gz" )
                quality = P.snip( infile, ".csfasta.gz" ) + ".qual.gz"
                if not os.path.exists( quality ):
                    raise ValueRerror( "no quality file for %s" % infile )
                statement.append(  """gunzip < %(infile)s 
                                      > %(tmpdir_fastq)s/%(track)s.csfasta""" % locals() )
                statement.append(  """gunzip < %(quality)s 
                                      > %(tmpdir_fastq)s/%(track)s.qual""" % locals() )
                fastqfiles.append( ("%s/%s.csfasta" % (tmpdir_fastq, track),
                                    "%s/%s.qual" % (tmpdir_fastq, track) ) )
                self.datatype = "solid"

            elif infile.endswith( ".fastq.1.gz" ):

                track = P.snip( os.path.basename( infile ), ".gz" )
                track = P.snip( infile, ".fastq.1.gz" )
                infile2 = "%s.fastq.2.gz" % track
                if not os.path.exists( infile2 ):
                    raise ValueError("can not find paired ended file '%s' for '%s'" % (infile2, infile))
                statement.append( """gunzip < %(infile)s 
                                     | python %%(scriptsdir)s/fastq2fastq.py --change-format=sanger --guess-format=phred64 --log=%(outfile)s.log
                                     > %(tmpdir_fastq)s/%(track)s.1.fastq;
                                     gunzip < %(infile2)s 
                                     | python %%(scriptsdir)s/fastq2fastq.py --change-format=sanger --guess-format=phred64 --log=%(outfile)s.log
                                     > %(tmpdir_fastq)s/%(track)s.2.fastq
                                 """ % locals() )
                fastqfiles.append( ("%s/%s.1.fastq" % (tmpdir_fastq, track),
                                    "%s/%s.2.fastq" % (tmpdir_fastq, track) ) )
            else:
                raise NotImplementedError( "unknown file format %s" % infile )

        
        self.tmpdir_fastq = tmpdir_fastq

        assert len(fastqfiles) > 0, "no fastq files for mapping"

        return "; ".join( statement) + ";", fastqfiles

    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.
        '''
        return ""
    
    def postprocess( self, infiles, outfile ):
        '''collect output data and postprocess.'''
        return ""

    def cleanup( self, outfile ):
        '''clean up.'''
        statement = '''rm -rf %s;''' % (self.tmpdir_fastq)
        
        return statement

    def build( self, infiles, outfile ):
        '''run mapper.'''

        cmd_preprocess, mapfiles = self.preprocess( infiles, outfile )
        cmd_mapper = self.mapper( mapfiles, outfile )
        cmd_postprocess = self.postprocess( infiles, outfile )
        cmd_clean = self.cleanup( outfile )
        
        assert cmd_preprocess.strip().endswith(";")
        assert cmd_mapper.strip().endswith(";")
        if cmd_postprocess:
           assert cmd_postprocess.strip().endswith(";")
        if cmd_clean:
           assert cmd_clean.strip().endswith(";")

        statement = " checkpoint; ".join( (cmd_preprocess, 
                                           cmd_mapper,
                                           cmd_postprocess,
                                           cmd_clean ) )

        return statement


class fastqFilter( Mapper ):
    
    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.'''

        num_files = [ len( x ) for x in infiles ]
        
        if max(num_files) != min(num_files):
            raise ValueError("mixing single and paired-ended data not possible." )

        nfiles = max(num_files)

        tmpdir_fastq = self.tmpdir_fastq
        
        statement = []
        for f in infiles:
            if nfiles == 1:
                track = P.snip( os.path.basename( f[0] ), ".fastq" )
                statement.append('''fastx_collapser -v -i %(f[0])s -o %(tmpdir_fastq)s/%(track)s.rmdup.fastq > filtered_fastq/%(track)s_filter_stats.txt; ''' % locals() )
                statement.append('''fastq_quality_filter -q 20 -p 95 -i %(tmpdir_fastq)s/%(track)s.rmdup.fastq -o %(tmpdir_fastq)s/%(track)s.rmdup.qf.fastq >> filtered_fastq/%(track)s_filter_stats.txt; ''' % locals() )
                statement.append('''fastx_artifacts_filter -v -z -i %(tmpdir_fastq)s/%(track)s.rmdup.qf.fastq -o filtered_fastq/%(track)s_filt.fastq >> filtered_fastq/%(track)s_filter_stats.txt; ''' % locals() )

            elif nfiles == 2:
                track = P.snip( os.path.basename( f[0] ), ".1.fastq" )
                read1 = f[0]
                read2 = f[1]
                statement.append('''python /ifs/apps/bio/galaxy-dist/tools/fastq/fastq_paired_end_joiner.py %(read1)s illumina %(read2)s illumina %(tmpdir_fastq)s/%(track)s.merged.fastq; ''' % locals() )
                statement.append('''fastx_collapser -v -i %(tmpdir_fastq)s/%(track)s.merged.fastq -o %(tmpdir_fastq)s/%(track)s.merged.rmdup.fastq > filtered_fastq/%(track)s_filter_stats.txt; ''' % locals() )
                statement.append('''fastq_quality_filter -q 20 -p 95 -i %(tmpdir_fastq)s/%(track)s.merged.rmdup.fastq -o %(tmpdir_fastq)s/%(track)s.merged.rmdup.qf.fastq >> filtered_fastq/%(track)s_filter_stats.txt; ''' % locals() )
                statement.append('''fastx_artifacts_filter -v -z -i %(tmpdir_fastq)s/%(track)s.merged.rmdup.qf.fastq -o %(tmpdir_fastq)s/%(track)s.merged.rmdup.qf.af.fastq >> filtered_fastq/%(track)s_filter_stats.txt; ''' % locals() )
                statement.append('''python /ifs/apps/bio/galaxy-dist/tools/fastq/fastq_paired_end_splitter.py %(tmpdir_fastq)s/%(track)s.merged.rmdup.qf.af.fastq illumina %(tmpdir_fastq)s/%(track)s_filt.1.fastq %(tmpdir_fastq)s/%(track)s_filt.2.fastq ; ''' % locals() )
                statement.append('''gzip -c %(tmpdir_fastq)s/%(track)s_filt.1.fastq > filtered_fastq/%(track)s_filt.fastq.1.gz;''' % locals() )
                statement.append('''gzip -c %(tmpdir_fastq)s/%(track)s_filt.2.fastq > filtered_fastq/%(track)s_filt.fastq.2.gz;''' % locals() )
            else:
                raise ValueError( "unexpected number read files to map: %i " % nfiles )
        return " ".join( statement )



class fastqc( Mapper ):
    
    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.
        
        The output is created in exportdir
        '''
        
        statement = []
        for f in infiles:
            for i, x in enumerate(f):
                track = P.snip( os.path.basename( x ), ".fastq" )
                statement.append( '''fastqc --outdir=%%(exportdir)s/fastqc %(x)s >& %(outfile)s;''' % locals() )
        return " ".join( statement )


class bwa( Mapper ):
    
    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.'''

        num_files = [ len( x ) for x in infiles ]
        
        if max(num_files) != min(num_files):
            raise ValueError("mixing single and paired-ended data not possible." )
        
        nfiles = max(num_files)
        
        tmpdir_bwa = os.path.join( self.tmpdir_fastq + "bwa" )
        statement = [ "mkdir -p %s;" % tmpdir_bwa ]
        tmpdir_fastq = self.tmpdir_fastq

        if nfiles == 1:
            infiles = ",".join( [ x[0] for x in infiles ] )
            track = P.snip( os.path.basename( infiles ), ".fastq" )
            statement.append('''
            bwa aln %%(bwa_aln_options)s %%(bwa_index_dir)s/%%(genome)s %(infiles)s > %(tmpdir_bwa)s/%(track)s.sai 2>%(outfile)s.bwa.stderr; 
            bwa samse %%(bwa_index_dir)s/%%(genome)s %(tmpdir_bwa)s/%(track)s.sai %(infiles)s > %(tmpdir_bwa)s/%(track)s.sam 2>%(outfile)s.bwa.stderr;
            ''' % locals() )

        elif nfiles == 2:
            infiles1 = ",".join( [ x[0] for x in infiles ] )
            infiles2 = ",".join( [ x[1] for x in infiles ] )
            track = P.snip( os.path.basename( infiles1 ), ".1.fastq" )
            track1 = P.snip( os.path.basename( infiles1 ), ".fastq" )
            track2 = P.snip( os.path.basename( infiles2 ), ".fastq" )

            statement.append('''
            bwa aln %%(bwa_aln_options)s %%(bwa_index_dir)s/%%(genome)s %(infiles1)s > %(tmpdir_bwa)s/%(track1)s.sai 2>>%(outfile)s.bwa.stderr; checkpoint;
            bwa aln %%(bwa_aln_options)s %%(bwa_index_dir)s/%%(genome)s %(infiles2)s > %(tmpdir_bwa)s/%(track2)s.sai 2>>%(outfile)s.bwa.stderr; checkpoint;
            bwa sampe %%(bwa_sampe_options)s %%(bwa_index_dir)s/%%(genome)s %(tmpdir_bwa)s/%(track1)s.sai %(tmpdir_bwa)s/%(track2)s.sai %(infiles1)s %(infiles2)s > %(tmpdir_bwa)s/%(track)s.sam 2>>%(outfile)s.bwa.stderr;
            ''' % locals() )
        else:
            raise ValueError( "unexpected number read files to map: %i " % nfiles )

        self.tmpdir_bwa = tmpdir_bwa

        return " ".join( statement )
    
    def postprocess( self, infiles, outfile ):
        '''collect output data and postprocess.'''
        
        track = P.snip( os.path.basename(outfile), ".bam" )
        tmpdir_bwa = self.tmpdir_bwa

        statement = '''
            samtools view -buS %(tmpdir_bwa)s/%(track)s.sam | samtools sort -o - - | samtools rmdup - %(outfile)s 2>>%(outfile)s.bwa.stderr; 
            samtools index %(outfile)s;''' % locals()

        return statement


class Tophat( Mapper ):
    
    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.
        '''

        num_files = [ len( x ) for x in infiles ]
        
        if max(num_files) != min(num_files):
            raise ValueError("mixing single and paired-ended data not possible." )
        
        nfiles = max(num_files)
        
        tmpdir_tophat = os.path.join( self.tmpdir_fastq + "tophat" )
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        data_options = []
        if self.datatype == "solid":
            data_options.append( "--quals --integer-quals --color" )
            index_file = "%(bowtie_index_dir)s/%(genome)s_cs"
        else:
            index_file = "%(bowtie_index_dir)s/%(genome)s"

        data_options = " ".join( data_options )

        if nfiles == 1:
            infiles = ",".join( [ x[0] for x in infiles ] )
            statement = '''
            tophat --output-dir %(tmpdir_tophat)s
                   --num-threads %%(tophat_threads)i
                   %(data_options)s
                   %%(tophat_options)s
                   %(index_file)s
                   %(infiles)s 
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()

        elif nfiles == 2:
            # this section works both for paired-ended fastq files
            # and single-end color space mapping (separate quality file)
            infiles1 = ",".join( [ x[0] for x in infiles ] )
            infiles2 = ",".join( [ x[1] for x in infiles ] )

            statement = '''
            tophat --output-dir %(tmpdir_tophat)s
                   --mate-inner-dist %%(tophat_mate_inner_dist)i
                   --num-threads %%(tophat_threads)i
                   %(data_options)s
                   %%(tophat_options)s
                   %(index_file)s
                   %(infiles1)s %(infiles2)s 
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()
        else:
            raise ValueError( "unexpected number reads to map: %i " % nfiles )

        self.tmpdir_tophat = tmpdir_tophat

        return statement
    
    def postprocess( self, infiles, outfile ):
        '''collect output data and postprocess.'''
        
        track = P.snip( outfile, ".bam" )
        tmpdir_tophat = self.tmpdir_tophat

        statement = '''
            mv %(tmpdir_tophat)s/accepted_hits.bam %(outfile)s; 
            gzip < %(tmpdir_tophat)s/junctions.bed > %(track)s.junctions.bed.gz; 
            mv %(tmpdir_tophat)s/logs %(outfile)s.logs;
            samtools index %(outfile)s;
            ''' % locals()

        return statement

class Bowtie( Mapper ):
    '''map with bowtie against genome.'''

    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.

        .. note:: a filter on bamfiles removes any /1 and /2
            markers from reads. The reason is that these
            markers are removed for paired-end data, but
            not for single-end data and will cause
            problems using read name lookup.
        '''

        num_files = [ len( x ) for x in infiles ]
        
        if max(num_files) != min(num_files):
            raise ValueError("mixing single and paired-ended data not possible." )

        nfiles = max(num_files)

        # transpose files
        infiles = zip( *infiles )

        # add options specific to data type
        data_options = []
        if self.datatype == "solid":
            data_options.append( "--quals --integer-quals --color" )
#            data_options.append( "-f -C" )
            if nfiles == 2:
                # single end,
                # second file will colors (unpaired data)
                data_options.append( "--quals %s" % ",".join( infiles[1] ) )
                nfiles -= 1
            elif nfiles == 4:
                data_options.append( "-Q1 %s -Q2 %s" % (",".join(infiles[2], infiles[3])) )
                nfiles -= 2
            else:
                raise ValueError( "unexpected number of files" )
            index_file = "%(bowtie_index_dir)s/%(genome)s_cs"
        else:
            index_file = "%(bowtie_index_dir)s/%(genome)s"

        data_options = " ".join( data_options )

        data_options = " ".join( data_options )
        tmpdir_fastq = self.tmpdir_fastq

        if nfiles == 1:
            infiles = ",".join( infiles[0])
            statement = '''
                bowtie --quiet --sam
                       --threads %%(bowtie_threads)i
                       %(data_options)s
                       %%(bowtie_options)s
                       %(index_file)s
                       %(infiles)s
                       2>%(outfile)s.log
               | awk -v OFS="\\t" '{sub(/\/[12]$/,"",$1);print}'
               | samtools import %%(reffile)s - %(tmpdir_fastq)s/out.bam 1>&2 2>> %(outfile)s.log;
            ''' % locals()

        elif nfiles == 2:
            infiles1 = ",".join( infiles[0] )
            infiles2 = ",".join( infiles[1] )

            statement = '''
                bowtie --quiet --sam
                       --threads %%(bowtie_threads)i
                       %(data_options)s
                       %%(bowtie_options)s
                       %(index_prefix)s
                       -1 %(infiles1)s -2 %(infiles2)s 
                       2>%(outfile)s.log
               | samtools import %%(reffile)s - %(tmpdir_fastq)s/out.bam 1>&2 2>> %(outfile)s.log;
            ''' % locals()            
        else:
            raise ValueError( "unexpected number reads to map: %i " % nfiles )

        return statement

    def postprocess( self, infiles, outfile ):
        '''collect output data and postprocess.'''
        
        track = P.snip( outfile, ".bam" )
        tmpdir_fastq = self.tmpdir_fastq

        statement = '''
             samtools sort %(tmpdir_fastq)s/out.bam %(track)s;
             samtools index %(outfile)s;
             ''' % locals()

        return statement


class BowtieTranscripts( Mapper ):
    '''map with bowtie against transcripts.'''

    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.

        .. note:: a filter on bamfiles removes any /1 and /2
            markers from reads. The reason is that these
            markers are removed for paired-end data, but
            not for single-end data and will cause
            problems using read name lookup.
        '''

        num_files = [ len( x ) for x in infiles ]
        
        if max(num_files) != min(num_files):
            raise ValueError("mixing single and paired-ended data not possible." )

        nfiles = max(num_files)

        # transpose files
        infiles = zip( *infiles )

        # add options specific to data type
        data_options = []
        if self.datatype == "solid":
            data_options.append( "-f -C" )
            if nfiles == 2:
                # single end,
                # second file will colors (unpaired data)
                data_options.append( "--quals %s" % ",".join( infiles[1] ) )
                nfiles -= 1
            elif nfiles == 4:
                data_options.append( "-Q1 %s -Q2 %s" % (",".join(infiles[2], infiles[3])) )
                nfiles -= 2
            else:
                raise ValueError( "unexpected number of files" )
            index_prefix = "%(prefix)s_cs"
        else:
            index_prefix = "%(prefix)s"

        data_options = " ".join( data_options )
        tmpdir_fastq = self.tmpdir_fastq

        if nfiles == 1:
            infiles = ",".join( infiles[0])
            statement = '''
                bowtie --quiet --sam
                       --threads %%(bowtie_threads)i
                       %(data_options)s
                       %%(bowtie_options)s
                       %(index_prefix)s
                       %(infiles)s
                       2>%(outfile)s.log
               | awk -v OFS="\\t" '{sub(/\/[12]$/,"",$1);print}'
               | samtools import %%(reffile)s - %(tmpdir_fastq)s/out.bam 1>&2 2>> %(outfile)s.log;
            ''' % locals()

        elif nfiles == 2:
            infiles1 = ",".join( infiles[0] )
            infiles2 = ",".join( infiles[1] )

            statement = '''
                bowtie --quiet --sam
                       --threads %%(bowtie_threads)i
                       %(data_options)s
                       %%(bowtie_options)s
                       %(index_prefix)s
                       -1 %(infiles1)s -2 %(infiles2)s 
                       2>%(outfile)s.log
               | samtools import %%(reffile)s - %(tmpdir_fastq)s/out.bam 1>&2 2>> %(outfile)s.log;
            ''' % locals()            
        else:
            raise ValueError( "unexpected number reads to map: %i " % nfiles )

        return statement

    def postprocess( self, infiles, outfile ):
        '''collect output data and postprocess.'''
        
        track = P.snip( outfile, ".bam" )
        tmpdir_fastq = self.tmpdir_fastq

        statement = '''
             samtools sort %(tmpdir_fastq)s/out.bam %(track)s;
             samtools index %(outfile)s;
             ''' % locals()

        return statement

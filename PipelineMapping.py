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
PipelineMapping.py - Utility functions for mapping short reads
==============================================================

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
   * bowtie mapping against transcriptome, genome and junctions
   * bwa against genome
   * stampy against genome

It implements:
   * .sra: paired-end and single-end
   * .fastq: paired-end and single-end
   * .csfasta: colour-space, single-end

Code
----

'''

import os, sys, shutil, glob, collections, re
import Pipeline as P
import Experiment as E
import IOTools
import Fastq
import pysam

SequenceInformation = collections.namedtuple( "SequenceInformation",
                                              """paired_end 
                                                 filename_first
                                                 filename_second
                                                 readlength_first
                                                 readlength_second
                                                 is_colour""" )

def getReadLengthFromFastq( filename ):
    '''return readlength from a fasta/fastq file.

    Only the first read is inspected. If there are
    different read lengths in the file, though luck.

    '''

    with IOTools.openFile( filename ) as infile:
        record = iterate( infile ).next()
        readlength = len(record.seq)
        return readlength

def getReadLengthFromBamfile( filename ):
    '''return readlength from a bam file.

    Only the first read is inspected. If there are
    different read lengths in the file, though luck.
    '''

    samfile = pysam.Samfile( filename, "rb" )
    record = samfile.fetch().next()
    readlength = record.rlen
    samfile.close()
    return readlength

def getSequencingInformation( track ):
    '''glean sequencing information from *track*.'''

    colour = False
    if os.path.exists( "%s.fastq.gz" % track): 
        first_pair = "%s.fastq.gz" % track
        second_pair = None
    elif os.path.exists( "%s.fastq.1.gz" % track): 
        first_pair = "%s.fastq.1.gz" % track
        second_pair = "%s.fastq.2.gz" % track
    elif os.path.exists( "%s.csfasta.gz" % track):
        first_pair = "%s.csfasta.gz" % track
        second_pair = None
        colour = True

    second_length = None
    if second_pair:
        if not os.path.exists( second_pair ):
            raise IOError("could not find second pair %s for %s" % (second_pair, first_pair) )
        second_length = getReadLength( second_pair )
    
    return SequenceInformation._make( (second_pair != None,
                                       first_pair, second_pair,
                                       getReadLength( first_pair ),
                                       second_length,
                                       colour ) )

class Mapper( object ):
    '''map reads.

    preprocesses the input data, calls mapper and post-process the output data.
    
    All in a single statement to be send to the cluster.
    '''
    
    datatype = "fastq"

    # set to True if you want to preserve colour space files.
    # By default, they are converted to fastq.
    preserve_colourspace = False

    # compress fastq files with gzip
    compress = False

    # convert to sanger quality scores
    convert = False

    def __init__(self, executable = None):
        if executable:
            self.executable = executable

    def quoteFile( self, filename ):
        '''add uncompression for compressed files.
        and programs that expect uncompressed files.

        .. warn::
            This will only work if the downstream programs read the
            file only once.
        '''
        if filename.endswith( ".gz" ) and not self.compress:
            return "<( gunzip < %s )" % filename 
        else:
            return filename

    def preprocess( self, infiles, outfile ):
        '''build preprocessing statement

        Build a command line statement that extracts/converts 
        various input formats to fastq formatted files.

        Mapping qualities are changed to solexa format.

        returns the statement and the fastq files to map.
        '''

        assert len(infiles) > 0, "no input files for mapping"

        tmpdir_fastq = P.getTempDir()

        # create temporary directory again for nodes
        statement = [ "mkdir -p %s" % tmpdir_fastq ]
        fastqfiles = []

        # get track by extension of outfile
        track = os.path.splitext( os.path.basename( outfile ) )[0]

        if self.compress:
            compress_cmd = "| gzip"
            extension = ".gz"
        else:
            compress_cmd = ""
            extension = ""

        for infile in infiles:

            if infile.endswith( ".export.txt.gz"):
                # single end illumina export
                statement.append( """gunzip < %(infile)s 
                     | awk '$11 != "QC" || $10 ~ /(\d+):(\d+):(\d+)/ \
                        { if ($1 != "") 
                             { readname=sprintf( "%%%%s_%%%%s:%%%%s:%%%%s:%%%%s:%%%%s", $1,$2,$3,$4,$5,$6);}
                        else { readname=sprintf( "%%%%s:%%%%s:%%%%s:%%%%s:%%%%s", $1,$3,$4,$5,$6); }
                       printf("@%%%%s\\n%%%%s\\n+\\n%%%%s\\n",readname,$9,$10);}'
                     %(compress_cmd)s
                     > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" % locals() )
                fastqfiles.append( ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension ),) )
            elif infile.endswith( ".fa.gz" ):
                statement.append( '''gunzip < %(infile)s > %(tmpdir_fastq)s/%(track)s.fa''' % locals() )
                fastqfiles.append( ("%s/%s.fa" % (tmpdir_fastq, track ),) )
                self.datatype = "fasta"
                
            elif infile.endswith( ".sra"):
                # sneak preview to determine if paired end or single end
                outdir = P.getTempDir()
                P.execute( "fastq-dump --gzip -X 1000 --outdir %(outdir)s %(infile)s" % locals() )
                f = glob.glob( os.path.join( outdir, "*.fastq.gz" ) )
                if len(f) == 3:
                    f = glob.glob( os.path.join( outdir, "*_[12].fastq.gz" ) )
                E.info("sra file contains the following files: %s" % f )
                shutil.rmtree( outdir )
                fastqfiles.append( [ "%s/%s" % (tmpdir_fastq, os.path.basename( x )) for x in sorted(f) ] )
                statement.append( "fastq-dump --gzip --outdir %(tmpdir_fastq)s %(infile)s" % locals() )
                
            elif infile.endswith( ".fastq.gz" ):
                format = Fastq.guessFormat( IOTools.openFile( infile, "r"), raises = False)
                if 'sanger' not in format and self.convert:
                    statement.append(  """gunzip < %(infile)s 
                                      | python %%(scriptsdir)s/fastq2fastq.py --change-format=sanger --guess-format=phred64 --log=%(outfile)s.log
                                      %(compress_cmd)s
                                      > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" % locals() )
                    fastqfiles.append( ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),) )
                else:
                    E.debug( "%s: assuming quality score format %s" % (infile, format ) ) 
                    fastqfiles.append( (infile, ) )

            elif infile.endswith( ".csfasta.gz" ):
                # single end SOLiD data
                if self.preserve_colourspace:
                    quality = P.snip( infile, ".csfasta.gz" ) + ".qual.gz"
                    if not os.path.exists( quality ):
                        raise ValueError( "no quality file for %s" % infile )
                    statement.append(  """gunzip < %(infile)s 
                                          > %(tmpdir_fastq)s/%(track)s.csfasta%(extension)s""" % locals() )
                    statement.append(  """gunzip < %(quality)s 
                                          > %(tmpdir_fastq)s/%(track)s.qual%(extension)s""" % locals() )
                    fastqfiles.append( ("%s/%s.csfasta%s" % (tmpdir_fastq, track, extension ),
                                        "%s/%s.qual%s" % (tmpdir_fastq, track, extension) ) )
                    self.datatype = "solid"
                else:
                    quality = P.snip( infile, ".csfasta.gz" ) + ".qual.gz"

                    statement.append( """solid2fastq <(gunzip < %(infile)s) <(gunzip < %(quality)s)
                                      %(compress_cmd)s
                                      > %(tmpdir_fastq)s/%(track)s.fastq%(extension)""" % locals() )
                    fastqfiles.append( ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),) )

            elif infile.endswith( ".csfasta.F3.gz" ):
                # paired end SOLiD data
                if self.preserve_colourspace:
                    bn = P.snip( infile, ".csfasta.F3.gz" )
                    # order is important - mirrors tophat reads followed by quals
                    f = []
                    for suffix in ("csfasta.F3", "csfasta.F5", "qual.F3", "qual.F5" ):
                        fn = "%(bn)s.%(suffix)s" % locals()
                        if not os.path.exists( fn + ".gz"): raise ValueError( "expected file %s.gz missing" % fn )
                        statement.append( """gunzip < %(fn)s.gz
                                          %(compress_cmd)s
                                          > %(tmpdir_fastq)s/%(track)s.%(suffix)s%(extension)s""" % locals() )
                        f.append( "%(tmpdir_fastq)s/%(track)s.%(suffix)s%(extension)s" % locals() )
                    fastqfiles.append( f )
                    self.datatype = "solid"
                else:
                    quality = P.snip( infile, ".csfasta.gz" ) + ".qual.gz"

                    statement.append( """solid2fastq <(gunzip < %(infile)s) <(gunzip < %(quality)s)
                                      %(compress_cmd)s
                                      > %(tmpdir_fastq)s/%(track)s.fastq%(extension)s""" % locals() )
                    fastqfiles.append( ("%s/%s.fastq%s" % (tmpdir_fastq, track, extension),) )
                

            elif infile.endswith( ".fastq.1.gz" ):

                bn = P.snip( infile, ".fastq.1.gz" )
                infile2 = "%s.fastq.2.gz" % bn
                if not os.path.exists( infile2 ):
                    raise ValueError("can not find paired ended file '%s' for '%s'" % (infile2, infile))
                
                format = Fastq.guessFormat( IOTools.openFile( infile ), raises = False )
                if 'sanger' not in format:
                    statement.append( """gunzip < %(infile)s 
                                     | python %%(scriptsdir)s/fastq2fastq.py --change-format=sanger --guess-format=phred64 --log=%(outfile)s.log
                                     %(compress_cmd)s
                                     > %(tmpdir_fastq)s/%(track)s.1.fastq%(extension)s;
                                     gunzip < %(infile2)s 
                                     | python %%(scriptsdir)s/fastq2fastq.py --change-format=sanger --guess-format=phred64 --log=%(outfile)s.log
                                     %(compress_cmd)s
                                     > %(tmpdir_fastq)s/%(track)s.2.fastq%(extension)s
                                 """ % locals() )
                    fastqfiles.append( ("%s/%s.1.fastq%s" % (tmpdir_fastq, track, extension),
                                        "%s/%s.2.fastq%s" % (tmpdir_fastq, track, extension) ) )

                else:
                    E.debug( "%s: assuming quality score format %s" % (infile, format ) ) 
                    fastqfiles.append( (infile, 
                                        infile2, ) )
                    
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
        #statement = '''rm -rf %s;''' % (self.tmpdir_fastq)
        statement = ""

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
    
class FastQc( Mapper ):
    '''run fastqc to test read quality.'''

    compress = True

    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.
        
        The output is created in exportdir
        '''
        
        statement = []
        for f in infiles:
            for i, x in enumerate(f):
                track = os.path.basename(  re.sub(".fastq.*", "", x) )
                statement.append( '''fastqc --outdir=%%(exportdir)s/fastqc %(x)s >& %(outfile)s;''' % locals() )
        return " ".join( statement )

class Counter( Mapper ):
    '''count number of reads in fastq files.'''

    compress = True
    
    def mapper( self, infiles, outfile ):
        '''count number of reads by counting number of lines
        in fastq files.
        '''
        
        statement = []
        for f in infiles:
            x = " ".join( f )
            statement.append( '''zcat %(x)s | awk '{n+=1;} END {printf("nreads\\t%%%%i\\n",n/4);}' >& %(outfile)s;''' % locals() )
        return " ".join( statement )

class BWA( Mapper ):
    '''run bwa to map reads against genome.

    * colour space not implemented
    '''

    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.'''

        num_files = [ len( x ) for x in infiles ]
        
        if max(num_files) != min(num_files):
            raise ValueError("mixing single and paired-ended data not possible." )
        
        nfiles = max(num_files)
        
        tmpdir = os.path.join( self.tmpdir_fastq + "bwa" )
        statement = [ "mkdir -p %s;" % tmpdir ]
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        # note: not fully implemented
        data_options = ["%(bwa_align_options)s"]
        if self.datatype == "solid":
            data_options.append( "-c" )
            index_prefix = "%(bwa_index_dir)s/%(genome)s_cs"
        elif self.datatype == "fasta":
            index_prefix = "%(bwa_index_dir)s/%(genome)s"
        else:
            index_prefix = "%(bwa_index_dir)s/%(genome)s"

        data_options = " ".join( data_options )

        tmpdir_fastq = self.tmpdir_fastq

        track = P.snip( os.path.basename( outfile ), ".bam" )

        if nfiles == 1:
            infiles = ",".join( [ self.quoteFile(x[0]) for x in infiles ] )
            
            statement.append('''
            bwa aln %%(bwa_aln_options)s -t %%(bwa_threads)i %(index_prefix)s %(infiles)s 
            > %(tmpdir)s/%(track)s.sai 2>>%(outfile)s.bwa.log; 
            bwa samse %%(bwa_index_dir)s/%%(genome)s %(tmpdir)s/%(track)s.sai %(infiles)s 
            > %(tmpdir)s/%(track)s.sam 2>>%(outfile)s.bwa.log;
            ''' % locals() )

        elif nfiles == 2:
            track1 = track + ".1"
            track2 = track + ".2"
            infiles1 = ",".join( [ self.quoteFile(x[0]) for x in infiles ] )
            infiles2 = ",".join( [ self.quoteFile(x[1]) for x in infiles ] )

            statement.append('''
            bwa aln %%(bwa_aln_options)s -t %%(bwa_threads)i %(index_prefix)s %(infiles1)s 
            > %(tmpdir)s/%(track1)s.sai 2>>%(outfile)s.bwa.log; checkpoint;
            bwa aln %%(bwa_aln_options)s -t %%(bwa_threads)i %(index_prefix)s %(infiles2)s 
            > %(tmpdir)s/%(track2)s.sai 2>>%(outfile)s.bwa.log; checkpoint;
            bwa sampe %%(bwa_sampe_options)s %(index_prefix)s
                      %(tmpdir)s/%(track1)s.sai 
                      %(tmpdir)s/%(track2)s.sai 
                      %(infiles1)s 
                      %(infiles2)s 
            > %(tmpdir)s/%(track)s.sam 2>>%(outfile)s.bwa.log;
            ''' % locals() )
        else:
            raise ValueError( "unexpected number read files to map: %i " % nfiles )

        self.tmpdir = tmpdir

        return " ".join( statement )
    
    def postprocess( self, infiles, outfile ):
        '''collect output data and postprocess.'''
        
        track = P.snip( os.path.basename(outfile), ".bam" )
        outf = P.snip( outfile, ".bam" )
        tmpdir = self.tmpdir
        
        statement = '''
            samtools view -buS %(tmpdir)s/%(track)s.sam | samtools sort - %(outf)s 2>>%(outfile)s.bwa.log; 
            samtools index %(outfile)s;''' % locals()

        return statement

class Stampy( BWA ):
    '''map reads against genome using STAMPY.
    '''

    # compress fastq files with gzip
    compress = True

    executable = "stampy.py"

    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.'''
        
        num_files = [ len( x ) for x in infiles ]
        
        if max(num_files) != min(num_files):
            raise ValueError("mixing single and paired-ended data not possible." )
        
        nfiles = max(num_files)
        executable = self.executable

        tmpdir = os.path.join( self.tmpdir_fastq + "stampy" )
        statement = [ "mkdir -p %s;" % tmpdir ]
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        # note: not fully implemented
        data_options = ["%(stampy_align_options)s"]
        if self.datatype == "solid":
            data_options.append( "-c" )
            bwa_index_prefix = "%(bwa_index_dir)s/%(genome)s_cs"
        elif self.datatype == "fasta":
            bwa_index_prefix = "%(bwa_index_dir)s/%(genome)s"
        else:
            bwa_index_prefix = "%(bwa_index_dir)s/%(genome)s"

        track = P.snip( os.path.basename( outfile ), ".bam" )

        if nfiles == 1:
            infiles = ",".join( [ self.quoteFile(x[0]) for x in infiles ] )

            statement.append('''
            %(executable)s -v 3 -g %%(stampy_index_dir)s/%%(genome)s -h %%(stampy_index_dir)s/%%(genome)s 
                      --bwaoptions="-q10 %(bwa_index_prefix)s" 
                      -M %(infiles)s 
            > %(tmpdir)s/%(track)s.sam 2>%(outfile)s.log;
            ''' % locals() )

        elif nfiles == 2:
            track1 = track + ".1"
            track2 = track + ".2"
            infiles1 = ",".join( [ self.quoteFile(x[0]) for x in infiles ] )
            infiles2 = ",".join( [ self.quoteFile(x[1]) for x in infiles ] )

            statement.append('''
            %(executable)s -v 3 -g %%(stampy_index_dir)s/%%(genome)s -h %%(stampy_index_dir)s/%%(genome)s 
                      --bwaoptions="-q10 %(bwa_index_prefix)s" 
                      -M %(infiles1)s %(infiles2)s 
            > %(tmpdir)s/%(track)s.sam 2>%(outfile)s.log;
            ''' % locals() )
        else:
            raise ValueError( "unexpected number of read files to map: %i " % nfiles )

        self.tmpdir = tmpdir

        return " ".join( statement )
    
class Tophat( Mapper ):

    # tophat can map colour space files directly
    preserve_colourspace = True

    # newer versions of tophat can work of compressed files
    compress = True
    
    executable = "tophat"

    def __init__(self, *args, **kwargs ):
        Mapper.__init__(self, *args, **kwargs)

    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.
        '''

        executable = self.executable

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
            index_prefix = "%(bowtie_index_dir)s/%(genome)s_cs"
        else:
            index_prefix = "%(bowtie_index_dir)s/%(genome)s"

        data_options = " ".join( data_options )

        if nfiles == 1:
            infiles = ",".join( [ x[0] for x in infiles ] )
            statement = '''
            %(executable)s --output-dir %(tmpdir_tophat)s
                   --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %(index_prefix)s
                   %(infiles)s 
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()

        elif nfiles == 2:
            # this section works both for paired-ended fastq files
            # and single-end color space mapping (separate quality file)
            infiles1 = ",".join( [ x[0] for x in infiles ] )
            infiles2 = ",".join( [ x[1] for x in infiles ] )

            statement = '''
            %(executable)s --output-dir %(tmpdir_tophat)s
                   --mate-inner-dist %%(tophat_mate_inner_dist)i
                    --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %(index_prefix)s
                   %(infiles1)s %(infiles2)s 
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()
        elif nfiles == 4:
            # this section works both for paired-ended fastq files
            # in color space mapping (separate quality file)
            # reads1 reads2 qual1 qual2
            infiles1 = ",".join( [ x[0] for x in infiles ] )
            infiles2 = ",".join( [ x[1] for x in infiles ] )
            infiles3 = ",".join( [ x[2] for x in infiles ] )
            infiles4 = ",".join( [ x[3] for x in infiles ] )

            statement = '''%(executable)s
                   --output-dir %(tmpdir_tophat)s
                   --mate-inner-dist %%(tophat_mate_inner_dist)i
                   --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %(index_prefix)s
                   %(infiles1)s %(infiles2)s 
                   %(infiles3)s %(infiles4)s 
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

class TopHat_fusion( Mapper ):
    
    # tophat can map colour space files directly
    preserve_colourspace = True
    
    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.
        '''

        num_files = [ len( x ) for x in infiles ]
        
        if max(num_files) != min(num_files):
            raise ValueError("mixing single and paired-ended data not possible." )
        
        nfiles = max(num_files)
        
        tmpdir_tophat = os.path.join(  self.tmpdir_fastq  + "tophat" )
        tmpdir_fastq = self.tmpdir_fastq

        # add options specific to data type
        data_options = []
        if self.datatype == "solid":
            data_options.append( "--quals --integer-quals --color" )
            index_prefix = "%(bowtie_index_dir)s/%(genome)s_cs"
        else:
            index_prefix = "%(bowtie_index_dir)s/%(genome)s"

        data_options = " ".join( data_options )

        if nfiles == 1:
            infiles = ",".join( [ x[0] for x in infiles ] )
            statement = '''
            module load tophatfusion;
            tophat-fusion --output-dir %(tmpdir_tophat)s
                   --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %%(tophatfusion_options)s
                   %(index_prefix)s
                   %(infiles)s 
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()

        elif nfiles == 2:
            # this section works both for paired-ended fastq files
            # and single-end color space mapping (separate quality file)
            infiles1 = ",".join( [ x[0] for x in infiles ] )
            infiles2 = ",".join( [ x[1] for x in infiles ] )

            statement = '''
            module load tophatfusion;
            tophat-fusion --output-dir %(tmpdir_tophat)s
                   --mate-inner-dist %%(tophat_mate_inner_dist)i
                    --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                  %%(tophat_options)s
                  %%(tophatfusion_options)s
                   %(index_prefix)s
                   %(infiles1)s %(infiles2)s 
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()
        elif nfiles == 4:
            # this section works both for paired-ended fastq files
            # in color space mapping (separate quality file)
            # reads1 reads2 qual1 qual2
            infiles1 = ",".join( [ x[0] for x in infiles ] )
            infiles2 = ",".join( [ x[1] for x in infiles ] )
            infiles3 = ",".join( [ x[2] for x in infiles ] )
            infiles4 = ",".join( [ x[3] for x in infiles ] )

            statement = '''
            module load bio/tophatfusion;
            tophat-fusion --output-dir %(tmpdir_tophat)s
                   --mate-inner-dist %%(tophat_mate_inner_dist)i
                   --num-threads %%(tophat_threads)i
                   --library-type %%(tophat_library_type)s
                   %(data_options)s
                   %%(tophat_options)s
                   %%(tophatfusion_options)s
                   %(index_prefix)s
                   %(infiles1)s %(infiles2)s 
                   %(infiles3)s %(infiles4)s 
                   >> %(outfile)s.log 2>&1 ;
            ''' % locals()


        else:
            raise ValueError( "unexpected number reads to map: %i " % nfiles )

        self.tmpdir_tophat = tmpdir_tophat

        return statement
    
    def postprocess( self, infiles, outfile ):
        '''collect output data and postprocess.'''
        
        track = P.snip( outfile, "/accepted_hits.sam" )
        tmpdir_tophat = self.tmpdir_tophat
        
        if not os.path.exists('%s' % track):
            os.mkdir('%s' % track)

        statement = '''
            mv -f %(tmpdir_tophat)s/* %(track)s/;  
            ''' % locals()

        return statement

class Bowtie( Mapper ):
    '''map with bowtie against genome.'''

    # bowtie can map colour space files directly
    preserve_colourspace = True

    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.

        .. note:: a filter on bamfiles removes any /1 and /2
            markers from reads. The reason is that these
            markers are removed for paired-end data, but
            not for single-end data and will cause
            problems using read name lookup.
        '''

        executable = self.executable

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
            if nfiles == 2:
                # single end,
                # second file will colors (unpaired data)
                data_options.append( "--quals %s" % ",".join( [ self.quoteFile(x) for x in infiles[1] ] ) ) 
                nfiles -= 1
            elif nfiles == 4:
                raise NotImplementeError()
                data_options.append( "-Q1 %s -Q2 %s" % (",".join(infiles[2], infiles[3])) )
                nfiles -= 2
            else:
                raise ValueError( "unexpected number of files" )
            index_prefix = "%(bowtie_index_dir)s/%(genome)s_cs"
        elif self.datatype == "fasta":
            data_options.append( "-f" )
            index_prefix = "%(bowtie_index_dir)s/%(genome)s"
        else:
            index_prefix = "%(bowtie_index_dir)s/%(genome)s"

        data_options = " ".join( data_options )

        tmpdir_fastq = self.tmpdir_fastq

        if nfiles == 1:
            infiles = ",".join( [ self.quoteFile(x) for x in infiles[0] ] )
            statement = '''
                %(executable)s --quiet --sam
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
            infiles1 = ",".join( [ self.quoteFile( x ) for x in infiles[0] ] )
            infiles2 = ",".join( [ self.quoteFile( x ) for x in infiles[1] ] )

            statement = '''
                %(executable)s --quiet --sam
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

    # bowtie can map colour space files directly
    preserve_colourspace = True

    compress = True

    executable = "bowtie"

    def mapper( self, infiles, outfile ):
        '''build mapping statement on infiles.

        .. note:: a filter on bamfiles removes any /1 and /2
            markers from reads. The reason is that these
            markers are removed for paired-end data, but
            not for single-end data and will cause
            problems using read name lookup.
        '''
        executable = self.executable

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
                data_options.append( "-Q1 <( zcat %s ) -Q2 <( zcat %s)" % (",".join(infiles[2], infiles[3])) )
                nfiles -= 2
            else:
                raise ValueError( "unexpected number of files" )
            index_prefix = "%(prefix)s_cs"
        else:
            index_prefix = "%(prefix)s"

        data_options = " ".join( data_options )
        tmpdir_fastq = self.tmpdir_fastq

        if nfiles == 1:
            infiles = ",".join( ["<(zcat %s)" % x for x in infiles[0] ] )
            statement = '''
                %(executable)s --quiet --sam
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
            infiles1 = ",".join( ["<(zcat %s)" % x for x in infiles[0] ] )
            infiles2 = ",".join( ["<(zcat %s)" % x for x in infiles[1] ] )

            statement = '''
                %(executable)s --quiet --sam
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

class BowtieJunctions( BowtieTranscripts ):
    '''map with bowtie against junctions.

    In post-processing, reads are mapped from junction coordinates
    to genomic coordinates.
    '''

    def postprocess( self, infiles, outfile ):
        '''collect output data and postprocess.'''
        
        track = P.snip( outfile, ".bam" )
        tmpdir_fastq = self.tmpdir_fastq

        statement = '''
             cat %(tmpdir_fastq)s/out.bam
             | python %%(scriptsdir)s/bam2bam.py --set-nh --log=%(outfile)s.log
             | python %%(scriptsdir)s/rnaseq_junction_bam2bam.py --contig-sizes=%%(contigsfile)s --log=%(outfile)s.log
             | samtools sort - %(track)s;
             checkpoint;
             samtools index %(outfile)s;
             ''' % locals()

        return statement



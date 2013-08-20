#!/usr/bin/env python

import sys
import os
import re
import glob
import shutil
import tempfile
import argparse
import subprocess

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()
    
def timenow():
    """return current time as a string
    """
    return time.strftime('%d/%m/%Y %H:%M:%S', time.localtime(time.time()))


class CGATBase():
    """
    simple base class with some utilities for Picard
    adapted and merged with Kelly Vincent's code april 2011 Ross
    lots of changes...
    """
    
    def __init__(self, opts=None,arg0=None):
        """ common stuff needed at init for a picard tool
        """
        pass

    def baseName(self,name=None):
        return os.path.splitext(os.path.basename(name))[0]

    def setLogging(self,logfname="picard_wrapper.log"):
        """setup a logger
        """
        logging.basicConfig(level=logging.INFO,
                    filename=logfname,
                    filemode='a')


    def readLarge(self,fname=None):
        """ read a potentially huge file.
        """
        try:
            # get stderr, allowing for case where it's very large
            tmp = open( fname, 'rb' )
            s = ''
            buffsize = 1048576
            try:
                while True:
                    more = tmp.read( buffsize )
                    if len(more) > 0:
                        s += more
                    else:
                        break
            except OverflowError:
                pass
            tmp.close()
        except Exception, e:
            stop_err( 'Read Large Exception : %s' % str( e ) )   
        return s
    
    def runStatement(self,cl=None,output_dir=None):
        """ construct and run a command line
        we have galaxy's temp path as opt.temp_dir so don't really need isolation
        sometimes stdout is needed as the output - ugly hacks to deal with potentially vast artifacts
        """
        assert cl <> None, 'PicardBase runCL needs a command line as cl'
        process = subprocess.Popen(cl, shell=True )
        rval = process.wait()
    
    def runPic(self, jar, cl):
        """
        cl should be everything after the jar file name in the command
        """
        runme = ['java -Xmx%s' % self.opts.maxjheap]
        runme.append(" -Djava.io.tmpdir='%s' " % self.opts.tmpdir)
        runme.append('-jar %s' % jar)
        runme += cl
        s,stdouts,rval = self.runCL(cl=runme, output_dir=self.opts.outdir)
        return stdouts,rval

    def samToBam(self,infile=None,outdir=None):
        """
        use samtools view to convert sam to bam
        """
        fd,tempbam = tempfile.mkstemp(dir=outdir,suffix='rgutilsTemp.bam')
        cl = ['samtools view -h -b -S -o ',tempbam,infile]
        tlog,stdouts,rval = self.runCL(cl,outdir)
        return tlog,tempbam,rval

    def sortSam(self, infile=None,outfile=None,outdir=None):
        """
        """
        print '## sortSam got infile=%s,outfile=%s,outdir=%s' % (infile,outfile,outdir)
        cl = ['samtools sort',infile,outfile]
        tlog,stdouts,rval = self.runCL(cl,outdir)
        return tlog

    def cleanup(self):
        for fname in self.delme:
            try:
                os.unlink(fname)
            except:
                pass

def __main__():

    # use argparse to ignore unknown options
    parser = argparse.ArgumentParser()

    parser.add_argument( "--version", action="version", version="%(prog)s")
    parser.add_argument( "--wrapper-command", dest="command", type=str)
    parser.add_argument( "--wrapper-bam-file", dest="bam_file", type=str)
    parser.add_argument( "--wrapper-bam-option", dest="bam_option", type=str)
    parser.add_argument( "--wrapper-bai-file", dest="bai_file", type=str)
    parser.add_argument( "--wrapper-dry-run", dest="dry_run", action="store_true")
    parser.add_argument( "--wrapper-html-dir", dest="html_dir", type=str)
    parser.add_argument( "--wrapper-html-file", dest="html_file", type=str)
        
    options, unknown = parser.parse_known_args( )

    cgat = CGATBase( options )

    option_map = []

    if options.bai_file or options.bam_file:
        if not (options.bai_file and options.bam_file):
            raise ValueError("wrapper called with bam or bai file, but not both")
        
        if not options.bam_option:
            options.bam_option = "bam-file"

        tmp_fd, tmp_name = tempfile.mkstemp()
        tmp_bam_name = '%s.bam' % tmp_name
        tmp_bai_name = '%s.bai' % tmp_bam_name
        os.symlink( options.bam_file, tmp_bam_name )
        os.symlink( options.bai_file, tmp_bai_name )
        if options.bam_option.startswith("--"):
            # long option
            option_map.append( "%s=%s" % (options.bam_option, tmp_bam_name) )
        else:
            # short option
            option_map.append( "%s %s" % (options.bam_option, tmp_bam_name) )

    if options.html_dir:
        os.mkdir( options.html_dir )
        option_map.append( "%s=%s/%%s" % ( "--output-filename-pattern", options.html_dir ) )

    statement = "python " + " ".join( [options.command] + unknown + option_map )

    if options.dry_run:
        sys.stdout.write( statement + "\n" )
        return

    else:
        cgat.runStatement( statement )
    
    if options.bai_file:
        os.unlink( tmp_bam_name )
        os.unlink( tmp_bai_name )
    
    if options.html_file:
        with open( options.html_file, "w" ) as outf:
            outf.write('<h1>%s - Output</h1>' % os.path.basename( options.wrapper_command))
            for fn in glob.glob( os.path.join( options.html_dir, "*.*" ) ):
                dirname, basename = os.path.split(fn)
                outf.write('''<li><a href="%s">%s</a></li>\n''' % (basename, basename))
                
if __name__=="__main__": __main__()


import sys, os, re, subprocess, optparse, tempfile, collections, itertools, random

from multiprocessing import Process
from threading import Thread as Process

import Experiment as E
import Stats

USAGE = """benchmark the relative speed of nodes on the cluster.

execute in a networked directory available on all nodes.
"""

CODE_CPU="""
import Experiment as E
import operator

E.Start()

def fib(n):
   if n == 0 or n == 1:
      return n
   else:
      return fib(n-1) + fib(n-2)

for i in range(%(iterations)s):
   a,b = (i, fib(i))

E.Stop()
"""

class Host: pass

from analyze_logfiles import LogFileData

def getSGEHosts():
    """get sun grid engine hosts."""
    try:
        x = subprocess.Popen( "qhost", 
                              shell=True,
                              stdout=subprocess.PIPE,
                              )

    except OSError, e:
        E.fail( "Execution of %s failed" % cmd)
        raise
    
    hosts = []
    for line in x.stdout:
        if line.startswith("HOSTNAME"): continue
        if line.startswith("global"): continue
        if line.startswith("-"): continue
        x = Host()
        x.mHost, x.mArch, x.mNCpu, x.mLoad, x.mMemtot, x.mMemuse, x.mSwapto, x.mSwapus = re.split( "\s+", line[:-1] )
        if x.mNCpu == "-" or x.mArch == "-": continue
        try: x.mNCpu = int(x.mNCpu )
        except ValueError: pass

        hosts.append( x )
    return hosts

def writeResults( outfile, results ):
    fields = ( "wall", "user", "sys", "cuser", "csys", "nchunks" )

    outfile.write( "host\t%s\n" % "\t".join( ["%s_%s" % (x, y) for x,y in itertools.product( fields, Stats.Summary().getHeaders() ) ] ) )

    hosts = results.keys()
    hosts.sort()

    for host in hosts:
        result = results[host]
        outfile.write( "%s" % host )
        for f in fields:
            d = [ y.__getitem__(f) for y in result ]
            outfile.write( "\t%s" % Stats.Summary( d ) )            
        outfile.write( "\n" )            

def main():

    parser = optparse.OptionParser( version = "%prog version: $Id: benchmark_cluster.py 2794 2009-09-16 15:27:49Z andreas $", usage = USAGE)

    parser.add_option( "-s", "--suite", dest="suite", type="choice",
                       choices=("CPU",),
                       help="suite to use [default=%default]." )

    parser.set_defaults(
        iterations_cpu = 36,
        iterations_sample = 5,
        saturation = 2,
        host_block = 10,
        suite = "CPU",
        cluster_prefix = """qrsh -q benchmark_jobs.q@%(host)s -cwd -v BASH_ENV=~/.bashrc -now n 'python %(script)s | grep -e "job started" -e "job finished"'""",
        )

    (options, args) = E.Start( parser )

    hosts = [ x for x in getSGEHosts() if x.mNCpu != "-" ]

    E.info( "found %i active hosts" % len(hosts) )

    script = os.path.abspath( os.path.join( "benchmarkscript" ) )

    if os.path.exists( script ):
        raise OSError( "file %s already exists - not overwritten" % script )

    outfile = open(script,"w")
    if options.suite == "CPU":
        iterations = options.iterations_cpu
        outfile.write( CODE_CPU % locals() )
    outfile.close()

    results = collections.defaultdict( list )

    def run( script, host, options ):
        """run python script on node in the cluster.
        """

        command = options.cluster_prefix % locals()

        try:
            child = subprocess.Popen( command, 
                                       shell=True,
                                       stdout=subprocess.PIPE,
                                       close_fds = True )
            (stdout, stderr) = child.communicate()
        except OSError, e:
            E.warn( "execution of %s failed" % cmd)
            raise

        data = LogFileData()
        for line in stdout.split("\n"): 
            data.add( line )
        results[host].append( data )

    # do on host-by-host basis, as otherwise 
    # we run out of ports. Several hosts are run
    # as one block. The order of hosts is randomized
    # at each step.
    idx = 0
    for n in range(0,options.iterations_sample):
        myhosts = list( hosts[:] )
        random.shuffle( myhosts )

        for h in range(0, len(myhosts), options.host_block) :
            for host in (myhosts[x] for x in range(h, min(len(myhosts), h+options.host_block))):
                processes = [] 
                for cpu in range(0, int( host.mNCpu * options.saturation) ):
                    p = Process( target = run, args = ( script, host.mHost, options ) )
                    processes.append( p )
                    p.start()
                E.info( "iteration %i: host %s: started %i threads" % (n, host.mHost, len(processes)))
            for p in processes: p.join()

            filename = "results_%i.table" % idx
            E.info( "iteration %i: output goes to %s" % (n, filename ))
            idx += 1
            outfile = open( filename, "w")
            writeResults( outfile, results )
            outfile.close()

    writeResults( options.stdout, results )

if __name__ == "__main__":
    sys.exit(main())

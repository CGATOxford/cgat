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
benchmark_cluster.py - script for testing CPU speed in SGE cluster
=====================================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Benchmark CPU speed of nodes in the cluster.

The script sends off a script to compute fibonacci numbers to
all nodes in the cluster. 

Testing is done in multiple iterations and on host-by-host basis.
Several hosts are run and the order of hosts is randomized at each
step.

More jobs than slots are send to each host.

The output is a collection of tab separated tables with benchmarking
data for each iteration separately and a combined analysis.

Usage
-----

Example::

   python benchmark_cluster.py 

Type::

   python benchmark_cluster.py --help

for command line help.

Documentation
-------------

Command line options
---------------------

'''
import sys
import os
import re
import subprocess
import optparse
import tempfile
import collections
import itertools
import random

from multiprocessing import Process
from threading import Thread as Process

import CGAT.Experiment as E
import CGAT.Stats as Stats
import CGAT.Logfile as Logfile

CODE_CPU="""
import sys
sys.path.append( "/ifs/devel/andreas/cgat" )
import CGAT.Experiment as E
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

CODE_SCRIPT="""
#! /bin/bash

python %(dir)s/benchmarkscript.py
"""

HOST = collections.namedtuple("HOST", "host arch ncpu load memtot memuse swapto swapus")
class Host: pass

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
        x = HOST( *re.split( "\s+", line[:-1] ) )
        if x.ncpu == "-" or x.arch == "-": continue

        try: x = x._replace( ncpu = int(x.ncpu ),  )
        except ValueError: pass

        assert x.memtot.endswith("G")
        x = x._replace( memtot = float(x.memtot[:-1] ) )

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

    parser = optparse.OptionParser( version = "%prog version: $Id: benchmark_cluster.py 2794 2009-09-16 15:27:49Z andreas $", usage = globals()["__doc__"])

    parser.add_option( "-s", "--suite", dest="suite", type="choice",
                       choices=("CPU", "hosts"),
                       help="suite to use [default=%default]." )

    parser.set_defaults(
        iterations_cpu = 36,
        iterations_sample = 5,
        saturation = 2,
        host_block = 10,
        suite = "CPU",
        cluster_prefix = """qrsh -q %(queue)s@%(host)s -cwd -now n 'python %(script)s | grep -e "job started" -e "job finished"'""",
        )

    (options, args) = E.Start( parser, add_cluster_options = True )

    hosts = [ x for x in getSGEHosts() if x.ncpu != "-" ]
    E.info( "found %i active hosts" % len(hosts) )

    curdir = os.path.abspath( os.getcwd() )

    outfile = options.stdout

    if options.suite == "hosts":
        outfile.write( "%s\n" % "\t".join(HOST._fields))
        for host in hosts:
            outfile.write( "%s\n" % "\t".join(map(str, host )))
    elif suite == "CPU":

        pythonscript = os.path.join( curdir, "benchmarkscript.py" ) 

        if os.path.exists( pythonscript ):
            raise OSError( "file %s already exists - not overwritten" % pythonscript )

        with open(pythonscript,"w") as outfile:
            if options.suite == "CPU":
                iterations = options.iterations_cpu
                outfile.write( CODE_CPU % locals() )

        results = collections.defaultdict( list )

        def run( script, host, options ):
            """run python script on node in the cluster.
            """

            queue = options.cluster_queue

            command = options.cluster_prefix % locals()
            E.debug( "started: %s" % command )

            try:
                child = subprocess.Popen( command, 
                                          shell=True,
                                          stdout=subprocess.PIPE,
                                          stdin=None,
                                          close_fds = True )
                (stdout, stderr) = child.communicate()
            except OSError, e:
                E.warn( "execution of %s failed" % cmd)
                raise

            E.debug( "finished: %s" % command )

            data = Logfile.LogFileData()
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
                    for cpu in range(0, int( host.ncpu * options.saturation) ):
                        p = Process( target = run, args = ( pythonscript, host.host, options ) )
                        processes.append( p )
                        p.start()
                    E.info( "iteration %i: host %s: started %i threads" % (n, host.host, len(processes)))
                for p in processes: p.join()

            filename = "results_%i.table" % idx
            E.info( "iteration %i: output goes to %s" % (n, filename ))
            idx += 1
            outfile = open( filename, "w")
            writeResults( outfile, results )
            outfile.close()

        writeResults( options.stdout, results )

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

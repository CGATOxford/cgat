#!/bin/env python

'''
ruffus_profile.py - analyze ruffus logfile
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script examines the ruffus logfile and collates
summary information.

**-t/--time** choice
   Report times either as ``milliseconds`` or ``seconds``.

.. note::

   All times are wall clock times.

'''

import sys, os, imp, cStringIO, re, types, glob, optparse, shutil, datetime, logging
import collections

import Experiment as E
import IOTools

class Counter(object):
    '''
    
    This class stores the calls per source as calls can be made by
    different sources in any order.
    '''

    def __init__(self):
        self._durations = collections.defaultdict( list )
        self._started = collections.defaultdict( int )
        self._calls = collections.defaultdict( int )

    def add( self, started, dt, source = None ):
        if self._started[source] != 0 and not started:
            self._durations[source].append( dt - self._started[source] )            
            self._started[source] = 0
        elif self._started[source] == 0 and started:
            self._calls[source] += 1
            self._started[source] = dt
        else:
            raise ValueError("""inconsistent time points for %s, has_started=%s, is_started=%s.
Possibly two pipelines have been running concurrently.
""" % (source,
       self._started[source], started))

    def reset( self, source = None):
        '''reset last event.'''
        self._started[source] = None
        self._calls[source] = 0

    def getDuration( self ): 
        if self._durations:
            x = None
            for source, durations in self._durations.iteritems():
                if x == None: x = durations[0] 
                for y in durations[1:]: x += y
            return x
        else:
            return datetime.timedelta()
        
    def getCalls( self ): 
        return sum( self._calls.values() )

    def getRunning( self ):
        '''get numbers of tasks unfinished or still running.'''
        return len( [x for x,y in self._started.iteritems() if y != 0 ] )

    duration = property(getDuration)
    calls = property(getCalls)
    running = property(getRunning)

def main( argv = sys.argv ):

    parser = optparse.OptionParser( version = "%prog version: $Id$", 
                                    usage = globals()["__doc__"])
    
    parser.add_option( "-l", "--logfile", dest="logfile", type="string",
                       help="name of logfile [default=%default]" )

    parser.add_option( "-t", "--time", dest="time", type="choice",
                       choices=("seconds", "milliseconds" ),
                       help="time to show [default=%default]" )

    parser.add_option( "-f", "--filter", dest="filter", type="choice",
                       choices=("unfinished", "running", "completed", "all" ),
                       help="apply filter to output [default=%default]" )

    parser.add_option( "-i", "--ignore-errors", dest="ignore_errors", action="store_true",
                       help="ignore errors [default=%default]" )

    parser.set_defaults( sections = [],
                         logfile = "pipeline.log",
                         filter = "all",
                         time = "seconds" )

    (options, args) = E.Start( parser, argv )

    rx = re.compile("^[0-9]+" )

    if options.sections:
        profile_sections = options.sections
    else:
        profile_sections = ("task", "job" )

    counts = {}
    for section in profile_sections:
        counts[section] = collections.defaultdict( Counter )

    rootpath = os.path.abspath(".")

    infile = IOTools.openFile( options.logfile )

    for line in infile:
        if not rx.match( line ): continue
        data = line[:-1].split()
        if len(data) < 5: continue
        date, time, level, source = data[:4]

        if re.search( "output generated by", line ): 
            E.info( "resetting counts at line=%s" % line[:-1] )
            for section in profile_sections:
                counts[section] = collections.defaultdict( Counter )
            continue

        if not re.match( "task\.", source ): continue

        dt = datetime.datetime.strptime( " ".join( (date, time) ), "%Y-%m-%d %H:%M:%S,%f")

        msg = "".join( data[4:] )
        
        started_task, completed_task, started_job, completed_job = None, None, None, None

        if re.search( "task.log_at_level.\d+Task=(\S+)", msg ):
            checked_task = re.search( "task.log_at_level.\d+Task=(\S+)", msg ).groups()[0]
        elif re.search( "Job=\[(\S+)->(\S+)\]Missingfile[s]*\[(\S+)\]", msg):
            started_infiles, started_job, missing = re.search( "Job=\[(\S+)->(\S+)\]Missingfile[s]*\[(\S+)\]", msg).groups()
        elif re.search( "Taskentersqueue=(\S+)", msg):
            started_task = re.search( "Taskentersqueue=(\S+)", msg).groups()[0]
        elif re.search( "Job=\[(\S+)->(\S+)\]completed", msg ):
            completed_infiles, completed_job = re.search( "Job=\[(\S+)->(\S+)\]completed", msg ).groups()
        elif re.search( "CompletedTask=(\S+)", msg ):
            completed_task = re.search( "CompletedTask=(\S+)", msg ).groups()[0]
        else:
            continue

        try:
            if started_task:
                counts["task"][started_task].add( True, dt, started_task )
            elif completed_task:
                counts["task"][completed_task].add( False, dt, completed_task )
            elif started_job:
                counts["job"][started_job].add( True, dt, started_job )
            elif completed_job:
                counts["job"][completed_job].add( False, dt, completed_job )
            else:
                raise ValueError( "unknown action")
        except ValueError, msg:
            if not options.ignore_errors:
                raise 

    if options.time == "milliseconds":
        f = lambda d: d.seconds + d.microseconds / 1000
    elif options.time == "seconds":
        f = lambda d: d.seconds + d.microseconds / 1000000

    for section in profile_sections:
        options.stdout.write( "\t".join( ("section", "object", "ncalls", "duration", "percall", "running") ) + "\n" )

        running = []
        for objct, c in counts[section].iteritems():
            
            # apply filters
            if options.filter in ("unfinished", "running") and c.running == 0: 
                continue

            d = f(c.duration)
            if c.calls > 0:
                percall = "%6.3f" %( d / float(c.calls))
            else:
                percall = "na"


            options.stdout.write( "\t".join( \
                    (map( str, \
                              (section, objct, 
                               c.calls,
                               d,
                               percall,
                               c.running,
                               )))) + "\n" )
            
            running.extend( [x for x,y in c._started.iteritems() if y != 0 ] )
            
        options.stdout.write( "#//\n\n" )

        if running:        
            options.stdout.write( "# running %ss\n" % section )
            options.stdout.write( "\n".join( map(str, running )) + "\n" )
            options.stdout.write( "#//\n\n" )

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

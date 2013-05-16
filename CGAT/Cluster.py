################################################################################
#   Gene prediction pipeline 
#
#   $Id: Cluster.py 2784 2009-09-10 11:41:14Z andreas $
#
#   Copyright (C) 2004 Andreas Heger
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
"""
Cluster.py - module for running a job in parallel on the cluster
================================================================
"""

import sys, os, subprocess

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

    def __str__(self):
        return str(self.message)

class ClusterError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

class Cluster:

    def __init__(self, options):
        
        try:
            self.mPriority = options.cluster_priority
        except AttributeError:
            self.mPriority = None

        try:
            self.mQueue = options.cluster_queue
        except AttributeError:
            self.mCluster = None

        try:
            self.mNumJobs = options.cluster_num_jobs
        except AttributeError:
            self.mNumJobs = None
            
    def addJob( self, statement ):
        
        self.mJobs.append( statement )

    def runJobs( self ):
        """run jobs.

        This uses forking which is not without problems:

        * too many open files.

        Thus, running in chunks.
        """
        
        self.mResults = []

        if self.mNumJobs:
            chunks = range( 0, len(self.mJobs), self.mNumJobs)
            chunks.append( len(self.mJobs) )
        else:
            chunks = [0, len(self.mJobs) ]

        for x in range(len(chunks)-1):
            
            pids = []
            
            for job in self.mJobs[chunks[x]:chunks[x+1]]:

                pid = os.fork()

                if not pid:
                    os.nice( 20 )
                    ## TODO: deal with return values
                    out, err = self.__runJob( job )
                    self.mResults.append( (out, err) )
                    break
                else:
                    pids.append(pid)
                
            for pid in pids:
                os.waitpid(pid,0)

        return self.mResults

    def runJob( self, job ):
        """submit job to queue.

        This subroutine should only exit if the job completed.
        """
        s = subprocess.Popen( self.getStatement(job),
                              shell = True,
                              stdin = subprocess.PIPE,
                              stdout = subprocess.PIPE,
                              stderr = subprocess.PIPE,
                              cwd = os.curdir,
                              close_fds = True)                              

        (out, err) = s.communicate()

        if s.returncode != 0:
            raise ClusterError ("Error while running job %s: %s\n%s\n" % (job, out, err))

        return out, err

    def getStatement( self, job):
        return job

    def setPriority( self, priority ):
        self.mPriority = priority

    def setQueue( self, queue):
        self.mQueue = queue
    
class ClusterSGE(Cluster):

    def __init__(self, *args):
        
        Cluster.__init__(self, *args)

        if not self.mQueue:
            self.mQueue = "medium_jobs.q"
        if not self.mPriority:
            self.mPriority = -10

        self.mExportedVariables = ( "PYTHONPATH", "LD_LIBRARY_PATH", "PATH", "BLASTMAT", "DIALIGN2_DIR", "LANG", "WISECONFIGDIR" )

        self.mVariables = []
        for x in self.mExportedVariables:
            try:
                self.mVariables.append( ("-v %s='%s'" % (x, os.environ[x]) ) )
            except KeyError:
                pass
        self.mVariables = " ".join( self.mVariables )

        self.mOptions ="-l arch=lx24-x86"
        self.mName = os.path.basename(sys.argv[0])
        if not self.mName: self.mName = "generic"

    def getStatement(self, job):

        return "qrsh -now n -cwd -q %s -p %i -N %s %s %s %s " % ( self.mQueue,
                                                                  self.mPriority,
                                                                  self.mName,
                                                                  self.mVariables,
                                                                  self.mOptions,
                                                                  job )
    
        

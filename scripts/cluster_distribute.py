################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $
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
cluster_distribute.py - distribute files to cluster nodes
=========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script distributes one or more files to the local
scratch directory on all nodes in the cluster.

The script uses ``rsync`` to only copy files that are newer.

Usage
-----

Example::

   python cluster_distribute.py -collection=blast /net/cpp-group/tools/polyphen-2.0.18/nrdb/uniref100.*.{pin,psd,psi,phr,psq}

The above command mirrors uniprot blast indexed databases
into the directory :file:`/scratch/blast` on the nodes.

Type::

   python cluster_distribute.py --help

for command line help.

.. note::

   The directory needs to be cleaned up as disk space on the
   nodes is limited.

.. todo::
   
   Currently files are copied to all nodes. This is potentially
   wasteful if jobs will only be executed on a few nodes.

Documentation
-------------

Code
----

'''

import os
import sys
import re
import optparse

import CGAT.Experiment as E

def getNodes( nodes = None):
    '''hack - allow ranges, ...'''
    if nodes == None or len(nodes) == 0:
        # return [ "node%02i" % x for x in range(1,15) + range(20,34) + range(48,62) ]
        return [ "cgat%03i" % x for x in range( 1,10) + range(101,117) ] + \
            [ "cgat150", "cgatsmp1", "andromeda", "gandalf", "saruman" ]
    return nodes

def main( argv = None ):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv: argv = sys.argv

    # setup command line parser
    parser = E.OptionParser( version = "%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $", 
                                    usage = globals()["__doc__"] )

    parser.add_option("-s", "--scratchdir", dest="scratchdir", type="string",
                      help="the scratch directory on the nodes [default=%default]."  )

    parser.add_option("-c", "--collection", dest="collection", type="string",
                      help="files will be put into collection. This is a directory that will be"
                      " created just below the scratch directory [default=%default]."  )

    parser.set_defaults(
        scratchdir = "/scratch",
        collection = "",
        nodes = [],
        )

    ## add common options (-h/--help, ...) and parse command line 
    (options, args) = E.Start( parser, argv = argv )

    if len(args) == 0:
        raise ValueError("please specify a collection of files/directories that should be mirrored." )

    targetdir = os.path.join( options.scratchdir, options.collection )

    nodes = getNodes( options.nodes )

    E.info( "copying to %s on nodes %s" % (targetdir, ",".join(nodes) ))

    ninput, noutput, nskipped = 0, 0, 0

    filenames = " ".join(args)

    for node in nodes:
        E.info( "copying to node %s" % node )
        ninput += 1
        statement = '''
               ssh %(node)s mkdir %(targetdir)s >& /dev/null;
               rsync --progress -az %(filenames)s %(node)s:%(targetdir)s
        ''' % locals()
        E.run( statement )
        noutput += 1

    E.info( "ninput=%i, noutput=%i, nskipped=%i" % (ninput, noutput,nskipped) )


    ## write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit( main( sys.argv) )

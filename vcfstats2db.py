#! /bin/env python
################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id$
#
#   Copyright (C) 2011 David Sims
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
vcfstats_sqlite.py - reformat output of vcf-stats for database loading
====================================

:Author: David Sims
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

create a csv separated file for loading into a database from 
output of vcf-stats utility in vcf-tools package.

Usage
-----

Example::

   python vcfstats_sqlite.py [files] > [outfile]

Type::

   python vcfstats_sqlite.py --help

for command line help.


Code
----

'''
import os, sys, string, re, time, optparse, tempfile, subprocess, types
import Experiment as E
import csv, CSV
import IOTools

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id: vcfstats_sqlite.py 0001 2011-04-13 davids $", usage = globals()["__doc__"])

    (options, args) = E.Start( parser )

    options.filenames = args

    if len(options.filenames) < 1:
        print USAGE, "no vcf-stats files specified/found."
        sys.exit(1)

    E.info( "Parsing %i file(s)" % len(options.filenames) )

    for filename in options.filenames:

        prefix = os.path.basename( filename )
        trackname = prefix.replace(".vcfstats","")

        if os.path.exists(filename):
            lines = [ x for x in IOTools.openFile(filename, "r").readlines() if x.find("..") and x.strip()]
        else:
            lines = []

        if len(lines) == 0:
            print USAGE, "Empty vcf-stats file found: $(filename)s"
            sys.exit(1)
        else:
            E.info( "File %i contains %i lines" % (len(options.filenames),len(lines)) )
            vcf_stats = dict(track=trackname)
            for line in lines:
                if line.find("..") > 0: 
                    fields = line.split("..")
                    key = fields[0].strip().replace(" ","_").replace("/","_").replace(">","-")
                    val = fields[1].strip()
                    if key != "1":
                         vcf_stats[key] = val 
            
            # Write header (for first file only)
            sep=""
            if filename == options.filenames[0]:
                for k in vcf_stats.iterkeys():
                    options.stdout.write("%s%s" % (sep,k))
                    sep="\t"
                options.stdout.write("\n")

            # Write data
            sep=""
            for k in vcf_stats.iterkeys():
                options.stdout.write("%s%s" % (sep,vcf_stats[k]))
                sep="\t"
            options.stdout.write("\n")

    E.Stop()
    sys.exit(0)


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
======================================================================

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
import os
import sys
import string
import re
import time
import optparse
import tempfile
import subprocess
import types
import CGAT.Experiment as E
import csv
import CGAT.CSV as CSV
import CGAT.IOTools as IOTools

if __name__ == '__main__':

    parser = E.OptionParser( version = "%prog version: $Id: vcfstats_sqlite.py 0001 2011-04-13 davids $", usage = globals()["__doc__"])

    (options, args) = E.Start( parser )

    options.filenames = args

    if len(options.filenames) < 1:
        options.stdout.write("# Error: no vcf-stats files specified/found.")
        sys.exit(1)

    E.info( "Parsing %i file(s)" % len(options.filenames) )

    # set up output files
    vcf_file = open('vcfstats.txt', 'w')
    indel_file = open('indelstats.txt', 'w')
    snp_file = open('snpstats.txt', 'w')
    shared_file = open('sharedstats.txt', 'w')

    for fileno,filename in enumerate(options.filenames):

        prefix = os.path.basename( filename )
        trackname = prefix.replace(".vcfstats","")

        if os.path.exists(filename):
            lines = [ x for x in IOTools.openFile(filename, "r").readlines()]
        else:
            lines = []

        if len(lines) == 0:
            options.stdout.write("# Error: empty vcf-stats file found: $(filename)s")
            sys.exit(1)
        else:
            E.info( "File %i contains %i lines" % (fileno,len(lines)) )
            vcf_stats = dict(track=trackname)
            snp_stats = dict(track=trackname)
            indel_stats = dict()
            shared_stats = dict() 
            all_vars = False
            indels = False
            snps = False
            shared = False 
            for i, line in enumerate(lines):
                line = line.strip()
                if line.find("'all'") > -1: 
                    all_vars = True
                    E.info("Found 'all'")
                    continue

                if all_vars:
                    if line.find("=>") > -1: 
                        fields = line.split("=>")
                        key = fields[0].strip().replace("'","").replace(">","_")
                        val = fields[1].strip().replace(",","")
                    else: 
                        key = "NA"
                        val = "NA"
                    if key=="indel" and val=="{":
                        indels = True
                        E.info("Found 'indels'")
                        continue
                    elif key=="snp" and val=="{":
                        snps = True
                        E.info("Found 'SNPs'")
                        continue
                    elif key=="shared" and val=="{":
                        shared = True
                        E.info("Found 'Shared'")
                        continue

                    if indels:
                        if line.find("}") > -1:
                           indels = False
                           E.info("Processed 'indels'")
                           continue
                        else:
                            indel_stats[key] = val
                    elif snps:
                        if line.find("}") > -1:
                            snps = False
                            E.info("Processed 'SNPs'")
                            continue
                        else:
                            snp_stats[key] = val
                    elif shared:
                        if line.find("}") > -1:
                            shared = False
                            E.info("Processed 'Shared'")
                            continue
                        else:
                            shared_stats[key] = val
                    elif key != "NA":
                        vcf_stats[key] = val 
            
            # Ensure all keys are present
            allkeys = ["nalt_1","nalt_2","nalt_3","nalt_4","nalt_5","track","count","snp_count","indel_count"]
            for k in allkeys:
                if vcf_stats.has_key(k):
                    continue
                else:
                    vcf_stats[k] = "0"

            # Write header (for first file only)
            if filename == options.filenames[0]:

                # Ensure keys are sorted
                srt = vcf_stats.keys()
                srt.sort()
                sep=""
                for k in srt:
                    vcf_file.write("%s%s" % (sep,k))
                    sep="\t"
                vcf_file.write("\n")

                indel_file.write("track\tindel_length\tindel_count\n")
                shared_file.write("track\tno_samples\tvar_count\n")

                sep=""
                for k in snp_stats.iterkeys():
                    snp_file.write("%s%s" % (sep,k))
                    sep="\t"
                snp_file.write("\n")

            # Write data
            sep=""
            srt = vcf_stats.keys()
            srt.sort()
            for k in srt:
                vcf_file.write("%s%s" % (sep,vcf_stats[k]))
                sep="\t"
            vcf_file.write("\n")


            # Check all indel lengths are covered
            r = range(-20,20,1)
            for i in r:
                if indel_stats.has_key(str(i)):
                    continue
                else:
                    indel_stats[i]="0"
            for k in indel_stats.iterkeys():
                indel_file.write("%s\t%s\t%s\n" % (trackname,k,indel_stats[k]))
 
            for k in shared_stats.iterkeys():
                shared_file.write("%s\t%s\t%s\n" % (trackname,k,shared_stats[k]))

            sep=""
            for k in snp_stats.iterkeys():
                snp_file.write("%s%s" % (sep,snp_stats[k]))
                sep="\t"
            snp_file.write("\n")

    # close files
    vcf_file.close()
    indel_file.close()
    snp_file.close()

    E.Stop()
    sys.exit(0)


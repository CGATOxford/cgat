##########################################################################
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
##########################################################################
'''
Sra.py - Methods for dealing with short read archive files
==========================================================

Utility functions for dealing with :term:`SRA` formatted files from
the Short Read Archive.

Requirements:
* fastq-dump >= 2.1.7

Code
----

'''
import os
import glob
import tempfile
import shutil
import CGAT.Experiment as E
import CGAT.Fastq as Fastq
import CGAT.IOTools as IOTools


def peek(sra, outdir=None):
    """return the full file names for all files which will be extracted

    Parameters
    ----------

    outdir : path
        perform extraction in outdir. If outdir is None, the extraction
        will take place in a temporary directory, which will be deleted
        afterwards.

    Returns
    -------
    files : list
        A list of fastq formatted files that are contained in the archive.
    format : string
        The quality score format in the :term:`fastq` formatted files.

    """
    
    if outdir is None:
        workdir = tempfile.mkdtemp()
    else:
        workdir = outdir

    # --split-files creates files called prefix_#.fastq.gz,
    # where # is the read number.
    # If file cotains paired end data:
    # output = prefix_1.fastq.gz, prefix_2.fastq.gz
    #    *special case: unpaired reads in a paired end --> prefix.fastq.gz
    #    *special case: if paired reads are stored in a single read,
    #                   fastq-dump will split. There might be a joining
    #                   sequence. The output would thus be:
    #                   prefix_1.fastq.gz, prefix_2.fastq.gz, prefix_3.fastq.gz
    #                   You want files 1 and 3.

    E.run("""fastq-dump --split-files --gzip -X 1000
                 --outdir %(workdir)s %(sra)s""" % locals())
    f = sorted(glob.glob(os.path.join(workdir, "*.fastq.gz")))
    ff = [os.path.basename(x) for x in f]

    if len(f) == 1:
        # sra file contains one read: output = prefix.fastq.gz
        pass

    elif len(f) == 2:
        # sra file contains read pairs:
        # output = prefix_1.fastq.gz, prefix_2.fastq.gz
        assert ff[0].endswith(
            "_1.fastq.gz") and ff[1].endswith("_2.fastq.gz")

    elif len(f) == 3:
        if ff[2].endswith("_3.fastq.gz"):
            f = glob.glob(os.path.join(workdir, "*_[13].fastq.gz"))
        else:
            f = glob.glob(os.path.join(workdir, "*_[13].fastq.gz"))

    # check format of fastqs in .sra
    fastq_format = Fastq.guessFormat(IOTools.openFile(f[0], "r"), raises=False)

    if outdir is None:
        shutil.rmtree(workdir)

    return f, fastq_format


def extract(sra, outdir):
    """return statement for extracting the SRA file in `outdir`."""

    statement = """fastq-dump --split-files --gzip --outdir
                 %(outdir)s %(sra)s""" % locals()

    return statement

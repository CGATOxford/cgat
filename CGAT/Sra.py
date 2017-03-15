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
import itertools
import CGAT.Experiment as E
import CGAT.Fastq as Fastq
import CGAT.IOTools as IOTools
from future.moves.urllib.request import urlopen


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
    fastq_datatype = Fastq.guessDataType(
        IOTools.openFile(f[0], "r"), raises=True)

    if outdir is None:
        shutil.rmtree(workdir)

    return f, fastq_format, fastq_datatype


def extract(sra, outdir, tool="fastq-dump"):
    """return statement for extracting the SRA file in `outdir`.
    possible tools are fastq-dump and abi-dump. Use abi-dump for colorspace"""

    if tool == "fastq-dump":
        tool += " --split-files"

    statement = """%(tool)s --gzip --outdir %(outdir)s %(sra)s""" % locals()

    return statement


def prefetch(sra):
    """Use prefetch from the SRA toolkit to download the local cache"""

    statement = """prefetch %(sra)s -a "$ASCP_BIN_PATH|$ASCP_KEY_PATH" """ % locals()
    return statement


def clean_cache(sra):
    """Remove the specified SRA file from the cache."""

    statement = """rm `srapath %(sra)s*`""" % locals()
    return statement

    
def fetch_ENA(dl_path, outdir, protocol="ascp"):
    """Fetch fastq from ENA given accession"""

    if protocol == "ascp":
        statement = """ascp -QT -l %%(aspera_bandwidth)s -i $ASCP_KEY_PATH
                       era-fasp@fasp.sra.ebi.ac.uk:/%(dl_path)s %(outdir)s """ % locals()
    elif protocol == "http":
        fn = os.path.basename(dl_path)
        outfile = os.path.join(outdir, fn)
        statement = "wget -O %(outfile)s ftp://ftp.sra.ebi.ac.uk:/%(dl_path)s" % locals()

    return statement


def fetch_ENA_files(accession):
    """Get the names of the files matching the ENA accession"""
    
    # Get the paths to the files
    url = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run&fields=fastq_ftp" \
          % accession

    try:
        paths = urlopen(url).readlines()[1:]
    except:
        E.debug("couldn't access %s" % url)
        raise

    paths = list(itertools.chain.from_iterable(
        [p.strip().split(";") for p in paths]))
    filenames = [os.path.basename(x) for x in paths]
    dl_paths = [os.path.join(*x.split("/")[1:]) for x in paths]
    return filenames, dl_paths


def fetch_TCGA_fastq(acc, filename, token=None, outdir="."):
    """Get Fastq file from TCGA repository. Because of the nature of the
    TCGA repository it assumes certain things:

            1) That data is paired-end fastq
            2) That the files end in _1.fastq or _2.fastq
            """

    statement = []

    if token:
        token = "-t %s" % token
    else:
        token = ""

    statement.append("""gdc-client download
                         %(token)s
                         --no-related-files
                         --no-annotations
                         -d %(outdir)s
                         %(acc)s""")

    tar_file = os.path.join(outdir, acc, filename)

    statement.append("tar -xvf %(tar_file)s -C %(outdir)s")

    out_name = os.path.join(outdir, acc)

    statement.append("gzip %(outdir)s/*.fastq")
    statement.append("rm -r %(out_name)s")
    statement.append("rename %(outdir)s/*_1.fastq.gz %(outdir)s/%(acc)s_1.fastq.gz %(outdir)s/*")
    statement.append("rename %(outdir)s/*_2.fastq.gz %(outdir)s/%(acc)s_2.fastq.gz %(outdir)s/*")
    statement = "; checkpoint;".join(statement)

    return statement % locals()


def fetch_TCGA_BAM(acc, token, outdir=".", filter_bed=None):
    """Get BAM file from TCGA repository based on UUID. Will return
    statement and path/filename of downloaded file. A bed file may be
    provided to filter to remove contigs not present in the
    reference genome"""

    statement = []

    if token:
        token = "-t %s" % token
    else:
        token = ""

    statement.append("""gdc-client download
                        %(token)s
                        -d %(outdir)s
                        %(acc)s """)

    from_prefix = os.path.join(outdir, acc, "*")
    to_prefix = os.path.join(outdir, acc, acc)

    if filter_bed:
        statement.append("""samtools view -hb
                                     %(from_prefix)s.bam
                                     -L %(filter_bed)s >
                                     %(to_prefix)s.bam;
                            checkpoint;
                            samtools index %(to_prefix)s.bam""")
    else:
        statement.append('''rename %(from_prefix)s.bam
                                   %(to_prefix)s.bam
                                   %(from_prefix)s''')
    
    statement = "; checkpoint;".join(statement)
    return statement % locals(), to_prefix + ".bam"


def process_remote_BAM(infile, token=None, outdir=".", filter_bed=None):
    """generate statement from .remote file"""

    statement = []
    files = []
    for line in open(infile):
        f = line.strip().split()
    
        if f[0] == "TCGA" or f[0] == "GDC":
            s, path = fetch_TCGA_BAM(f[1], token, outdir,
                                     filter_bed=filter_bed)
            statement.append(s)
            files.append(path)
        else:
            raise ValueError(
                "Repository %s not implimented for BAM files" % f[0])

    return "; checkpoint;".join(statement), files

'''
cgat_script_template.py
=============================================

:Tags: Python

Purpose
-------

Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import re
from future.moves.urllib.parse import urlencode
from future.moves.urllib.request import urlopen
import xml.etree.ElementTree as ET

import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    params = {'db': 'gds',
              'term': 'transcriptome[All Fields] AND "Homo sapiens"[Organism] AND high throughput sequencing[Platform Technology Type]',
              'retmax': 5000,
              'usehistory': 'y',
              }

    params = urlencode(params)
    query_filter = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    query_retrieve = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
    query_fetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    query_summary = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

    data = urlopen(query_filter, params)
    etree = ET.parse(data)
    root = etree.getroot()

    webenv = root.find("WebEnv").text
    query_key = root.find("QueryKey").text

    uids = [x.text for x in root.findall("*/Id")]

    pubmedlist = []

    params = {'dbfrom': 'gds',
              'db': 'pubmed',
              }

    params = urlencode(params)
    # necessary to preserve 1to1 links
    params += "&" + "&".join(["id=%s" % x for x in uids])
    data = urlopen(query_retrieve, params)

    etree = ET.parse(data)
    root = etree.getroot()

    map_uid2pmid = {}
    for linkset in root.findall("LinkSet"):
        uid = linkset.find("*/Id").text
        try:
            pmid = linkset.find("./LinkSetDb/Link/Id").text
        except AttributeError:
            pmid = None
        map_uid2pmid[uid] = pmid

    params = {'db': 'gds',
              'id': ",".join(uids)}

    params = urlencode(params)
    data = urlopen(query_fetch, params).read()

    map_uid2accession = {}
    map_uid2description = {}
    map_pmid2accession = {}

    for block in data.split("\n\n"):
        uid = re.search("ID: (\d+)", block).groups()[0]
        accession = re.search("Accession: (\S+)", block).groups()[0]
        map_uid2accession[uid] = accession
        lines = block.split("\n")
        description = lines[0]
        map_uid2description[uid] = description
        map_pmid2accession[map_uid2pmid[uid]] = accession.encode(
            "ascii", "ignore")

    url_pmid = "http://www.ncbi.nlm.nih.gov/pubmed/%s"
    url_geo = "http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s"

    outfile_uid = IOTools.openFile("uids.tsv", "w")
    outfile_uid.write("uid\tdescription\tpmid\taccession\n")
    for uid in uids:
        outfile_uid.write("\t".join(map(str,
                                        (uid,
                                         map_uid2description[uid],
                                         url_pmid % map_uid2pmid.get(uid, ""),
                                         url_geo % map_uid2accession.get(uid, "")))) + "\n")

    outfile_pmid = IOTools.openFile("pmid.tsv", "w")
    outfile_pmid.write("pmid\tyear\tjournal\ttitle\tabstract\tgeo\n")

    E.info("retrieving pubmed records")
    # output by pubmed id
    for pmid in list(map_uid2pmid.values()):
        if pmid is None:
            continue
        print(pmid)
        # retrieve record

        params = {'db': 'pubmed',
                  'id': pmid,
                  'retmode': 'xml'}

        params = urlencode(params)
        data = urlopen(query_fetch, params)

        etree = ET.parse(data)
        root = etree.getroot()
        article = root.find("PubmedArticle")
        assert article is not None
        journal = article.find("*//Journal")
        assert journal is not None
        year = journal.find("./JournalIssue/PubDate/Year").text
        journal_title = journal.find("Title").text

        title = article.find("*//ArticleTitle").text.encode("ascii", "ignore")
        try:
            abstract = article.find(
                "*//AbstractText").text.encode("ascii", "ignore")
        except AttributeError:
            abstract = ""

        outfile_pmid.write("\t".join(map(str,
                                         (url_pmid % pmid,
                                          year,
                                          journal_title,
                                          title,
                                          abstract,
                                          url_geo % map_pmid2accession.get(pmid, "")))) + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

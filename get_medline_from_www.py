################################################################################
#   Gene prediction pipeline 
#
#   $Id: get_medline_from_www.py 2782 2009-09-10 11:40:29Z andreas $
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
import os, sys, string, re, optparse, time, datetime

"""read Medline entries from Web
"""

import Experiment

from Bio import Medline
from Bio import PubMed

parser = optparse.OptionParser( version = "%prog version: $Id: get_medline_from_www.py 2782 2009-09-10 11:40:29Z andreas $")

if __name__ == "__main__":

    parser.add_option("-p", "--pmid", dest="pmid", type="string",
                      help="search for pmid."  )

    parser.add_option("-f", "--format", dest="format", type="string",
                      help="output format."  )

    parser.add_option("-c", "--clipboard", dest="clipboard",
                      help="output format." , action = "store_true"  )

    parser.add_option("-l", "--library", dest="library",
                      help="library to add medline entries to." , metavar="FILE"  )

    parser.set_defaults(
        pmid = None,
        format= "oo",
        clipboard = False,
        library = None,
        )

    (options, args) = Experiment.Start( parser )

    if options.pmid:
        ids = PubMed.search_for(options.pmid)

    outlines = []

    medline_parser = Medline.RecordParser()
    medline_dict = PubMed.Dictionary(parser = medline_parser)
            
    for id in ids:

        this_record = medline_dict[id]

        year = this_record.publication_date.split(" ")[0]

        last_names = map( lambda x: x.split(" ")[0], this_record.authors)

        if options.format == "oo":
            
            if len(last_names) > 2:
                a = last_names[0] + " et al."
            elif len(last_names) == 2:
                a = "%s & %s" % (last_names[0], last_names[1])
            elif len(last_names) >=1:
                a = last_names[0]
            else:
                a = "unknown_%s" % id

            brief = a + " (%s)" % str(year)
            outlines.append( brief )

            if options.loglevel >= 1:
                print "retrieved: %s '%s'" % (brief, this_record.title)
                
            outlines.append('\'%s\'' % this_record.title)
            
            outlines.append("%s" % " ".join(this_record.abstract.split("\n")))

            outlines.append("%s" % ", ".join(this_record.authors))
            outlines.append("%s %s" % (this_record.source, this_record.pubmed_id))
            
            outlines.append("http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=%s" % this_record.pubmed_id)

        else:
            raise "unknown format"
        
    if options.clipboard:
        os.system('dcop klipper klipper setClipboardContents \"%s\"' % '\n'.join(outlines))
    else:
        print "\n".join(outlines)

    if options.library:

        nids = []
        if os.path.exists( options.library ):
            pmids = map(lambda x: x[:-1].split(" ")[1], filter(lambda x: re.match( "PMID", x), open(options.library).readlines()))
            for id in ids:
                if id in pmids: continue
                nids.append(id)
            outfile = open(options.library, "a")                            
        else:
            outfile = open(options.library, "w")            

        f = lambda id, x: outfile.write( "\n%s\n" % (str(x)))
        PubMed.download_many( nids, f )
        outfile.close()

        print "added %i out of %i entries to library %s" % (len(nids), len(ids), options.library)
        
    Experiment.Stop()



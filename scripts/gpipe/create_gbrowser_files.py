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
gpipe/create_gbrowser_files.py - 
======================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python gpipe/create_gbrowser_files.py --help

Type::

   python gpipe/create_gbrowser_files.py --help

for command line help.

Documentation
-------------

Code
----

'''
import os
import sys
import string
import re
import random
import optparse
import CGAT.Experiment as E
import pgdb

"""python gpipe/create_gbrowser_files.py.

Script to create configuration files for gbrowser.

Default options will create the files in the devel section.
"""

a = """"""


def GetContigNames(dbhandle, tablename, max_size=0, order_by_rand=False):

    genome_lengths = {}

    statement = "SELECT DISTINCT(sbjct_token) FROM %s" % tablename

    if order_by_rand:
        statement += " ORDER BY sbjct_token "

    cc = dbhandle.cursor()
    try:
        cc.execute(statement)
        result = map(lambda x: x[0], cc.fetchall())
    except pgdb.DatabaseError, msg:
        print "# query failed with message", msg
        result = None

    if order_by_rand:
        random.shuffle(result)

    if max_size:
        return result[:max_size]
    else:
        return result

TEMPLATE_HEADER =\
    """
[GENERAL]
description   = %(description)s
db_adaptor    = Bio::DB::GFF
db_args       = -adaptor dbi::mysql
	        -dsn %(database)s:localhost
user          = nobody
pass          =

aggregators = match
		processed_transcript
		coding

plugins     = BatchDumper FastaDumper GFFDumper

# list of tracks to turn on by default
default features = Transcripts

# The class of the objects used to establish the reference coordinates.
reference class  =  Contig

# examples to show in the introduction
examples = %(examples)s

# "automatic" classes to try when an unqualified identifier is given
automatic classes = Gene Transcript Contig Chromosome

### HTML TO INSERT AT VARIOUS STRATEGIC LOCATIONS ###
# inside the <head></head> section
head = 

# at the top...
header = <div class="banner">
              <div id="title">OPTIC - Relase %(release)s</div>
        <div class="logo_mrc" id="logo_mrc.gif" >
        <a href="http://www.mrc.ac.uk"><img src="%(url_server)s/images/logo_mrc.gif" id="logo_mrc" alt="MRC logo" /></a></div>
      </div>
      [<a href="%(url_server)s/clades/%(clade)s">OPTIC clade server</a>] [<a href="%(url_server)s/clades/%(clade)s/helpGBrowser">Help</a>]

# a footer
footer = 

# Various places where you can insert your own HTML -- see configuration docs
html1 =
html2 = 
html3 = 
html4 = 
html5 = 
html6 = 

# what image widths to offer
image widths  = 450 640 800 1024

# default width of detailed view (pixels)
default width = 800

# Web site configuration info
stylesheet  = /gbrowse/gbrowse.css
buttons     = /gbrowse/images/buttons
tmpimages   = /gbrowse/tmp

# max and default segment sizes for detailed view
max segment     = 500000
default segment = 50000

# zoom levels
zoom levels    = 100 200 1000 2000 5000 10000 20000 40000 50000 100000 200000 300000 400000 500000

# colors of the overview, detailed map and key
overview bgcolor = lightgrey
detailed bgcolor = lightgoldenrodyellow
key bgcolor      = beige

label density = 25
bump density  = 100

########################
# Default glyph settings
########################

[TRACK DEFAULTS]
glyph         = generic
height        = 10
bgColor       = lightgrey
fgcolor       = black
font2color    = blue
label density = 25
bump density  = 100
# where to link to when user clicks in detailed view
link          = %(url_gene)s$name
"""

TEMPLATE_TRACKS_DEVEL =\
    """
################## TRACK CONFIGURATION ####################
# the remainder of the sections configure individual tracks
###########################################################

[Transcripts]
feature      = processed_transcript:gpipe_clean
glyph        = processed_transcript
bgcolor      = peachpuff
description  = 1
key          = Cleaned set of predictions
citation     = Transcripts are predictions after filtering. They are grouped into
	       genes if they share identical exons.

[GeneStructureTranscripts]
feature      = mRNA:gpipe_clean
glyph        = segments
key          = Cleaned set of predictions: Gene structure conservation
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = sub {
               my $feature = shift;
               return "peachpuff" unless $feature->has_tag('Status');
               my %%map_status2colour = ( 'CG' => "blue" ,
                                         'PG' => "lightblue",
                                         'SG' => "cyan",
                                         'CP' => "red",
                                         'PP' => "orange",
                                         'SP' => "orange",
                                         'RP' => "red" );
               my @c = $feature->get_tag_values('Status');
               if (defined $map_status2colour{$c[0]}) {
                   return $map_status2colour{$c[0]};
               } else {
                   return "lightgrey";
               }
               }
citation     = Colour bars indicating the quality of a particular prediction.
	       Genes are blue colours, pseudogenes are red colours and fragments
	       are grey colours.<P>
	       The detailed colour coding is:
	       <OL><LI>CG: blue</LI>
	           <LI>PG: lightblue</LI>
	           <LI>SG: cyan</LI>
	           <LI>CP: red</LI>
	           <LI>PP/SP: orange</LI>
	           <LI>RP: red</LI>
	       </OL>

[CdsIdentityTranscripts]
feature      = processed_transcript:gpipe_clean
glyph        = graded_segments
description  = 0
key          = Cleaned set of predictions: Cds identity
decorate_introns = 1
label        = 0
min_score    = 0
max_score    = 100
bgcolor      = red
fgcolor      = black
vary_fg      = 1
height       = 10
citation     = Percent identity of predicted coding exons to the
	       the query. The more intense the color, the higher
	       the percent identity. 

#################################################################
[Filtered]
feature      = processed_transcript:gpipe_filtered
glyph        = processed_transcript
bgcolor      = lightblue
description  = 1
key          = Filtered set of predictions
citation     = Predictions after filtering (removing redundant and conflicting predictions). Predictions are 
	       grouped into genes if they share identical exons.

[GeneStructureFiltered]
feature      = mRNA:gpipe_filtered
glyph        = segments
key          = Filtered set of predictions: Gene structure conservation
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = sub {
               my $feature = shift;
               return "peachpuff" unless $feature->has_tag('Status');
               my %%map_status2colour = ( 'CG' => "blue" ,
                                         'PG' => "lightblue",
                                         'SG' => "cyan",
                                         'CP' => "red",
                                         'PP' => "orange",
                                         'SP' => "orange",
                                         'RP' => "red" );
               my @c = $feature->get_tag_values('Status');
               if (defined $map_status2colour{$c[0]}) {
                   return $map_status2colour{$c[0]};
               } else {
                   return "lightgrey";
               }
               }
citation     = Colour bars indicating the quality of a particular prediction.
	       Genes are blue colours, pseudogenes are red colours and fragments
	       are grey colours.<P>
	       The detailed colour coding is:
	       <OL><LI>CG: blue</LI>
	           <LI>PG: lightblue</LI>
	           <LI>SG: cyan</LI>
	           <LI>CP: red</LI>
	           <LI>PP/SP: orange</LI>
	           <LI>RP: red</LI>
	       </OL>

[CdsIdentityFiltered]
feature      = processed_transcript:gpipe_filtered
glyph        = graded_segments
description  = 0
key          = Filtered set of Transcripts: Cds identity
decorate_introns = 1
label        = 0
min_score    = 0
max_score    = 100
bgcolor      = red
fgcolor      = black
vary_fg      = 1
height       = 10
citation     = Percent identity of predicted coding exons to the
	       the query. The more intense the color, the higher
	       the percent identity. 

###########################################################
[Predictions]
feature      = processed_transcript:gpipe_raw
glyph        = processed_transcript 
bgcolor      = #CCCC99
description  = 1
key          = Full set of predictions
citation     = The <B>Predictions</B> track contains all transcripts
	     that were predicted by the pipeline. This is the raw data before
	     filtering. For a cleaner subset, look at the <B>Transcripts</B> 
	     track. The description contains the following information from
	     left to right:
	     <OL><LI>Query identitifier</LI>
	         <LI>Coverage of query by prediction</LI>
		 <LI>Percent identity between query and prediction</LI>
		 <LI>Number of frameshifts in prediction</LI>	
	         <LI>Number of stopcodons in prediction</LI>	
             </OL>
	     
[GeneStructurePrediction]
feature      = mRNA:gpipe_raw
glyph        = segments
key          = Full set of predictions: Gene structure conservation 
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = sub {
	       my $feature = shift;
	       return "green" unless $feature->has_tag('Status');
	       my %%map_status2colour = ( 'CG' => "blue" ,
				         'PG' => "lightblue",
					 'SG' => "cyan",
					 'CP' => "red",
					 'PP' => "orange",
					 'SP' => "orange",
				         'RP' => "red" );
	       my @c = $feature->get_tag_values('Status');
               if (defined $map_status2colour{$c[0]}) {
	           return $map_status2colour{$c[0]};
               } else {
                   return "lightgrey";
               }
	       }
citation     = Colour bars indicating the quality of a particular prediction.
	       Genes are blue colours, pseudogenes are red colours and fragments
	       are grey colours.<P>
	       The detailed colour coding is:
	       <OL><LI>CG: blue</LI>
	           <LI>PG: lightblue</LI>
	           <LI>SG: cyan</LI>
	           <LI>CP: red</LI>
	           <LI>PP/SP: orange</LI>
	           <LI>RP: red</LI>
	       </OL>

[CdsIdentityPrediction]
feature      = processed_transcript:gpipe_raw
glyph        = graded_segments 
description  = 0
key          = Full set of predictions: Cds identity
decorate_introns = 1
label        = 0
min_score    = 0
max_score    = 100 
bgcolor      = red 
fgcolor      = black 
vary_fg      = 1
height       = 10
citation     = Percent identity of predicted coding exons to the
	       the query. The more intense the color, the higher
	       the percent identity. 

#######################################################################################
[DNA]
glyph          = dna
global feature = 1
height         = 40
do_gc          = 1
fgcolor        = red
axis_color     = blue
strand         = both
key            = DNA/GC Content

[CDS]
feature      = coding
glyph        = cds
key          = Frame usage

[Translation]
glyph          = translation
global feature = 1
height         = 40
fgcolor        = purple
start_codons   = 0
stop_codons    = 1
translation    = 6frame
key            = 6-frame translation
"""

TEMPLATE_TRACKS_RELEASE =\
    """
################## TRACK CONFIGURATION ####################
# the remainder of the sections configure individual tracks
###########################################################

[Transcripts]
feature      = processed_transcript:gpipe_clean
glyph        = processed_transcript
bgcolor      = peachpuff
description  = 1
key          = Cleaned set of predictions
citation     = Transcripts are predictions after filtering. They are grouped into
	       genes if they share identical exons.

[GeneStructureTranscripts]
feature      = mRNA:gpipe_clean
glyph        = segments
key          = Cleaned set of predictions: Gene structure conservation
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = sub {
               my $feature = shift;
               return "peachpuff" unless $feature->has_tag('Status');
               my %%map_status2colour = ( 'CG' => "blue" ,
                                         'PG' => "lightblue",
                                         'SG' => "cyan",
                                         'CP' => "red",
                                         'PP' => "orange",
                                         'SP' => "orange",
                                         'RP' => "red" );
               my @c = $feature->get_tag_values('Status');
               if (defined $map_status2colour{$c[0]}) {
                   return $map_status2colour{$c[0]};
               } else {
                   return "lightgrey";
               }
               }
citation     = Colour bars indicating the quality of a particular prediction.
	       Genes are blue colours, pseudogenes are red colours and fragments
	       are grey colours.<P>
	       The detailed colour coding is:
	       <OL><LI>CG: blue</LI>
	           <LI>PG: lightblue</LI>
	           <LI>SG: cyan</LI>
	           <LI>CP: red</LI>
	           <LI>PP/SP: orange</LI>
	           <LI>RP: red</LI>
	       </OL>

[CdsIdentityTranscripts]
feature      = processed_transcript:gpipe_clean
glyph        = graded_segments
description  = 0
key          = Cleaned set of predictions: Cds identity
decorate_introns = 1
label        = 0
min_score    = 0
max_score    = 100
bgcolor      = red
fgcolor      = black
vary_fg      = 1
height       = 10
citation     = Percent identity of predicted coding exons to the
	       the query. The more intense the color, the higher
	       the percent identity. 

#################################################################
[Filtered]
feature      = processed_transcript:gpipe_filtered
glyph        = processed_transcript
bgcolor      = lightblue
description  = 1
key          = Filtered set of predictions
citation     = Predictions after filtering (removing redundant and conflicting predictions). Predictions are 
	       grouped into genes if they share identical exons.

[GeneStructureFiltered]
feature      = mRNA:gpipe_filtered
glyph        = segments
key          = Filtered set of predictions: Gene structure conservation
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = sub {
               my $feature = shift;
               return "peachpuff" unless $feature->has_tag('Status');
               my %%map_status2colour = ( 'CG' => "blue" ,
                                         'PG' => "lightblue",
                                         'SG' => "cyan",
                                         'CP' => "red",
                                         'PP' => "orange",
                                         'SP' => "orange",
                                         'RP' => "red" );
               my @c = $feature->get_tag_values('Status');
               if (defined $map_status2colour{$c[0]}) {
                   return $map_status2colour{$c[0]};
               } else {
                   return "lightgrey";
               }
               }
citation     = Colour bars indicating the quality of a particular prediction.
	       Genes are blue colours, pseudogenes are red colours and fragments
	       are grey colours.<P>
	       The detailed colour coding is:
	       <OL><LI>CG: blue</LI>
	           <LI>PG: lightblue</LI>
	           <LI>SG: cyan</LI>
	           <LI>CP: red</LI>
	           <LI>PP/SP: orange</LI>
	           <LI>RP: red</LI>
	       </OL>

[CdsIdentityFiltered]
feature      = processed_transcript:gpipe_filtered
glyph        = graded_segments
description  = 0
key          = Filtered set of Transcripts: Cds identity
decorate_introns = 1
label        = 0
min_score    = 0
max_score    = 100
bgcolor      = red
fgcolor      = black
vary_fg      = 1
height       = 10
citation     = Percent identity of predicted coding exons to the
	       the query. The more intense the color, the higher
	       the percent identity. 

###########################################################
[Predictions]
feature      = processed_transcript:gpipe_raw
glyph        = processed_transcript 
bgcolor      = #CCCC99
description  = 1
key          = Full set of predictions
citation     = The <B>Predictions</B> track contains all transcripts
	     that were predicted by the pipeline. This is the raw data before
	     filtering. For a cleaner subset, look at the <B>Transcripts</B> 
	     track. The description contains the following information from
	     left to right:
	     <OL><LI>Query identitifier</LI>
	         <LI>Coverage of query by prediction</LI>
		 <LI>Percent identity between query and prediction</LI>
		 <LI>Number of frameshifts in prediction</LI>	
	         <LI>Number of stopcodons in prediction</LI>	
             </OL>
	     
[GeneStructurePrediction]
feature      = mRNA:gpipe_raw
glyph        = segments
key          = Full set of predictions: Gene structure conservation 
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = sub {
	       my $feature = shift;
	       return "green" unless $feature->has_tag('Status');
	       my %%map_status2colour = ( 'CG' => "blue" ,
				         'PG' => "lightblue",
					 'SG' => "cyan",
					 'CP' => "red",
					 'PP' => "orange",
					 'SP' => "orange",
				         'RP' => "red" );
	       my @c = $feature->get_tag_values('Status');
               if (defined $map_status2colour{$c[0]}) {
	           return $map_status2colour{$c[0]};
               } else {
                   return "lightgrey";
               }
	       }
citation     = Colour bars indicating the quality of a particular prediction.
	       Genes are blue colours, pseudogenes are red colours and fragments
	       are grey colours.<P>
	       The detailed colour coding is:
	       <OL><LI>CG: blue</LI>
	           <LI>PG: lightblue</LI>
	           <LI>SG: cyan</LI>
	           <LI>CP: red</LI>
	           <LI>PP/SP: orange</LI>
	           <LI>RP: red</LI>
	       </OL>

[CdsIdentityPrediction]
feature      = processed_transcript:gpipe_raw
glyph        = graded_segments 
description  = 0
key          = Full set of predictions: Cds identity
decorate_introns = 1
label        = 0
min_score    = 0
max_score    = 100 
bgcolor      = red 
fgcolor      = black 
vary_fg      = 1
height       = 10
citation     = Percent identity of predicted coding exons to the
	       the query. The more intense the color, the higher
	       the percent identity. 

#######################################################################################
[DNA]
glyph          = dna
global feature = 1
height         = 40
do_gc          = 1
fgcolor        = red
axis_color     = blue
strand         = both
key            = DNA/GC Content

[CDS]
feature      = coding
glyph        = cds
key          = Frame usage

[Translation]
glyph          = translation
global feature = 1
height         = 40
fgcolor        = purple
start_codons   = 0
stop_codons    = 1
translation    = 6frame
key            = 6-frame translation
"""

TEMPLATE_REPEATS =\
    """
[CompTE]
feature      = mobile_element:CompTE
glyph        = segments
key          = Repeats: CompTE predicted transposable elements
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = red
citation     = please contact Casey Bergman for more information.

[RepeatRunner]
feature      = repeat_region:RepeatRunner
glyph        = segments
key          = Repeats: RepeatRunner predicted repeats.
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = green
citation     = please contact Casey Bergman for more information.

[RepeatRunnerPiler]
feature      = repeat_region:RepeatRunnerPiler
glyph        = segments
key          = Repeats: RepeatRunner predicted repeats (Piler).
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = blue
citation     = please contact Casey Bergman for more information.
"""

TEMPLATE_TRNAS =\
    """
[TRNA]
feature = processed_transcript:TRNA_ATT
glyph        = segments
key          = tRNAs predicted by Casey Bergman
description  = 0
label        = 0
height       = 2
bgcolor      = black
fgcolor      = blue
citation     = please contact Casey Bergman for more information.
"""

TEMPLATE_ORTHOLOGY =\
    """
[%s]
feature      = ortho_%s
glyph        = span
description  = 1
key          = Orthologs in %s
decorate_introns = 1
label        = 1
bgcolor      = red 
fgcolor      = black
link         = ../%s?name=$name
citation     = Orthologs between %s and %s. The code presents the gene/transcript degeneracy.
                For example, 1m-mm is a one-to-many gene and many-to-many transcript orthology
                relationship.
"""

# name, database, DroSpeGe, Eisen

FLY_DATA = (
    ("Dyakuba", "dyak_vs_dmel7", "dyak", "yakuba"),
    ("Dmelanogaster", "dmel_vs_dmel4", "dmel", "melanogaster"),
    ("Dpseudoobscura", "dpse_vs_dmel8", "dpse", "pseudoobscura"),
    ("Dvirilis", "dvir_vs_dmel6", "dvir", "virilis"),
    ("Dmojavensis", "dmoj_vs_dmel6", "dmoj", "mojavensis"),
    ("Dananassae", "dana_vs_dmel6", "dana", "ananassae"),
    ("Dgrimshawi", "dgri_vs_dmel5", "dgri", "grimshawi"),
    ("Dsimulans", "dsim_vs_dmel6", "dsim", "simulans"),
    ("Derecta", "dere_vs_dmel6", "dere", "erecta"),
    ("Dsechellia", "dsec_vs_dmel3", "dsec", "sechellia"),
    ("Dpersimilis", "dper_vs_dmel3", "dper", "persimilis"),
    ("Dwillistoni", "dwil_vs_dmel1", "dwil", "willistoni"),
)

WORM_DATA = (
    ("Celegans",  "ww_cele1", "", ""),
    ("Cbriggsae", "ww_cbri2", "", ""),
    ("Cfour",     "ww_cfou2", "", ""),
    ("Cremanei",  "ww_crem2", "", ""),
)


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gpipe/create_gbrowser_files.py 2781 2009-09-10 11:33:14Z andreas $")

    parser.add_option("-t", "--target-dir", dest="target", type="string",
                      help="target directory.", metavar="FILE")
    parser.add_option("-d", "--database", dest="database", type="string",
                      help="database.", metavar="FILE")
    parser.add_option("-u", "--url", dest="", type="string",
                      help="gr.", metavar="FILE")
    parser.add_option("-r", "--release", dest="release", type="string",
                      help="release [%default].")
    parser.add_option("-s", "--style", dest="style", type="choice",
                      choices=("devel", "release"))
    parser.add_option("-c", "--contrib", dest="contrib", type="string",
                      help="additional tracks with contributed data [repeats]")
    parser.add_option("-p", "--project", dest="project", type="choice",
                      choices=("flies", "worms", "mammals"),
                      help="project to add")

    parser.set_defaults(
        connection="db:andreas",
        target="/var/www/conf/gbrowse-devel.conf",
        database="gbrowser_devel_%s",
        release="1v5",
        style="devel",
        contrib="",
        url_server="http://genserv.anat.ox.ac.uk",
        url_gene="http://genserv.anat.ox.ac.uk%s/clades/%s/queryForGene?name=",
    )

    (options, args) = E.Start(parser, add_database_options=True)

    if options.contrib:
        options.contrib = options.contrib.lower().split(",")

    dbhandle = pgdb.connect(options.psql_connection)

    if options.style == "release":
        style_code = ""
    else:
        style_code = "/devel"

    # selection of data and contig prefix
    use_prefix = False
    if options.project == "flies":
        data = FLY_DATA
        use_prefix = True
        if options.contrib:
            contrib = options.contrib
        else:
            contrib = ["repeats", ]
        options.url_gene = options.url_gene % (style_code, "flies")
        clade = "flies"

    elif options.project == "worms":
        data = WORM_DATA
        contrib = None
        options.url_gene = options.url_gene % (style_code, "worms")
        clade = "worms"

    for name, schema, map2drospege, map2eisen in data:

        filename = "%s/%s.conf" % (options.target, name)

        if options.loglevel >= 1:
            print "setting up %s" % filename
            sys.stdout.flush()

        outfile = open(filename, "w")

        database = options.database % schema

        contig_names = GetContigNames(dbhandle, "%s.%s" % (schema, "predictions"),
                                      max_size=30,
                                      order_by_rand=True)

        if contig_names:
            # patch to add species prefix
            if use_prefix:
                prefix = re.sub("vs_.*", "", schema)
                if not re.match(prefix, contig_names[0]):
                    contig_names = map(lambda x: prefix + x, contig_names)

            contig_names.sort()

            contig_names = " ".join(contig_names)
        else:
            contig_names = ""

        params = {
            "release": options.release,
            "description": name,
            "database": database,
            "examples": contig_names,
            "map2drospege": map2drospege,
            "map2eisen": map2eisen,
            "url_gene": options.url_gene,
            "url_server": options.url_server,
            "clade": clade,
        }

        outfile.write(TEMPLATE_HEADER % params)
        if options.style == "devel":
            outfile.write(TEMPLATE_TRACKS_DEVEL % params)
        elif options.style == "release":
            outfile.write(TEMPLATE_TRACKS_RELEASE % params)

        for name2, schema2, map2dropsege2, map2eisen2 in data:

            if schema == schema2:
                continue

            outfile.write(TEMPLATE_ORTHOLOGY % (
                "Orthology%s" % name2,
                schema2,
                name2,
                name2,
                name,
                name2))

            outfile.write("\n")

        if contrib:
            for con in options.contrib:
                if con == "repeats":
                    outfile.write(TEMPLATE_REPEATS)
                elif con == "trnas":
                    outfile.write(TEMPLATE_TRNAS)

        outfile.close()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

##############################################################################
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
###############################################################################
'''
cgat_script_template.py
=============================================

:Author: Mike Morgan
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Example::

   python controller_script.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import optparse
import CGAT.Experiment as E
import interface


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-t", "--test", dest="test", type="string",
                      help="supply help")

    parser.add_option("--input-bam-file", dest="bam_file", type="string",
                      help="input bam file, where applicable")

    parser.add_option("--config-file", dest="config", type="string",
                      help="location of config file if not"
                      " in current directory")

    parser.add_option("--input-gtf-file", dest="gtf_file", type="string",
                      help="input gtf/gff file, where applicable")

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    parser.set_defaults(config=None)

    # find the config file and parse into a dict if not specified
    if options.config is None:
        dir_contents = os.listdir(os.getcwd())
        config_file = [x for x in dir_contents if re.search("config.ini", x)]
        E.info("parsing config file %s" % config_file)
        params = interface.parse_config(config_file)
    else:
        E.info("parsing config file %s" % options.config)
        params = interface.parse_config(options.config)

    # instantiate interface object, then assign parsed config file
    # and set genome view from config file
    inter_inst = interface.Interface()
    inter_inst.set_params(params)
    PARAMS = inter_inst.get_params()
    inter_inst.set_genome_view()

    E.info("genome view is %s:%s-%s" % (PARAMS['genome_contig'],
                                        PARAMS['genome_start'],
                                        PARAMS['genome_end']))

    # parse data and intervals then store in interface object
    E.info("getting data from %s" % options.bam_file)
    inter_inst.set_data_vector(options.bam_file,
                               data_format=PARAMS['data_format'])
    E.info("getting ranges from %s" % options.gtf_file)
    inter_inst.set_intervals_vector(options.gtf_file,
                                    file_format=PARAMS['intervals_type'])

    # count number of track instances required
    n_tracks = int(PARAMS['graphs_n']) + int(PARAMS['intervals_n'])
    inter_inst.set_n_tracks(n_tracks)

    # collate interval types
    track_types = []
    track_types.append(PARAMS['graphs_type'])
    track_types.append(PARAMS['genes_models'])

    # instantiate Layout
    figure = inter_inst.makeLayout(inter_inst.get_genome_view())

    E.info("plotting tracks")
    for n in range(n_tracks):
        track = track_types[n]
        if track == "transcripts":
            track_type = "ranges"
            data_format = PARAMS['intervals_type']
            inter_inst.makeTrack(track_type=track_type,
                                 genome_view=inter_inst.get_genome_view(),
                                 vector=inter_inst.get_ranges(data_format),
                                 layout=figure,
                                 priority=int(PARAMS['genes_priority']))

        elif track == "histogram":
            data_format = PARAMS['data_format']
            track_type = "data"
            inter_inst.makeTrack(track_type=track_type,
                                 genome_view=inter_inst.get_genome_view(),
                                 vector=inter_inst.get_data(data_format),
                                 layout=figure,
                                 priority=int(PARAMS['graphs_priority']))

    figure.generatePlot()

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

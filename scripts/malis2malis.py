'''
malis2malis.py - extract subsets of malis
=========================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

extract sequences from a set of multiple alignments and distribute into new directories.

Files can be split according to an input file given to the option --filename-components.
This file is a tab separated table containing the following fields::

   1       sequence        a sequence identifier
   2       input_id        The id under which the alignment is found.
   3       component_id    The component_id which under which the directory shall be stored.
                           This is optional, if no third column is specified, the input_id
                           is used.

Usage
-----

Example::

   python malis2malis.py --help

Type::

   python malis2malis.py --help

for command line help.

Command line options
--------------------

'''
import os
import sys
import re
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Malis as Malis

USAGE = """python %s [OPTIONS]

""" % sys.argv[0]


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id$",
        usage=globals()["__doc__"])

    parser.add_option("-d", "--pattern-output",
                      dest="pattern_output",
                      type="string",
                      help="filename pattern for output multiple "
                      "alignment files.")

    parser.add_option("-f", "--filter-tsv-file",
                      dest="filename_filter",
                      type="string",
                      help="filename with strings to filter by.")

    parser.add_option("--list-filter", dest="list_filter", type="string",
                      help="list of strings to filter by.")

    parser.set_defaults(
        pattern_output="%s.mali",
        methods="",
        parameters="",
        filename_filter=None,
        list_filter=None,
    )

    Malis.addOptions(parser)

    (options, args) = E.Start(parser)

    options.methods = options.methods.split(",")
    options.parameters = options.parameters.split(",")

    if not options.pattern_mali:
        raise "Please specifiy a pattern to find the malis using --pattern-mali"

    ####################################################################
    ####################################################################
    ####################################################################
    # Read components
    ####################################################################
    map_seq_id2component, map_component2seq_id, map_component2input_id = \
        Malis.readComponents(options)

    ####################################################################
    ####################################################################
    ####################################################################
    # Read filtering information
    ####################################################################
    if options.filename_filter:
        id_filter, nerrors = IOTools.ReadList(
            open(options.filename_filter, "r"))
        if options.loglevel >= 1:
            options.stdlog.write(
                "# read %i identifiers to filter each multiple alignment with.\n" %
                len(id_filter))
            options.stdlog.flush()
    elif options.list_filter:
        id_filter = options.list_filter.split(",")
    else:
        id_filter = None

    ####################################################################
    ####################################################################
    ####################################################################
    # Read regions to mask
    ####################################################################
    map_component2masks = Malis.readMasks(options, map_component2input_id)

    ####################################################################
    ####################################################################
    ####################################################################
    # Read regions to extract
    ####################################################################
    map_component2extracts = Malis.readExtracts(
        options, map_component2input_id)

    ####################################################################
    ####################################################################
    ####################################################################
    # Read regions to annotate
    ####################################################################
    map_component2annotations = Malis.readAnnotations(
        options, map_component2input_id)

    ####################################################################
    ####################################################################
    ####################################################################
    # Prepare for run
    ####################################################################
    component_ids = map_component2seq_id.keys()
    component_ids.sort()

    if options.loglevel >= 1:
        options.stdlog.write(
            "# %i component ids to start with.\n" % (len(component_ids)))

    component_ids, map_sample2reference = Malis.selectComponents(
        component_ids,
        map_component2seq_id,
        map_component2input_id,
        id_filter,
        options)

    if options.test:
        component_ids = component_ids[:options.test]

    if options.loglevel >= 1:
        options.stdlog.write(
            "# %i component ids selected for output.\n" % (len(component_ids)))

    ninput = 0
    noutput = 0
    nskipped = 0
    nskipped_length = 0

    for component_id in component_ids:

        ninput += 1

        if options.loglevel >= 3:
            options.stdlog.write(
                "# processing component %s\n" % (component_id))

        mali = Malis.getMali(component_id,
                             map_component2seq_id,
                             map_component2input_id,
                             id_filter,
                             options)

        if mali is None:
            E.warn("empty mali returned for component %s" % (component_id))
            nskipped += 1
            continue

        if mali.getNumColumns() == 0:
            E.warn("skipping output of empty alignment for component %s" %
                   (component_id))
            nskipped += 1
            continue

        mali.setName(str(component_id))

        ###############################################################
        # add annotations
        if map_component2annotations is not None:
            Malis.annotateAlignment(mali,
                                    map_component2annotations,
                                    options)

        ###############################################################
        # mask the alignment
        Malis.maskAlignment(mali,
                            map_component2masks,
                            map_component2extracts,
                            map_sample2reference,
                            options)

        if mali.getNumColumns() < options.minimum_mali_length:
            nskipped_length += 1
            if options.loglevel >= 1:
                options.stdlog.write("# component %s: skipped, because length %i less than threshold.\n" % (
                    component_id, mali.getNumColumns()))
            continue

        ###############################################################
        # prepare the mali for output
        if "%s" not in options.pattern_output:
            append = True
        else:
            append = False

        output_filename = re.sub("%s", component_id, options.pattern_output)
        input_id = map_component2input_id[component_id]

        if options.loglevel >= 2:
            options.stdlog.write("# component %s: input from %s, goes to %s\n" % (
                component_id, input_id, output_filename))

        dirname = os.path.dirname(output_filename)

        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname)

        if not os.path.exists(output_filename):
            mali.writeToFile(
                open(output_filename, "w"), format=options.output_format)
            noutput += 1
        else:
            if append:
                mali.writeToFile(
                    open(output_filename, "a"), format=options.output_format)
                noutput += 1
            else:
                if options.loglevel >= 1:
                    options.stdlog.write("# skipping because output for component %s already exists: %s\n" % (
                        component_id, output_filename))
                nskipped += 1

        # if we only sample, stop if you have reached
        # the desired number
        if options.sample and noutput == options.sample:
            break

    E.info("ninput=%i, noutput=%i, nskipped=%i, nskipped_length=%i" %
           (ninput, noutput, nskipped, nskipped_length))

    E.Stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))

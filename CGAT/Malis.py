'''Malis.py
==============

Utility functions for working with multiple multiple alignments.
'''

import os
import sys
import re
import random
import collections

import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Intervals as Intervals
import CGAT.Mali as Mali


def addOptions(parser):
    """add options to parser concerning components and masks."""

    parser.add_option("-b", "--pattern-component", dest="pattern_component", type="string",
                      help="how to extract identifier of mali from component name.")

    parser.add_option("-c", "--filename-components", dest="filename_components", type="string",
                      help="filename of components to choose for each multiple alignment.")

    parser.add_option("-e", "--pattern-filter", dest="pattern_filter", type="string",
                      help="pattern of how to extract filter from identifier.")

    parser.add_option("-i", "--input-format", dest="input_format", type="choice",
                      choices=("plain", "fasta", "clustal", "phylip"),
                      help="input format of multiple alignment")

    parser.add_option("-o", "--output-format", dest="output_format", type="choice",
                      choices=(
                          "fasta", "codeml", "phylip", "profile", "stockholm"),
                      help="output format of multiple alignment")

    parser.add_option("--filename-mask-regions", dest="filename_mask_regions", type="string",
                      help="""mask regions given file. Input format is component_id\tstart\tend with coordinates in
the multiple alignment.""" )

    parser.add_option("--filename-extract-regions", dest="filename_extract_regions", type="string",
                      help="""extract regions given file. Input format is component_id\tstart\tend with coordinates in
the multiple alignment. The alignment is restricted to the given coordinates by masking everything else.""" )

    parser.add_option("--filename-annotate-regions", dest="filename_annotate_regions", type="string",
                      help="""annotate regions given file. Input format is component_id\tstart\tend\tlabel with coordinates in
the multiple alignment. The label should be a single letter. If the file is empty, the full alignment is output with label N.""")

    parser.add_option("--use-input-id", dest="use_input_id", action="store_true",
                      help="""use input-id instead of component_id in --filename-X-regions""" )

    parser.add_option("-a", "--pattern-mali", dest="pattern_mali", type="string",
                      help="filename pattern for input multiple alignment files.")

    parser.add_option("--test", dest="test", type="int",
                      help="test output with first # components.")

    parser.add_option("--skip-doubles", dest="skip_doubles", action="store_true",
                      help="skip doubles, takes first entry encountered.")

    parser.add_option("--ignore-missing", dest="ignore_missing", action="store_true",
                      help="ignore missing alignments [%default].")

    parser.add_option("--remove-all-gaps", dest="remove_all_gaps", type=int,
                      help="remove fully masked or gapped positions in a multiple alignment. The integer parameter supplies the frame. Use 3 for codon alignments.")

    parser.add_option("--remove-any-gaps", dest="remove_any_gaps", type=int,
                      help="remove positions in a multiple alignment containing at least one gap. The integer parameter supplies the frame. Use 3 for codon alignments.")

    parser.add_option("--sample", dest="sample", type="int",
                      help="sample # components.")

    parser.add_option("--sample-method", dest="sample_method", type="choice",
                      choices=(
                          "simple-without-replacement", "length-without-replacement", ),
                      help="sampling method: without-replacement=sample without replacement. length-sample: sample components of same length (without replacement).")

    parser.add_option("--filename-sample-reference", dest="filename_sample_reference", type="string",
                      help="filename to with reference components. The components should be part of the file supplied by --filename-components.")

    parser.add_option("--minimum-mali-length", dest="minimum_mali_length", type="int",
                      help="minimum multiple alignment length.")

    parser.set_defaults(
        pattern_mali=None,
        pattern_component="(^\S+)",
        pattern_filter="(^\S+)",
        input_format="fasta",
        output_format="fasta",
        test=0,
        skip_doubles=False,
        remove_all_gaps=None,
        remove_any_gaps=None,
        sample=0,
        filename_sample_reference=None,
        sample_method="simple-without-replacement",
        filename_components=None,
        filename_mask_regions=None,
        filename_annotated_regions=None,
        filename_extract_regions=None,
        minimum_mali_length=0,
        use_input_id=False,
        ignore_missing=False,
    )


def readComponents(options):
    """read components from filename supplied in the options.
    """
    if options.filename_components:
        map_seq_id2component =\
            IOTools.ReadMap(open(options.filename_components, "r"),
                            columns="all",
                            both_directions=False)

        map_component2seq_id = {}
        map_component2input_id = {}
        for key, val in map_seq_id2component.items():
            if isinstance(val, str):
                input_id = val
                output_id = val
            elif isinstance(val, str):
                if len(val) == 2:
                    input_id = val[0]
                    output_id = val[1]
                else:
                    input_id = val[0]
                    output_id = val[0]
            else:
                raise ValueError("error in reading %s: %s->%s" %
                                 (options.filename_components, key, val))

            if output_id not in map_component2seq_id:
                map_component2seq_id[output_id] = []
            map_component2seq_id[output_id].append(key)
            map_component2input_id[output_id] = input_id
        return (map_seq_id2component,
                map_component2seq_id,
                map_component2input_id)
    else:
        return None, None, None


def readMasks(options, map_component2input_id):
    """read masking information from filename supplied in the options.
    """

    if options.use_input_id:
        map_id2component = IOTools.getInvertedDictionary(
            map_component2input_id)
    else:
        map_id2component = {}
        for component in map_component2input_id.keys():
            map_id2component[component] = (component,)

    if options.filename_mask_regions:
        map_component2masks = collections.defaultdict(list)
        if not os.path.exists(options.filename_mask_regions):
            options.stdlog.write(
                "# could not find %s - ignored \n" % options.filename_mask_regions)
        else:
            for line in open(options.filename_mask_regions, "r"):
                if line[0] == "#":
                    continue
                id, start, end = line[:-1].split("\t")
                try:
                    start, end = int(start), int(end)
                except ValueError:
                    continue
                if id not in map_id2component:
                    continue
                for x in map_id2component[id]:
                    map_component2masks[x].append((start, end))
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# read masks for %i malis.\n" % len(map_component2masks))
                options.stdlog.flush()
    else:
        map_component2masks = None

    return map_component2masks


def readExtracts(options, map_component2input_id):
    """read extract information from filename supplied in the options.
    """
    if not options.filename_extract_regions:
        return None

    if options.use_input_id:
        map_id2component = IOTools.getInvertedDictionary(
            map_component2input_id)
    else:
        map_id2component = {}
        for component in map_component2input_id.keys():
            map_id2component[component] = (component,)

    map_component2extracts = collections.defaultdict(list)
    if not os.path.exists(options.filename_extract_regions):
        options.stdlog.write(
            "# could not find %s - ignored \n" % options.filename_extract_regions)
    else:
        for line in open(options.filename_extract_regions, "r"):
            if line[0] == "#":
                continue
            id, start, end = line[:-1].split("\t")
            start, end = int(start), int(end)
            if id not in map_id2component:
                continue
            for x in map_id2component[id]:
                map_component2extracts[x].append((start, end))
        if options.loglevel >= 1:
            options.stdlog.write(
                "# read extracts for %i malis.\n" % len(map_component2extracts))
            options.stdlog.flush()

    return map_component2extracts


def readAnnotations(options, map_component2input_id):
    """read annotation information from filename supplied in the options.
    """
    if not options.filename_annotate_regions:
        return None

    if options.use_input_id:
        map_id2component = IOTools.getInvertedDictionary(
            map_component2input_id)
    else:
        map_id2component = {}
        for component in map_component2input_id.keys():
            map_id2component[component] = (component,)

    map_component2annotations = collections.defaultdict(list)
    if not os.path.exists(options.filename_annotate_regions):
        options.stdlog.write(
            "# could not find %s - ignored \n" % options.filename_annotate_regions)
    else:
        for line in open(options.filename_annotate_regions, "r"):
            if line[0] == "#":
                continue
            try:
                id, start, end, label = line[:-1].split("\t")
            except ValueError:
                raise ValueError("parsing error in line %s\n" % (line[:-1]))

            start, end = int(start), int(end)
            if id not in map_id2component:
                continue
            for x in map_id2component[id]:
                map_component2annotations[x].append((start, end, label))

        if options.loglevel >= 1:
            options.stdlog.write(
                "# read annotations for %i malis.\n" % len(map_component2annotations))
            options.stdlog.flush()

    return map_component2annotations


def maskAlignment(mali,
                  map_component2masks,
                  map_component2extracts,
                  map_sample2reference,
                  options):
    """mask an alignment.

    If map_sample2reference is given, coordinates in references are used to
    mask residues in sample.
    """

    id = mali.getName()

    if options.loglevel >= 5:
        options.stdout.write("# multiple alignment %s before masking:\n" % id)
        mali.writeToFile(sys.stdout)

    def getMasks(id, map1, map_sample2reference):
        if map_sample2reference and id in map_sample2reference:
            xid = map_sample2reference[id]
            if options.loglevel >= 1:
                options.stdlog.write(
                    "# using mapped coordinates from %s for %s\n" % (xid, id))
        else:
            xid = id
        if xid in map1:
            return map1[xid]
        else:
            return []

    ##########################################################################
    # mask alignment
    if map_component2masks:
        masks = getMasks(id, map_component2masks, map_sample2reference)
        if masks:
            for start, end in masks:

                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# component: %s: masking region due to mask %i-%i\n" % (id, start, end))

                mali.maskColumns(range(start, min(end, mali.getWidth())))

    ##########################################################################
    # extract regions from an alignment by masking everything else
    if map_component2extracts:
        masks = getMasks(id, map_component2extracts, map_sample2reference)

        if masks:
            other_masks = Intervals.complementIntervals(
                masks, 0, mali.getWidth())
            for start, end in other_masks:
                if options.loglevel >= 1:
                    options.stdlog.write(
                        "# component: %s: masking region due to extract %i-%i\n" % (id, start, end))

                mali.maskColumns(range(start, min(end, mali.getWidth())))

    if options.loglevel >= 5:
        options.stdout.write("# multiple alignment after masking:\n")
        mali.writeToFile(sys.stdout)

    if mali.getAlphabet() == "aa":
        gap_chars = "Xx-"
    else:
        gap_chars = "XxNn-"

    if options.remove_all_gaps:

        width_before = mali.getNumColumns()

        mali.removePattern(
            match_function=lambda x: x in gap_chars,
            minimum_matches=mali.getNumSequences(),
            search_frame=1,
            delete_frame=options.remove_all_gaps)

        width_after = mali.getNumColumns()

        if options.loglevel >= 1:
            options.stdlog.write("# component: %s: removed %i fully gapped/masked columns, old size=%i, new size=%i\n" %
                                 (id, width_before - width_after, width_before, width_after))

    if options.remove_any_gaps:

        width_before = mali.getNumColumns()

        mali.removePattern(
            match_function=lambda x: x in gap_chars,
            minimum_matches=1,
            search_frame=1,
            delete_frame=options.remove_any_gaps)

        width_after = mali.getNumColumns()

        if options.loglevel >= 1:
            options.stdlog.write("# component: %s: removed %i columns containing at least one gap/mask, old size=%i, new size=%i\n" %
                                 (id, width_before - width_after, width_before, width_after))

    if options.loglevel >= 5:
        options.stdout.write("# multiple alignment after cleaning:\n")
        mali.writeToFile(sys.stdout)


def annotateAlignment(mali,
                      map_component2annotations,
                      options):
    """annotate an alignment.

    Do not use this function if options.clean_mali is True and 
    regions have been extracted. Both might result in a change
    of the length of the multiple alignment and thus the annotated
    coordinates would be incorrect.
    """

    id = mali.getName()

    if map_component2annotations is None:
        return

    max_length = mali.getNumColumns()

    if len(map_component2annotations) == 0:
        mali.addAnnotation("STATE", "N" * max_length)
        return

    if id not in map_component2annotations:

        if options.loglevel >= 1:
            options.stdlog.write("# no annotations found for %s.\n" % id)
        return

    max_length = mali.getNumColumns()
    annotation = ["-"] * max_length

    for start, end, label in map_component2annotations[id]:
        end = min(end, max_length)
        annotation[start:end] = [label] * (end - start)

    mali.addAnnotation("STATE", "".join(annotation))

master_mali = None


def getMali(component_id,
            map_component2seq_id,
            map_component2input_id,
            id_filter,
            options):

    global master_mali

    rx_component = re.compile(options.pattern_component)

    mali = Mali.Mali()

    nsubstitutions = len(re.findall("%s", options.pattern_mali))

    input_id = rx_component.search(component_id).groups()[0]
    input_id = map_component2input_id[input_id]

    if nsubstitutions == 0:

        if master_mali is None:

            master_mali = Mali.Mali()

            E.debug("retrieving multiple alignment from file %s" %
                    (options.pattern_mali))

            master_mali.readFromFile(
                open(options.pattern_mali, "r"), format=options.input_format)

        for s in map_component2seq_id[component_id]:

            if options.pattern_filter and id_filter:
                f = re.search(options.pattern_filter, s).groups()[0]

                if f not in id_filter:
                    E.debug("removing %s from %s: not in filter" %
                            (f, component_id))
                    continue

            if options.output_format == "codeml":
                if len(master_mali[s]) % 3 != 0:
                    raise ValueError(
                        "length of sequence %s is not a multiple of 3: %i" % (s, len(master_mali[s])))

            if s in mali:
                if options.skip_doubles:
                    E.warn("skipped double entry %s in component %s" %
                           (s, component_id))
                    return None
                else:
                    raise ValueError(
                        "duplicate entry %s in component %s" % (s, component_id))

            mali.addEntry(master_mali.getEntry(s))

    else:

        input_filename = options.pattern_mali % tuple(
            [input_id] * nsubstitutions)

        E.debug("retrieving multiple alignment for component %s from file %s" %
                (component_id, input_filename))

        if not os.path.exists(input_filename):
            if options.ignore_missing:
                E.warn("alignment %s not found" % input_filename)
                return None
            else:
                raise OSError("alignment %s not found" % input_filename)

        mali.readFromFile(
            open(input_filename, "r"), format=options.input_format)

        # get identifiers (and make a copy)
        s = tuple(mali.getIdentifiers())
        for ss in s:
            if options.pattern_filter and id_filter:
                f = re.search(options.pattern_filter, ss).groups()[0]

                if f not in id_filter:
                    mali.deleteEntry(ss)
                    if options.loglevel >= 5:
                        options.stdlog.write(
                            "# removing %s from %s: not in filter.\n" % (ss, component_id))
                    continue

            if ss not in map_component2seq_id[component_id]:
                if options.loglevel >= 5:
                    options.stdlog.write(
                        "# removing %s from %s: not in component list.\n" % (ss, component_id))
                mali.deleteEntry(ss)
            else:
                if options.output_format == "codeml":
                    if len(mali[ss]) % 3 != 0:
                        raise "length of sequence %s is not a multiple of 3: %i" % (
                            ss, len(mali[ss]))

    mali.setName(component_id)

    return mali


def selectComponents(component_ids,
                     map_component2seq_id,
                     map_component2input_id,
                     id_filter,
                     options):
    """select a set of components from component_ids.
    """

    map_sample2reference = {}

    if options.sample:

        if options.sample_method == "simple-without-replacement":
            random.shuffle(component_ids)

        elif options.sample_method == "length-without-replacement":

            map_component_id2length = {}

            for component_id in component_ids:

                mali = getMali(component_id,
                               map_component2seq_id,
                               map_component2input_id,
                               id_filter,
                               options)

                if not mali:
                    continue

                map_component_id2length[component_id] = mali.getWidth()

            reference_ids, nerrors = IOTools.ReadList(
                open(options.filename_sample_reference, "r"))

            # do not sample from the reference set
            sampled_components = set(reference_ids)

            new_component_ids = []

            ninput, noutput, nskipped = 0, 0, 0
            # now go through reference set
            for ref_id in reference_ids:

                ninput += 1
                if ref_id not in map_component_id2length:
                    if options.loglevel >= 1:
                        options.stdlog.write(
                            "# reference component %s not found.\n" % (str(ref_id)))
                    nskipped += 1
                    continue

                ref_length = map_component_id2length[ref_id]
                ref_length_min = ref_length - 50
                ref_length_max = ref_length + 50

                # find all components with a length similar to ref_length
                # excluding previously sampled ones.
                test_components = filter(lambda x: ref_length_min < map_component_id2length[
                                         x] < ref_length_max, component_ids)
                test_components = list(
                    set(test_components).difference(sampled_components))

                if len(test_components) == 0:
                    if options.loglevel >= 1:
                        options.stdlog.write("# reference components %s: skipped - no others with equivalent length around %i found." %
                                             (ref_id, ref_length))
                    nskipped += 1
                    continue

                random.shuffle(test_components)

                component_id = test_components[0]
                sampled_components.add(component_id)

                map_sample2reference[component_id] = ref_id

                if options.loglevel >= 1:
                    options.stdlog.write("# reference component mapping: %s\t%s\t%i\t%i\t%i\n" % (
                        ref_id,
                        component_id,
                        ref_length,
                        map_component_id2length[component_id],
                        len(test_components)))

                new_component_ids.append(component_id)
                noutput += 1

            if options.loglevel >= 1:
                options.stdlog.write("# sampling results: ninput=%i, noutput=%i, nskipped=%i\n" % (
                    ninput, noutput, nskipped))

            component_ids = new_component_ids
            options.sample = len(new_component_ids)

    return component_ids, map_sample2reference

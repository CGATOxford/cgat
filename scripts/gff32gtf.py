'''gff32gtf.py - various methods for converting gff3 files to gtf
====================================================

:Author: Ian Sudbery
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Provide a range of methods for converting GFF3 formated files to valid GTF 
format files.

Background
----------

While the various flavours of GFF format are supposedly backward
compatible, this is broken by GTF2.2 and GFF3. GTF requires the
presence of gene_id and transcript_id fields for each record. This not
so for GFF3. Further key,value tags in the attributes fields of GTF
are " " delimited, but are "=" delimited in GFF.

Conversion is non-trivial. GFF3 records are hierachical. To find the
gene_id and transcript_id one must traverse the hierarchy to the
correct point. Futher records can have multiple parents.

                                               -> Exon
While the standard structure is Gene -> mRNA -|       ,
                                               -> CDS

this is not manditory, and it is possible the conversion will want to
be done in a different way.

Usage
-----

Example::

   python gff32gtf.py --method=[METHOD] [options]

Their are several ways in which the conversion can be done:

hierachical
+++++++++++

By default this script will read in the entire GFF3 file, and then for
each entry traverse the hierarchy until an object of type GENE_TYPE
("gene" by default") or an object with no parent is found. This
becomes the "gene_id". Any object of TRANSCRIPT_TYPE encountered on
the way is set as the transcript_id. If not such object is encountered
then the object directly below the gene object is used as the
trancript_id. Objects that belong to multipe transcripts or genes are
duplicated.

This method requires ID and Parent fields to be present.

Because this method reads the whole file in, it uses the most memory, although
see --read-twice and --by-chrom for tricks that might help.

set-field
+++++++++

The gene_id and transcript_id fields are set to the  value of a provided field.
Records that don't have these fields are discarded. By default:

transcript_id=ID
gene_id=Parent

set-pattern
+++++++++++

As above, but the fieldnames are set by a string format involving the
fields of the record.

set-none
++++++++

transcript_id and gene_id are set to None.

Command line options
--------------------

'''

import sys

import CGAT.Experiment as E
import CGAT.GFF3 as GFF3
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools


def search_hierarchy(ID, hierarchy, options):
    '''Returns a three element tuple of lists.

        * The first two lists are the gene_ids and transcript_ids that
        * are associated with specified IDs.  The third is a list of
        * possible transcript_ids - that is trancript_ids that are one
        * level below where the gene id came from.

        All three lists are guarenteed to be the same length, but both
        the transcript lists could contain None values where no
        transcript_id has been found.

        Works by calling it self recursively, not efficient, but does
        deal with the problem of cicular references: the recursion
        limit will be quickly reached.

        Can also raise ValueError if no feature of type
        options.gene_type is found and options.missing_gene is false

    '''

    gene_id = []
    transcript_id = []
    possible_transcript_id = []

    entry = hierarchy[ID]

    if entry['type'] == options.gene_type:
        gene_id.append(hierarchy[ID]['gene_id'])

        if not entry['type'] == options.transcript_type:
            transcript_id = [None]
            possible_transcript_id = [None]
        else:
            transcript_id = [entry['transcript_id']]
            possible_transcript_id = [None]

        return (gene_id, transcript_id, possible_transcript_id)

    for parent in entry['Parent']:

        new_gene_id, new_transcript_id, new_possible_transcript_id = search_hierarchy(
            parent, hierarchy, options)

        gene_id.extend(new_gene_id)
        transcript_id.extend(new_transcript_id)
        possible_transcript_id.extend(new_possible_transcript_id)

    if options.missing_gene:
        possible_transcript_id = [
            entry['transcript_id'] if x is None else x for x in possible_transcript_id]

    if len(gene_id) == 0 and options.missing_gene:
        gene_id = [entry['gene_id']]
        transcript_id = [None]
        possible_transcript_id = [None]
    elif len(gene_id) == 0 and not options.missing_gene:
        raise ValueError(
            "Root found without finding an object of type %s" % options.gene_type)

    if entry['type'] == options.transcript_type:
        transcript_id = [
            entry['transcript_id'] if x is None else x for x in transcript_id]

    assert len(gene_id) == len(transcript_id) and len(
        transcript_id) == len(possible_transcript_id)
    assert len(gene_id) > 0

    return gene_id, transcript_id, possible_transcript_id


def convert_hierarchy(first_gffs, second_gffs, options):
    ''' Converts GFF to GTF by parsing the hierarchy.
    First parses :param:first_gffs to build the hierarchy then iterates over second_gffs
    using a call to the recursive function search_hierarchy to identify gene_ids and transcript_ids.

    If multiple gene and transcript_ids are found outputs a record for each combination.

    If no definitive transcript_id is found and options.missing_gene is True, it will use the 
    possible_transcript_id as transcript_id, which is the ID one level below the entry used as gene_id.
    If this is also None (that is there was only on level), sets transcript_id to gene_id.

    Might raise ValueError if options.missing_gene is false and either no gene or no transcript_id
    was found for an entry.

    Might raise RuntimeError if the recursion limit was reached because the input contains circular
    references. '''

    hierarchy = {}

    for gff in first_gffs:

        if not(options.parent == "Parent"):
            if options.parent in gff.asDict():
                gff['Parent'] = gff[options.parent].split(",")
            else:
                gff['Parent'] = []

        hierarchy[gff['ID']] = {
            "type": gff.feature,
            "Parent": gff.asDict().get("Parent", []),
            "gene_id": gff.attributes.get(
                options.gene_field_or_pattern, gff['ID']),
            "transcript_id": gff.attributes.get(
                options.transcript_field_or_pattern, gff['ID'])}

    for gff in second_gffs:

        if options.discard and (
                (options.missing_gene and options.parent not in gff) or (
                gff.feature in (options.gene_type, options.transcript_type))):

            continue

        gene_ids, transcript_ids, poss_transcript_ids = search_hierarchy(
            gff['ID'], hierarchy, options)

        assert len(gene_ids) > 0 and len(transcript_ids) > 0

        if options.missing_gene:

            transcript_ids = [poss if found is None else found
                              for found, poss in
                              zip(transcript_ids, poss_transcript_ids)]

            transcript_ids = [gid if found is None else found
                              for found, gid in
                              zip(transcript_ids, gene_ids)]

        elif None in transcript_ids:
            raise ValueError("failed to find transcript id for %s" % gff['ID'])

        for gene_id, transcript_id in zip(gene_ids, transcript_ids):

            gff.gene_id = gene_id
            gff.transcript_id = transcript_id

            gtf_entry = GTF.Entry()
            gtf_entry.copy(gff)
            if "Parent" in gtf_entry:
                gtf_entry['Parent'] = ",".join(gtf_entry['Parent'])

            options.stdout.write(str(gtf_entry) + "\n")


def convert_set(gffs, gene_pattern, transcript_pattern, options):
    ''' creates the gene_id and transcript_id fields from a string format pattern using
    fields of the gff. '''

    for gff in gffs:

        gff.gene_id = str(gene_pattern) % gff.asDict()
        gff.transcript_id = str(gene_pattern) % gff.asDict()

        gtf_entry = GTF.Entry()

        gtf_entry.copy(gff)
        if "Parent" in gtf_entry:
            gtf_entry['Parent'] = ",".join(gtf_entry['Parent'])

        options.stdout.write(str(gtf_entry) + "\n")


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-m", "--method", dest="method", type="choice", action="store",
                      choices=(
                          "hierarchy", "set-field", "set-pattern", "set-none"),
                      help="Method to use for conversion")

    parser.add_option("-g", "--gene-type", dest="gene_type", type="string",
                      help="feature type to get gene_id from if possible [%default]")

    parser.add_option("-t", "--transcript-type", dest="transcript_type", type="string",
                      help="feature type to get transcript_id from if possible [%default]")

    parser.add_option("-d", "--no-discard", dest="discard", action="store_false",
                      help="Do not discard feature types specified by GENE_TYPE and TRANSCRIPT_TYPE")

    parser.add_option("--gene-id", dest="gene_field_or_pattern", type="string",
                      help="Either field or pattern for the gene_id [%default]")

    parser.add_option("--transcript-id", dest="transcript_field_or_pattern", type="string",
                      help="Either field or pattern for the transcript_id [%default]")

    parser.add_option("--parent-field", dest="parent", type="string",
                      help="field that specifies the parent relationship. Currently only"
                      "if left as Parent will features with multiple parents be parsed"
                      "correctly""")

    parser.add_option("--read-twice", dest="read_twice", action="store_true",
                      help="Instead of holding the whole file in memory, read once for parsing the "
                      "hierarchy, and then again for actaully doing the conversion. Means a real file "
                      "and not a pipe must be provided.""")

    parser.add_option("--by-chrom", dest="by_chrom", action="store_true",
                      help="Parse input file one choromosome at a time. Reduces memory usage, "
                      "but input must be sorted by chromosome and features may not split accross "
                      " multiple chromosomes""")

    parser.add_option("--fail-missing-gene", dest="missing_gene", action="store_false",
                      help="Fail if no feature of type GENE_TYPE is found instead of using "
                      "defaulting to highest object in hierarchy""")

    parser.set_defaults(
        method="hierarchy",
        gene_type="gene",
        transcript_type="mRNA",
        discard=True,
        gene_field_or_pattern="ID",
        transcript_field_or_pattern="ID",
        read_twice=False,
        by_chrom=False,
        missing_gene=True,
        parent="Parent"
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    gffs = GFF3.flat_file_iterator(options.stdin)

    if options.by_chrom:
        gffs = GFF3.chrom_iterator(gffs)
    else:
        gffs = [gffs]

    # running early so that fails early if configuration is wrong
    if options.read_twice:
        # Will throw IOError if options.stdin is not a normal file
        second_gff = GFF3.flat_file_iterator(
            IOTools.openFile(options.stdin.name))

        if options.by_chrom:
            second_gff = GFF3.chrom_iterator(second_gff)
        else:
            second_gff = iter([second_gff])
    else:
        second_gff = None

    for chunk in gffs:

        if options.read_twice:
            second_gff_chunk = second_gff.next()
        else:
            chunk = list(chunk)
            second_gff_chunk = chunk

        if options.method == "hierarchy":

            convert_hierarchy(chunk, second_gff_chunk, options)
        elif options.method == "set-field":
            gene_id_pattern = "%%(%s)s" % options.gene_field_or_pattern
            transcript_id_pattern = "%%(%s)s" % options.transcript_field_or_pattern
            convert_set(chunk, gene_id_pattern, transcript_id_pattern, options)
        elif options.method == "set-pattern":
            convert_set(chunk, options.gene_field_or_pattern,
                        options.transcript_field_or_pattern, options)
        elif options.method == "set-none":
            convert_set(chunk, None, None, options)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

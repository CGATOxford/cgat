'''
genes2ensembl.py - map gene identifiers to ENSEMBL identifiers
==============================================================

:Author:
:Tags: Python

Purpose
-------

Map identifiers/gene names into ENSEMBL identifiers.

Options
-------

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import collections
import CGAT.Experiment as E
import CGAT.Biomart as Biomart


def buildIdentifierMap(query_species):
    columns = ('ensembl_gene_id',
               'entrezgene',
               'hgnc_id',
               'hgnc_symbol')

    data = Biomart.biomart_iterator(
        columns,
        dataset="%s_gene_ensembl" % query_species)

    map_identifiers = collections.defaultdict(set)
    for row in data:
        ensid = row['ensembl_gene_id']
        for column in columns[1:]:
            xid = str(row[column])
            if xid == "NA" or xid == "":
                continue

            if ensid.startswith('LRG_'):
                continue

            map_identifiers[xid].add((ensid, column))

    # convert to lists
    map_identifiers = dict([(x, list(y)) for x, y in list(map_identifiers.items())])

    return map_identifiers


def buildOrthologyMap(query_species,
                      target_species,
                      filter_type="ortholog_one2one"):
    '''build map of genes in query species to
    those in target species.'''

    columns = ('ensembl_gene_id',
               '%s_homolog_ensembl_gene' % target_species,
               '%s_homolog_orthology_type' % target_species,
               '%s_homolog_orthology_confidence' % target_species)

    data = Biomart.biomart_iterator(
        columns,
        dataset="%s_gene_ensembl" % query_species)

    map_query2target = dict([(
        x['ensembl_gene_id'],
        x['%s_homolog_ensembl_gene' % target_species]) for x in data
        if x['%s_homolog_orthology_type' % target_species] == filter_type])

    return map_query2target


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-q", "--query-species",
        dest="query_species", type="string",
        help="query species [%default]")

    parser.add_option(
        "-t", "--target-species",
        dest="target_species", type="string",
        help="target species [%default]")

    parser.add_option(
        "-c", "--column",
        dest="column", type="string",
        help="column number or name with gene name to map [%default]")

    parser.set_defaults(
        query_species="hsapiens",
        target_species=None,
        column=1)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.target_species:
        map_orthologs = buildOrthologyMap(
            options.query_species,
            options.target_species)
    else:
        map_orthologs = None

    E.info("orthology map: %i -> %i (%i unique)" %
           (len(map_orthologs),
            len(list(map_orthologs.values())),
            len(set(map_orthologs.values()))))

    map_identifiers = buildIdentifierMap(
        options.query_species)

    E.info("identifier map: %i -> %i" %
           (len(map_identifiers),
            len(list(map_identifiers.values()))))

    first = True
    outfile = options.stdout
    c = E.Counter()

    for line in options.stdin:
        if line.startswith("#"):
            continue
        data = line[:-1].split("\t")

        if first:
            try:
                column = data.index(options.column)
                data[column] = "gene_id"
            except ValueError:
                column = int(options.column) - 1
            outfile.write("\t".join(data) + "\n")
            first = False

        orig_id = data[column]
        gene_id = data[column].upper()
        c.input += 1

        if gene_id in map_identifiers:
            m = map_identifiers[gene_id]
            if len(m) > 1:
                c.skipped_multiple_identifiers += 1
                E.warn("skipped: %s - multiple identifiers: %s" %
                       (gene_id, m))
                continue

            new_id, method = m[0]
            c[method] += 1
            gene_id = new_id
        else:
            c.skipped_no_identifiers += 1
            continue

        if map_orthologs:
            if gene_id in map_orthologs:
                c.has_ortholog += 1
                gene_id = map_orthologs[gene_id]
            else:
                c.skipped_no_orthologs += 1
                E.warn("skipped: no ortholog %s->%s" % (orig_id, gene_id))
                continue

        c.output += 1
        data[column] = gene_id

        outfile.write("\t".join(data) + "\n")

    E.info(c)

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

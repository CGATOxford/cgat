'''vcf2vcf.py - manipulate vcf files
=================================

Purpose
-------

Manipulate vcf-formatted files.


Usage
-----

Example::

   cat in.vcf | python vcf2vcf.py - --reorder=alphabetical > sorted.vcf

This command generates a sorted vcf with the sample columns in in.vcf
in alphabetical order.

Type::

   python vcf2vcf.py --help

for command line usage.

Methods
-------

This script provides the following methods:

re-order
   reorder sample columns in vcf formatted file according to a given sort order

Documentation
-------------

This is a tool for manipulating vcf-formatted files.  The following
options are available:

+-----------+-------------------------+
|--reorder  |reorders sample columns  |
+-----------+-------------------------+

Reorder
^^^^^^^

To sort sample columns into alphabetical order::

   cat example.vcf | python vcf2vcf.py - --reorder=alphabetical

This will sort the columns in the example.vcf into the order "#CHROM,
POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE_A, SAMPLE_B,
SAMPLE_C, ..."

To specify a non-alphabetical order::

   cat example.vcf | python vcf2vcf.py - --reorder=SAMPLE_C,SAMPLE_A,SAMPLE_B,...

This will sort the columns in the example.vcf into the order "#CHROM,
POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, | SAMPLE_C, SAMPLE_A,
SAMPLE_B, ..."

Command line options
--------------------

'''

import sys
import CGAT.Experiment as E
import CGAT.VCF as VCF


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "--method", dest="methods", type="choice", action="append",
        choices=("reorder-samples",),
        help="method to apply [default=%default]")

    parser.add_option(
        "--sort-order", dest="sort_order",
        help="sort order for sample names. Give column names as "
        "comma-separated list or specify ``alphabetical`` "
        "[default=%default]")

    parser.set_defaults(
        methods=[],
        sort_order="alphabetical",
    )

    (options, args) = E.Start(parser, add_pipe_options=True)

    if not options.methods:
        raise ValueError("no method specified")

    infile = VCF.VCFFile(options.stdin)

    sort_order = False
    if "reorder-samples" in options.methods:
        if options.sort_order:
            sort_order = options.sort_order.split(",")
            if "alphabetical" in sort_order:
                sort_order = sorted(infile.samples)

    infile.writeHeader(options.stdout, order=sort_order)

    for vcf in infile:
        if sort_order:
            vcf.order = sort_order
        options.stdout.write(str(vcf) + "\n")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

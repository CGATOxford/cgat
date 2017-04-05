'''
substitute_tokens.py - substitute fields in tables
==================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script substitutes the elements within columns within an table.


Usage
-----

Example::

   python substitute_tokens.py --help

Type::

   python substitute_tokens.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import os
import CGAT.IOTools as IOTools
import CGAT.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if not argv:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id: cgat_script_template.py 2871 2010-03-03 10:20:44Z andreas $",
                            usage=globals()["__doc__"])

    parser.add_option("-c", "--create", dest="create", type="string",
                      help="create substitution list [%default] ")

    parser.add_option("-r", "--regex-token", dest="regex_token", type="string",
                      help="regular expression for tokens (has to create one pair of brackets) [%default] ")

    parser.add_option("-p", "--pattern-sub", dest="pattern_sub", type="string",
                      help="pattern for substitution [%default] ")

    parser.add_option("-a", "--map-tsv-file", dest="apply", type="string",
                      help="apply substitution list [%default] ")

    parser.add_option("-x", "--extended", dest="extended", action="store_true",
                      help="replace not just with second column in map, but all columns. [%default] ")

    parser.add_option("-i", "--invert", dest="invert", action="store_true",
                      help="pairs of substitution patterns is to,from [%default] ")

    parser.add_option("-m", "--multiple", dest="multiple", action="store_true",
                      help="do multiple substitutions per row [%default] ")

    parser.add_option("-e", "--echo", dest="echo", action="store_true",
                      help="echo susbstituted column [%default] ")

    parser.add_option("-k", "--keep", dest="keep", action="store_true",
                      help="keep column that is substituted [%default] ")

    parser.add_option("-f", "--method=filter --filter-method", dest="filter", action="store_true",
                      help="remove lines not matching [%default] ")

    parser.add_option("-y", "--reverse-filter", dest="reverse_filter", action="store_true",
                      help="remove lines matching [%default] ")

    parser.add_option("-n", "--inplace", dest="inplace", action="store_true",
                      help="do inplace subsitutions of all files on command line [%default] ")

    parser.add_option("-b", "--backup", dest="backup", action="store_true",
                      help="keep backup (with ending .bak) [%default] ")

    parser.add_option("--keep-header", dest="keep_header", action="store_true",
                      help="do not apply transformation to header [%default] ")

    parser.add_option("-o", "--columns-token", dest="columns_token", type="string",
                      help="substitute tokens in columns [%default] ")

    parser.add_option("-s", "--select-rows", dest="regex_rows", type="string",
                      help="regular expression for rows to use. [%default] ")

    parser.set_defaults(create=None,
                        regex_token=None,
                        pattern_sub="%s",
                        apply=None,
                        invert=False,
                        multiple=False,
                        columns_token=None,
                        filter=None,
                        reverse_filter=None,
                        inplace=False,
                        backup=False,
                        regex_rows=None,
                        echo=False,
                        extended=False,
                        keep=False,
                        keep_header=False)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.regex_token:
        options.regex_token = re.compile(options.regex_token)

    if options.regex_rows:
        options.regex_rows = re.compile(options.regex_rows)

    if options.columns_token:
        if options.columns_token != "all":
            options.columns_token = [int(x) - 1 for x in options.columns_token.split(",")]

    file_id = 0

    keys = {}

    if options.apply:
        infile = IOTools.openFile(options.apply, "r")
        for line in infile:
            if line[0] == "#":
                continue

            d = line[:-1].split("\t")
            try:
                a, b = d[:2]
            except ValueError:
                print("# invalid map skipped in line: %s" % line)
                continue

            if options.invert:
                a, b = b, a
                if options.extended:
                    b = "\t".join(d[0] + d[2:])
            else:
                if options.extended:
                    b = "\t".join(d[1:])

            if a not in keys:
                keys[a] = []

            if options.keep:
                b = a + "\t" + b

            keys[a].append(b)

    files = args

    if not options.inplace and len(args) == 0:
        files = ["-"]

    for file in files:

        close_infile = False
        close_outfile = False
        if file == "-":
            infile = sys.stdin
            outfile = sys.stdout
        else:
            if options.inplace:
                os.rename(file, file + ".bak")
                infile = IOTools.openFile(file + ".bak", "r")
                outfile = IOTools.openFile(file, "w")
                close_infile = True
                close_outfile = True
            else:
                infile = IOTools.openFile(file, "r")
                outfile = sys.stdout
                close_infile = True

        first = True

        for line in infile:
            if line[0] == "#":
                outfile.write(line)
                continue

            if first:
                first = False
                if options.keep_header:
                    outfile.write(line)
                    continue

            if options.regex_rows:
                if options.regex_rows.search(line):
                    outfile.write(line)
                    continue

            new_lines = []
            if options.regex_token:
                r = options.regex_token.search(line[:-1])
                while r:
                    key = r.group(1)
                    if key not in keys:
                        if options.create:
                            keys[key] = [options.pattern_sub % str(len(keys))]
                        else:
                            new_lines.append(line[:-1])
                            break

                    for k in keys[key]:
                        new_lines.append(
                            line[:r.start(1)] + k + line[r.end(1):-1])

                    if options.multiple:
                        r = options.regex_token.search(line[r.end(1):-1])
                    else:
                        break
                else:
                    if not options.filter:
                        new_lines.append(line[:-1])

            elif options.columns_token:
                data = line[:-1].split("\t")
                if options.columns_token == "all":
                    columns = list(range(len(data)))
                else:
                    columns = options.columns_token
                keep = not options.reverse_filter
                first_multiple = True
                for c in columns:
                    k = data[c]
                    if k in keys:
                        if len(keys[k]) > 1:
                            if not first_multiple:
                                raise "warning: could not substitute multiple keys for %s in multiple columns in line: %s" % (
                                    k, line)
                            first_multiple = False
                        for v in keys[k]:
                            if options.echo:
                                data.append(data[c])
                            # multiple substitutions: write data now
                            data[c] = v
                            if keep:
                                new_lines.append("\t".join(data))
                        keep = False
                    else:
                        if options.create:
                            keys[k] = [options.pattern_sub % str(len(keys))]
                            data[c] = keys[k][0]
                        elif options.filter:
                            keep = False
                        elif options.reverse_filter:
                            keep = True
                if keep:
                    new_lines.append("\t".join(data))

            elif options.apply:
                for key in keys:
                    for k in keys[key]:
                        line = line.replace(key, k)
                new_lines.append(line[:-1])

            if new_lines:
                outfile.write("\n".join(new_lines) + "\n")

        if options.create:
            create_file = IOTools.openFile(options.create, "w")
            for key in keys:
                for k in keys[key]:
                    create_file.write("%s\t%s\n" % (key, str(k)))
            create_file.close()

        if close_outfile:
            outfile.close()
        if close_infile:
            infile.close()

        if options.inplace and not options.backup:
            os.remove(file + ".bak")

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

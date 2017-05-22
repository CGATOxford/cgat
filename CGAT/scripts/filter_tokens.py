'''
filter_tokens.py - select columns from a table 
==============================================

:Tags: Python

Purpose
-------

.. todo::
   
   describe purpose of the script.

Usage
-----

Example::

   python filter_tokens.py --help

Type::

   python filter_tokens.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import CGAT.Experiment as E

USAGE = """python %s [token1 [token2 [...]]] < stdin > stdout

filter lines where pattern matches

""" % sys.argv[0]


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: filter_tokens.py 2782 2009-09-10 11:40:29Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-r", "--regex-token", dest="regex_token", type="string",
                      help="regular expression for tokens (has to create one pair of brackets).")

    parser.add_option("-f", "-a", "--map-tsv-file", "--filename-tokens", dest="filename_tokens", type="string",
                      help="filename with tokens.")

    parser.add_option("-o", "--tokens-tsv-file", "--tokens-tsv-file", dest="tokens", type="string",
                      help="',' separated list of tokens [default=%default].")

    parser.add_option("-c", "--columns", dest="columns", type="string",
                      help="comma separated list of columns.")

    parser.add_option("-i", "--invert-match", dest="invert_match", action="store_true",
                      help="only keep non-matching entries.")

    parser.add_option("-t", "--no-titles", dest="titles", action="store_false",
                      help="no titles.")

    parser.set_defaults(
        regex_token="^(\S+)",
        filename_tokens=None,
        columns=None,
        invert=False,
        titles=True,
        tokens=None,
    )

    (options, args) = E.Start(parser)

    if options.columns:
        options.columns = [int(x) - 1 for x in options.columns.split(",")]

    if len(args) == 0 and not options.filename_tokens:
        print(USAGE, "please specify tokens on command line or supply file with tokens.")
        sys.exit(1)

    regex_token = re.compile(options.regex_token)

    file_id = 0

    ##########################################
    # build list of keys to keep
    keys = {}

    for a in args:
        keys[a] = 1

    if options.filename_tokens:
        infile = IOTools.openFile(options.filename_tokens, "r")
        for line in infile:
            if line[0] == "#":
                continue
            a = line[:-1].split("\t")[0]
            keys[a] = 1

    if options.tokens:
        for a in options.tokens.split(","):
            keys[a] = 1

    ##########################################

    ninput = 0
    nkept = 0

    for line in sys.stdin:

        if line[0] == "#":
            continue
        if line[0] == ">":
            print(line[:-1])
            continue

        ninput += 1

        # skip first line if there are titles
        if ninput == 1 and options.titles:
            print(line[:-1])
            continue

        found = False

        if options.columns:
            data = line[:-1].split("\t")
            for c in options.columns:
                try:
                    if data[c] in keys:
                        found = True
                        break
                except IndexError:
                    continue

        elif regex_token:
            r = regex_token.search(line[:-1])
            if r:
                if r.groups()[0] in keys:
                    found = True

        if options.invert_match:
            found = not found

        if found:
            nkept += 1
            print(line[:-1])

    if options.loglevel >= 1:
        options.stdlog.write(
            "# input=%i, kept=%i, discarded=%i\n" % (ninput, nkept, ninput - nkept))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

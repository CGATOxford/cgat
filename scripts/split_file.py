'''
split_file.py - split a file into parts
=======================================

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

   python split_file.py --help

Type::

   python split_file.py --help

for command line help.

Command line options
--------------------

'''
import sys
import re
import string
import os
import getopt
import CGAT.Experiment as E

USAGE = """python %s < stdin > stdout

split a file into chunks.

OPTIONS:
-h, --help                      print this message.
-v, --verbose=                  loglevel.
-r, --split-regex               split at regular expression
-a, --after                     split after match
-s, --skip                      do not echo match
-p, --pattern-output            pattern of output files (has to contain %s)
-c, --column=                   split according to column
-m, --map=                      split according to map
-d, --dry-run                   echo files that would be created,
                                but do not create any.
-e, --header-names                    add header to each file
-r, --remove-key                remove key column
-append                         append data to existing files.
--pattern-identifier            if given, use this pattern to extract
                                id from column.
--chunk-size                    Number of matching records in each output file
--version                       output version information
""" % (sys.argv[0], "s")


def CreateOpen(file, mode="w", dry_run=False, header=None):
    """open file. Check first, if directory exists.
    """

    if dry_run:
        print "# opening file %s" % file
        return open("/dev/null", mode)

    if mode in ("w", "a"):
        dirname = os.path.dirname(file)
        if dirname and not os.path.exists(dirname):
            os.makedirs(dirname)

    if os.path.exists(file):
        existed = True
    else:
        existed = False

    f = open(file, mode)

    if header and not existed:
        f.write(header + "\n")

    return f


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    param_long_options = [
        "verbose=", "help", "split-regex=", "after", "pattern-output=", "skip",
        "column=", "map=", "dry-run",
        "header", "remove-key", "append", "pattern-identifier=", "version",
        "chunk-size="]

    param_short_options = "v:hr:ap:sc:dek"

    param_loglevel = 1
    param_split_at_regex = None
    param_after = None
    param_skip = None
    param_pattern_output = "%s.chunk"
    param_split_column = None
    param_filename_map = None
    param_dry_run = False
    param_header = False
    param_remove_key = False
    param_append = "w"
    param_pattern_identifier = None
    param_chunk_size = 1

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      param_short_options,
                                      param_long_options)

    except getopt.error, msg:
        print USAGE, msg
        sys.exit(1)

    for o, a in optlist:
        if o in ("-v", "--verbose"):
            param_loglevel = int(a)
        elif o in ("--version", ):
            print "version="
            sys.exit(0)
        elif o in ("-h", "--help"):
            print USAGE
            sys.exit(0)
        elif o in ("-r", "--split-regex"):
            param_split_at_regex = re.compile(a)
        elif o in ("-a", "--after"):
            param_after = 1
        elif o in ("-s", "--skip"):
            param_skip = 1
        elif o in ("-p", "--pattern-output"):
            param_pattern_output = a
        elif o in ("-c", "--column"):
            param_split_column = int(a) - 1
        elif o in ("-m", "--map"):
            param_filename_map = a
        elif o in ("-d", "--dry-run"):
            param_dry_run = True
        elif o in ("-e", "--header-names"):
            param_header = True
        elif o in ("-r", "--remove-key"):
            param_remove_key = True
        elif o == "--append":
            param_append = "a"
        elif o == "--pattern-identifier":
            param_pattern_identifier = re.compile(a)
        elif o == "--chunk-size":
            param_chunk_size = int(a)


    print E.GetHeader()
    print E.GetParams()

    mymap = {}
    if param_filename_map:
        infile = open(param_filename_map, "r")
        for line in infile:
            if line[0] == "#":
                continue
            data = line[:-1].split("\t")[:2]
            mymap[data[0]] = data[1]

    filenames = set()
    found = set()
    ninput, noutput = 0, 0

    if param_split_column is not None:

        header = None
        files = {}
        for line in sys.stdin:

            if line[0] == "#":
                continue

            ninput += 1

            if param_header:
                if not header:
                    header = line[:-1]
                    continue
            else:
                header = None

            data = line[:-1].split("\t")

            try:
                key = data[param_split_column]
            except ValueError:
                continue

            if param_pattern_identifier:
                key = param_pattern_identifier.search(key).groups()[0]

            if mymap:
                if key in mymap:
                    key = mymap[key]
                else:
                    continue

            found.add(key)

            filename = re.sub("%s", key, param_pattern_output)
            filenames.add(filename)

            if filename not in files:

                # reset if too many files are open
                if len(files) > 1000:
                    if param_loglevel >= 1:
                        print "# resetting all files."
                        sys.stdout.flush()

                    for f in files.values():
                        f.close()
                    files = {}

                files[filename] = CreateOpen(
                    filename, "a", param_dry_run, header)

            if param_remove_key:
                del data[param_split_column]
                files[filename].write(string.join(data, "\t") + "\n")
            else:
                files[filename].write(line)

            noutput += 1

        for f in files.values():
            f.close()

    else:
        file_id = 0

        filename = re.sub("%s", str(file_id), param_pattern_output)
        outfile = CreateOpen(filename, param_append, param_dry_run)
        nlines = 0

        header = param_header
        split = 0

        for line in sys.stdin:

            if param_split_at_regex and param_split_at_regex.search(line[:-1]):
                split += 1

            if split == param_chunk_size:
                if param_after:
                    nlines += 1
                    outfile.write(line)
                if nlines > 0:
                    outfile.close()
                    file_id += 1
                    filename = re.sub("%s", str(file_id), param_pattern_output)
                    outfile = CreateOpen(
                        filename, param_append, param_dry_run, header)
                    filenames.add(filename)
                    split = 0

                nlines = 0
                if param_after or param_skip:
                    continue

            outfile.write(line)
            nlines += 1

        outfile.close()

    if param_loglevel >= 1:
        sys.stdout.write(
            "# ninput=%i, noutput=%i, nfound=%i, nnotfound=%i, nfiles=%i\n" % (
                ninput,
                noutput,
                len(found),
                len(set(mymap).difference(found)),
                len(filenames)))

    print E.GetFooter()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

'''
cgat_scan_email.py - read email folders in a directory
=======================================================

:Author:
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

This script reads all emails and returns the date
of the last email send to a project directory.

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

import os
import glob
import re
import sys
import CGAT.Experiment as E

# See: http://stackoverflow.com/questions/2301789/read-a-file-in-reverse-order-using-python


def reverse_readline(filename, buf_size=8192):
    """a generator that returns the lines of a file in reverse order"""
    with open(filename) as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        total_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(total_size, offset + buf_size)
            fh.seek(-offset, os.SEEK_END)
            buffer = fh.read(min(remaining_size, buf_size))
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # the first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # if the previous chunk starts right from the beginning of line
                # do not concact the segment to the last line of new chunk
                # instead, yield the segment first
                if buffer[-1] is not '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                yield lines[index]
        yield segment


def generate_message(lines):
    """yield a message block."""

    block = []
    for line in lines:
        if not line:
            continue
        block.append(line)
        if line.startswith("From - "):
            yield reversed(block)
            block = []


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-g", "--glob", dest="glob", type="string",
                      help="glob of emails to scan")

    parser.add_option("-m", "--multiple-projects", dest="multiple",
                      type="string",
                      action="append",
                      help="project folder with multiple projects")

    parser.set_defaults(
        glob="/ifs/home/andreas/mail/Local Folders/Projects.sbd/*",
        regex="CGAT-PROJECT-(\d+)@JISCMAIL.AC.UK",
        max_emails=100,
        multiple=[],
        )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    if options.multiple:
        options.multiple = dict([x.split("=") for x in options.multiple])

    filenames = glob.glob(options.glob)

    rx_project = re.compile(options.regex, re.I)
    max_emails = options.max_emails
    outfile = options.stdout
    outfile.write("project_id\treceived\tfilename\tscanned\n")

    for filename in filenames:
        if filename.endswith(".msf"):
            continue

        E.debug("working on %s" % filename)

        basename = os.path.basename(filename)
        if basename in options.multiple:
            expected = set([x.strip() for x in
                            options.multiple[basename].split(":")])
        else:
            expected = set([])

        found = {}

        received = None
        project_id = None

        for n, message in enumerate(
                generate_message(
                    reverse_readline(filename))):

            for line in message:
                if line.startswith("From - "):
                    received = line[7:]
                    continue
                x = rx_project.search(line)
                if x is not None:
                    project_id = x.groups()[0]
                    break
                if line.startswith("MIME-Version"):
                    break

            if project_id is not None:
                found[project_id] = received

            if expected is not None:
                if project_id in expected:
                    expected.remove(project_id)

            if n > max_emails:
                break

            if found and len(expected) == 0:
                break

        for project_id, received in found.items():
            line = "\t".join(map(str, (
                project_id,
                received,
                filename,
                n)))
            # remove ^M characters in fields
            line = re.sub("\r", "", line)
            outfile.write(line + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

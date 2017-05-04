'''
convert_time2seconds.py -
=============================================

:Tags: Python

Purpose
-------

convert output of the time command into seconds:

input: for example something like

1h33m47.32s

Usage
-----

Example::

   python cgat_script_template.py --help

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import re
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

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    for line in sys.stdin:

        if line[0] == "#":
            continue

        data = line[:-1].split("\t")

        new_fields = []

        for field in data:

            x = re.match("^(\d+)h(\d+)m([\d.]+)s", field)
            if not x:
                x = re.match("^(\d+)m([\d.]+)s", field)
                if x:
                    minutes, seconds = x.groups()
                else:
                    new_fields.append(field)
                    continue
                hours = "0"
            else:
                hours, minutes, seconds = x.groups()

            new_fields.append("%.2f" % (
                float(hours) * 3600.0 + float(minutes) * 60.0 + float(seconds)))

        options.stdout.write(("\t".join(new_fields)) + "\n")

    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

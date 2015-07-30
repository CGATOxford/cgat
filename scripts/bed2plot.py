'''
bed.plot.py - create genomic snapshots using the IGV Viewer
===========================================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Create genomic plots in a set of intervals using
the IGV snapshot mechanism.

The script can use a running instance of IGV identified
by host and port. Alternatively, it can start IGV and load
a pre-built session.

Usage
-----

Example::

   python bed2plot.py < in.bed

Type::

   python script_template.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import re
import signal

import CGAT.IGV as IGV
import CGAT.Bed as Bed

import CGAT.Experiment as E


def main(argv=sys.argv):

    # setup command line parser
    parser = E.OptionParser(
        version="%prog version: $Id$", usage=globals()["__doc__"])

    parser.add_option("-n", "--new-instance", dest="new_instance",
                      action="store_true",
                      help="create a new IGV instance [%default]")

    parser.add_option("-s", "--session", dest="session",
                      type="string",
                      help="load session before creating plots "
                      "[%default]")

    parser.add_option("-d", "--snapshot-dir", dest="snapshotdir",
                      type="string",
                      help="directory to save snapshots in [%default]")

    parser.add_option("-f", "--format", dest="format", type="choice",
                      choices=("png", "eps", "svg"),
                      help="output file format [%default]")

    parser.add_option("-o", "--host", dest="host", type="string",
                      help="host that IGV is running on [%default]")

    parser.add_option("-p", "--port", dest="port", type="int",
                      help="port that IGV listens at [%default]")

    parser.add_option("-e", "--extend", dest="extend", type="int",
                      help="extend each interval by a number of bases "
                      "[%default]")

    parser.add_option("-x", "--expand", dest="expand", type="float",
                      help="expand each region by a certain factor "
                      "[%default]")

    parser.add_option("--session-only", dest="session_only",
                      action="store_true",
                      help="plot session after opening, "
                      "ignore intervals "
                      "[%default]")

    parser.add_option("--keep", dest="keep_open",
                      action="store_true",
                      help="keep a newly created IGV session open "
                      "[%default]")

    parser.set_defaults(
        command="igv.sh",
        host='127.0.0.1',
        port=61111,
        snapshotdir=os.getcwd(),
        extend=0,
        format="png",
        expand=1.0,
        session=None,
        session_only=False,
        keep_open=False,
    )

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv, add_output_options=True)

    igv_process = None
    if options.new_instance:
        E.info("starting new IGV process")
        igv_process = IGV.startIGV(command=options.command,
                                   port=options.port)
        E.info("new IGV process started")

    E.info("connection to process on %s:%s" % (options.host, options.port))
    E.info("saving images in %s" % options.snapshotdir)
    igv = IGV.IGV(host=options.host,
                  port=options.port,
                  snapshot_dir=os.path.abspath(options.snapshotdir))

    if options.session:
        E.info('loading session from %s' % options.session)
        igv.load(options.session)
        E.info('loaded session')

    if options.session_only:
        E.info('plotting session only ignoring any intervals')
        fn = "%s.%s" % (os.path.basename(options.session), options.format)
        E.info("writing snapshot to '%s'" %
               os.path.join(options.snapshotdir, fn))
        igv.save(fn)

    else:
        c = E.Counter()
        for bed in Bed.iterator(options.stdin):

            c.input += 1

            # IGV can not deal with white-space in filenames
            name = re.sub("\s", "_", bed.name)

            E.info("going to %s:%i-%i for %s" %
                   (bed.contig, bed.start, bed.end, name))

            start, end = bed.start, bed.end
            extend = options.extend
            if options.expand:
                d = end - start
                extend = max(extend, (options.expand * d - d) // 2)

            start -= extend
            end += extend

            igv.go("%s:%i-%i" % (bed.contig, start, end))

            fn = E.getOutputFile("%s.%s" % (name, options.format))
            E.info("writing snapshot to '%s'" % fn)
            igv.save(fn)

            c.snapshots += 1

        E.info(c)

    if igv_process is not None and not options.keep_open:
        E.info('shutting down IGV')
        igv_process.send_signal(signal.SIGKILL)

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

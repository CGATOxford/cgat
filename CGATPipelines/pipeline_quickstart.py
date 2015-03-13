
'''
pipeline_quickstart.py - setup a new pipeline
=============================================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

Usage
-----

Example::

   python pipeline_quickstart.py --set-name=chipseq

Type::

   python pipeline_quickstart.py --help

for command line help.

Documentation
-------------

Code
----

'''
import sys
import re
import os
import shutil
import CGAT.Experiment as E


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-d", "--dest", dest="destination", type="string",
                      help="destination directory.")

    parser.add_option(
        "-n", "--set-name", dest="name", type="string",
        help="name of this pipeline. 'pipeline_' will be prefixed.")

    parser.add_option(
        "-f", "--force-output", dest="force", action="store_true",
        help="overwrite existing files.")

    parser.add_option(
        "-t", "--pipeline-type", dest="pipeline_type", type="choice",
        choices=("full", "minimal"),
        help="type of pipeline to output. "
        "full=a complete pipeline for the CGAT environment "
        "minimum=minimum pipeline "
        "[%default]")

    parser.set_defaults(
        destination=".",
        name=None,
        force=False,
        pipeline_type="full",
    )

    (options, args) = E.Start(parser)

    if not options.name:
        raise ValueError("please provide a pipeline name")

    reportdir = os.path.abspath("src/pipeline_docs/pipeline_%s" % options.name)
    confdir = os.path.abspath("src/pipeline_%s" % (options.name))

    destination_dir = options.destination

    # create directories
    for d in ("", "src", "report",
              "src/pipeline_docs",
              "src/pipeline_%s" % options.name,
              reportdir,
              "%s/_templates" % reportdir,
              "%s/pipeline" % reportdir,
              "%s/trackers" % reportdir):

        dd = os.path.join(destination_dir, d)
        if not os.path.exists(dd):
            os.makedirs(dd)

    # copy files
    # replaces all instances of template with options.name within
    # filenames and inside files.
    rx_file = re.compile("template")
    rx_type = re.compile("_%s" % options.pipeline_type)
    rx_template = re.compile("@template@")
    rx_reportdir = re.compile("@reportdir@")

    srcdir = os.path.dirname(__file__)

    def copy(src, dst, name):

        # remove "template" and the pipeline type from file/directory
        # names.
        fn_dest = os.path.join(
            destination_dir,
            dst,
            rx_type.sub("", rx_file.sub(name, src)))

        fn_src = os.path.join(srcdir,
                              "pipeline_template_data", src)

        E.debug("fn_src=%s, fn_dest=%s, src=%s, dest=%s" %
                (fn_src, fn_dest, src, dst))

        if os.path.exists(fn_dest) and not options.force:
            raise OSError(
                "file %s already exists - not overwriting." % fn_dest)

        outfile = open(fn_dest, "w")
        infile = open(fn_src)
        for line in infile:
            outfile.write(rx_reportdir.sub(reportdir,
                                           rx_template.sub(name, line)))

        outfile.close()
        infile.close()

    def copytree(src, dst, name):

        fn_dest = os.path.join(destination_dir, dst, rx_file.sub(name, src))
        fn_src = os.path.join(srcdir, "pipeline_template_data", src)

        if os.path.exists(fn_dest) and not options.force:
            raise OSError(
                "file %s already exists - not overwriting." % fn_dest)

        shutil.copytree(fn_src, fn_dest)

    for f in ("conf.py",
              "pipeline.ini"):
        copy(f, 'src/pipeline_%s' % options.name, name=options.name)

    # copy the script
    copy("pipeline_template_%s.py" % options.pipeline_type, 'src',
         name=options.name)

    # create links
    for src, dest in (("conf.py", "conf.py"),
                      ("pipeline.ini", "pipeline.ini")):
        d = os.path.join("report", dest)
        if os.path.exists(d) and options.force:
            os.unlink(d)
        os.symlink(os.path.join(confdir, src), d)

    for f in ("cgat_logo.png",):
        copy(f, "%s/_templates" % reportdir,
             name=options.name)

    for f in ("themes",):
        copytree(f, "src/pipeline_docs",
                 name=options.name)

    for f in ("contents.rst",
              "pipeline.rst",
              "__init__.py"):
        copy(f, reportdir,
             name=options.name)

    for f in ("Dummy.rst",
              "Methods.rst"):
        copy(f, "%s/pipeline" % reportdir,
             name=options.name)

    for f in ("TemplateReport.py", ):
        copy(f, "%s/trackers" % reportdir,
             name=options.name)

    absdest = os.path.abspath(destination_dir)

    name = options.name

    print """
Welcome to your new %(name)s CGAT pipeline.

All files have been successfully copied to `%(destination_dir)s`. In
order to start the pipeline, go to `%(destination_dir)s/report`

   cd %(destination_dir)s/report

You can start the pipeline by typing:

   python ../src/pipeline_%(name)s.py -v 5 -p 5 make full

To build the report, type:

   python ../src/pipeline_%(name)s.py -v 5 -p 5 make build_report

The report will be in file:/%(absdest)s/report/report/html/index.html.

The source code for the pipeline is in %(destination_dir)s/src.

""" % locals()

    E.Stop()

if __name__ == "__main__":
    sys.exit(main())

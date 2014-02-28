import os
import sys
import re
import types
import itertools
import glob
import collections

from SphinxReport.Tracker import *
from SphinxReport.Utils import PARAMS as P
from SphinxReport.odict import OrderedDict as odict

from SphinxReport.ResultBlock import ResultBlock, ResultBlocks
from SphinxReport import Utils

###################################################################
###################################################################
# parameterization

EXPORTDIR = P['testing_exportdir']
DATADIR = P['testing_datadir']
DATABASE = P['testing_backend']

###########################################################################


class TestingTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

LogFileSummary = collections.namedtuple(
    "LogFileSummary", "info debug warning error")


def summarizeLogFile(filename):
    '''summarize a SphinxReport logfile.'''

    info, debug, warning, error = 0, 0, 0, 0
    with open(filename) as f:
        for line in f:
            words = line.split()
            if len(words) < 3:
                continue
            w = words[2]
            if w == "INFO":
                info += 1
            elif w == "DEBUG":
                debug += 1
            elif w == "WARNING":
                warning += 1
            elif w == "ERROR":
                error += 1

    return LogFileSummary._make((info, debug, warning, error))

##############################################################
##############################################################
##############################################################


class PipelineStatus(Status):

    tracks = [x[:-4] for x in glob.glob("*.dir")]

    def testCompletion(self, track):
        '''check if pipeline completed successfully.
        '''

        lines = open(track + ".log").readlines()
        started = "not started"

        if len(lines) < 1:
            return 'FAIL', started

        x = re.search("# job started at ([^-]*) on", lines[1])
        if x:
            started = x.groups()[1]

        x = re.search(
            "# job finished in (\d+) seconds at ([^-]*) -- ", lines[-1])
        if not x:
            return 'FAIL', started
        else:
            return 'PASS', x.groups()[1]

    def testReport(self, track):
        '''check if report completed successfully.
        '''

        lines = open(track + ".log").readlines()
        started = "not started"

        if len(lines) < 1:
            return 'FAIL', started

        x = re.search("# job started at ([^-]*) on", lines[1])
        if x:
            started = x.groups()[1]

        x = re.search(
            "# job finished in (\d+) seconds at ([^-]*) -- ", lines[-1])
        if not x:
            return 'FAIL', started
        else:
            return 'PASS', x.groups()[1]


class ReportTable(TestingTracker):

    tracks = [x[:-4] for x in glob.glob("*.dir")]

    def __call__(self, track):

        try:
            logfileresult = summarizeLogFile(
                os.path.join(track + ".dir", "sphinxreport.log"))
        except IOError:
            return

        report_file = os.path.join(track + ".dir", "report.log")
        fn = os.path.abspath(
            os.path.join(track + ".dir", "report", "html", "pipeline.html"))

        r = odict()
        r["link"] = "`%(track)s <%(fn)s>`_" % locals()
        r["error"] = logfileresult.error
        r["warning"] = logfileresult.warning
        r["info"] = logfileresult.info
        r["debug"] = logfileresult.debug
        return r


class XReportTable(TestingTracker):

    tracks = [x[:-4] for x in glob.glob("*.dir")]

    def __call__(self, track):

        logfileresult = summarizeLogFile(
            os.path.join(track + ".dir", "sphinxreport.log"))
        report_file = os.path.join(track + ".dir", "report.log")

        toc_text = []
        link_text = []

        fn = os.path.join(track + ".dir", "report", "html", "pipeline.html")
        toc_text.append("* %(track)s_" % locals())
        link_text.append(".. _%(track)s: %(fn)s" % locals())

        toc_text = "\n".join(toc_text)
        link_text = "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

        return odict((("text", rst_text),))


class X:

    def __call__(self, track, slice=None):

        edir = EXPORTDIR

        toc_text = []
        link_text = []

        filenames = sorted([x.asFile() for x in TRACKS])

        for fn in filenames:
            if PE == "True":
                fn1 = fn + ".1"
                fn2 = fn + ".2"
                toc_text.append("* %(fn1)s_" % locals())
                toc_text.append("* %(fn2)s_" % locals())
                link_text.append(
                    ".. _%(fn1)s: %(edir)s/fastqc/%(fn1)s_fastqc/fastqc_report.html" % locals())
                link_text.append(
                    ".. _%(fn2)s: %(edir)s/fastqc/%(fn2)s_fastqc/fastqc_report.html" % locals())
            else:
                toc_text.append("* %(fn)s_" % locals())
                link_text.append(
                    ".. _%(fn)s: %(edir)s/fastqc/%(fn)s_fastqc/fastqc_report.html" % locals())

        toc_text = "\n".join(toc_text)
        link_text = "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

        return odict((("text", rst_text),))

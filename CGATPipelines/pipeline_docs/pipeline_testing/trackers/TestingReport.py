import os
import re
import glob
import collections
from collections import OrderedDict as odict

from CGATReport.Tracker import *
from CGATReport.Utils import PARAMS as P

###################################################################
###################################################################
# parameterization

EXPORTDIR = P.get('testing_exportdir', P.get('exportdir', 'export'))
DATADIR = P.get('testing_datadir', P.get('datadir', '.'))
DATABASE = P.get('testing_backend', P.get('sql_backend', 'sqlite:///./csvdb'))


class TestingTracker(TrackerSQL):

    '''Define convenience tracks for plots'''

    def __init__(self, *args, **kwargs):
        TrackerSQL.__init__(self, *args, backend=DATABASE, **kwargs)

LogFileSummary = collections.namedtuple(
    "LogFileSummary", "info debug warning error")


def summarizeLogFile(filename):
    '''summarize a CGATReport logfile.'''

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

        lines = open(track + ".report").readlines()
        started = [x for x in lines if x.startswith("# job started")]
        finished = [x for x in lines if x.startswith("# job finished")]
        error = [x for x in lines if "ERROR" in lines]

        if len(started) == 0:
            return 'WARN', 'never started'

        if len(finished) == 0:
            return 'FAIL', 'started, but never finished'

        if error:
            return 'FAIL', 'report caused errors'
        else:
            return 'PASS', 'report completed'


class ReportTable(TestingTracker):

    tracks = [x[:-4] for x in glob.glob("*.dir")]

    def __call__(self, track):

        try:
            logfileresult = summarizeLogFile(
                os.path.join(track + ".dir", "cgatreport.log"))
        except IOError:
            return

        report_file = os.path.join(track + ".dir", "report.log")
        fn = os.path.abspath(
            os.path.join(track + ".dir", "report", "html", "contents.html"))

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
            os.path.join(track + ".dir", "cgatreport.log"))
        report_file = os.path.join(track + ".dir", "report.log")

        toc_text = []
        link_text = []

        fn = os.path.join(track + ".dir", "report", "html", "contents.html")
        toc_text.append("* %(track)s_" % locals())
        link_text.append(".. _%(track)s: %(fn)s" % locals())

        toc_text = "\n".join(toc_text)
        link_text = "\n".join(link_text)

        rst_text = '''
%(toc_text)s

%(link_text)s
''' % locals()

        return odict((("text", rst_text),))


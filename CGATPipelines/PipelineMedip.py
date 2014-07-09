'''
PipelineMedip.py - tasks associated with MedipSeq analysis
==========================================================

'''

import re
import os
import collections
import sqlite3

import CGAT.Experiment as E
import CGAT.Pipeline as P

import CGAT.Database as Database
import CGAT.IOTools as IOTools

from rpy2.robjects import r as R
import rpy2.robjects as ro

PARAMS = {}


def buildDMRStats(tables, method, outfile):
    '''build dmr summary statistics.

    Creates some diagnostic plots in

    <exportdir>/<method> directory.

    Tables should be labeled <tileset>_<design>_<method>.

    '''

    dbhandle = sqlite3.connect(PARAMS["database"])

    def togeneset(tablename):
        return re.match("([^_]+)_", tablename).groups()[0]

    keys_status = "OK", "NOTEST", "FAIL", "NOCALL"

    outf = IOTools.openFile(outfile, "w")
    outf.write("\t".join(("tileset", "design", "track1", "track2", "tested",
                          "\t".join(["status_%s" % x for x in keys_status]),
                          "significant",
                          "up", "down",
                          "twofold",
                          "twofold_up", "twofold_down",
                          )) + "\n")

    all_tables = set(Database.getTables(dbhandle))
    outdir = os.path.join(PARAMS["exportdir"], "diff_methylation")

    for tablename in tables:

        prefix = P.snip(tablename, "_%s" % method)
        tileset, design = prefix.split("_")

        def toDict(vals, l=2):
            return collections.defaultdict(int, [(tuple(x[:l]), x[l]) for x in vals])

        E.info("collecting data from %s" % tablename)

        tested = toDict(Database.executewait(dbhandle,
                                             """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                GROUP BY treatment_name,control_name""" % locals() ).fetchall() )
        status = toDict(Database.executewait(dbhandle,
                                             """SELECT treatment_name, control_name, status, COUNT(*) FROM %(tablename)s 
                                GROUP BY treatment_name,control_name,status""" % locals() ).fetchall(), 3 )
        signif = toDict(Database.executewait(dbhandle,
                                             """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE significant
                                GROUP BY treatment_name,control_name""" % locals() ).fetchall() )
        fold2 = toDict(Database.executewait(dbhandle,
                                            """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE (l2fold >= 1 or l2fold <= -1) AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )

        up = toDict(Database.executewait(dbhandle,
                                         """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE l2fold > 0 AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )

        down = toDict(Database.executewait(dbhandle,
                                           """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE l2fold < 0 AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )

        fold2up = toDict(Database.executewait(dbhandle,
                                              """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE l2fold > 1 AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )

        fold2down = toDict(Database.executewait(dbhandle,
                                                """SELECT treatment_name, control_name, COUNT(*) FROM %(tablename)s 
                                WHERE l2fold < -1 AND significant
                                GROUP BY treatment_name,control_name,significant""" % locals() ).fetchall() )

        groups = tested.keys()

        for treatment_name, control_name in groups:
            k = (treatment_name, control_name)
            outf.write("\t".join(map(str, (
                tileset,
                design,
                treatment_name,
                control_name,
                tested[k],
                "\t".join([str(status[(treatment_name, control_name, x)])
                          for x in keys_status]),
                signif[(k)],
                up[k], down[k],
                fold2[k],
                fold2up[k], fold2down[k]))) + "\n")

        ###########################################
        ###########################################
        ###########################################
        # plot length versus P-Value
        data = Database.executewait(dbhandle,
                                    '''SELECT end - start, pvalue 
                             FROM %(tablename)s
                             WHERE significant''' % locals() ).fetchall()

        # require at least 10 datapoints - otherwise smooth scatter fails
        if len(data) > 10:
            data = zip(*data)

            pngfile = "%(outdir)s/%(tileset)s_%(design)s_%(method)s_pvalue_vs_length.png" % locals()
            R.png(pngfile)
            R.smoothScatter(R.log10(ro.FloatVector(data[0])),
                            R.log10(ro.FloatVector(data[1])),
                            xlab='log10( length )',
                            ylab='log10( pvalue )',
                            log="x", pch=20, cex=.1)

            R['dev.off']()

    outf.close()

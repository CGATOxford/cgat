from ChipseqReport import *
import glob


class BioProspectorMotifs(ChipseqTracker):

    pattern = "(.*)_bioprospector$"

    def __call__(self, track, slice=None):

        resultsdir = os.path.join(EXPORTDIR, "bioprospector")

        if not os.path.exists(resultsdir):
            return []

        data = odict()

        files = glob.glob(os.path.join(resultsdir, "%s_*" % track))

        # only top 10% of intervals are used
        ninput = self.getValue(
            "SELECT COUNT(*) FROM %(track)s_intervals" % locals()) / 10

        for file in sorted(files):
            motif = re.match(".*_(\d+).png", file).groups()[0]
            nintervals = self.getValue(
                "SELECT COUNT(DISTINCT id) FROM %(track)s_bioprospector WHERE motif = '%(motif)s'" % locals())
            nmotifs = self.getValue(
                "SELECT COUNT(*) FROM %(track)s_bioprospector WHERE motif = '%(motif)s'" % locals())

            # the / at image is necessary - don't understand why.
            # the path is correct in rst output.
            data[motif] = odict((
                ("ninput", ninput),
                ("nmotifs", nmotifs),
                ("nintervals", nintervals),
                ("nintervals / %", "%5.2f" % (100.0 * nintervals / ninput)),
                ("image", '''.. image:: /%s''' % file ),
            ))

        return data

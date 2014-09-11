from GeneSetsReport import *

import CGAT.Pipeline as P
import CGAT.IOTools as IOTools

PARAMS = P.peekParameters(
    ".",
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True)


class HypergeometricStatus(Status):
    '''status information for annotations.
    '''
    pattern = 'hypergeometric_(\S+)_summary$'

    slices = ('SignificantResults',)

    def testSignificantResults(self, track):
        '''
        PASS: Genes have been found in foreground and significant results
              exist.

        WARN: Genes have been found in foreground, but no significant results
              exist.

        FAIL: No genes in foreground sets.

        The value indicates the number of significant results.
        '''

        significant, nforeground = self.getFirstRow(
            """SELECT SUM(significant), SUM(nforeground_mapped) FROM
            hypergeometric_%(track)s_summary""")

        if significant > 0:
            return "PASS", significant
        elif nforeground > 0:
            return "WARN", significant
        else:
            return "FAIL", 0

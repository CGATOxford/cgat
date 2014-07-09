from AnnotationReport import *

import CGAT.Pipeline as P
import CGAT.IOTools as IOTools

PARAMS = P.peekParameters(
    ".",
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True)


class AnnotationStatus(Status):
    '''status information for annotations.


    '''
    tracks = [x for x, y in PARAMS.items() if str(y).endswith(
        (".bed.gz", ".gtf.gz", ".gff.gz",
         ".tsv.gz", ".tsv"))]

    slices = ('AnnotationIsPresent',)

    def testAnnotationIsPresent(self, track):
        '''
        PASS: File exists and is not empty

        FAIL: File exists and is empty (no data except comments)

        NA: File does not exist. This might indicate an error
            or simply that the annotation has not been computed.

        The value indicates the number of lines in the file.
        '''
        
        fn = PARAMS[track]
        if not os.path.exists(fn):
            return ('NA', 0)

        nlines = IOTools.getNumLines(fn)
        if nlines > 0:
            return ('PASS', nlines)
        else:
            return ('FAIL', nlines)

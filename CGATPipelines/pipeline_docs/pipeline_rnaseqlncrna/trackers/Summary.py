from CGATReport.Tracker import *

import CGAT.IOTools as IOTools
import CGAT.GTF as GTF
import numpy as np
import scipy.stats
import collections
import CGATPipelines.PipelineLncRNA as PipelineLncRNA

#################################################
#################################################
#################################################


class Summary(Tracker):

    def getTracks(self):

        return ["lncrna", "lncrna_filtered", "lncrna_final"]

    def getFilename(self, track):

        return track + ".gtf.gz"

    def __call__(self, track, slice=None):

        return odict((("single_exon", PipelineLncRNA.CounterSingleExonGenes(os.path.join("gtfs", track) + ".gtf.gz").count()), ("multi_exon", PipelineLncRNA.CounterMultiExonGenes(os.path.join("gtfs", track) + ".gtf.gz").count())))

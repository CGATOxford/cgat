import pandas as pd
from ReadqcReport import *


class PerExperimentSequenceQuality(TrackerSQL):
    def __call__(self, track):
        statement = ("SELECT * FROM experiment_per_sequence_quality")
        df = self.getDataFrame(statement)
        df = pd.melt(df, id_vars=["Quality", ])

        return df

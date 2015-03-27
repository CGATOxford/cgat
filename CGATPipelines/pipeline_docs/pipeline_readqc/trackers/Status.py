import os
import collections
from CGATReport.Tracker import Status


class FastqcStatus(Status):

    '''status information for mapping stage.'''

    def __init__(self, *args, **kwargs):
        Status.__init__(self, *args, **kwargs)

        data = self.getAll("SELECT * FROM status_summary")
        # reorganize into dictionary with keys [track][quality]
        keys = data.keys()
        vv = collections.defaultdict(dict)
        for row in zip(*[data[x] for x in keys]):
            d = dict(zip(keys, row))
            for s in keys:
                vv[os.path.basename(d['filename'])][s.lower()] = d[s]
        self.data = vv
        self.tracks = vv.keys()

        self.cache = False

    def testBasicStatistics(self, track):
        '''basic statistics.
        '''
        return self.data[track].get(
            'basic_statistics', "NA"), 0

    def testPerTileSequenceQuality(self, track):
        '''per tile sequence quality.
        
        requires fastqc >= 0.11
        '''
        return self.data[track].get(
            'per_base_sequence_quality', "NA"), 0

    def testPerBaseSequenceQuality(self, track):
        '''per base sequence quality
        '''
        return self.data[track].get(
            'per_base_sequence_quality', "NA"), 0

    def testPerSequenceQualityScores(self, track):
        '''per base sequence quality
        '''
        return self.data[track].get(
            'per_sequence_quality_scores', "NA"), 0

    def testPerBaseSequenceContent(self, track):
        '''per base sequence quality
        '''
        return self.data[track].get(
            'per_base_sequence_content', "NA"), 0

    def testPerBaseGCContent(self, track):
        '''per base sequence quality

        not available in fastqc >= 0.11
        '''
        return self.data[track].get(
            'per_base_gc_content', "NA"), 0

    def testPerSequenceGCContent(self, track):
        '''per base sequence quality
        '''
        return self.data[track].get(
            'per_sequence_gc_content', "NA"), 0

    def testPerBaseNContent(self, track):
        '''per base sequence quality
        '''
        return self.data[track].get(
            'per_base_n_content', "NA"), 0

    def testSequenceLengthDistribution(self, track):
        '''Sequence length distribution
        '''
        return self.data[track].get(
            'sequence_length_distribution', "NA"), 0

    def testSequenceDuplicationLevels(self, track):
        '''Sequence duplication levels
        '''
        return self.data[track].get(
            'sequence_duplication_levels', "NA"), 0

    def testOverrepresentedSequences(self, track):
        '''Overrepresented sequences
        '''
        return self.data[track].get(
            'overrepresented_sequences', "NA"), 0

    def testAdapterContent(self, track):
        '''Adapter content

        requires fastqc >= 0.11
        '''
        return self.data[track].get(
            'adapter_content', "NA"), 0

    def testKmerContent(self, track):
        '''Kmer content
        '''
        return self.data[track].get(
            'kmer_content', "NA"), 0

import CGAT.IOTools as IOTools
import CGAT.Pipeline as Pipeline
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import ConfigParser
import collections
import pysam
from bx.bbi.bigwig_file import BigWigFile
import numpy as np
import os
import tracks
import layout_code

'''
TODO:
* write unit tests

'''


def parse_config(config_file, allow_missing=True):
    '''
    Parse config file into a params dictionary.  Use Pipeline.py
    config file parsing code
    '''
    # use ConfigParser to parse config file into
    # an collections.orderedDict (default data structure)
    config = ConfigParser.ConfigParser(allow_no_value="allow_missing")

    open_config = IOTools.openFile(config_file, "r")
    config.readfp(open_config)
    params = Pipeline.configToDictionary(config)
    open_config.close()

    return params


class Interface(object):
    '''
    CGATViewer interface object.  Parses input files into a params
    dictionary and instantiates Layout and Track objects.

    Interface parameters can be defined explicitly or in a config
    file.

    Data from input files are parsed and passed by value to each Track
    instance as required as a numpy array of values.

    Interval files are likewise parsed using GTF.py or BED.py
    and passed as an array of tuples with (start, end).
    '''

    def __init__(self):
        self.params = None
        self.genome_view = None
        self.n_tracks = None
        self.data_container = collections.OrderedDict()
        self.intervals_container = collections.OrderedDict()

    def set_n_tracks(self, n_tracks):
        self.n_tracks = n_tracks

    def get_n_tracks(self):
        return self.n_tracks

    def set_genome_view(self):
        self.genome_view = self.parsed_genome_view(self.params)

    def get_genome_view(self):
        return self.genome_view

    def set_params(self, params):
        self.params = params

    def get_params(self):
        return self.params

    def parsed_genome_view(self, params):
        '''
        Returns the genome view as a tuple:
        (contig, start, end)
        '''

        contig = params['genome_contig']
        start = int(params['genome_start'])
        end = int(params['genome_end'])
        parsed_view = (contig, start, end)

        return parsed_view

    def get_data(self, key):
        return self.data_container[key]

    def get_ranges(self, key):
        return self.intervals_container[key]

    def set_data_vector(self, file_path, data_format):
        '''
        return a numpy ndarray of data values over
        the genome-view interval.
        Currently handles BAM and BigWig format
        TODO: allow data from BEDn  and wiggle formats
        '''

        parsed_view = self.get_genome_view()

        if data_format == "BAM":
            # use pysam to parse BAM file
            # index reads over a region, rather than iterating
            # check for .bai file to prevent need to iterate over file

            samfile = pysam.Samfile(file_path, "rb")

            try:
                assert os.path.exists(file_path + ".bai")
            except AssertionError:
                print "No index file for this BAM file. "
                print "Creating %s.bai. This may take some time" % file_path
                os.system("samtools index %s" % file_path)

            piled = samfile.pileup(reference=parsed_view[0],
                                   start=parsed_view[1],
                                   end=parsed_view[2])

            scores = np.array(val for val in piled)
            data_type = self.get_params()['data_format']
            try:
                self.data_container[data_type].append(scores)
            except KeyError:
                self.data_container[data_type] = scores

        elif data_format == "BigWig":
            # use biopyton to pull out data with indexed file
            file_handle = IOTools.openFile(file_path, "rb")
            bw_index = BigWigFile(file_handle)
            scores = bw_index.get_as_array(parsed_view[0],
                                           parsed_view[1],
                                           parsed_view[2])
            scores[np.isnan(scores)] = 0
            file_handle.close()

            self.data_container[self.get_genome_view()] = scores

        else:
            raise ValueError("data format not recognised\n"
                             "please input data as either BAM, "
                             "wiggle or BigWig")
            self.data_container[self.get_genome_view()] = scores

    def col_iter(iterable):
        '''
        iteratre over a pysam.IteratorColumn object
        return number of reads at each nt
        From bam2.wiggle.py
        '''

        start = None
        end = 0
        n = None
        for t in iterable:
            if (t.pos - end) > 1 or n != t.n:
                if start:
                    yield n
                start = t.pos
                end = t.pos
                n = t.n
            end = t.pos
        yield n

    def set_intervals_vector(self, file_path, file_format):
        '''
        Parse a .gtf, .gff or .bed file into a vector of
        start, end co-ordinates within the genome view specified
        '''

        # index the genome view in the intervals file and iterate
        # over each interval adding to the intervals container
        # currently only store contig, start and end co-ordinates

        file_handle = IOTools.openFile(file_path, "rb")
        interval_set = set()
        file_format = file_format.lower()
        if file_format == "gtf" or file_format == "gff":
            gtf_idx = GTF.readAndIndex(GTF.iterator(file_handle))
            intersects = gtf_idx.get(self.get_genome_view()[0],
                                     self.get_genome_view()[1],
                                     self.get_genome_view()[2])
            # GTF.readAndIndex returns and IndexedGenome object, when indexed
            # into this returns a nested containment list object which is an
            # iterator of tuples containing the start and end coordinates and a
            # pysam.TabProxies object this last object contains the GTF entries
            # that intersect the genome view
            # TODO: incorporate gene/transcript/feature ids into plotting

            for interval in intersects:
                entry = interval[-1]
                start = entry.start
                end = entry.end
                coord = (start, end)
                interval_set.add(coord)

        elif file_format == "bed":
            # Bed.readAndIndex returns a dictionary of NCLSimple objects
            # with contigs as keys and nested containment list object as values
            # use NCL.find(start, end) to pull out genome view region iterator
            # each object in iterator is a tuple(start, end, index)
            bed_idx = Bed.readAndIndex(file_handle)
            ncl = bed_idx[self.get_genome_view()[0]]
            overlaps = ncl.find(self.get_genome_view()[1],
                                self.get_genome_view()[2])
            for interval in overlaps:
                coord = (interval[0], interval[1])
                interval_set.add(coord)

        interval_gen = [x for x in interval_set]

        try:
            inter_type = self.get_params()['intervals_type']
            self.intervals_container[inter_type].append(interval_gen)
        except KeyError:
            self.intervals_container[inter_type] = interval_gen

    def makeTrack(self, track_type, genome_view, vector, layout, priority):
        '''
        instantiate a track.  Values passed are the track type,
        either 'ranges' or 'data'.
        The vector is either an array of values or a container of tuples
        with interval start and end positions.
        '''

        # instantiate track
        if track_type == "ranges":
            ttype = "r"
        elif track_type == "data":
            ttype = "b"
        else:
            raise ValueError("Unsupported track type")

        inst = tracks.regions(ttype=ttype,
                              trackinfo=vector,
                              L=layout)

        inst.set_Priority(new_priority=priority)

    def makeLayout(self, genome_view):
        '''
        Instantiate the Layout class
        '''

        # add in aesthetics parameters from the config file
        layout = layout_code.layout(region_start=genome_view[1],
                                    region_finish=genome_view[2],
                                    title=self.get_params()['title'])

        return layout

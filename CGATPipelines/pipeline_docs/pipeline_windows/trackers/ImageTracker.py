from CGATReport.Tracker import Tracker
import glob
from collections import OrderedDict as odict


class SortedTrackerImages(Tracker):
    '''
    Collect image files, sort them lexicographically
    by filename and arrange them in a gallery.
    '''

    def __init__(self, *args, **kwargs):
        Tracker.__init__(self, *args, **kwargs)
        if "glob" not in kwargs:
            raise ValueError("TrackerImages requires a:glob: parameter")
        self.glob = kwargs["glob"]

    def getTracks(self, subset=None):
        names = glob.glob(self.glob)
        names.sort()
        return names

    def __call__(self, track, **kwargs):
        """return a data structure for track:param: track and slice:slice:"""

        return odict((('name', track), ('filename', track)))

    def sort_names(self, track_dict, **kwargs):
        '''sort globbed images by name'''

        track_items = track_dict.items()
        sort_track = sorted(track_items, key=lambda t: t[0])
        sort_dict = odict(sort_track)

        return sort_dict


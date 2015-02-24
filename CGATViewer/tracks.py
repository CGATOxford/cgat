import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import layout_code


# L = layout_code.layout(0,100)


class regions(object):
    def __init__(self,
                 ttype="t",
                 trackinfo=None,
                 L=None,
                 name="default",
                 priority="default"):
        '''
        Specify the start and end of the region to be viewed.
        ttypes (for now) :

        name = name of track e.g. one of these:
        ["gene_model","read_possitions",
        "read_count", "splice_varient", "default"]

        priority = number set by user in config file to determine ordering of
        graphs on page - default will be sorted in accordance to possition in
        list of prioritys by layout class
        ttype = "r": genomic ranges with start and end positions
        trackinfo is a list of tuples representing the ranges
        e.g.
        r = regions(ttype="r", trackinfo=[(0,10), (5,30), (50,80)])

        ttype = "b": bar chart - list of heights each for one base
        trackinfo is a list of ints representing bar heights (1 per position)
        e.g.
        r = regions(ttype="b", trackinfo=[2, 3, 0, 4, 5, 6, 0, 2 ,3, 6])
        L is a layout object
        '''
        # start of the region to view
        self.start = L.getStart()
        # end of the region to view
        self.end = L.getFinish()
        self.name = name
        self.priority = priority
        self.ttype = ttype  # type of chart
        self.trackinfo = trackinfo
        L.addTracks(self)

    def set_Priority(self, new_priority):
        self.priority = new_priority

    def get_Priority(self):
        return self.priority

    def get_Name(self):
        return self.name

    def draw_regions(self):
        ''' makes a collection of matplotlib patches
        patches are rectangles representing regions '''
        ttype = self.ttype
        start = self.start
        end = self.end
        trackinfo = self.trackinfo
        if ttype == "b":  # generates pos_pairs for barchart
            pos_pairs = zip(range(start, end), range(start+1, end+1))
            heights = trackinfo
        elif ttype == "r":  # generates heights for range chart
            pos_pairs = trackinfo
            # MM: pos_pairs is passed as a generator object by interface
            # this is to limit memory usage for very large plots
            # will just use a list for now
            heights = [1] * len(pos_pairs)
        patches = []  # list of rectangles
        i = 0
        rows = []  # list of lists representing rows in tracks chart where 0
        # is an empty position and 1 an occupied position
        self.new_row(rows)
        rn = 0
        for pair in pos_pairs:
            sstart = pair[0]  # start of range
            eend = pair[1]  # end of range
            height = heights[i]  # height of range
            if ttype == "r":
                rn = self.check_rowocc(rows, (sstart, eend))
                # finds a row in the chart with space for this rectangle
                if rn is False:  # if there is no free space
                    self.new_row(rows)  # add a new row and put it in there
                    rn = len(rows) - 1
                self.occupy_row(rows, rn, (sstart, eend))
                # prevent more rectangles being added in the same spot
            r = mpatches.Rectangle([sstart, rn], (eend-sstart), height)
            # make a rectangle with corners at:
            # (sstart,0), (eend, 0), (sstart, height), (eend,height)
            if ttype == "r":
                max_height = len(rows)
                #  ranges chart height depends on number of overlapping ranges
            elif ttype == "b":
                max_height = max(heights)
                #  bar chart height depends on tallest bar
            self.max_height = max_height
            patches.append(r)  # store the patch in the patch list
            i += 1
        p = PatchCollection(patches, match_original=True)
        self.patches = p
        # convert rectangles to a patch collection
        return p

    def check_rowocc(self, rows, positions):
        '''
        Used to build charts showing genomic ranges.
        Tests if any row in the chart is free to display a rectangle from
        positions[0] to positions[1]
        '''
        sites = range(positions[0], positions[1]+1)
        rn = 0
        for row in rows:
            x = 0
            for site in sites:
                if row[site] == 1:
                    x += 1
            if x == 0:
                return rn  # returns the row number of the free row
            rn += 1
        return False  # returns False if there is no space

    def new_row(self, rows):
        '''
        Used to build charts showing genomic ranges.
        Adds a new row of 0s to "rows"
        '''
        rowlen = self.end - self.start + 1
        row = [0] * rowlen
        rows.append(row)

    def occupy_row(self, rows, rn, positions):
        '''
        Used to build charts showing genomic ranges.
        Replaces 0s with 1s in "rows" row rn from positions[0] to positions[1]
        represents a new rectangle which has been added
        '''
        sites = range(positions[0], positions[1]+1)
        for site in sites:
            rows[rn][site] = 1

    def get_plotheight(self):
        '''
        returns the highest y value needed to display the plot
        '''
        ttype = self.ttype
        if ttype == "r":
            max_height = self.get_coverage()
        elif ttype == "b":
            max_height = max(self.trackinfo)
        return max_height

    def get_coverage(self):
        ''' counts how many tracks cover each site.  returns the maximum'''
        start = self.start
        end = self.end
        pos_pairs = self.trackinfo
        coverage = []
        for p in range(start, end+1):
            count = 0
            for pair in pos_pairs:
                if p >= pair[0] and p <= pair[1]:
                    count += 1
            coverage.append(count)
        return max(coverage)

    def get_patches(self):
        '''
        returns a matplotlib PatchCollection object representing the data
        '''
        return self.patches

    def get_start(self):
        '''
        returns the start position of the genome region to be displayed
        '''
        return self.start

    def get_end(self):
        '''
        returns the end position of the genome region to be displayed
        '''
        return self.end

    def get_type(self):
        '''
        returns the chart type requested
        currently "r" is genomic ragnes, "b" is bar chart
        '''
        if self.ttype == "r":
            return "genomic ranges"
        elif self.ttype == "b":
            return "bar chart"

# #####Testing Only########
# r = regions(ttype="r", trackinfo=[(0,25), (10, 20), (20,60), (30, 60), (50,80), (90,100)])
# r.draw_regions()
# fig, ax = plt.subplots(figsize=(10, r.get_plotheight()))  # initiate figure and axis
# plt.axis([0, 100, 0, r.get_plotheight()])
# ax.add_collection(r.get_patches())  # put rectangles on axes

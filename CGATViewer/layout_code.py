
import matplotlib.pyplot as plt

'''
Need:
    region plotting to determine x axis - region_start, region_finish
    types graph to plot i.e list of tracks
    number of graphs to determine number of subplot
    size of page to plot on
    data plotting
    title
'''


class layout(object):
    """ class to co-ordinate plotting of sevaral tracks on one layout,
    as defined by user interface

    instantiates layout plot for particular region as speicifed by user
    interface

    list_of_tracks = list of track objects (i.e graph types) that need
    to be ploted
    title = title of layout, default = region layout depicts
    region_start = int representing start of x axis for each plot
    region_finish = int representing end of x axis for each plot
    page_size = size of layout"""

    def __init__(self,
                 region_start,
                 region_finish,
                 title=0,
                 page_size="a4",
                 default_graph_priority=["gene_model",
                                         "read_possitions",
                                         "read_count",
                                         "splice_varient",
                                         "default"]):

        if title == 0:
            start_str = str(region_start)
            finish_str = str(region_finish)
            self.title = "Region " + start_str + "to " + finish_str
        else:
            self.title = title

        # create empty list to collect track objects in
        self.list_of_track_objects = []

        # region over which to create plot
        self.region_start = region_start
        self.region_finish = region_finish

        # know page size
        self.page_size = page_size

        # default priorities for graph subplots using ordered list
        self.default_graph_priority = default_graph_priority

    def getStart(self):
        return self.region_start

    def getFinish(self):
        return self.region_finish

    def getTracks(self):
        return self.list_of_track_objects

    def getTitle(self):
        return self.title

    def addTracks(self, new_track_obj):
        """ adds track object to to list_of_tracks to be ploted,
        called by track plotter """
        if new_track_obj not in self.list_of_track_objects:
            self.list_of_track_objects.append(new_track_obj)
        else:
            pass  # Should i raise an error here?

        return self.list_of_track_objects

    def generatePlot(self):
        """tell plotter how much plot has of page?
        define order of tracks
        return - write plot to file """

        # number of tracks to put on subplot
        number_of_tracks = len(self.list_of_track_objects)

        # create plot with shared x axis
        # figure = plt.subplot(number_of_tracks,1,1, sharex=True)

        fig, axarr = plt.subplots(number_of_tracks, sharex=True)
        axarr[0].set_title(self.getTitle())
        for number in range(number_of_tracks):
            axarr[number].set_xlim([self.region_start,
                                    self.region_finish])

        # empty lists to hold priority order
        graphs_with_predef_priority = []
        graphs_no_priority = []

        # need to define order of tracks by position?
        # need to store user-requested priority as a track.attribute..

        for track in self.list_of_track_objects:
            """graphs with predetermined priority will be ploted first,
            then those with default priorities will be plotted second"""

            if track.get_Priority != "default":
                # find prederemined set of priorites
                graphs_with_predef_priority.append(track)
            else:
                # if no priority set
                gpriority = self.default_graph_priority.index[track.get_Name()]
                track.set_Priority(gpriority)
                graphs_no_priority.append(track)

        # sort list of proritised tracks for plotting onto subplot
        # NEED TO TEST THIS?
        sorted_priority_list = sorted(graphs_with_predef_priority,
                                      key=lambda track: track.get_Priority)
        sorted_nopriority_list = sorted(graphs_no_priority,
                                        key=lambda track: track.get_Priority)

        # append tracks into final ordered list which
        # dictates order they need to be plotted in
        final_sorted_figure_list = sorted_priority_list
        final_sorted_figure_list.extend(sorted_nopriority_list)

        # tracker for number of tracks that have been ploted
        subplot_no = 0

        # plot subplots
        for track in sorted_priority_list:
            # ax = ax[subplot_no,1].plot(number_of_tracks,1,subplot_no)
            axarr[subplot_no].set_ylim([0, track.get_plotheight()])
            # I think i need to call the draw otherwise they don't get made?
            axarr[subplot_no].add_collection(track.draw_regions())
            axarr[subplot_no].set_ylabel(track.get_Name())
            subplot_no += 1

        return fig

        """
        while number < self.number_of_tracks:
            for track in self.list_of_track_objects:
                # add subplot to Figure
                # this needs to be updated somehow with number so subplot
                # in right possition
                self.figure.add_subplot(number_of_tracks,1,number)
                number += 1
        """

    def showPlot(self):
        """ show plot on screen or write to different files"""
        # return plt.show(self.figure, block=False)
        # block = lets return to command line when figures are open
        # in non-interative mode
        return plt.show(self.generatePlot(), block=False)

    def save_plot(self,
                  file_name,
                  dpi=None,
                  facecolor="w",
                  edgecolor="w",
                  orientation="portrait",
                  papertype=None,
                  format=None):
        """ file_name= string containing file path to save file to"""

        return plt.savefig(file_name,
                           dpi,
                           facecolor,
                           edgecolor,
                           orientation,
                           papertype,
                           format)

"""
    def removeTrack(self, unwanted_track_objects):
    # removes track from subplot
    # returns subplot without track

        for track in unwanted_track_objects:
            if track in self.list_of_track_objects:
                self.list_of_track_objects.remove(track)

        return self.list_of_track_objects
# test
# 1. create track objects

# A) CREATE HISTOGRAm

# B) create SCATTER PLOT

# CALL TO LAYout to create subplot field
'''Things that need to be tested'''

# That all the tracks are for the right regions have same start-end regions
# That order of priority works
"""

'''gff2plot.py - plot genomic data
===============================

:Tags: Python

Purpose
-------

This script creates plots from genome wide data supplied in :term:`gff` format.

Usage
-----

The data is given in gff files denoting windows. The score field
in the gff file denotes the value of a window.

Multiple tracks can be displayed in the same plot or in several plots
underneath each other.

Data can be displayed as a histogram with various styles or
as a heat-map.

The layout of the plot is given by a config file. The config
file is arranged in plots that contain one more data tracks.

Here is an example of a track file::

    ################################################
    ## Start of track file
    ################################################
    ## tracks
    [paralogs-ds]
    filename=paralogs_ds_median_score.gff
    style=matrix

    [orthologs-ds]
    filename=orthologs_ds_median_score.gff
    style=matrix

    [variable-ds]
    filename=variable_ds_median_score.gff
    style=matrix

    [median<intronic>]
    filename=regions_intronic_median_score.gff

    [<intronic>]
    filename=regions_intronic_mean_score.gff

    ## add multi-track tracks
    [median ds]
    tracks=paralogs-ds, orthologs-ds, variable-ds
    text=dS computed from paralog trees using a
            single omega for each ortholog group.
            For paralogs, the distance is between
            C. elegans and C. briggsae. For
            orthologs, the distance is for the
            C. elegans terminal lineage, only.


    [intron-lengths]
    tracks=median<intronic>,<intronic>

    # add vertical lines to the graph
    [vlines]
    filename=boxes.gff
    color=k
    linewidth=3

    # global options for the figure
    [figure]
    # size of a figure in inches
    figsize=12,12

    # create a plot for the legend
    [legend]
    # size of the legend in inches
    figsize=12,12

    ################################################
    ## End of track file
    ################################################

Example::

   python gff2plot.py --help

Type::

   python gff2plot.py --help

for command line help.

Command line options
--------------------

'''
import sys
import os
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

import matplotlib
import pylab
import matplotlib.ticker
import numpy

import CGAT.Experiment as E
import CGAT.IOTools as IOTools


def formatGenomicCoordinate(value, pos=None):

    decimals = 0
    exps = ['', 'k', 'm', 'g', 't']
    suffix = ""

    value = int(value)
    format = "%%.%if%%s%s" % (int(decimals), suffix)

    for exp in range(len(exps) - 1, -1, -1):

        if value < 1000 ** (exp):
            continue
        else:
            return format % (float(value) / 1000 ** (exp), exps[exp])
            break
    else:
        return str(value)


def enterParams(o, params):
    """enter params from object o if they exists.

    Supply a tuple of name, conversion function for non-string
    options.
    """
    r = {}
    for p in params:
        if isinstance(p, tuple):
            p, f = p
        else:
            f = str
        if hasattr(o, p):
            r[p] = f(getattr(o, p))
    return r


def normalizeValuesByWindows(data, window_size=None):
    """normalize to constant window size.

    If no window size is given, the smallest window is used. Windows
    smaller than windows-size are ignored.
    """

    if not options.window_size:
        window_size = min([x[1] - x[0] for x in data])

    new_values = []

    for start, end, value in data:
        if end - start < window_size:
            continue
        start = start - start % window_size
        for z in range(start, end, window_size):
            new_values.append((z, value))

    new_values.sort()

    # interpolate values for the same windows with average
    xvals = []
    yvals = []
    last_x = None
    values = []
    for x, value in new_values:

        if last_x != x:
            if last_x is not None:
                xvals.append(last_x)
                yvals.append(numpy.mean(values))
            values = []

        last_x = x
        values.append(value)

    return xvals, yvals


def addPlot(ax, track, contig, nplotted,
            nsubplotted=None,
            nsubplots=None,
            y_min=None,
            y_max=None):
    """add a track to an axes object.
    """

    if contig not in track.mData:
        return

    # step function
    # datapoint is in window average
    xvals = [(x[1] + x[0]) / 2.0 for x in track.mData[contig]]
    yvals = [x[2] for x in track.mData[contig]]
    l = len(xvals)

    if nsubplots:
        plotnum = nsubplotted
    else:
        plotnum = nplotted

    if track.style == "matrix":

        # unequal window sizes confuse the image. Therefore, normalize
        # xvals and yvals to a constant image size
        matrix = pylab.array(yvals)
        matrix.shape = (1, len(yvals))

        # make sure that the extent of the image and the plot coincide by
        # using extent
        if nsubplots is not None:
            y_width = float(y_max - y_min) / nsubplots
            extent = (min(xvals), max(xvals), y_min + y_width *
                      nsubplotted, y_min + y_width * (nsubplotted + 1))
        else:
            extent = (min(xvals), max(xvals), min(yvals), max(yvals))
        ax.imshow(matrix,
                  cmap=track.color_scheme,
                  extent=extent,
                  interpolation="nearest")

        symbol = options.symbols[plotnum % len(options.symbols)]
        plot = ax.plot(xvals, yvals, symbol, lw=2)
    else:
        symbol = options.symbols[plotnum % len(options.symbols)]
        plot = ax.plot(xvals, yvals, symbol)

    return plot


def plotContig(contig, tracks, options, plot_legend=False,
               extra_features=None):
    """plot data for contig."""

    if extra_features and "figure" in extra_features:
        figure = pylab.figure(**enterParams(extra_features['figure'],
                                            (("figsize", lambda x: list(map(int, x.split(",")))),
                                             ("dpi", int),
                                             "facecolor",
                                             "edgecolor")))
    else:
        figure = pylab.figure()

    if plot_legend:
        if extra_features and "legend" in extra_features:
            legend = pylab.figure(**enterParams(extra_features['legend'],
                                                (("figsize", lambda x: list(map(int, x.split(",")))),
                                                 ("dpi", int),
                                                 "facecolor",
                                                 "edgecolor")))
        else:
            legend = pylab.figure()
        lx = legend.add_axes((0.1, 0.1, 0.9, 0.9))
        lx.set_title("Legend")
        lx.set_axis_off()
    else:
        legend = None

    axes = []

    ywidth = 0.8 / float(len(tracks))
    yoffset = 0.1

    axprops = {}
    ayprops = {}

    min_x, max_x = 1000000000, 0

    for track in tracks:
        if track.mData:
            min_x = min(min_x, min([(x[0]) for x in track.mData[contig]]))
            max_x = max(max_x, max([(x[1]) for x in track.mData[contig]]))

    # make sure that we use the same view for all axes
    axprops['xlim'] = (min_x, max_x)

    nplotted = 0
    for track in tracks:

        labels, plots = [], []

        ax = figure.add_axes(
            (0.1, track.mYOffset, 0.8, track.mYWidth), **axprops)

        if 'sharex' not in axprops:
            ax.xaxis.set_major_formatter(
                matplotlib.ticker.FuncFormatter(formatGenomicCoordinate))
            ax.set_xlabel("Genomic position / Mb")
            axprops['sharex'] = ax
        else:
            pylab.setp(ax.get_xticklabels(), visible=False)

        ax.set_ylabel(track.mTitle, **ayprops)

        if track.mSubTracks:

            # compute maximum extent of y-axis in all of subtracks
            first = True
            for tt in track.mSubTracks:
                if first:
                    min_y = min([x[2] for x in tt.mData[contig]])
                    max_y = max([x[2] for x in tt.mData[contig]])
                    first = False
                else:
                    min_y = min(
                        min_y, min([x[2] for x in tt.mData[contig]]))
                    max_y = max(
                        max_y, max([x[2] for x in tt.mData[contig]]))

            nsubplotted = 0
            for tt in track.mSubTracks:
                plot = addPlot(ax, tt, contig, nplotted,
                               nsubplotted, len(track.mSubTracks),
                               min_y, max_y)
                nsubplotted += 1
                plots.append(plot)
                if hasattr(tt, "legend"):
                    labels.append(tt.legend)
                else:
                    labels.append(tt.mTitle)

        else:
            min_y = min([x[2] for x in track.mData[contig]])
            max_y = max([x[2] for x in track.mData[contig]])

            if options.global_colours:
                n_for_colour = nplotted
            else:
                n_for_colour = 0

            plot = addPlot(ax, track, contig, n_for_colour)
            plots.append(plot)
            if hasattr(track, "legend"):
                lables.append(track.legend)
            else:
                labels.append(track.mTitle)

        # reduce number of ticks by 2
        old_ticks = ax.get_yticks()
        step_size = (old_ticks[1] - old_ticks[0]) * 2
        new_ticks = list(pylab.arange(old_ticks[0], old_ticks[-1], step_size))
        ax.set_yticks(new_ticks)

        if nplotted % 2 == 0:
            ax.yaxis.set_ticks_position("right")
            ax.yaxis.set_label_position("right")
        else:
            ax.yaxis.set_ticks_position("left")
            ax.yaxis.set_label_position("left")

        # deal with extra_features
        if extra_features:
            for key, config in list(extra_features.items()):
                if key == "vlines":
                    if contig not in config.mData:
                        continue
                    lines = []
                    for start, end, value in config.mData[contig]:
                        lines.append(start)
                        lines.append(end)
                    ax.vlines(
                        lines, min_y, max_y, **enterParams(config, ("colour", "linewidth")))

        nplotted += 1

        if legend:
            lx = legend.add_axes((0.1, track.mYOffset, 0.8, track.mYWidth))
            pylab.setp(lx.get_xticklabels(), visible=False)
            lx.set_xticks([])
            lx.set_yticks([])
            lx.text(0.4, 0.5, track.mTitle)
            lx.legend(plots, labels, 'center left')
            if hasattr(track, "text"):
                lx.text(0.6, 0.2, track.text, size="smaller",
                        clip_on=True)

    ax.set_title(contig)
    # has to be set at the end, otherwise re-scaled?
    ax.set_xlim(min_x, max_x)

    return figure, legend


def readData(infile):
    """read data from infile."""
    dd = {}
    for line in infile:
        if line[0] == "#":
            continue
        d = line[:-1].split("\t")
        contig, start, end, score = d[0], int(d[3]), int(d[4]), float(d[5])
        # if contig != "I": continue
        if contig not in dd:
            dd[contig] = []
        dd[contig].append((start, end, score))
    return dd


class Track:

    def __init__(self, title,
                 priority=0,
                 data=None,
                 subtracks=None,
                 config=None):
        self.mTitle = title
        self.mData = data
        self.mSubTracks = subtracks
        self.height = 1.0
        self.priority = priority
        self.style = "default"
        self.color_scheme = None

        if config:
            for key, value in config.items(section):
                setattr(self, key, value)


def layoutTracks(tracks):
    """layout tracks."""

    # combine subtracks - these are ordered by appearance
    for track in list(tracks.values()):
        if track.mSubTracks:
            s = []
            for subtrack in track.mSubTracks:
                s.append(tracks[subtrack])
                del tracks[subtrack]
            track.mSubTracks = s

    # sort by priority and compute heights
    total_width = sum([x.height for x in list(tracks.values())])

    # convert to list - this is the order in which the tracks will
    # be output
    tracks = list(tracks.values())
    tracks.sort(lambda x, y: cmp(x.priority, y.priority))

    yoffset = 0.1
    ymax_width = 0.8

    for track in tracks:
        track.mYOffset = yoffset
        width = ymax_width * track.height / total_width
        track.mYWidth = width
        yoffset += width

    return tracks


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    parser = E.OptionParser(
        version="%prog version: $Id: gff2plot.py 2781 2009-09-10 11:33:14Z andreas $", usage=globals()["__doc__"])

    parser.add_option("-f", "--file", dest="filenames", type="string",
                      help="files[s] to take data from,stdin = -.")
    parser.add_option("", "--symbols", dest="symbols", type="string",
                      help="symbols to use for each histogram [steps|...].")
    parser.add_option("--slide-show", dest="slide_show", type="choice",
                      choices=("first", "all", "sequence"),
                      help="do a slide show - otherwise, write image to file.")
    parser.add_option("--config", dest="filename_config", type="string",
                      help="filename of track configuration file.")
    parser.add_option("--dpi", dest="dpi", type="int",
                      help="dpi for hardcopy output.")
    parser.add_option("--window-size", dest="window_size", type="int",
                      help="window-size.")
    parser.add_option("--output-filename-pattern", dest="output_pattern_image", type="string",
                      help="output pattern for images. Should contain a '%(contig)s' pattern .")
    parser.add_option("--global-colours", dest="global_colours", action="store_true",
                      help="cycle through colours for all tracks.")

    parser.set_defaults(
        filenames=None,
        symbols="k-,b-,r-,c-,m-,y-,g-",
        output_pattern_image="%(contig)s.png",
        slide_show=None,
        window_size=None,
        filename_config=None,
        dpi=None,
        global_colours=False,
    )

    (options, args) = E.Start(parser)
    options.symbols = options.symbols.split(",")

    #--------------------------------------------------------
    # collect all the data
    # list of data per source and contig
    tracks = {}
    extra_features = {}

    if options.filenames:
        options.filenames = options.filenames.split(",")

        if len(args) > 0:
            options.filenames = args

    if options.filenames:

        for filename in options.filenames:

            if filename == "-":
                infile = sys.stdin
            else:
                infile = IOTools.openFile(filename)

            data = readData(infile)

            if filename != "-":
                infile.close()

            track[filename] = Track(title=filename, data=data)

    elif options.filename_config:
        # get track information from config file
        config = configparser.ConfigParser()
        config.read(os.path.expanduser(options.filename_config))

        # first extract special sections
        for section in config.sections():
            if section == "vlines":
                infile = IOTools.openFile(config.get(section, "filename"), "r")
                data = readData(infile)
                infile.close()
                extra_features[section] = Track(title=section,
                                                data=data,
                                                config=config)
                config.remove_section(section)
            elif section in ("figure", "legend"):
                extra_features[section] = Track(title=section,
                                                data=None,
                                                config=config)
                config.remove_section(section)
        n = 0
        for section in config.sections():

            if config.has_option(section, "filename"):
                infile = IOTools.openFile(config.get(section, "filename"), "r")
                data = readData(infile)
                infile.close()

                tracks[section] = Track(title=section,
                                        data=data,
                                        priority=n,
                                        config=config)

            elif config.has_option(section, "tracks"):
                subtracks = config.get(section, "tracks")
                subtracks = [x.strip() for x in subtracks.split(",")]

                tracks[section] = Track(title=section,
                                        data=None,
                                        config=config,
                                        priority=n,
                                        subtracks=subtracks)
            n += 1

    # compile set of all contigs
    contigs = set()
    for track in list(tracks.values()):
        if track.mData:
            contigs = contigs.union(list(track.mData.keys()))

    # re-arrange tracks and subtracks
    tracks = layoutTracks(tracks)

    nplots = 0
    figures = []
    legend = None
    for contig in contigs:
        figure, l = plotContig(contig, tracks, options,
                               plot_legend=legend is None,
                               extra_features=extra_features)
        figures.append(figure)
        if l:
            legend = l

    if options.slide_show:
        if options.slide_show == "first":
            pylab.show()
        elif options.slide_show == "all":
            pylab.show()
        elif options.slide_show == "sequence":
            pylab.show()
    else:

        extra_args = {}
        if options.dpi:
            extra_args['dpi'] = options.dpi

        for contig, figure in zip(contigs, figures):
            params = {'contig': contig}
            filename = options.output_pattern_image % params
            E.info("# creating image: %s" % filename)
            figure.savefig(os.path.expanduser(filename), **extra_args)
        if legend:
            params = {'contig': "legend"}
            filename = options.output_pattern_image % params
            E.info("creating image: %s" % filename)
            legend.savefig(os.path.expanduser(filename), **extra_args)

    E.info("ninput=%i, ncontigs=%i, nplots=%i" %
           (len(tracks), nplots, len(contigs)))

    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

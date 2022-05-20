import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from genome_browser._util import (
    ax_off,
    despine,
    disjoint_bins,
    ticklabels_to_percent,
    ticklabels_to_thousands_sep)

__all__ = [
    'GraphTrack',
    'IntervalTrack',
    'Track']


class Track(object):
    def __init__(self, name="", height_ratio=1, step=500):
        self.annotate = True
        self.height_ratio = height_ratio
        self.name = name
        self.step = step

        self.ALPHA = 0.9
        self.PADDING = 0.3

    def highlight_interval(self, interval):
        ...


class FeatureTrack(Track):
    def __init__(self, name="", height_ratio=1, step=500):
        Track.__init__(self, name, height_ratio=height_ratio, step=step)
        self.features = []

    def add_feature(self, feature):
        self.features.append(feature)


class GraphTrack(Track):
    def __init__(self, name="", height_ratio=1, is_proportional=False, step=500):
        Track.__init__(self, name, height_ratio=height_ratio, step=step)
        self.graphs = []
        self.is_proportional = is_proportional

        self.RESOLUTION = np.int(1e5)

    @property
    def xlimits(self):
        left, right = zip(*[(min(graph.x), max(graph.x)) for graph in self.graphs])
        return min(left), max(right)

    @property
    def is_empty(self):
        return len(self.graphs) == 0

    def add_graph(self, graph):
        self.graphs.append(graph)

    def plot(self, ax=None):
        ax = ax or plt.gca()

        for spine in ('top', 'right'):
            ax.spines[spine].set_visible(False)

        # Consider order of plotted graphs to increase z-order appropriately
        # This ensure that each graph neatly overlays over the last and is not
        # stacked in the same z-layer.
        for i, graph in enumerate(self.graphs):
            x, y = graph.x, graph.y

            if graph.fmt == 'interpolate':
                from scipy.interpolate import Akima1DInterpolator

                x_new = np.linspace(min(x), max(x), self.RESOLUTION)
                x, y = x_new, Akima1DInterpolator(x, y)(x_new)

            # Fill under x, y values with the same color as the plot. Otherwise
            # plot just the line.
            if graph.fill:
                ax.fill_between(
                    x, [0] * len(x), y,
                    interpolate=True,
                    lw=None,
                    edgecolor=None,
                    facecolor=graph.color,
                    antialiased=True,
                    alpha=graph.alpha,
                    zorder=2.7 + i)
            else:
                ax.plot(x, y, color=graph.color, alpha=graph.alpha, zorder=2.7 + i)

        if max(ax.get_ylim()) > 999:
            ax = ticklabels_to_thousands_sep(ax, 'y')

        # If this is proportional data then convert y-axis to percents.
        if self.is_proportional is True:
            ax = ticklabels_to_percent(ax, 'y')

        for k, spine in ax.spines.items():
            spine.set_zorder(2.7 + i + 1)

        ax.yaxis.grid(True, color='0.8', ls='-', zorder=0)
        ax.set_ylim(0, ax.get_ylim()[1])

        return ax


class IntervalTrack(Track):
    def __init__(self, name="", height_ratio=1, step=500):
        Track.__init__(self, name, height_ratio=height_ratio, step=step)
        self.intervals = []
        self.pullback = 1

    @property
    def interval_pretty_sort(self):
        """Function to sort the intervals by three attributes in priority order
        of start position, length, and strand orientation. First intervals are
        sorted by in ascending order by how soon they appear in the figure.
        Then, intervals are sorted in descending order with regard to their
        length. This ensures longer left-aligning intervals will appear earlier
        and lower on the stacking diagram. Finally strandness is sorted for to
        help keep order beyond placements.

        """
        return sorted(
            self.intervals,
            key=lambda _: (_.start, -len(_), _.strand))

    @property
    def is_empty(self):
        return not self.intervals

    @property
    def xlimits(self):
        """Return a tuple of the smallest and largest break in all features."""
        left, right = zip(*[(interval.start, interval.end) for interval in self.intervals])
        return min(left), max(right)

    def add_interval(self, interval):
        self.intervals.append(interval)

    def from_bam_file():
        ...

    def from_bed_file():
        ...

    def plot(self, ax=None):
        ax = ax or plt.gca()

        # Height of the features is some percent less than one.
        height = 1 / (self.PADDING + 1)

        # The non-overlapping disjoint intervals are computed using the
        # logic which priortizes first position of interval, then length.
        levels = list(disjoint_bins(self.interval_pretty_sort))

        # pull_back is the distance in units to pull back corners for
        # directional intervals
        pull_back = self.pullback * 0.005 * abs(np.subtract(*self.xlimits))

        for interval, level in zip(self.interval_pretty_sort, levels):

            width = interval.end - interval.start

            # pull_back is used on either the start or end of the interval
            # depending on the strand, if the pull_back is greater than the
            # width of the interval, then just pull back the entire width.
            start_taper = min(pull_back if interval.strand == '-' else 0, width)
            end_taper = min(pull_back if interval.strand == '+' else 0, width)

            # The polygon is simply a rectangle with two variable midpoints at
            # the middle of the left and right sides which act as anchors.
            # The four corners can be 'pulled back' (either left or right) to
            # simulate a directional rectangle.

            ax.add_patch(mpatches.Polygon(
                [[interval.start + start_taper, level],
                 [interval.start, level + height / 2],
                 [interval.start + start_taper, level + height],

                 [interval.start + width - end_taper, level + height],
                 [interval.start + width, level + height / 2],
                 [interval.start + width - end_taper, level]],
                lw=0,
                closed=True,
                color=interval.get('color', '0.2'),
                alpha=self.ALPHA))

        # For features, remove y-axis by default.
        ax = despine(ax_off(ax, axis='y'))

        # Adjust y-limits to include padding, scales with the number of levels.
        ax.set_ylim(
            (0 - self.PADDING) * (max(levels) + 1) / 2,
            max(levels) + 1)
        ax.set_xlim(*self.xlimits)

        return ax

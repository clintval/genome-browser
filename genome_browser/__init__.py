import collections as collections

from operator import itemgetter as itemgetter

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as patches
import numpy as np

from itertools import chain as chain
from more_itertools import flatten as flatten
from scipy.interpolate import Akima1DInterpolator as interpolate

ALPHA = 0.9
ASPECT = 2.8
HSPACE = 0.05
PADDING = 0.3
RESOLUTION = 1e5

PALETTE = {
    'gene': '#00A388',
    'misc_recomb': '#FF4A32',
    'repl_origin': '#97D54C',
    'nice_set': ['#E74C3C', '#3498DB', '#238b45']
}


def ax_off(ax, axis='x'):
    getattr(ax, 'get_{}axis'.format(axis))().set_visible(False)
    return ax


def cleanup_chart_junk(ax):
    ax.grid(False)
    ax.patch.set_facecolor('white')
    plt.gcf().set_facecolor('white')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for tic in ax.yaxis.get_major_ticks() + ax.xaxis.get_major_ticks():
        tic.tick2On = False
    return ax


def despine(ax):
    for spine in ['top', 'left', 'bottom', 'right']:
        ax.spines[spine].set_visible(False)
    return ax


def ticklabels_to_percent(ax, axis='y'):
    getattr(ax, '{}axis'.format(axis)).set_major_formatter(
        mticker.FuncFormatter(lambda s, position: '{:0.2%}'.format(s)))
    return ax


def ticklabels_to_thousands_sep(ax, axis='y'):
    getattr(ax, '{}axis'.format(axis)).set_major_formatter(
        mticker.FuncFormatter(lambda s, position: '{:,}'.format(int(s))))
    return ax


def disjoint_bins(intervals):
    bins = collections.defaultdict(list)

    for interval in intervals:
        start, end = interval
        interval = set(range(start, end + 1))

        level = None
        for bin, features in sorted(bins.items()):
            if not any([interval.intersection(f) for f in features]):
                level = bin
                break

        if not bins:
            level = 0
        elif level is None:
            level = max(bins) + 1
        bins[level].append(interval)

        yield level


def impute_zeros(x, y):
    imputed_y = np.zeros(np.ptp(x) + 1)
    np.put(imputed_y,
           ind=np.array(x) - min(x),
           v=np.array(y),
           mode='clip')
    return imputed_y.tolist()


class GenomeDiagram(object):
    def __init__(self, name=None):
        self.annotation = None
        self.name = name
        self.tracks = []

    @property
    def height_ratios(self):
        return [track.height_ratio for track in self.tracks]

    @property
    def xlim(self):
        flat = list(chain.from_iterable(
                    [t.xlim for t in self.tracks if t._empty is False]))
        return min(flat), max(flat)

    def add_track(self, track, ifempty=True):
        if ifempty is False and track._empty is True:
            return None
        else:
            self.tracks.append(track)

    def draw(self):
        # Construct a one-column figure with rows that follow the hight_ratios
        # defined in each of the tracks defined. At the moment only a generic
        # horizontal space parameter (HSPACE) is supplied, but spacing can be
        # tweaked post-draw.
        fig, axes = plt.subplots(
            nrows=len(self.tracks),
            ncols=1,
            figsize=(20, len(self.tracks) * 0.75 * ASPECT),
            gridspec_kw={
                'height_ratios': self.height_ratios,
                'hspace': HSPACE,
                'wspace': 0.0})

        # If this figure has a name, make it the title.
        if self.name is not None:
            fig.suptitle(self.name, x=0.5, y=0.94, fontsize=24)

        # If we are only plotting one track then put it in an iterable for zip.
        if len(self.tracks) == 1:
            axes = np.array([axes])

        for i, (ax, track) in enumerate(zip(axes.flatten(), self.tracks)):
            if track._empty is True:
                ax.axis('off')
                continue

            # Cleanup each track before plotting.
            ax = track._plot(cleanup_chart_junk(ax))

            # Set each ax to the limits of the greatest data range.
            ax.set_xlim(*self.xlim)

            # If this is not the last ax then turn off the x-axis. If not,
            # plot xticks at the step interval defined.
            if i != len(self.tracks) - 1:
                ax = ax_off(ax, axis='x')
            else:
                ax.spines['bottom'].set_visible(True)
                # Set the positions every step size after the 0-defined min.
                ax.set_xticks(range(*map(int, ax.get_xlim()), track.step))
                # Set the labels as numbers increasing by step, after 0.
                ax.set_xticklabels(
                    range(0, int(abs(np.subtract(*ax.get_xlim()))),
                          track.step))
                # Figure annotations will be applied to the last ax in offset
                # coordinates in a lightgray text. Clipping is ignored as the
                # text clearly cips with the axes outboard frame.
                if self.annotation is not None:
                    ax.annotate(xy=(1, 0),
                                xycoords='axes fraction',
                                s=self.annotation,
                                xytext=(0, -60),
                                textcoords='offset points',
                                va='bottom',
                                ha='right',
                                color='0.6',
                                clip_on=False)

            # Provide simple logic for plotting a track annotation in the
            # top left of each track. The position will remain constant as its
            # defined using proportional values of the ax's x and y limits.
            if track.annotate is True and track.name is not None:
                ax.annotate(track.name,
                            xy=(ax.get_xlim()[0] +
                                abs(np.subtract(*ax.get_xlim())) / 100,
                                ax.get_ylim()[1] / 1.01),
                            va='top', ha='left', annotation_clip=False)
        return fig, axes


class Track(object):
    def __init__(self, name=None, height_ratio=1):
        self._empty = True
        self.annotate = True
        self.height_ratio = height_ratio
        self.name = name
        self.step = 500

    def highlight_interval(self):
        ...


class Feature(Track):
    def __init__(self, name=None, height_ratio=1):
        Track.__init__(self, name, height_ratio=height_ratio)
        self.features = []
        self.pullback = 1

    @property
    def _sorted_features(self):
        """
        Priority sort, in ascending order, features based on position, length,
        and strand.
        """
        return sorted(self.features, key=itemgetter(0, 1, 2), reverse=False)

    @property
    def _intervals(self):
        """Return a list of all open-ended intervals of all features."""
        positions, widths, *_ = zip(*self._sorted_features)
        return [(pos, pos + width) for pos, width in zip(positions, widths)]

    @property
    def xlim(self):
        """Return a tuple of the smallest and largest break in all features."""
        edges = tuple(chain.from_iterable(self._intervals))
        return min(edges), max(edges)

    def add_feature(self, feature):
        self._empty = False
        self.features.append(feature)

    def _plot(self, ax=None):
        if ax is None:
            ax = plt.gca()

        # Height of the features is some percent less than one.
        height = 1 / (PADDING + 1)

        # The non-overlapping disjoint intervals are computed using the
        # logic which priortizes first position of interval, then length.
        levels = list(disjoint_bins(self._intervals))

        # pull_back is the distance in units to pull back corners for
        # directional intervals
        pull_back = self.pullback * 0.005 * abs(np.subtract(*self.xlim))

        for (position, width, strand, color), level in zip(
                self._sorted_features, levels):

            # pull_back is used on either the start or end of the interval
            # depending on the strand, if the pull_back is greater than the
            # width of the interval, then just pull back the entire width.
            start_taper = min(pull_back if strand == '-' else 0, width)
            end_taper = min(pull_back if strand == '+' else 0, width)

            # The polygon is simply a rectangle with two variable midpoints at
            # the middle of the left and right sides which act as anchors.
            # The four corners can be 'pulled back' (either left or right) to
            # simulate a directional rectangle.
            ax.add_patch(patches.Polygon(
                [[position + start_taper, level],
                 [position, level + height / 2],
                 [position + start_taper, level + height],

                 [position + width - end_taper, level + height],
                 [position + width, level + height / 2],
                 [position + width - end_taper, level]],
                lw=0, closed=True, color=color, alpha=ALPHA))

        # For features, remove y-axis by default.
        ax = despine(ax_off(ax, axis='y'))

        # Adjust y-limits to include padding, scales with the number of levels.
        ax.set_ylim((0 - PADDING) * (max(levels) + 1) / 2, max(levels) + 1)
        ax.set_xlim(*self.xlim)
        return ax


class Graph(Track):
    def __init__(self, name=None, height_ratio=1, is_proportional=False):
        Track.__init__(self, name, height_ratio=height_ratio)
        self._graphs = []
        self.is_proportional = is_proportional

    @property
    def xlim(self):
        independents, *_ = zip(*self._graphs)
        flat = list(chain.from_iterable(
                    [x for x in independents]))
        return min(flat), max(flat)

    def new_graph(self, x, y, color='0.1', alpha=1, fmt='interpolate',
                  fill=False):
        self._empty = False
        self._graphs.append((x, y, color, fmt, fill, alpha))

    def _plot(self, ax=None):
        if ax is None:
            ax = plt.gca()

        # Consider order of plotted graphs to increase z-order appropriately
        # This ensure that each graph neatly overlays over the last and is not
        # stacked in the same z-layer.
        for i, (x, y, color, fmt, fill, alpha) in enumerate(self._graphs):
            if fmt == 'interpolate':
                x_new = np.linspace(min(x), max(x), RESOLUTION)
                x, y = x_new, interpolate(x, y)(x_new)

            # Fill under x, y values with the same color as the plot. Otherwise
            # plot just the line.
            if fill is True:
                ax.fill_between(
                    x, [0] * len(x), y,
                    interpolate=True,
                    lw=None,
                    edgecolor=None,
                    facecolor=color,
                    antialiased=True,
                    alpha=alpha,
                    zorder=2.7 + i)
            else:
                ax.plot(x, y, color=color, alpha=alpha, zorder=2.7 + i)

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

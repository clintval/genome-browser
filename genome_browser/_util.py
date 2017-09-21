import collections

import matplotlib.ticker as mticker
import numpy as np


def ax_off(ax, axis='x'):
    getattr(ax, 'get_{}axis'.format(axis))().set_visible(False)
    return ax


def impute_zeros(x, y):
    imputed_y = np.zeros(np.ptp(x) + 1)
    np.put(
        imputed_y,
        ind=np.array(x) - min(x),
        v=np.array(y),
        mode='clip')
    return imputed_y.tolist()


def despine(ax):
    for spine in ['top', 'left', 'bottom', 'right']:
        ax.spines[spine].set_visible(False)
    return ax


def disjoint_bins(intervals):
    bins = collections.defaultdict(list)

    for interval in intervals:

        if isinstance(interval, type(interval)):  ###TODO
            start, end = interval.start, interval.end
        elif isinstance(interval, tuple):
            start, end = interval
        else:
            raise ValueError(f'Interval {repr(interval)} not of class Interval or tuple')

        level = None
        sites = set(range(start, end + 1))

        for bin, features in sorted(bins.items()):
            if not any([sites.intersection(f) for f in features]):
                level = bin
                break

        if not bins:
            level = 0
        elif level is None:
            level = max(bins) + 1

        bins[level].append(sites)

        yield level


def ticklabels_to_percent(ax, axis='y'):
    getattr(ax, '{}axis'.format(axis)).set_major_formatter(
        mticker.FuncFormatter(lambda s, position: '{:0.2%}'.format(s)))
    return ax


def ticklabels_to_thousands_sep(ax, axis='y'):
    getattr(ax, '{}axis'.format(axis)).set_major_formatter(
        mticker.FuncFormatter(lambda s, position: '{:,}'.format(int(s))))
    return ax

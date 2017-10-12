from collections import defaultdict
from tempfile import NamedTemporaryFile

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


def despine(ax, top=True, left=True, bottom=True, right=True):
    for spine, on in zip(
        ('top', 'left', 'bottom', 'right'), (top, left, bottom, right)
    ):
        ax.spines[spine].set_visible(not on)
    return ax


def disjoint_bins(intervals):
    """Given an interable of interval features provide the most efficient
    vertical packing of them such that none overlap. Horizontal packing is
    determined solely by the order in which they are presented. For instance,
    an iterable sorted by length will prioritize placing long intervals at
    lower levels. Returned are the levels each interval must exist on to ensure
    there are no overlaps.

    Parameters
    ----------
    intervals : iterable
        Can be a list of (float, float) tuples or a list of interval classes
        with 'start' and 'stop' attributes such as a BedTool interval or a VCF
        Record class. Tuple-style intervals and class intervals can be mixed.

    Returns
    -------
    levels : generator
        The levels corresponding to each interval in `iterable` such that none
        overlap.

    """
    bins = defaultdict(list)

    for interval in intervals:
        if hasattr(interval, 'start') and hasattr(interval, 'end'):
            start, end = interval.start, interval.end
        else:
            start, end = interval

        loci = set(range(start, end + 1))

        level = None
        for bin, features in sorted(bins.items()):
            if not any((loci.intersection(feature) for feature in features)):
                level = bin
                break

        if not bins:
            level = 0
        elif level is None:
            level = max(bins) + 1
        bins[level].append(loci)

        yield level


def string_as_temporary_file(string):
    handle = NamedTemporaryFile(mode='w+', delete=False)
    handle.write(string)
    handle.close()
    return handle


def ticklabels_to_percent(ax, axis='y'):
    getattr(ax, '{}axis'.format(axis)).set_major_formatter(
        mticker.FuncFormatter(lambda s, position: '{:0.2%}'.format(s)))
    return ax


def ticklabels_to_thousands_sep(ax, axis='y'):
    getattr(ax, '{}axis'.format(axis)).set_major_formatter(
        mticker.FuncFormatter(lambda s, position: '{:,}'.format(int(s))))
    return ax

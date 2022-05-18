from typing import List
import itertools
import matplotlib.pyplot as plt
import numpy as np

from genome_browser._util import ax_off

__all__ = [
    'GenomeDiagram']


class GenomeDiagram(object):
    def __init__(self, targets: List[tuple], name: str = None):
        self.targets = [targets] if not isinstance(targets, (list, tuple)) else targets
        self.name = name

        self.annotation = None
        self.tracks = []

        self.ASPECT = 2.8
        self.HSPACE = 0.05

    @property
    def bottom_track(self):
        return self.tracks[-1] if len(self.tracks) > 0 else None

    @property
    def height_ratios(self):
        return [track.height_ratio for track in self.tracks]

    @property
    def xlimits(self):
        all_xlimits = [track.xlimits for track in self.tracks if not track.is_empty]
        flat = list(itertools.chain.from_iterable(all_xlimits))

        return min(flat), max(flat)

    def add_track(self, track, add_if_empty=True):
        if track.is_empty and not add_if_empty:
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
            ncols=len(self.targets),
            figsize=(20, len(self.tracks) * 0.75 * self.ASPECT),
            gridspec_kw={
                'height_ratios': self.height_ratios,
                'hspace': self.HSPACE,
                'wspace': 0.0})

        if self.name is not None:
            fig.suptitle(self.name, x=0.5, y=0.94, fontsize=24)

        axes = np.array([axes]) if len(self.tracks) == 1 else axes

        for i, (ax, track) in enumerate(zip(axes.flatten(), self.tracks)):
            if track.is_empty:
                ax.axis('off')
                continue

            ax = track.plot(ax)

            # Set each ax to the limits of the greatest data range.
            ax.set_xlim(*self.xlimits)

            # If this is not the last ax then turn off the x-axis. If not,
            # plot xticks at the step interval defined.
            if i != len(self.tracks) - 1:
                ax = ax_off(ax, axis='x')
            else:
                ax.spines['bottom'].set_visible(True)

                # Set the positions every step size after the 0-defined min.
                ax.set_xticks(range(*map(int, ax.get_xlim()), track.step))

                # Set the labels as numbers increasing by step, after 0.
                ax.set_xticklabels(range(
                        0,
                        int(abs(np.subtract(*ax.get_xlim()))),
                        track.step))

                # Figure annotations will be applied to the last ax in offset
                # coordinates in a lightgray text. Clipping is ignored as the
                # text clearly cips with the axes outboard frame.
                if self.annotation is not None:
                    ax.annotate(
                        xy=(1, 0),
                        xycoords='axes fraction',
                        text=self.annotation,
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
                ax.annotate(
                    text=track.name,
                    xy=(ax.get_xlim()[0] + abs(np.subtract(*ax.get_xlim())) / 100,
                        ax.get_ylim()[1] / 1.01),
                    va='top',
                    ha='left',
                    annotation_clip=False)

        return fig, axes

    def __repr__(self):
        return f'GenomeDiagram(targets={repr(self.targets[0])}...)'

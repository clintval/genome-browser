from matplotlib import pyplot as plt
import numpy as np

import genome_browser as gb

n = 105  # Length of random genomic interval.

# Plot a density of random data, interpolated and filled.
track1 = gb.GraphTrack("Random Density")
track1.add_graph(
    gb.Graph(x=np.arange(n), y=np.abs(np.random.randn(n)), fmt="interpolate", fill=True)
)

# Plot 9 random interval features (random start, length, orientation, and color).
track2 = gb.GraphTrack("Random Intervals", height_ratio=0.4)
genomic_intervals = []
for _ in range(9):
    # Feature must follow iterable as (position, width, strand, color)
    genomic_intervals.append(
        gb.Interval(
            chrom="chr3",
            start=np.random.randint(0, 20),
            end=np.random.randint(20, n - 15),
            strand=np.random.choice(["+", "-", "."]),
            metadata={"color": np.random.choice(["#E74C3C", "#3498DB", "0.2"])},
        )
    )

g = gb.GenomeDiagram(genomic_intervals)

g.add_track(track1)
g.add_track(track2)

# Annotate the figure with interval specific metadata. Will always appear in lower-right
g.annotation = "{}:{:,}-{:,}".format("chr3", 20000, 812383)

fig, axes = g.draw()
plt.show()

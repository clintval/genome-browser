## Genome Diagram

A lightweight Python 3 graphing API for genomic graphs over genomic intervals.

```python
from matplotlib import pyplot as plt
import numpy as np

import genome_browser as gb

n = 105  # Length of random genomic interval.

g = gb.GenomeDiagram()

# Plot a density of random data, interpolated and filled.
graph_track = gb.GraphTrack("Random Density")
graph_track.add_graph(
    gb.Graph(x=np.arange(n), y=np.abs(np.random.randn(n)), fmt="interpolate", fill=True)
)
g.add_track(graph_track)

# Plot 9 random interval features (random start, length, orientation, and color).
interval_track = gb.IntervalTrack("Random Intervals", height_ratio=0.4, step=50)
for _ in range(9):
    # Feature must follow iterable as (position, width, strand, color)
    interval = gb.Interval(
        chrom="chr3",
        start=np.random.randint(0, n),
        end=np.random.randint(0, n),
        strand=np.random.choice(["+", "-", "."]),
        color=np.random.choice(["#E74C3C", "#3498DB", "0.2"]),
    )
    interval_track.add_interval(interval)


g.add_track(interval_track)

# Annotate the figure with interval specific metadata. Will always appear in lower-right
g.annotation = "{}:{:,}-{:,}".format("chr3", 20000, 812383)

fig, axes = g.draw()
plt.show()
```

![test_interval](https://raw.githubusercontent.com/clintval/genome-browser/master/img/gb_test.png "Test Interval")

"""
Plot GPU vs CPU t-SNE benchmark results.
"""

import pandas as pd
from plotnine import (
    ggplot, aes, geom_line, geom_point,
    scale_x_log10, scale_y_log10, scale_color_manual,
    labs, theme
)
from theme_nxn import theme_nxn, get_nxn_palette

# Benchmark data
data = pd.DataFrame({
    "cells": [50_000, 250_000, 500_000, 1_058_909] * 2,
    "time_seconds": [4.9, 16.1, 40.2, 129.0, 25.5, 79.4, 159.9, 332.0],
    "method": ["GPU (serverless, Modal T4)"] * 4 + ["openTSNE (local, 12 cores)"] * 4,
})

# Format y-axis as minutes:seconds
def format_time(seconds):
    mins = int(seconds // 60)
    secs = int(seconds % 60)
    if mins > 0:
        return f"{mins}m {secs:02d}s"
    return f"{secs}s"

# Create breaks for y-axis (log scale)
y_breaks = [5, 10, 30, 60, 120, 300]
y_labels = [format_time(s) for s in y_breaks]

# Create breaks for x-axis (log scale)
x_breaks = [50_000, 100_000, 250_000, 500_000, 1_000_000]
x_labels = ["50k", "100k", "250k", "500k", "1M"]

palette = get_nxn_palette()

p = (
    ggplot(data, aes(x="cells", y="time_seconds", color="method"))
    + geom_line(size=1)
    + geom_point(size=3)
    + scale_x_log10(breaks=x_breaks, labels=x_labels)
    + scale_y_log10(breaks=y_breaks, labels=y_labels)
    + scale_color_manual(values=[palette[0], palette[1]])
    + labs(
        x="Dataset size",
        y="Runtime",
        color="Method",
        title="TSNE performance",
    )
    + theme_nxn(base_size=16)
    + theme(figure_size=(8, 5))
)

p.save("benchmark_plot.png", dpi=150)
print("Saved benchmark_plot.png")

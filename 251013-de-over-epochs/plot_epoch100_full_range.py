#!/usr/bin/env python3
"""
Create volcano plot for epoch 100 with full auto-scaled axes.
"""

import pandas as pd
import numpy as np
from plotnine import *
from theme_nxn import theme_nxn, get_nxn_palette

# Load DE results
print("Loading DE results...")
de_long = pd.read_parquet('results/de_results.parquet')

# Filter to epoch 100
de_100 = de_long[de_long['epoch'] == 100].copy()

print(f"Epoch 100 data shape: {de_100.shape}")
print(f"Cell types: {de_100['celltype'].unique()}")

# Get all cell types
celltypes = sorted(de_100['celltype'].unique())

# Create volcano plot for each cell type
for celltype in celltypes:
    print(f"Plotting {celltype}...")

    df = de_100[de_100['celltype'] == celltype].copy()

    # Create the plot
    p = (
        ggplot(df, aes(x='lfc_mean', y='-log10_pval'))
        + geom_point(alpha=0.3, size=0.5, color=get_nxn_palette()[0])
        + labs(
            title=f'Epoch 100 - {celltype}',
            x='Log fold change (Inflamed vs NonInvolved)',
            y='-log10(p-value)'
        )
        + theme_nxn()
        + theme(
            figure_size=(10, 8)
        )
    )

    # Save
    filename = f'epoch_100_fullrange_{celltype.replace("~", "_").replace("[", "_").replace("]", "_").replace("*", "")}.png'
    p.save(filename, dpi=150)
    print(f"  Saved: {filename}")

# Also create a faceted version with all cell types
print("\nCreating faceted plot with all cell types...")

p_all = (
    ggplot(de_100, aes(x='lfc_mean', y='-log10_pval'))
    + geom_point(alpha=0.3, size=0.3, color=get_nxn_palette()[0])
    + facet_wrap('~celltype', ncol=3, scales='free')
    + labs(
        title='Epoch 100 - All Cell Types (Auto-scaled)',
        x='Log fold change (Inflamed vs NonInvolved)',
        y='-log10(p-value)'
    )
    + theme_nxn()
    + theme(
        figure_size=(14, 10),
        strip_text=element_text(size=9)
    )
)

p_all.save('epoch_100_fullrange_all_celltypes.png', dpi=150)
print("  Saved: epoch_100_fullrange_all_celltypes.png")

print("\nDone!")

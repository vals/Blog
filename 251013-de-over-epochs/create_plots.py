#!/usr/bin/env python3
"""
Create volcano plots from saved DE results with fixed coordinate limits.

This script loads pre-computed DE results and generates plots with
consistent axes for smooth animations.
"""

import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from plotnine import *
from theme_nxn import theme_nxn, get_nxn_palette

# Paths
RESULTS_DIR = Path("./results")
PLOT_DIR = Path("./plots")
PLOT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("Creating Volcano Plots from Saved Results")
print("="*80)

# Load results
print("\nLoading saved results...")
de_long = pd.read_parquet(RESULTS_DIR / 'de_results.parquet')

with open(RESULTS_DIR / 'metadata.pkl', 'rb') as f:
    metadata = pickle.load(f)

valid_celltypes = metadata['valid_celltypes']
epochs = sorted(de_long['epoch'].unique())

print(f"Loaded DE results for {len(epochs)} epochs")
print(f"Cell types: {valid_celltypes}")
print(f"Total rows: {len(de_long):,}")


# Calculate global coordinate limits
print("\nCalculating global coordinate limits...")

# Calculate -log10(p) if not already in dataframe
if '-log10_pval' not in de_long.columns:
    de_long['-log10_pval'] = -np.log10(de_long['proba_not_de'] + 1e-10)

# Get all values for limit calculation
all_log2fc = de_long['lfc_mean'].values
all_neglog10p = de_long['-log10_pval'].values

# Set limits with some padding
log2fc_lim = (-2.5, 2.5)

neglog10p_max = np.percentile(all_neglog10p, 99.5)
neglog10p_lim = (0, neglog10p_max * 1.1)

print(f"Log2FC limits: {log2fc_lim}")
print(f"-log10(p) limits: {neglog10p_lim}")


# Create volcano plots
def create_volcano_plot(de_data_epoch, epoch, save_path):
    """Create multi-facet volcano plot with fixed axes."""

    # Filter data for this epoch
    plot_df = de_data_epoch[de_data_epoch['epoch'] == epoch].copy()

    if len(plot_df) == 0:
        print(f"No data for epoch {epoch}")
        return

    # Prepare plotting dataframe
    plot_df = plot_df.rename(columns={'lfc_mean': 'log2FC', '-log10_pval': 'neg_log10_pval'})

    # Create plot with fixed limits
    p = (
        ggplot(plot_df, aes(x='log2FC', y='neg_log10_pval')) +
        geom_point(alpha=0.3, size=0.5, color=get_nxn_palette()[0]) +
        geom_hline(yintercept=-np.log10(0.05), linetype='dashed', color='red', alpha=0.5) +
        geom_vline(xintercept=0, linetype='dashed', color='gray', alpha=0.5) +
        facet_wrap('~celltype', ncol=3) +
        xlim(log2fc_lim[0], log2fc_lim[1]) +
        ylim(neglog10p_lim[0], neglog10p_lim[1]) +
        labs(
            x='Log2 Fold Change (Inflamed vs NonInvolved)',
            y='-log10(p-value)',
            title=f'Differential Expression - Epoch {epoch}'
        ) +
        theme_nxn() +
        theme(
            figure_size=(14, 12),
            strip_text=element_text(size=10, face='bold')
        )
    )

    p.save(save_path, dpi=150, verbose=False)
    return save_path


# Generate all plots
print("\nGenerating volcano plots...")
print(f"Creating plots for {len(epochs)} epochs...")

for i, epoch in enumerate(epochs, 1):
    save_path = PLOT_DIR / f"volcano_epoch_{epoch:03d}.png"

    create_volcano_plot(de_long, epoch, save_path)

    if i % 5 == 0:
        print(f"  Created {i}/{len(epochs)} plots")

print(f"  Created {len(epochs)}/{len(epochs)} plots ✓")


# Create training curves
print("\nCreating training curves...")

history_df = pd.read_csv(RESULTS_DIR / 'training_history.csv')

history_long = pd.melt(
    history_df,
    id_vars=['epoch'],
    value_vars=['train_loss', 'val_loss'],
    var_name='metric',
    value_name='loss'
)

p_train = (
    ggplot(history_long, aes(x='epoch', y='loss', color='metric')) +
    geom_line(size=1) +
    scale_color_manual(
        values=[get_nxn_palette()[0], get_nxn_palette()[1]],
        labels=['Training Loss', 'Validation Loss']
    ) +
    labs(
        x='Epoch',
        y='ELBO Loss',
        title='SCVI Training Curves',
        color='Metric'
    ) +
    theme_nxn() +
    theme(
        figure_size=(10, 6),
        legend_position='right'
    )
)

train_curve_path = PLOT_DIR / "training_curves.png"
p_train.save(train_curve_path, dpi=150, verbose=False)
print(f"Saved training curves: {train_curve_path}")


print("\n" + "="*80)
print("Plotting complete!")
print("="*80)
print(f"\nPlots directory: {PLOT_DIR}")
print(f"- {len(epochs)} volcano plots")
print(f"- training_curves.png")
print(f"\nCoordinate limits used:")
print(f"- Log2FC: [{log2fc_lim[0]:.2f}, {log2fc_lim[1]:.2f}]")
print(f"- -log10(p): [0, {neglog10p_lim[1]:.2f}]")
print(f"\nTo create animation, run:")
print(f"  ./create_animation.sh")

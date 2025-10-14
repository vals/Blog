#!/usr/bin/env python3
"""
Check the actual ranges of log fold change values.
"""

import pandas as pd
import numpy as np

# Load DE results
print("Loading DE results...")
de_long = pd.read_parquet('results/de_results.parquet')

print(f"Data shape: {de_long.shape}")
print(f"Epochs: {de_long['epoch'].min()} to {de_long['epoch'].max()}")

# Look at all epochs
print("\n" + "="*80)
print("LOG FOLD CHANGE RANGES (ALL EPOCHS)")
print("="*80)

print(f"\nOverall statistics:")
print(f"  Min: {de_long['lfc_mean'].min():.2f}")
print(f"  Max: {de_long['lfc_mean'].max():.2f}")
print(f"  Mean: {de_long['lfc_mean'].mean():.2f}")
print(f"  Median: {de_long['lfc_mean'].median():.2f}")
print(f"  Std: {de_long['lfc_mean'].std():.2f}")

print(f"\nPercentiles:")
for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
    val = de_long['lfc_mean'].quantile(p/100)
    print(f"  {p:2d}%: {val:6.2f}")

# By epoch
print("\n" + "="*80)
print("LFC RANGES BY EPOCH (selected epochs)")
print("="*80)

epochs_to_show = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
epoch_summary = []

for epoch in epochs_to_show:
    de_epoch = de_long[de_long['epoch'] == epoch]
    epoch_summary.append({
        'epoch': epoch,
        'min': de_epoch['lfc_mean'].min(),
        'max': de_epoch['lfc_mean'].max(),
        'mean': de_epoch['lfc_mean'].mean(),
        'std': de_epoch['lfc_mean'].std(),
        'p01': de_epoch['lfc_mean'].quantile(0.01),
        'p99': de_epoch['lfc_mean'].quantile(0.99)
    })

epoch_df = pd.DataFrame(epoch_summary)
print(epoch_df.to_string(index=False))

# Late epochs only
print("\n" + "="*80)
print("LFC RANGES FOR LATE EPOCHS (50-100)")
print("="*80)

de_late = de_long[de_long['epoch'] >= 50]

print(f"\nOverall statistics (epochs 50-100):")
print(f"  Min: {de_late['lfc_mean'].min():.2f}")
print(f"  Max: {de_late['lfc_mean'].max():.2f}")
print(f"  Mean: {de_late['lfc_mean'].mean():.2f}")
print(f"  Median: {de_late['lfc_mean'].median():.2f}")
print(f"  Std: {de_late['lfc_mean'].std():.2f}")

print(f"\nPercentiles (epochs 50-100):")
for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
    val = de_late['lfc_mean'].quantile(p/100)
    print(f"  {p:2d}%: {val:6.2f}")

# By cell type
print("\n" + "="*80)
print("LFC RANGES BY CELL TYPE (epochs 50-100)")
print("="*80)

ct_summary = de_late.groupby('celltype')['lfc_mean'].agg([
    ('min', 'min'),
    ('p01', lambda x: x.quantile(0.01)),
    ('p25', lambda x: x.quantile(0.25)),
    ('median', 'median'),
    ('p75', lambda x: x.quantile(0.75)),
    ('p99', lambda x: x.quantile(0.99)),
    ('max', 'max'),
    ('std', 'std')
]).round(2)

print(ct_summary)

# Distribution of absolute LFC
print("\n" + "="*80)
print("DISTRIBUTION OF ABSOLUTE LFC (epochs 50-100)")
print("="*80)

de_late['abs_lfc'] = abs(de_late['lfc_mean'])

print("\nHow many genes fall in different LFC bins?")
bins = [0, 0.5, 1, 1.5, 2, 3, 5, 10, 100]
for i in range(len(bins)-1):
    count = ((de_late['abs_lfc'] >= bins[i]) & (de_late['abs_lfc'] < bins[i+1])).sum()
    pct = count / len(de_late) * 100
    print(f"  |LFC| in [{bins[i]:4.1f}, {bins[i+1]:4.1f}): {count:8d} ({pct:5.2f}%)")

# Extreme LFC genes
print("\n" + "="*80)
print("GENES WITH LARGEST LOG FOLD CHANGES (epochs 50-100)")
print("="*80)

# Average across late epochs
de_late_avg = de_late.groupby(['celltype', 'gene'])['lfc_mean'].mean().reset_index()
de_late_avg.columns = ['celltype', 'gene', 'avg_lfc']

print("\nTop 20 upregulated genes (largest positive LFC):")
top_up = de_late_avg.nlargest(20, 'avg_lfc')[['celltype', 'gene', 'avg_lfc']]
print(top_up.to_string(index=False))

print("\nTop 20 downregulated genes (largest negative LFC):")
top_down = de_late_avg.nsmallest(20, 'avg_lfc')[['celltype', 'gene', 'avg_lfc']]
print(top_down.to_string(index=False))

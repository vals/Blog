#!/usr/bin/env python3
"""
Analyze stability of log fold changes across epochs.
"""

import pandas as pd
import numpy as np

# Load DE results
print("Loading DE results...")
de_long = pd.read_parquet('results/de_results.parquet')

print(f"Data shape: {de_long.shape}")
print(f"Epochs: {de_long['epoch'].min()} to {de_long['epoch'].max()}")
print(f"Cell types: {de_long['celltype'].nunique()}")

# Focus on epochs >= 50
de_late = de_long[de_long['epoch'] >= 50].copy()

print(f"\nLate epochs (>=50) data shape: {de_late.shape}")

# Calculate variability metrics for each gene in each cell type
print("\nCalculating log fold change variability after epoch 50...")

variability = de_late.groupby(['celltype', 'gene'])['lfc_mean'].agg([
    ('mean_lfc', 'mean'),
    ('std_lfc', 'std'),
    ('min_lfc', 'min'),
    ('max_lfc', 'max'),
    ('range_lfc', lambda x: x.max() - x.min()),
    ('cv', lambda x: x.std() / abs(x.mean()) if x.mean() != 0 else np.nan)
]).reset_index()

# Summary statistics
print("\n" + "="*80)
print("SUMMARY: Log Fold Change Variability (Epochs 50-100)")
print("="*80)

print("\nStandard deviation of LFC across epochs 50-100:")
print(f"  Mean: {variability['std_lfc'].mean():.4f}")
print(f"  Median: {variability['std_lfc'].median():.4f}")
print(f"  95th percentile: {variability['std_lfc'].quantile(0.95):.4f}")
print(f"  Max: {variability['std_lfc'].max():.4f}")

print("\nRange of LFC (max - min) across epochs 50-100:")
print(f"  Mean: {variability['range_lfc'].mean():.4f}")
print(f"  Median: {variability['range_lfc'].median():.4f}")
print(f"  95th percentile: {variability['range_lfc'].quantile(0.95):.4f}")
print(f"  Max: {variability['range_lfc'].max():.4f}")

# By cell type
print("\n" + "="*80)
print("VARIABILITY BY CELL TYPE")
print("="*80)

ct_summary = variability.groupby('celltype').agg({
    'std_lfc': ['mean', 'median', 'max'],
    'range_lfc': ['mean', 'median', 'max']
}).round(4)

print(ct_summary)

# Compare to early epochs (0-20)
print("\n" + "="*80)
print("COMPARISON: Early (0-20) vs Late (50-100) Epochs")
print("="*80)

de_early = de_long[de_long['epoch'] <= 20].copy()

variability_early = de_early.groupby(['celltype', 'gene'])['lfc_mean'].agg([
    ('std_lfc', 'std'),
    ('range_lfc', lambda x: x.max() - x.min())
]).reset_index()

print("\nEarly epochs (0-20):")
print(f"  Mean std of LFC: {variability_early['std_lfc'].mean():.4f}")
print(f"  Mean range of LFC: {variability_early['range_lfc'].mean():.4f}")

print("\nLate epochs (50-100):")
print(f"  Mean std of LFC: {variability['std_lfc'].mean():.4f}")
print(f"  Mean range of LFC: {variability['range_lfc'].mean():.4f}")

reduction = (1 - variability['std_lfc'].mean() / variability_early['std_lfc'].mean()) * 100
print(f"\nReduction in variability: {reduction:.1f}%")

# Find most and least stable genes
print("\n" + "="*80)
print("MOST STABLE GENES (lowest variability after epoch 50)")
print("="*80)
most_stable = variability.nsmallest(10, 'std_lfc')[['celltype', 'gene', 'mean_lfc', 'std_lfc', 'range_lfc']]
print(most_stable.to_string(index=False))

print("\n" + "="*80)
print("LEAST STABLE GENES (highest variability after epoch 50)")
print("="*80)
least_stable = variability.nlargest(10, 'std_lfc')[['celltype', 'gene', 'mean_lfc', 'std_lfc', 'range_lfc']]
print(least_stable.to_string(index=False))

# Genes with large LFC changes
print("\n" + "="*80)
print("HIGHLY VARIABLE DE GENES (|mean_lfc| > 1 and std > 0.1)")
print("="*80)
variable_de = variability[
    (abs(variability['mean_lfc']) > 1) &
    (variability['std_lfc'] > 0.1)
].sort_values('std_lfc', ascending=False)

print(f"\nFound {len(variable_de)} genes with |LFC| > 1 and std > 0.1")
if len(variable_de) > 0:
    print("\nTop 10:")
    print(variable_de.head(10)[['celltype', 'gene', 'mean_lfc', 'std_lfc', 'range_lfc']].to_string(index=False))
else:
    print("All strongly DE genes are stable after epoch 50!")

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)
if variability['std_lfc'].mean() < 0.01:
    print("✓ Log fold changes are very stable after epoch 50")
    print(f"  Average std: {variability['std_lfc'].mean():.4f}")
elif variability['std_lfc'].mean() < 0.05:
    print("✓ Log fold changes are reasonably stable after epoch 50")
    print(f"  Average std: {variability['std_lfc'].mean():.4f}")
else:
    print("⚠ Log fold changes show notable variability after epoch 50")
    print(f"  Average std: {variability['std_lfc'].mean():.4f}")

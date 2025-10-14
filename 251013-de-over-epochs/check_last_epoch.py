#!/usr/bin/env python3
"""
Check log fold changes in the last epoch (100).
"""

import pandas as pd
import numpy as np

# Load DE results
print("Loading DE results...")
de_long = pd.read_parquet('results/de_results.parquet')

# Filter to epoch 100
de_last = de_long[de_long['epoch'] == 100].copy()

print(f"Epoch 100 data shape: {de_last.shape}")

print("\n" + "="*80)
print("LOG FOLD CHANGE RANGES - EPOCH 100")
print("="*80)

print(f"\nOverall statistics:")
print(f"  Min: {de_last['lfc_mean'].min():.2f}")
print(f"  Max: {de_last['lfc_mean'].max():.2f}")
print(f"  Mean: {de_last['lfc_mean'].mean():.2f}")
print(f"  Median: {de_last['lfc_mean'].median():.2f}")
print(f"  Std: {de_last['lfc_mean'].std():.2f}")

print(f"\nPercentiles:")
for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
    val = de_last['lfc_mean'].quantile(p/100)
    print(f"  {p:2d}%: {val:6.2f}")

# By cell type
print("\n" + "="*80)
print("LFC RANGES BY CELL TYPE - EPOCH 100")
print("="*80)

ct_summary = de_last.groupby('celltype')['lfc_mean'].agg([
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

# Distribution
print("\n" + "="*80)
print("DISTRIBUTION OF ABSOLUTE LFC - EPOCH 100")
print("="*80)

de_last['abs_lfc'] = abs(de_last['lfc_mean'])

print("\nHow many genes fall in different LFC bins?")
bins = [0, 0.5, 1, 1.5, 2, 3, 5, 10, 100]
for i in range(len(bins)-1):
    count = ((de_last['abs_lfc'] >= bins[i]) & (de_last['abs_lfc'] < bins[i+1])).sum()
    pct = count / len(de_last) * 100
    print(f"  |LFC| in [{bins[i]:4.1f}, {bins[i+1]:4.1f}): {count:8d} ({pct:5.2f}%)")

# Top genes
print("\n" + "="*80)
print("GENES WITH LARGEST LOG FOLD CHANGES - EPOCH 100")
print("="*80)

print("\nTop 20 upregulated genes:")
top_up = de_last.nlargest(20, 'lfc_mean')[['celltype', 'gene', 'lfc_mean']]
print(top_up.to_string(index=False))

print("\nTop 20 downregulated genes:")
top_down = de_last.nsmallest(20, 'lfc_mean')[['celltype', 'gene', 'lfc_mean']]
print(top_down.to_string(index=False))

# Compare with epoch 50
print("\n" + "="*80)
print("COMPARISON: EPOCH 50 vs EPOCH 100")
print("="*80)

de_50 = de_long[de_long['epoch'] == 50].copy()

print(f"\nEpoch 50:")
print(f"  Range: {de_50['lfc_mean'].min():.2f} to {de_50['lfc_mean'].max():.2f}")
print(f"  Mean: {de_50['lfc_mean'].mean():.2f}")
print(f"  Std: {de_50['lfc_mean'].std():.2f}")

print(f"\nEpoch 100:")
print(f"  Range: {de_last['lfc_mean'].min():.2f} to {de_last['lfc_mean'].max():.2f}")
print(f"  Mean: {de_last['lfc_mean'].mean():.2f}")
print(f"  Std: {de_last['lfc_mean'].std():.2f}")

# Check genes with |LFC| > 2
print("\n" + "="*80)
print("STRONGLY DE GENES (|LFC| > 2) - EPOCH 100")
print("="*80)

strong_de = de_last[abs(de_last['lfc_mean']) > 2.0].copy()
print(f"\nFound {len(strong_de)} genes with |LFC| > 2")

if len(strong_de) > 0:
    print("\nBy cell type:")
    ct_counts = strong_de['celltype'].value_counts()
    print(ct_counts)

    print("\nAll genes with |LFC| > 2:")
    strong_de_sorted = strong_de.sort_values('abs_lfc', ascending=False)[['celltype', 'gene', 'lfc_mean']]
    print(strong_de_sorted.to_string(index=False))

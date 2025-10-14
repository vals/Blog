#!/usr/bin/env python3
"""
Check the distribution of p-values in epoch 100.
"""

import pandas as pd
import numpy as np

# Load DE results
print("Loading DE results...")
de_long = pd.read_parquet('results/de_results.parquet')

# Filter to epoch 100
de_100 = de_long[de_long['epoch'] == 100].copy()

print(f"Epoch 100 data shape: {de_100.shape}")

# The p-value column is 'proba_not_de' (probability of NOT being DE)
print("\n" + "="*80)
print("P-VALUE DISTRIBUTION (EPOCH 100)")
print("="*80)

print(f"\nColumn 'proba_not_de' statistics:")
print(f"  Min: {de_100['proba_not_de'].min():.2e}")
print(f"  Max: {de_100['proba_not_de'].max():.2e}")
print(f"  Mean: {de_100['proba_not_de'].mean():.2e}")
print(f"  Median: {de_100['proba_not_de'].median():.2e}")

print(f"\nPercentiles:")
for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
    val = de_100['proba_not_de'].quantile(p/100)
    print(f"  {p:2d}%: {val:.2e}")

# Count significant genes
print("\n" + "="*80)
print("SIGNIFICANCE THRESHOLDS")
print("="*80)

thresholds = [0.05, 0.01, 0.001, 0.0001, 1e-5, 1e-10, 1e-20]
for thresh in thresholds:
    count = (de_100['proba_not_de'] < thresh).sum()
    pct = count / len(de_100) * 100
    print(f"  p < {thresh:8.0e}: {count:7d} genes ({pct:5.1f}%)")

# Check how many are NOT significant
print("\n" + "="*80)
print("NON-SIGNIFICANT GENES (p >= 0.05)")
print("="*80)

non_sig = de_100[de_100['proba_not_de'] >= 0.05].copy()
print(f"\nTotal non-significant genes: {len(non_sig)} ({len(non_sig)/len(de_100)*100:.1f}%)")

if len(non_sig) > 0:
    print("\nBy cell type:")
    ct_counts = non_sig['celltype'].value_counts()
    print(ct_counts)

    print("\nSome examples of non-significant genes:")
    examples = non_sig.sample(min(20, len(non_sig)))[['celltype', 'gene', 'lfc_mean', 'proba_not_de', '-log10_pval']]
    print(examples.to_string(index=False))
else:
    print("\nAll genes are significant at p < 0.05!")

# Check the -log10_pval column
print("\n" + "="*80)
print("-LOG10(P-VALUE) DISTRIBUTION")
print("="*80)

print(f"\n-log10_pval statistics:")
print(f"  Min: {de_100['-log10_pval'].min():.2f}")
print(f"  Max: {de_100['-log10_pval'].max():.2f}")
print(f"  Mean: {de_100['-log10_pval'].mean():.2f}")
print(f"  Median: {de_100['-log10_pval'].median():.2f}")

# Note: -log10(0.05) = 1.30
print(f"\nNote: -log10(0.05) = {-np.log10(0.05):.2f}")
print(f"      So p < 0.05 corresponds to -log10(p) > 1.30")

count_sig = (de_100['-log10_pval'] > -np.log10(0.05)).sum()
print(f"\nGenes with -log10(p) > 1.30 (p < 0.05): {count_sig} ({count_sig/len(de_100)*100:.1f}%)")

# Distribution in bins
print("\n" + "="*80)
print("DISTRIBUTION BY -LOG10(P) BINS")
print("="*80)

bins = [0, 1.3, 2, 3, 5, 10, 20, 30, 1000]
bin_labels = ['<1.3 (p>0.05)', '1.3-2', '2-3', '3-5', '5-10', '10-20', '20-30', '>30']

for i in range(len(bins)-1):
    count = ((de_100['-log10_pval'] >= bins[i]) & (de_100['-log10_pval'] < bins[i+1])).sum()
    pct = count / len(de_100) * 100
    print(f"  {bin_labels[i]:15s}: {count:7d} ({pct:5.1f}%)")

#!/usr/bin/env python3
"""
Train SCVI model and save DE results for each epoch.

This script trains the model once and saves all DE results to disk,
allowing us to generate plots with consistent settings later.
"""

import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from pathlib import Path

# Set random seeds for reproducibility
np.random.seed(42)
scvi.settings.seed = 42

# File paths
DATA_PATH = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Single cell data/GSE235665/GSE235663/GSE235663.h5ad"
OUTPUT_DIR = Path("./results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("SCVI Training and DE Analysis")
print("="*80)

# Load data
print("\nLoading data...")
adata = sc.read_h5ad(DATA_PATH)
print(f"Data shape: {adata.shape}")

# Filter to paired UC patients
print("\nFiltering to paired UC patients...")
adata_paired = adata[adata.obs['condition'] == 'Paired'].copy()
print(f"Paired data shape: {adata_paired.shape}")

# Identify valid cell types
cell_counts = adata_paired.obs.groupby(['celltype_subset', 'group'], observed=True).size().unstack(fill_value=0)
tissue_cols = ['Inflamed', 'NonInvolved']
valid_celltypes = cell_counts[(cell_counts[tissue_cols] >= 50).all(axis=1)].index.tolist()
print(f"\nValid cell types ({len(valid_celltypes)}): {valid_celltypes}")

# Save metadata
metadata = {
    'valid_celltypes': valid_celltypes,
    'cell_counts': cell_counts,
    'n_cells_total': adata.n_obs,
    'n_cells_paired': adata_paired.n_obs,
    'n_genes': adata.n_vars
}
with open(OUTPUT_DIR / 'metadata.pkl', 'wb') as f:
    pickle.dump(metadata, f)
print(f"Saved metadata to {OUTPUT_DIR / 'metadata.pkl'}")


# Setup SCVI Model with default parameters
print("\n" + "="*80)
print("Setting up SCVI model (default parameters)...")
print("="*80)

scvi.model.SCVI.setup_anndata(adata)

# Use default parameters (don't specify n_layers or n_latent)
model = scvi.model.SCVI(adata, gene_likelihood="nb")

print(f"Model initialized with:")
print(f"  - Cells: {adata.n_obs}")
print(f"  - Genes: {adata.n_vars}")
print(f"  - Default architecture parameters")


# DE Analysis Function
def perform_de_analysis(model, adata_paired, valid_celltypes):
    """Perform DE analysis for each valid cell type."""
    de_results = {}

    for celltype in valid_celltypes:
        adata_ct = adata_paired[adata_paired.obs['celltype_subset'] == celltype].copy()

        try:
            de_df = model.differential_expression(
                adata=adata_ct,
                groupby='group',
                group1='Inflamed',
                group2='NonInvolved',
                mode='change'
            )

            de_df['-log10_pval'] = -np.log10(de_df['proba_not_de'] + 1e-10)
            de_results[celltype] = de_df

        except Exception as e:
            print(f"Error in DE analysis for {celltype}: {e}")
            continue

    return de_results


# Training Loop
print("\n" + "="*80)
print("Starting training (100 epochs)...")
print("="*80)

max_epochs = 100

# Storage for results
all_de_results = {}
training_history = {
    'epoch': [],
    'train_loss': [],
    'val_loss': []
}

# Initial DE analysis (epoch 0)
print(f"\nEpoch 0 - Performing DE analysis...")
de_results = perform_de_analysis(model, adata_paired, valid_celltypes)
all_de_results[0] = de_results

# Train epoch by epoch
for epoch in range(1, max_epochs + 1):
    model.train(
        max_epochs=1,
        train_size=0.9,
        check_val_every_n_epoch=1
    )

    train_loss = float(model.history['elbo_train'].iloc[-1].values[0])
    val_loss = float(model.history['elbo_validation'].iloc[-1].values[0]) if 'elbo_validation' in model.history else train_loss

    training_history['epoch'].append(epoch)
    training_history['train_loss'].append(train_loss)
    training_history['val_loss'].append(val_loss)

    # Perform DE analysis every epoch
    print(f"Epoch {epoch:3d}/{max_epochs} - Train: {train_loss:.2f}, Val: {val_loss:.2f} - DE...", end='')
    de_results = perform_de_analysis(model, adata_paired, valid_celltypes)
    all_de_results[epoch] = de_results
    print(" ✓")

print("\nTraining complete!")


# Save all results
print("\n" + "="*80)
print("Saving results...")
print("="*80)

# Convert DE results to long format DataFrame
print("Converting DE results to long format...")
de_records = []
for epoch, de_results in all_de_results.items():
    for celltype, de_df in de_results.items():
        df_temp = de_df.copy()
        df_temp['epoch'] = epoch
        df_temp['celltype'] = celltype
        df_temp['gene'] = df_temp.index
        de_records.append(df_temp)

de_long = pd.concat(de_records, ignore_index=True)

# Reorder columns for clarity
cols = ['epoch', 'celltype', 'gene'] + [c for c in de_long.columns if c not in ['epoch', 'celltype', 'gene']]
de_long = de_long[cols]

# Save as parquet
de_long.to_parquet(OUTPUT_DIR / 'de_results.parquet', index=False)
print(f"Saved DE results to {OUTPUT_DIR / 'de_results.parquet'}")
print(f"  - Format: Long format parquet")
print(f"  - Rows: {len(de_long):,}")
print(f"  - Columns: {list(de_long.columns)}")

# Save training history
history_df = pd.DataFrame(training_history)
history_df.to_csv(OUTPUT_DIR / 'training_history.csv', index=False)
print(f"Saved training history to {OUTPUT_DIR / 'training_history.csv'}")

# Save metadata
with open(OUTPUT_DIR / 'metadata.pkl', 'wb') as f:
    pickle.dump(metadata, f)
print(f"Saved metadata to {OUTPUT_DIR / 'metadata.pkl'}")

# Save model
model.save(OUTPUT_DIR / 'scvi_model', overwrite=True)
print(f"Saved model to {OUTPUT_DIR / 'scvi_model'}")

print("\n" + "="*80)
print("All results saved successfully!")
print("="*80)
print(f"\nResults directory: {OUTPUT_DIR}")
print(f"- de_results.parquet: DE results in long format")
print(f"- metadata.pkl: Cell type info and counts")
print(f"- training_history.csv: Training metrics")
print(f"- scvi_model/: Trained model")

#!/usr/bin/env python3
"""
SCVI DE Analysis over Training Epochs - TEST VERSION

Quick test with 10 epochs to verify everything works.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from pathlib import Path
from plotnine import *
from theme_nxn import theme_nxn, get_nxn_palette

# Set random seeds for reproducibility
np.random.seed(42)
scvi.settings.seed = 42

# File paths
DATA_PATH = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Single cell data/GSE235665/GSE235663/GSE235663.h5ad"
OUTPUT_DIR = Path("./output_test")
PLOT_DIR = OUTPUT_DIR / "volcano_plots"
PLOT_DIR.mkdir(parents=True, exist_ok=True)

print("Loading data...")
adata = sc.read_h5ad(DATA_PATH)

print(f"Data shape: {adata.shape}")

# Filter to paired UC patients
print("\nFiltering to paired UC patients...")
adata_paired = adata[adata.obs['condition'] == 'Paired'].copy()
print(f"Paired data shape: {adata_paired.shape}")

# Count cells per cell type and group
cell_counts = adata_paired.obs.groupby(['celltype_subset', 'group'], observed=True).size().unstack(fill_value=0)

# Identify valid cell types (≥50 cells in both tissues)
tissue_cols = ['Inflamed', 'NonInvolved']
valid_celltypes = cell_counts[(cell_counts[tissue_cols] >= 50).all(axis=1)].index.tolist()
print(f"\nValid cell types: {valid_celltypes}")
print(f"Number of valid cell types: {len(valid_celltypes)}")


# ============================================================================
# Setup SCVI Model on Full Dataset
# ============================================================================
print("\n" + "="*80)
print("Setting up SCVI model on full dataset...")
print("="*80)

scvi.model.SCVI.setup_anndata(adata)

model = scvi.model.SCVI(
    adata,
    n_layers=2,
    n_latent=30,
    gene_likelihood="nb"
)

print(f"Model initialized with {adata.n_obs} cells and {adata.n_vars} genes")


# ============================================================================
# Define DE Testing Function
# ============================================================================
def perform_de_analysis(model, adata_paired, valid_celltypes):
    """
    Perform differential expression analysis for each valid cell type.
    """
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
            de_df['celltype'] = celltype

            de_results[celltype] = de_df

        except Exception as e:
            print(f"Error in DE analysis for {celltype}: {e}")
            continue

    return de_results


# ============================================================================
# Define Volcano Plot Function
# ============================================================================
def create_volcano_plots(de_results, epoch, save_path):
    """
    Create multi-facet volcano plot for all cell types.
    """
    df_list = []
    for celltype, de_df in de_results.items():
        df_subset = de_df[['lfc_mean', '-log10_pval', 'celltype']].copy()
        df_subset = df_subset.reset_index()
        df_subset.columns = ['gene', 'log2FC', 'neg_log10_pval', 'celltype']
        df_list.append(df_subset)

    if not df_list:
        print(f"No DE results to plot for epoch {epoch}")
        return

    plot_df = pd.concat(df_list, ignore_index=True)

    p = (
        ggplot(plot_df, aes(x='log2FC', y='neg_log10_pval')) +
        geom_point(alpha=0.3, size=0.5, color=get_nxn_palette()[0]) +
        geom_hline(yintercept=-np.log10(0.05), linetype='dashed', color='red', alpha=0.5) +
        geom_vline(xintercept=0, linetype='dashed', color='gray', alpha=0.5) +
        facet_wrap('~celltype', ncol=3, scales='free') +
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
    print(f"Saved volcano plot: {save_path}")


# ============================================================================
# Training Loop - TEST VERSION (10 epochs)
# ============================================================================
print("\n" + "="*80)
print("Starting TEST training (10 epochs)...")
print("="*80)

max_epochs = 10
check_interval = 2

training_history = {
    'epoch': [],
    'train_loss': [],
    'val_loss': []
}

# Initial DE analysis (before training)
print(f"\nPerforming initial DE analysis (epoch 0)...")
de_results = perform_de_analysis(model, adata_paired, valid_celltypes)
save_path = PLOT_DIR / f"volcano_epoch_000.png"
create_volcano_plots(de_results, 0, save_path)

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

    print(f"Epoch {epoch}/{max_epochs} - Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}")

    if epoch % check_interval == 0 or epoch == max_epochs:
        print(f"  Performing DE analysis at epoch {epoch}...")
        de_results = perform_de_analysis(model, adata_paired, valid_celltypes)

        save_path = PLOT_DIR / f"volcano_epoch_{epoch:03d}.png"
        create_volcano_plots(de_results, epoch, save_path)

print("\nTest training complete!")


# ============================================================================
# Plot Training Curves
# ============================================================================
print("\n" + "="*80)
print("Plotting training curves...")
print("="*80)

history_df = pd.DataFrame(training_history)

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
        title='SCVI Training Curves (Test Run)',
        color='Metric'
    ) +
    theme_nxn() +
    theme(
        figure_size=(10, 6),
        legend_position='right'
    )
)

train_curve_path = OUTPUT_DIR / "training_curves.png"
p_train.save(train_curve_path, dpi=150, verbose=False)
print(f"Saved training curves: {train_curve_path}")

history_csv_path = OUTPUT_DIR / "training_history.csv"
history_df.to_csv(history_csv_path, index=False)
print(f"Saved training history: {history_csv_path}")


print("\n" + "="*80)
print("Test complete!")
print("="*80)
print(f"Volcano plots saved to: {PLOT_DIR}")
print(f"Training curves saved to: {train_curve_path}")

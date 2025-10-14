# SCVI DE Analysis Over Training Epochs - Project Plan

## Goal
Track how differential expression (DE) results change as an SCVI model trains on single-cell RNA-seq data.

## Data
**Source:** GSE235663.h5ad
**Location:** `/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Single cell data/GSE235665/GSE235663/GSE235663.h5ad`

**Dataset Details:**
- Total cells: 15,411
- Total genes: 14,203
- UC patients with paired Inflamed/NonInvolved tissue samples

## Analysis Approach

### Model Training
- **Model:** SCVI with default parameters (n_layers, n_latent, learning rate)
- **Training:** 100 epochs, one at a time
- **Data:** Train on entire dataset (all cells, all genes)
- **Train/Val Split:** 90/10 with validation tracked every epoch
- **Reproducibility:** Random seed set to 42

### DE Analysis
**Filtering:**
```python
adata_paired = adata[adata.obs['condition'] == 'Paired'].copy()
```

**Cell Type Selection:**
- Individual comparison for each celltype_subset
- Limit to cell types with ≥50 cells in both Inflamed and NonInvolved tissues

**Valid Cell Types (8 total):**
| Cell Type Subset | Inflamed | NonInvolved | Fold Change (I/NI) |
|-----------------|----------|-------------|-------------------|
| CD4~Stem-like | 270 | 111 | 2.43 |
| CD4~T[CM] | 225 | 299 | 0.75 |
| CD4~T[H]*17 | 59 | 207 | 0.29 |
| CD4~T[ITGA1] | 165 | 676 | 0.24 |
| CD4~T[REG] | 230 | 78 | 2.95 |
| CD8~T[C]*17 | 270 | 176 | 1.53 |
| CD8~T[GZMK] | 483 | 147 | 3.29 |
| CD8~T[RM] | 118 | 1,372 | 0.09 |

**Excluded (insufficient cells):**
- CD4~T[INT]: 274 Inflamed, 27 NonInvolved
- CD8~Stem-like: 97 Inflamed, 20 NonInvolved
- Cycling: 24 Inflamed, 9 NonInvolved

**DE Method:**
- SCVI differential_expression with mode='change'
- Comparison: Inflamed vs NonInvolved tissue
- Performed after every epoch (101 time points: epochs 0-100)

## Visualization

### Volcano Plots
- Multi-facet figure with 8 panels (one per cell type)
- Fixed coordinate limits across all epochs for smooth animation
- Limits calculated from 1st/99th percentiles with 10% padding
- Plotted with plotnine using custom theme_nxn styling

### Training Curves
- ELBO train and validation loss over epochs
- Separate curves for train vs validation

### Animation
- 101 frames (one per epoch)
- 24 fps for smooth, cinematic playback
- ~4.2 second duration
- Format: MP4 video

## Implementation

### Two-Step Workflow

**Step 1: Train and Save (`train_and_save_de.py`)**
- Train SCVI model for 100 epochs
- Perform DE analysis after every epoch
- Save results to parquet file (long format)
- Save training history (CSV)
- Save trained model

**Step 2: Generate Plots (`create_plots.py`)**
- Load saved DE results (parquet)
- Calculate global coordinate limits
- Generate 101 volcano plots with fixed axes
- Generate training curves
- All plots ready for animation

**Step 3: Create Animation (`create_animation.sh`)**
- Use ffmpeg to compile plots into MP4
- 24 fps for smooth playback
- High-quality H.264 encoding

### Output Structure
```
results/
├── de_results.parquet       # Long format: [epoch, celltype, gene, lfc_mean, ...]
├── metadata.pkl             # Cell type info, counts
├── training_history.csv     # Train/val loss per epoch
└── scvi_model/             # Saved SCVI model

plots/
├── volcano_epoch_000.png
├── volcano_epoch_001.png
├── ... (101 plots total)
├── training_curves.png
└── de_evolution.mp4         # Final animation
```

## Key Design Decisions

1. **Two-step workflow** - Train once, plot many times for rapid iteration
2. **Parquet format** - Portable, efficient, compatible with R/Python/Julia
3. **Fixed axes** - Essential for smooth animation across epochs
4. **Every epoch DE** - Maximum temporal resolution (101 time points)
5. **24 fps** - Standard cinematic frame rate for professional look
6. **Default parameters** - Use SCVI defaults for reproducibility
7. **Long format data** - Easy to filter, analyze, and visualize

## Benefits

- **Train once, plot many** - Adjust visualization without retraining (~1-2 hours)
- **Fast iteration** - Regenerate plots in seconds
- **Smooth animation** - Fixed axes and 24 fps ensure professional quality
- **Portable data** - Parquet files work across languages and tools
- **Reproducible** - Fixed seeds, default parameters, saved model
- **Flexible** - Easy to create subsets, change colors, adjust limits

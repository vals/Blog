# SCVI DE Analysis Over Training Epochs

This project tracks how differential expression (DE) results change as an SCVI model trains on single-cell RNA-seq data.

## Overview

The analysis:
1. Trains an SCVI model epoch-by-epoch on UC patient single-cell data
2. Saves DE results after each epoch (every 5 epochs)
3. Creates volcano plots with **fixed coordinate limits** for smooth animation
4. Generates training loss curves
5. Creates an MP4 animation showing DE evolution during training

## Workflow

### Step 1: Train Model and Save DE Results
```bash
python train_and_save_de.py
```

This trains the SCVI model for 100 epochs and saves:
- `results/de_results.parquet` - DE results in long format (portable, efficient)
- `results/metadata.pkl` - Cell type info and counts
- `results/training_history.csv` - Training metrics
- `results/scvi_model/` - Trained model

**Time:** ~1-2 hours for 100 epochs

### Step 2: Create Plots with Fixed Axes
```bash
python create_plots.py
```

This loads saved DE results and creates:
- `plots/volcano_epoch_*.png` - Volcano plots with **consistent axes**
- `plots/training_curves.png` - Training loss curves

The script calculates global coordinate limits from all epochs to ensure:
- Same x-axis (log2 fold change) across all plots
- Same y-axis (-log10 p-value) across all plots
- Smooth, professional-looking animation

### Step 3: Create Animation
```bash
brew install ffmpeg  # If not already installed
./create_animation.sh
```

Creates:
- `plots/de_evolution.mp4` - Smooth animation of DE results

## Files

**Main Scripts:**
- `train_and_save_de.py` - Train SCVI and save all DE results
- `create_plots.py` - Generate plots with fixed coordinate limits
- `create_animation.sh` - Create MP4 animation from plots

**Legacy Scripts (for reference):**
- `scvi_de_analysis.py` - Original all-in-one script
- `scvi_de_analysis_test.py` - Quick 10-epoch test version
- `check_history.py` - Utility to inspect model training history

## Data

**Input:** GSE235663.h5ad single-cell dataset

**Filtering:** UC patients with paired Inflamed/NonInvolved samples

**Valid cell types:** 8 cell type subsets with ≥50 cells in each tissue:
- CD4~Stem-like
- CD4~T[CM]
- CD4~T[H]*17
- CD4~T[ITGA1]
- CD4~T[REG]
- CD8~T[C]*17
- CD8~T[GZMK]
- CD8~T[RM]

## Model Details

- **Model:** SCVI (scvi-tools)
- **Parameters:** Default architecture (n_layers and n_latent)
- **Gene likelihood:** Negative binomial
- **Train/val split:** 90/10
- **Learning rate:** Default
- **Training:** 100 epochs total
- **DE Analysis:** Every 2 epochs (51 time points total)

## DE Analysis

- **Method:** SCVI differential_expression with mode='change'
- **Comparison:** Inflamed vs NonInvolved tissue
- **Visualization:** Multi-facet volcano plots (one panel per cell type)
- **Fixed axes:** Global coordinate limits computed from all results

## Output Structure

```
results/
├── de_results.parquet       # DE results in long format
├── metadata.pkl             # Cell type metadata
├── training_history.csv     # Training metrics
└── scvi_model/             # Saved SCVI model

plots/
├── volcano_epoch_000.png
├── volcano_epoch_002.png
├── volcano_epoch_004.png
├── ... (51 plots total)
├── training_curves.png
└── de_evolution.mp4         # Final animation
```

## Quick Start

For a complete analysis:

```bash
# 1. Train model and save results (~1-2 hours)
python train_and_save_de.py

# 2. Generate plots with fixed axes
python create_plots.py

# 3. Create animation (requires ffmpeg)
brew install ffmpeg
./create_animation.sh
```

## Why This Workflow?

The **two-step approach** (train → plot) offers several advantages:

1. **Train once, plot many times** - Adjust plotting parameters without retraining
2. **Fixed coordinate limits** - Ensures smooth, professional animations
3. **Faster iteration** - Regenerate plots in seconds instead of hours
4. **Flexibility** - Easy to customize plots, colors, or create subsets
5. **Portable format** - Parquet files can be read by R, Python, Julia, etc.

## Notes

- All random seeds set to 42 for reproducibility
- Volcano plots use plotnine with custom theme_nxn styling
- Animation runs at 5 fps (51 frames = ~10 second animation)
- Coordinate limits use 1st/99th percentiles with 10% padding
- DE analysis every 2 epochs for smooth animation

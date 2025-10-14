#!/usr/bin/env python3
"""Quick script to check what's in model.history after training."""

import numpy as np
import scanpy as sc
import scvi
from pathlib import Path

np.random.seed(42)
scvi.settings.seed = 42

DATA_PATH = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Single cell data/GSE235665/GSE235663/GSE235663.h5ad"

print("Loading data...")
adata = sc.read_h5ad(DATA_PATH)

print("Setting up SCVI...")
scvi.model.SCVI.setup_anndata(adata)

model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")

print("Training for 2 epochs...")
model.train(max_epochs=2)

print("\nModel history keys:")
print(model.history.keys())

print("\nModel history dataframe:")
print(model.history)

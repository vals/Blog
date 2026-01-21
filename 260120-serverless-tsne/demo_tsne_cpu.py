"""
Demo script timing CPU t-SNE (openTSNE) on a real dataset.
"""

import time
import numpy as np
import anndata as ad
from openTSNE import TSNE

H5AD_PATH = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Blog/250621 - Variational batch embeddings/scvi_output/variational_all_donors_20250618_055850.h5ad"


def main():
    # Load data
    print(f"Loading {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    X = np.asarray(adata.X, dtype=np.float32)
    print(f"Data shape: {X.shape}")

    # Full dataset
    print(f"\n--- Full dataset ({X.shape[0]:,} cells) ---")
    print("Using openTSNE with FFT, 12 cores")

    tsne = TSNE(
        n_components=2,
        perplexity=30,
        n_iter=1000,
        learning_rate=200,
        n_jobs=12,
        random_state=42,
    )

    start = time.perf_counter()
    coords = tsne.fit(X)
    full_time = time.perf_counter() - start
    print(f"t-SNE time: {full_time:.1f}s")
    print(f"Result shape: {coords.shape}")


if __name__ == "__main__":
    main()

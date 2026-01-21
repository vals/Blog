"""
Demo script timing GPU t-SNE on a real dataset.
"""

import time
import numpy as np
import anndata as ad

from gpu_embedding import gpu_tsne

H5AD_PATH = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Blog/250621 - Variational batch embeddings/scvi_output/variational_all_donors_20250618_055850.h5ad"


def main():
    # Load data
    print(f"Loading {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    X = np.asarray(adata.X, dtype=np.float32)
    print(f"Data shape: {X.shape}")

    # Warm up with small sample (first call has cold start overhead)
    print("\n--- Warmup (500 cells) ---")
    X_small = X[:500]
    start = time.perf_counter()
    _ = gpu_tsne(X_small, n_iter=250)
    warmup_time = time.perf_counter() - start
    print(f"Warmup time (includes cold start): {warmup_time:.1f}s")

    # Full dataset
    print(f"\n--- Full dataset ({X.shape[0]:,} cells) ---")
    start = time.perf_counter()
    coords = gpu_tsne(X)
    full_time = time.perf_counter() - start
    print(f"t-SNE time: {full_time:.1f}s")
    print(f"Result shape: {coords.shape}")

    # Save result
    adata.obsm["X_tsne_gpu"] = coords
    output_path = H5AD_PATH.replace(".h5ad", "_with_tsne.h5ad")
    adata.write_h5ad(output_path)
    print(f"\nSaved to: {output_path}")


if __name__ == "__main__":
    main()

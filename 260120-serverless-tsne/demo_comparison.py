"""
Compare GPU vs CPU t-SNE on 250k cells.
"""

import time
import numpy as np
import anndata as ad
from openTSNE import TSNE
from gpu_embedding import gpu_tsne

H5AD_PATH = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Blog/250621 - Variational batch embeddings/scvi_output/variational_all_donors_20250618_055850.h5ad"
N_SAMPLES = 250_000


def main():
    # Load and sample data
    print(f"Loading {H5AD_PATH}...")
    adata = ad.read_h5ad(H5AD_PATH)
    X_full = np.asarray(adata.X, dtype=np.float32)

    np.random.seed(42)
    idx = np.random.choice(X_full.shape[0], N_SAMPLES, replace=False)
    X = X_full[idx]
    print(f"Sampled {N_SAMPLES:,} cells, shape: {X.shape}")

    # Warm up GPU
    print("\n--- Warming up GPU ---")
    _ = gpu_tsne(X[:500], n_iter=250)
    print("Warm up complete")

    # GPU t-SNE
    print(f"\n--- GPU (Modal + RAPIDS cuML) ---")
    start = time.perf_counter()
    coords_gpu = gpu_tsne(X)
    gpu_time = time.perf_counter() - start
    print(f"Time: {gpu_time:.1f}s")

    # CPU t-SNE
    print(f"\n--- CPU (openTSNE, 12 cores) ---")
    tsne = TSNE(
        n_components=2,
        perplexity=30,
        n_iter=1000,
        learning_rate=200,
        n_jobs=12,
        random_state=42,
    )
    start = time.perf_counter()
    coords_cpu = tsne.fit(X)
    cpu_time = time.perf_counter() - start
    print(f"Time: {cpu_time:.1f}s")

    # Summary
    print(f"\n--- Summary ({N_SAMPLES:,} cells) ---")
    print(f"GPU: {gpu_time:.1f}s")
    print(f"CPU: {cpu_time:.1f}s")
    print(f"Speedup: {cpu_time / gpu_time:.1f}x")


if __name__ == "__main__":
    main()

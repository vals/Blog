# GPU vs CPU t-SNE Benchmark

Comparison of GPU-accelerated t-SNE (via Modal + RAPIDS cuML on T4) vs CPU t-SNE (openTSNE with FFT on 12 cores, Mac Mini M4 Pro).

## Dataset

- Source: `variational_all_donors_20250618_055850.h5ad`
- Full size: 1,058,909 cells × 10 features
- Data type: float32

## Results

| Cells | GPU (Modal T4) | CPU (12 cores) | Speedup |
|------:|---------------:|---------------:|--------:|
| 50,000 | 4.9s | 25.5s | 5.2x |
| 250,000 | 16.1s | 79.4s | 4.9x |
| 500,000 | 40.2s | 159.9s | 4.0x |
| 1,058,909 | 129s | 332s | 2.6x |

## Observations

1. **GPU speedup decreases with scale**: At 50k cells the GPU is 5.2x faster, but at 1M cells only 2.6x faster. This is due to network transfer overhead growing with data size (~40MB each way at 1M cells).

2. **Both scale roughly linearly**: CPU time roughly doubles when data doubles. GPU time increases faster due to network overhead.

3. **Practical threshold**: For datasets under ~100k cells, the GPU advantage is most pronounced (>5x). For million-cell datasets, the speedup is still meaningful (2-3x) but network transfer becomes a bottleneck.

4. **Cold start**: First GPU call includes ~10-15s cold start overhead (container spin-up). Subsequent calls within ~5 minutes avoid this.

## Cost

Modal T4 GPU: ~$0.59/hour

- 50k cells: ~$0.001
- 250k cells: ~$0.003
- 500k cells: ~$0.007
- 1M cells: ~$0.02

## Setup

- GPU: Modal serverless T4 with RAPIDS cuML (FFT method)
- CPU: Mac Mini M4 Pro, openTSNE with 12 threads
- t-SNE parameters: perplexity=30, n_iter=1000, learning_rate=200


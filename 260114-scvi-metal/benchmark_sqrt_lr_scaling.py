"""
Benchmark scVI training with sqrt learning rate scaling.
LR scales with square root of batch size: lr = base_lr * sqrt(batch_size / base_batch_size)
"""
import time
import math
import pandas as pd
import anndata as ad
import scvi

# Configuration
BATCH_SIZES = [128, 256, 512, 1024, 2048, 4096]
MAX_EPOCHS = 5
BASE_BATCH_SIZE = 128
BASE_LR = 1e-3  # scVI default

DATA_PATH = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Single cell data/OMIX002221/OMIX002221.h5ad"
OUTPUT_DIR = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Blog/260113 - d - SCVI joint embedding"

# Load data once
print("Loading data...")
adata = ad.read_h5ad(DATA_PATH)
n_cells = adata.n_obs
print(f"Data shape: {adata.shape} ({n_cells:,} cells)")

# Store results
all_results = []

for batch_size in BATCH_SIZES:
    # Sqrt LR scaling
    lr = BASE_LR * math.sqrt(batch_size / BASE_BATCH_SIZE)

    print(f"\n{'='*60}")
    print(f"Training with batch_size={batch_size}, lr={lr:.6f}")
    print('='*60)

    # Setup model fresh
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)

    # Train with timing
    start_time = time.time()
    model.train(
        max_epochs=MAX_EPOCHS,
        batch_size=batch_size,
        accelerator="mps",
        check_val_every_n_epoch=1,
        plan_kwargs={'lr': lr},
    )
    training_time = time.time() - start_time

    print(f"Training time: {training_time:.2f}s ({training_time/MAX_EPOCHS:.2f}s/epoch)")

    # Extract history
    history = model.history
    train_losses = history['train_loss']['train_loss'].values
    val_losses = history['validation_loss']['validation_loss'].values

    for epoch in range(MAX_EPOCHS):
        epoch_num = epoch + 1
        cells_processed = epoch_num * n_cells

        result = {
            'batch_size': batch_size,
            'lr': lr,
            'epoch': epoch_num,
            'cells_processed': cells_processed,
            'train_loss': train_losses[epoch],
            'validation_loss': val_losses[epoch],
            'total_training_time': training_time,
        }
        all_results.append(result)

# Create DataFrame and save
results_df = pd.DataFrame(all_results)
results_path = f"{OUTPUT_DIR}/benchmark_sqrt_lr_scaling_results.csv"
results_df.to_csv(results_path, index=False)
print(f"\nResults saved to: {results_path}")

# Print summary
print("\n" + "="*60)
print("TIMING & FINAL LOSS SUMMARY (with sqrt LR scaling)")
print("="*60)
for batch_size in BATCH_SIZES:
    subset = results_df[results_df['batch_size'] == batch_size]
    total_time = subset['total_training_time'].iloc[0]
    final_val_loss = subset['validation_loss'].iloc[-1]
    lr = subset['lr'].iloc[0]
    print(f"Batch {batch_size:4d} (lr={lr:.5f}): {total_time:6.2f}s, val_loss={final_val_loss:.1f}")

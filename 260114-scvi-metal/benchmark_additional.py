"""
Additional benchmark runs for batch sizes 128 and 256.
"""
import time
import pandas as pd
import anndata as ad
import scvi

# Configuration
BATCH_SIZES = [128, 256]
MAX_EPOCHS = 5
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
    print(f"\n{'='*60}")
    print(f"Training with batch_size={batch_size}")
    print('='*60)

    # Setup model (must be done fresh each time)
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)

    # Train with timing
    start_time = time.time()
    model.train(
        max_epochs=MAX_EPOCHS,
        batch_size=batch_size,
        accelerator="mps",
        check_val_every_n_epoch=1,
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
            'epoch': epoch_num,
            'cells_processed': cells_processed,
            'train_loss': train_losses[epoch],
            'validation_loss': val_losses[epoch],
            'total_training_time': training_time,
        }
        all_results.append(result)

# Create DataFrame with new results
new_results_df = pd.DataFrame(all_results)

# Load existing results and append
existing_results_path = f"{OUTPUT_DIR}/benchmark_results.csv"
existing_df = pd.read_csv(existing_results_path)
combined_df = pd.concat([existing_df, new_results_df], ignore_index=True)

# Save combined results
combined_df.to_csv(existing_results_path, index=False)
print(f"\nResults appended to: {existing_results_path}")

# Print timing summary for all batch sizes
print("\n" + "="*60)
print("TIMING SUMMARY (ALL BATCH SIZES)")
print("="*60)
timing = combined_df.groupby('batch_size')['total_training_time'].first().sort_index()
for bs, t in timing.items():
    print(f"Batch size {bs:4d}: {t:6.2f}s total, {t/MAX_EPOCHS:.2f}s/epoch")

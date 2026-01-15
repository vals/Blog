"""
Compare CPU vs MPS training with default settings.
"""
import time
import pandas as pd
import anndata as ad
import scvi

# Configuration
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

# Test both CPU and MPS with default batch size (128)
accelerators = [("cpu", "cpu"), ("mps", "mps")]

for accel_name, accel_value in accelerators:
    print(f"\n{'='*60}")
    print(f"Training with accelerator={accel_name} (default batch_size=128)")
    print('='*60)

    # Setup model fresh
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)

    # Train with timing (default batch_size=128)
    start_time = time.time()
    model.train(
        max_epochs=MAX_EPOCHS,
        accelerator=accel_value,
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
            'accelerator': accel_name,
            'epoch': epoch_num,
            'cells_processed': cells_processed,
            'train_loss': train_losses[epoch],
            'validation_loss': val_losses[epoch],
            'total_training_time': training_time,
        }
        all_results.append(result)

# Create DataFrame and save
results_df = pd.DataFrame(all_results)
results_path = f"{OUTPUT_DIR}/cpu_vs_mps_results.csv"
results_df.to_csv(results_path, index=False)
print(f"\nResults saved to: {results_path}")

# Print comparison
print("\n" + "="*60)
print("CPU vs MPS COMPARISON (default batch_size=128)")
print("="*60)
cpu_time = results_df[results_df['accelerator'] == 'cpu']['total_training_time'].iloc[0]
mps_time = results_df[results_df['accelerator'] == 'mps']['total_training_time'].iloc[0]

print(f"CPU:  {cpu_time:6.2f}s total, {cpu_time/MAX_EPOCHS:.2f}s/epoch")
print(f"MPS:  {mps_time:6.2f}s total, {mps_time/MAX_EPOCHS:.2f}s/epoch")
print(f"\nSpeedup: {cpu_time/mps_time:.2f}x faster with MPS")

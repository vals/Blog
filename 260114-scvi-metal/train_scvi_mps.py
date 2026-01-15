"""
Train scVI model with Apple Metal/MPS acceleration.
"""
import time
import anndata as ad
import scvi

# Load data
print("Loading data...")
data_path = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Single cell data/OMIX002221/OMIX002221.h5ad"
adata = ad.read_h5ad(data_path)
print(f"Data shape: {adata.shape}")
print(f"Observations: {adata.n_obs:,}")
print(f"Variables: {adata.n_vars:,}")

# Setup scVI
print("\nSetting up scVI model...")
scvi.model.SCVI.setup_anndata(adata)

# Create model
model = scvi.model.SCVI(adata)
print(f"Model architecture: {model.module}")

# Train with MPS and large batch size
# Large batch sizes work well with MPS
batch_size = 1024  # Large batch for MPS
max_epochs = 5

print(f"\nTraining configuration:")
print(f"  Device: MPS (Apple Metal)")
print(f"  Batch size: {batch_size}")
print(f"  Max epochs: {max_epochs}")

print("\nStarting training...")
start_time = time.time()

model.train(
    max_epochs=max_epochs,
    batch_size=batch_size,
    accelerator="mps",
)

end_time = time.time()
training_time = end_time - start_time

print(f"\nTraining completed!")
print(f"Total training time: {training_time:.2f} seconds ({training_time/60:.2f} minutes)")
print(f"Time per epoch: {training_time/max_epochs:.2f} seconds")

# Get latent representation
print("\nGetting latent representation...")
latent = model.get_latent_representation()
print(f"Latent shape: {latent.shape}")

# Save the model
model_path = "/Users/val/Library/CloudStorage/GoogleDrive-valentine.svensson@gmail.com/My Drive/Blog/260113 - d - SCVI joint embedding/scvi_model"
model.save(model_path, overwrite=True)
print(f"\nModel saved to: {model_path}")

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import torch
import os
import anndata
import scvi
import json
from datetime import datetime

def extract_variational_batch_embeddings():
    """Extract batch encodings from variational model using the dedicated batch encoder"""
    
    # Configuration - match the timestamp from extract_batch_embeddings_direct.py
    output_dir = "scvi_output"
    timestamp = "20250618_055850"  # Match the timestamp from embedding extraction
    
    print("Loading data and test set...")
    
    # Load full dataset
    adata = anndata.io.read_h5ad("data/2a99fd19-9a29-48c3-9d65-47467fd7cefe.h5ad")
    
    # Use raw counts from adata.raw.X (integer counts) instead of adata.X
    print("Using raw integer counts from adata.raw.X...")
    adata.X = adata.raw.X
    
    # Ensure donor_id is categorical for consistency
    adata.obs['donor_id'] = adata.obs['donor_id'].astype('category')
    
    # Load test set donors for labeling
    with open("test_set.txt", 'r') as f:
        test_donors = [line.strip() for line in f if line.strip()]
    
    print(f"Total donors: {adata.obs['donor_id'].nunique()}")
    print(f"Test donors: {len(test_donors)}")
    
    # Find variational model path
    model_dirs = [d for d in os.listdir(output_dir) if d.startswith('scvi_model_variational_')]
    if not model_dirs:
        print("No variational model found")
        return False
    
    latest_model_dir = sorted(model_dirs)[-1]
    variational_model_path = os.path.join(output_dir, latest_model_dir)
    
    print(f"Loading variational model from: {variational_model_path}")
    
    try:
        # Load the variational model with training data first (to match original setup)
        train_mask = ~adata.obs['donor_id'].isin(test_donors)
        adata_train = adata[train_mask].copy()
        adata_train.obs['donor_id'] = adata_train.obs['donor_id'].astype('category')
        
        variational_model = scvi.model.SCVI.load(variational_model_path, adata_train)
        print("Successfully loaded variational model")
        
        # Transfer to full dataset to handle unseen donors
        try:
            transferred_manager = variational_model.adata_manager.transfer_fields(adata, extend_categories=True)
        except TypeError:
            print("extend_categories parameter not available, trying alternative approach...")
            transferred_manager = variational_model.adata_manager.transfer_fields(adata)
        
        variational_model._adata_manager = transferred_manager
        variational_model._adata = adata
        print("Successfully transferred model to handle all donors (including held-out)")
        
        # Verify this is a variational batch representation model
        if not hasattr(variational_model.module, 'batch_encoder'):
            print("Error: This model does not have a batch encoder (not variational batch representation)")
            return False
        
        print("Confirmed: Model has batch encoder for variational batch representation")
        
        print("Creating pseudobulk profiles for all donors...")
        
        # Create pseudobulk representations for all donors (like VariationalBatchDataLoader does)
        all_donors = adata.obs['donor_id'].unique()
        pseudobulk_list = []
        donor_ids = []
        
        print(f"Processing {len(all_donors)} donors...")
        
        for i, donor in enumerate(all_donors):
            if i % 50 == 0:
                print(f"Processing donor {i+1}/{len(all_donors)}")
            
            donor_cells = adata[adata.obs['donor_id'] == donor]
            # Sum across cells for each gene (pseudobulk) - same as VariationalBatchDataLoader
            pseudobulk = np.array(donor_cells.X.sum(axis=0)).flatten()
            pseudobulk_list.append(pseudobulk)
            donor_ids.append(donor)
        
        pseudobulk_array = np.array(pseudobulk_list)
        print(f"Pseudobulk array shape: {pseudobulk_array.shape}")
        
        # Convert to torch tensor for batch encoder
        pseudobulk_tensor = torch.tensor(pseudobulk_array, dtype=torch.float32)
        
        # Move to same device as model
        device = next(variational_model.module.parameters()).device
        pseudobulk_tensor = pseudobulk_tensor.to(device)
        
        print(f"Pseudobulk tensor on device: {device}")
        
        # Extract batch encodings using the dedicated batch encoder
        print("Extracting variational batch encodings using batch encoder...")
        
        # Set model to evaluation mode
        variational_model.module.eval()
        
        with torch.no_grad():
            # Pass pseudobulk data through the batch encoder
            # This is the key: use the batch_encoder, not the main VAE encoder
            batch_qz, batch_z = variational_model.module.batch_encoder(pseudobulk_tensor)
            
            # Get the batch latent representations
            # batch_z is the sampled latent representation from the batch encoder
            batch_encodings_tensor = batch_z
            
            # Apply z_transformation if it exists (as in the generative method)
            if hasattr(variational_model.module.batch_encoder, 'z_transformation'):
                batch_encodings_tensor = variational_model.module.batch_encoder.z_transformation(batch_encodings_tensor)
        
        # Convert back to numpy
        batch_encodings = batch_encodings_tensor.cpu().numpy()
        print(f"Variational batch encodings shape: {batch_encodings.shape}")
        
        # Create labels for train vs test donors
        donor_split_labels = ['test' if donor in test_donors else 'train' for donor in donor_ids]
        is_held_out = [donor in test_donors for donor in donor_ids]
        
        n_train = sum(1 for label in donor_split_labels if label == 'train')
        n_test = sum(1 for label in donor_split_labels if label == 'test')
        
        print(f"Training batch encodings: {n_train}")
        print(f"Test (held-out) batch encodings: {n_test}")
        
        # Create DataFrame with all batch encodings
        variational_batch_df = pd.DataFrame(
            batch_encodings,
            index=donor_ids,
            columns=[f'variational_batch_{i}' for i in range(batch_encodings.shape[1])]
        )
        
        # Add metadata columns
        variational_batch_df['donor_split'] = donor_split_labels
        variational_batch_df['is_held_out'] = is_held_out
        
        print(f"Variational batch DataFrame shape: {variational_batch_df.shape}")
        
        # Save variational batch encodings (matching timestamp from embedding extraction)
        var_batch_file = f"{output_dir}/variational_batch_embeddings_{timestamp}.csv"
        variational_batch_df.to_csv(var_batch_file)
        print(f"Saved variational batch encodings to {var_batch_file}")
        
        # Also save as AnnData object for consistency
        adata_var_batch = anndata.AnnData(
            X=batch_encodings,
            obs=pd.DataFrame({
                'donor_id': donor_ids,
                'donor_split': donor_split_labels,
                'is_held_out': is_held_out
            }),
            var=pd.DataFrame(index=[f'variational_batch_{i}' for i in range(batch_encodings.shape[1])])
        )
        
        # Add comprehensive metadata
        adata_var_batch.uns['variational_batch_metadata'] = {
            'model_type': 'variational',
            'timestamp': timestamp,
            'model_path': variational_model_path,
            'encoding_dim': int(batch_encodings.shape[1]),
            'n_train_donors': int(n_train),
            'n_test_donors': int(n_test),
            'n_total_donors': int(len(donor_ids)),
            'extraction_method': 'batch_encoder_direct',
            'test_donors': test_donors
        }
        
        var_batch_h5ad_file = f"{output_dir}/variational_batch_embeddings_{timestamp}.h5ad"
        adata_var_batch.write(var_batch_h5ad_file)
        print(f"Saved variational batch encodings as AnnData to {var_batch_h5ad_file}")
        
        # Save separate training and test files for detailed analysis
        train_batch_df = variational_batch_df[variational_batch_df['donor_split'] == 'train'].copy()
        test_batch_df = variational_batch_df[variational_batch_df['donor_split'] == 'test'].copy()
        
        train_batch_file = f"{output_dir}/variational_batch_embeddings_train_{timestamp}.csv"
        test_batch_file = f"{output_dir}/variational_batch_embeddings_test_{timestamp}.csv"
        
        train_batch_df.to_csv(train_batch_file)
        test_batch_df.to_csv(test_batch_file)
        
        print(f"Saved training batch encodings to {train_batch_file}")
        print(f"Saved test batch encodings to {test_batch_file}")
        
        print(f"\n✅ Successfully extracted variational batch encodings!")
        print(f"   Combined shape: {batch_encodings.shape}")
        print(f"   Training donors: {n_train}")
        print(f"   Test (held-out) donors: {n_test}")
        print(f"   Files: {var_batch_file}, {var_batch_h5ad_file}")
        print(f"   Method: Batch encoder direct (CORRECT)")
        print(f"   Timestamp: {timestamp} (matched to embedding extraction)")
        
        return True
        
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = extract_variational_batch_embeddings()
    if success:
        print("\n✅ Variational batch encoding extraction completed successfully!")
    else:
        print("\n❌ Variational batch encoding extraction failed!")
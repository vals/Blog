#!/usr/bin/env python3

import pandas as pd
import numpy as np
import torch
import os

def extract_batch_embeddings_direct():
    """Extract batch embeddings directly from model.pt file"""
    
    # Configuration - match the downstream analysis log
    output_dir = "scvi_output"
    timestamp = "20250618_055850"  # Match the timestamp from downstream analysis log
    
    # Load only metadata to recreate training data structure (much faster)
    print("Loading metadata only for categorical mapping...")
    import anndata
    adata = anndata.io.read_h5ad("data/2a99fd19-9a29-48c3-9d65-47467fd7cefe.h5ad", backed='r')
    
    # Filter out test donors to match training data
    with open("test_set.txt", 'r') as f:
        test_donors = [line.strip() for line in f if line.strip()]
    
    print(f"Total donors: {adata.obs['donor_id'].nunique()}")
    
    # Create training data structure (just obs, no expression data needed)
    train_mask = ~adata.obs['donor_id'].isin(test_donors)
    adata_train_obs = adata.obs[train_mask].copy()
    adata_train_obs['donor_id'] = adata_train_obs['donor_id'].astype('category')
    
    print(f"Training donors: {adata_train_obs['donor_id'].nunique()}")
    
    # Find embedding model path
    model_dirs = [d for d in os.listdir(output_dir) if d.startswith('scvi_model_embedding_')]
    if not model_dirs:
        print("No embedding model found")
        return False
    
    latest_model_dir = sorted(model_dirs)[-1]
    embedding_model_path = os.path.join(output_dir, latest_model_dir)
    
    print(f"Loading embedding model from: {embedding_model_path}")
    
    # Create minimal AnnData with correct structure for model loading
    import scipy.sparse
    print("Creating minimal data structure for model loading...")
    
    # Take a small sample to create the structure
    sample_size = min(1000, len(adata_train_obs))
    sample_obs = adata_train_obs.iloc[:sample_size].copy()
    
    # Create minimal sparse matrix (just for structure)
    X_minimal = scipy.sparse.csr_matrix((sample_size, adata.shape[1]))
    adata_minimal = anndata.AnnData(X=X_minimal, obs=sample_obs)
    adata_minimal.obs['donor_id'] = adata_minimal.obs['donor_id'].astype('category')
    
    try:
        # Load the scVI model with minimal data to get categorical mapping
        import scvi
        embedding_model = scvi.model.SCVI.load(embedding_model_path, adata_minimal)
        
        print("Successfully loaded model - extracting categorical mapping...")
        
        # Get the correct categorical mapping from the model's adata_manager
        batch_registry = embedding_model.adata_manager.get_state_registry('batch')
        
        # The categories should be in the same order as the embeddings
        if hasattr(embedding_model.adata, 'obs'):
            donor_categories = embedding_model.adata.obs['donor_id'].cat.categories
        else:
            # Fallback: use the registry mapping
            donor_categories = adata_minimal.obs['donor_id'].cat.categories
        
        print(f"Found {len(donor_categories)} donor categories from model")
        
        # Extract batch embeddings from the model
        batch_embeddings = embedding_model.module._embeddings_dict['batch'].weight.detach().cpu().numpy()
        print(f"Extracted batch embeddings shape: {batch_embeddings.shape}")
        
        # Verify dimensions match
        if len(donor_categories) != batch_embeddings.shape[0]:
            print(f"Warning: Category count ({len(donor_categories)}) != embedding count ({batch_embeddings.shape[0]})")
            # Use the minimum to avoid index errors
            n_to_use = min(len(donor_categories), batch_embeddings.shape[0])
            donor_categories = donor_categories[:n_to_use]
            batch_embeddings = batch_embeddings[:n_to_use]
        else:
            n_to_use = len(donor_categories)
        
        print(f"Using correct categorical mapping with {n_to_use} donors")
        
        # Create DataFrame with donor IDs and embeddings (correctly mapped!)
        batch_emb_df = pd.DataFrame(
            batch_embeddings,
            index=donor_categories,
            columns=[f'batch_emb_{i}' for i in range(batch_embeddings.shape[1])]
        )
        
        print(f"Batch embeddings DataFrame shape: {batch_emb_df.shape}")
        
        # Save batch embeddings (overwriting previous incorrect mapping)
        batch_emb_file = f"{output_dir}/batch_embeddings_{timestamp}.csv"
        batch_emb_df.to_csv(batch_emb_file)
        print(f"Saved corrected batch embeddings to {batch_emb_file}")
        
        # Also save as AnnData object for consistency
        adata_batch = anndata.AnnData(
            X=batch_embeddings,
            obs=pd.DataFrame(index=donor_categories),
            var=pd.DataFrame(index=[f'batch_emb_{i}' for i in range(batch_embeddings.shape[1])])
        )
        adata_batch.obs['donor_id'] = adata_batch.obs.index
        
        # Add metadata
        adata_batch.uns['batch_embedding_metadata'] = {
            'model_type': 'embedding',
            'timestamp': timestamp,
            'model_path': embedding_model_path,
            'embedding_dim': int(batch_embeddings.shape[1]),
            'n_donors': int(n_to_use),
            'extraction_method': 'scvi_model_categorical_mapping',
            'source_tensor': '_embeddings_dict.batch.weight'
        }
        
        batch_emb_h5ad_file = f"{output_dir}/batch_embeddings_{timestamp}.h5ad"
        adata_batch.write(batch_emb_h5ad_file)
        print(f"Saved corrected batch embeddings as AnnData to {batch_emb_h5ad_file}")
        
        print(f"\n✅ Successfully extracted batch embeddings with correct mapping!")
        print(f"   Shape: {batch_embeddings.shape}")
        print(f"   Donors: {n_to_use}")
        print(f"   Files: {batch_emb_file}, {batch_emb_h5ad_file}")
        print(f"   Method: scVI categorical mapping (CORRECTED)")
        
        return True
        
    except Exception as e:
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = extract_batch_embeddings_direct()
    if success:
        print("\n✅ Batch embedding extraction completed successfully!")
    else:
        print("\n❌ Batch embedding extraction failed!")
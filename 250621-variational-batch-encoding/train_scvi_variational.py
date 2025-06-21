#!/usr/bin/env python3

import anndata
import pandas as pd
import scvi
import pickle
import json
from datetime import datetime
import argparse

def load_test_set_donors(test_set_file):
    """Load donor IDs to exclude from training"""
    with open(test_set_file, 'r') as f:
        test_donors = [line.strip() for line in f if line.strip()]
    return test_donors

def main():
    # Configuration
    input_file = "data/2a99fd19-9a29-48c3-9d65-47467fd7cefe.h5ad"
    test_set_file = "test_set.txt"
    output_dir = "scvi_output"
    batch_method = "variational"
    
    # Create output directory
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Loading data from {input_file}...")
    adata = anndata.io.read_h5ad(input_file)
    
    # Use raw counts from adata.raw.X (integer counts) instead of adata.X
    print("Using raw integer counts from adata.raw.X...")
    adata.X = adata.raw.X
    
    print(f"Original data shape: {adata.shape}")
    print(f"Unique donors: {adata.obs['donor_id'].nunique()}")
    
    # Load test set donors to exclude
    print(f"Loading test set donors from {test_set_file}...")
    test_donors = load_test_set_donors(test_set_file)
    print(f"Excluding {len(test_donors)} donors from training")
    
    # Filter out test set donors
    train_mask = ~adata.obs['donor_id'].isin(test_donors)
    adata_train = adata[train_mask].copy()
    
    print(f"Training data shape after filtering: {adata_train.shape}")
    print(f"Training donors: {adata_train.obs['donor_id'].nunique()}")
    
    # Setup scVI
    print("Setting up scVI model...")
    scvi.model.SCVI.setup_anndata(adata_train, batch_key='donor_id')
    
    # Create model with variational batch representation
    print(f"Creating model with {batch_method} batch representation...")
    model = scvi.model.SCVI(
        adata_train,
        n_layers=1,  # Changed from 2 to 1
        batch_representation='variational',
        batch_embedding_kwargs={"embedding_dim": 5},  # Changed from 2 to 5
        compute_pseudobulk=True,
        gene_likelihood='nb'
    )
    
    print("Starting training...")
    start_time = datetime.now()
    
    # Train model
    model.train(max_epochs=20, check_val_every_n_epoch=1)
    
    end_time = datetime.now()
    training_duration = end_time - start_time
    
    print(f"Training completed in {training_duration}")
    
    # Save training log
    training_log = {
        'batch_method': batch_method,
        'start_time': start_time.isoformat(),
        'end_time': end_time.isoformat(),
        'duration_seconds': training_duration.total_seconds(),
        'input_file': input_file,
        'test_set_file': test_set_file,
        'original_shape': adata.shape,
        'training_shape': adata_train.shape,
        'excluded_donors': test_donors,
        'training_donors': adata_train.obs['donor_id'].unique().tolist(),
        'model_params': {
            'n_layers': 1,
            'batch_representation': 'variational',
            'batch_embedding_dim': 5,
            'compute_pseudobulk': True,
            'gene_likelihood': 'nb',
            'max_epochs': 20
        },
        'note': 'Training completed successfully - history logging updated for current scVI version'
    }
    
    log_file = f"{output_dir}/training_log_{batch_method}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(log_file, 'w') as f:
        json.dump(training_log, f, indent=2)
    print(f"Training log saved to {log_file}")
    
    # Save model
    model_file = f"{output_dir}/scvi_model_{batch_method}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    model.save(model_file, overwrite=True)
    print(f"Model saved to {model_file}")
    
    print(f"Training script completed successfully using {batch_method} batch representation!")

if __name__ == "__main__":
    main()
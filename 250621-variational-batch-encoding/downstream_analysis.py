#!/usr/bin/env python3

import anndata
import pandas as pd
import scvi
import numpy as np
import json
from datetime import datetime
import os

def load_test_set_donors(test_set_file):
    """Load donor IDs to exclude from analysis"""
    with open(test_set_file, 'r') as f:
        test_donors = [line.strip() for line in f if line.strip()]
    return test_donors

def find_latest_model(output_dir, batch_method):
    """Find the most recent model for a given batch method"""
    model_dirs = [d for d in os.listdir(output_dir) if d.startswith(f'scvi_model_{batch_method}_')]
    if not model_dirs:
        raise FileNotFoundError(f"No model found for batch method: {batch_method}")
    # Sort by timestamp and get the latest
    latest_model = sorted(model_dirs)[-1]
    return os.path.join(output_dir, latest_model)

def main():
    # Configuration
    input_file = "data/2a99fd19-9a29-48c3-9d65-47467fd7cefe.h5ad"
    test_set_file = "test_set.txt"
    output_dir = "scvi_output"
    sample_size = 100000
    random_seed = 42
    
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
    print(f"Excluding {len(test_donors)} donors from analysis")
    
    # Filter out test set donors
    train_mask = ~adata.obs['donor_id'].isin(test_donors)
    adata_train = adata[train_mask].copy()
    
    # Ensure donor_id is categorical for consistency with training
    adata_train.obs['donor_id'] = adata_train.obs['donor_id'].astype('category')
    
    print(f"Training data shape after filtering: {adata_train.shape}")
    print(f"Training donors: {adata_train.obs['donor_id'].nunique()}")
    
    # Sample random cells
    np.random.seed(random_seed)
    if adata_train.shape[0] > sample_size:
        print(f"Sampling {sample_size} cells randomly...")
        sample_indices = np.random.choice(adata_train.shape[0], size=sample_size, replace=False)
        adata_sample = adata_train[sample_indices].copy()
    else:
        print(f"Using all {adata_train.shape[0]} cells (less than requested sample size)")
        adata_sample = adata_train.copy()
    
    print(f"Sample data shape: {adata_sample.shape}")
    print(f"Sample unique donors: {adata_sample.obs['donor_id'].nunique()}")
    
    # Ensure donor_id is categorical for consistency with training
    adata_sample.obs['donor_id'] = adata_sample.obs['donor_id'].astype('category')
    
    # Batch methods to process
    batch_methods = ['variational', 'embedding', 'one-hot']
    
    # Store all latent representations
    latent_results = {}
    
    for batch_method in batch_methods:
        print(f"\nProcessing {batch_method} model...")
        
        # Find latest model
        try:
            model_path = find_latest_model(output_dir, batch_method)
            print(f"Loading model from {model_path}")
            
            # Load the trained model with the new API to handle potential unseen categories
            # First load with training data to initialize adata_manager
            model = scvi.model.SCVI.load(model_path, adata_train)
            
            # Then transfer field mappings to sample data
            try:
                # Try with extend_categories=True parameter
                transferred_manager = model.adata_manager.transfer_fields(adata_sample, extend_categories=True)
            except TypeError:
                # If extend_categories is not a valid parameter, try without it
                print("extend_categories parameter not available, trying alternative approach...")
                transferred_manager = model.adata_manager.transfer_fields(adata_sample)
            
            model._adata_manager = transferred_manager
            model._adata = adata_sample
            
            # Get latent representation
            print(f"Generating latent representation for {batch_method}...")
            latent = model.get_latent_representation()
            
            # Evaluate training ELBO for fair comparison
            print(f"Computing training ELBO for {batch_method} model...")
            # Use same sample size as used for variational test evaluation 
            train_sample_size = 100000  # Standard sample size for ELBO comparison
            if len(adata_sample) < train_sample_size:
                train_sample_size = len(adata_sample)
            
            train_elbo = model.get_elbo(adata_sample)
            print(f"Training ELBO for {batch_method} (sample of {train_sample_size} cells): {train_elbo:.4f}")
            
            # Store results
            latent_results[batch_method] = {
                'latent': latent,
                'model_path': model_path,
                'shape': latent.shape,
                'train_elbo': float(train_elbo),
                'train_sample_size': train_sample_size
            }
            
            print(f"Latent representation shape for {batch_method}: {latent.shape}")
            
        except Exception as e:
            print(f"Error processing {batch_method} model: {str(e)}")
            continue
    
    # Save results
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Create AnnData object with latent representations
    for batch_method, result in latent_results.items():
        # Create new AnnData object with latent representation as X
        adata_latent = anndata.AnnData(
            X=result['latent'],
            obs=adata_sample.obs.copy(),
            var=pd.DataFrame(index=[f'latent_{i}' for i in range(result['latent'].shape[1])])
        )
        
        # Add metadata to uns
        adata_latent.uns['analysis_metadata'] = {
            'batch_method': batch_method,
            'timestamp': timestamp,
            'input_file': input_file,
            'test_set_file': test_set_file,
            'sample_size_requested': sample_size,
            'sample_size_actual': adata_sample.shape[0],
            'random_seed': random_seed,
            'model_path': result['model_path'],
            'latent_shape': list(result['shape'])  # Convert tuple to list for h5ad compatibility
        }
        
        # Save as h5ad file
        latent_file = f"{output_dir}/latent_{batch_method}_{timestamp}.h5ad"
        adata_latent.write(latent_file)
        print(f"Saved {batch_method} latent representation to {latent_file}")
    
    # Print training ELBO comparison summary
    if latent_results:
        print(f"\n=== Training ELBO Comparison ===")
        for method, result in latent_results.items():
            if 'train_elbo' in result:
                print(f"{method.capitalize():>12}: {result['train_elbo']:>8.4f}")
    
    # Save overall metadata
    metadata = {
        'timestamp': timestamp,
        'input_file': input_file,
        'test_set_file': test_set_file,
        'sample_size_requested': sample_size,
        'sample_size_actual': adata_sample.shape[0],
        'random_seed': random_seed,
        'excluded_donors': test_donors,
        'sample_donors': adata_sample.obs['donor_id'].unique().tolist(),
        'sample_shape': adata_sample.shape,
        'models_processed': {method: {'model_path': result['model_path'], 
                                    'latent_shape': list(result['shape']),
                                    'train_elbo': result.get('train_elbo'),
                                    'train_sample_size': result.get('train_sample_size')} 
                           for method, result in latent_results.items()}
    }
    
    metadata_file = f"{output_dir}/latent_metadata_{timestamp}.json"
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"Saved metadata to {metadata_file}")
    
    # Extract batch embeddings for embedding model
    if 'embedding' in latent_results:
        print(f"\nExtracting batch embeddings for embedding model...")
        try:
            # Load the embedding model with training data
            embedding_model_path = latent_results['embedding']['model_path']
            embedding_model = scvi.model.SCVI.load(embedding_model_path, adata_train)
            
            # Transfer to sample data
            try:
                emb_transferred_manager = embedding_model.adata_manager.transfer_fields(adata_sample, extend_categories=True)
            except TypeError:
                emb_transferred_manager = embedding_model.adata_manager.transfer_fields(adata_sample)
            
            embedding_model._adata_manager = emb_transferred_manager
            embedding_model._adata = adata_sample
            
            # Get unique training donors
            training_donors = adata_sample.obs['donor_id'].unique()
            print(f"Extracting batch embeddings for {len(training_donors)} training donors")
            
            # Extract batch embeddings from the model
            # The batch embeddings are stored in the model's batch embedding layer
            batch_embeddings = embedding_model.module.batch_embedding.embedding.weight.detach().cpu().numpy()
            
            # Create DataFrame with donor IDs and embeddings
            # Note: scVI uses categorical encoding, so we need to map back to donor IDs
            donor_mapping = adata_sample.obs['donor_id'].cat.categories
            
            batch_emb_df = pd.DataFrame(
                batch_embeddings,
                index=donor_mapping,
                columns=[f'batch_emb_{i}' for i in range(batch_embeddings.shape[1])]
            )
            
            # Filter to only training donors (should already be the case, but being explicit)
            batch_emb_df = batch_emb_df.loc[training_donors]
            
            print(f"Batch embeddings shape: {batch_emb_df.shape}")
            
            # Save batch embeddings
            batch_emb_file = f"{output_dir}/batch_embeddings_{timestamp}.csv"
            batch_emb_df.to_csv(batch_emb_file)
            print(f"Saved batch embeddings to {batch_emb_file}")
            
            # Also save as AnnData object for consistency
            adata_batch = anndata.AnnData(
                X=batch_embeddings,
                obs=pd.DataFrame(index=donor_mapping),
                var=pd.DataFrame(index=[f'batch_emb_{i}' for i in range(batch_embeddings.shape[1])])
            )
            adata_batch.obs['donor_id'] = adata_batch.obs.index
            
            # Add metadata
            adata_batch.uns['batch_embedding_metadata'] = {
                'model_type': 'embedding',
                'timestamp': timestamp,
                'model_path': embedding_model_path,
                'embedding_dim': int(batch_embeddings.shape[1]),
                'n_donors': int(len(training_donors))
            }
            
            batch_emb_h5ad_file = f"{output_dir}/batch_embeddings_{timestamp}.h5ad"
            adata_batch.write(batch_emb_h5ad_file)
            print(f"Saved batch embeddings as AnnData to {batch_emb_h5ad_file}")
            
        except Exception as e:
            print(f"Error extracting batch embeddings: {str(e)}")
    else:
        print("Embedding model not found, skipping batch embedding extraction")
    
    # Evaluate variational model on held-out test donors
    if 'variational' in latent_results:
        print(f"\nEvaluating variational model on held-out test donors...")
        try:
            # Load test set data (held-out donors)
            test_mask = adata.obs['donor_id'].isin(test_donors)
            adata_test = adata[test_mask].copy()
            
            print(f"Test data shape: {adata_test.shape}")
            print(f"Test donors: {adata_test.obs['donor_id'].nunique()}")
            
            # Ensure donor_id is categorical for consistency with training
            adata_test.obs['donor_id'] = adata_test.obs['donor_id'].astype('category')
            
            # Load the variational model using new API to handle unseen batch categories
            variational_model_path = latent_results['variational']['model_path']
            print(f"Loading model from {variational_model_path}")
            
            # First, try to load model with training data to get the original adata_manager
            # We need to use training data (without test donors) to initialize the model
            print("Loading model with training data to initialize adata_manager...")
            variational_model = scvi.model.SCVI.load(variational_model_path, adata_train)
            
            # Now transfer field mappings to handle unseen batch categories in test data
            print("Transferring field mappings for unseen batch categories...")
            try:
                # Try with extend_categories=True parameter
                transferred_manager = variational_model.adata_manager.transfer_fields(adata_test, extend_categories=True)
            except TypeError:
                # If extend_categories is not a valid parameter, try without it
                print("extend_categories parameter not available, trying alternative approach...")
                transferred_manager = variational_model.adata_manager.transfer_fields(adata_test)
            
            variational_model._adata_manager = transferred_manager
            variational_model._adata = adata_test
            print("Successfully configured model for test data with unseen donors")
            
            # Evaluate ELBO on test set
            print("Computing ELBO on held-out test donors...")
            test_elbo = variational_model.get_elbo(adata_test)
            
            print(f"Test set ELBO: {test_elbo:.4f}")
            
            # Also evaluate ELBO on training data for comparison
            print("Computing ELBO on training data for comparison...")
            # Switch model back to training data configuration
            train_transferred_manager = variational_model.adata_manager.transfer_fields(adata_train, extend_categories=True)
            variational_model._adata_manager = train_transferred_manager
            variational_model._adata = adata_train
            
            # Sample training data for faster ELBO computation (similar size to test set)
            train_sample_size = min(len(adata_test), len(adata_train))
            np.random.seed(random_seed)
            train_sample_indices = np.random.choice(len(adata_train), size=train_sample_size, replace=False)
            adata_train_sample = adata_train[train_sample_indices].copy()
            
            # Update model with training sample
            train_sample_manager = variational_model.adata_manager.transfer_fields(adata_train_sample, extend_categories=True)
            variational_model._adata_manager = train_sample_manager
            variational_model._adata = adata_train_sample
            
            train_elbo = variational_model.get_elbo(adata_train_sample)
            print(f"Training set ELBO (sample of {train_sample_size} cells): {train_elbo:.4f}")
            print(f"Test vs Training ELBO difference: {test_elbo - train_elbo:.4f}")
            
            # Switch back to test data for remaining operations
            test_transferred_manager = variational_model.adata_manager.transfer_fields(adata_test, extend_categories=True)
            variational_model._adata_manager = test_transferred_manager
            variational_model._adata = adata_test
            
            # Load training ELBO history for comparison
            # Extract full timestamp from model path (e.g., scvi_model_variational_20250617_073416 -> 20250617_073416)
            model_timestamp = variational_model_path.split('scvi_model_variational_')[-1]
            training_log_file = f"{output_dir}/training_log_variational_{model_timestamp}.json"
            
            try:
                with open(training_log_file, 'r') as f:
                    training_log = json.load(f)
                
                final_train_elbo = training_log['history']['elbo_train'][-1]
                final_val_elbo = training_log['history']['elbo_validation'][-1]
                
                print(f"Final training ELBO: {final_train_elbo:.4f}")
                print(f"Final validation ELBO: {final_val_elbo:.4f}")
                print(f"Test ELBO difference from training: {test_elbo - final_train_elbo:.4f}")
                print(f"Test ELBO difference from validation: {test_elbo - final_val_elbo:.4f}")
                
                # Save test evaluation results
                test_evaluation = {
                    'timestamp': timestamp,
                    'model_path': variational_model_path,
                    'test_elbo': float(test_elbo),
                    'current_train_elbo': float(train_elbo),
                    'current_train_vs_test_diff': float(train_elbo - test_elbo),
                    'final_train_elbo': final_train_elbo,
                    'final_validation_elbo': final_val_elbo,
                    'test_vs_train_diff': float(test_elbo - final_train_elbo),
                    'test_vs_val_diff': float(test_elbo - final_val_elbo),
                    'test_donors': test_donors,
                    'test_shape': adata_test.shape,
                    'train_sample_size': train_sample_size,
                    'n_test_donors': adata_test.obs['donor_id'].nunique()
                }
                
                test_eval_file = f"{output_dir}/test_evaluation_{timestamp}.json"
                with open(test_eval_file, 'w') as f:
                    json.dump(test_evaluation, f, indent=2)
                print(f"Saved test evaluation to {test_eval_file}")
                
            except FileNotFoundError:
                print(f"Warning: Could not find training log file {training_log_file}")
                
                # Save test ELBO and current training ELBO
                test_evaluation = {
                    'timestamp': timestamp,
                    'model_path': variational_model_path,
                    'test_elbo': float(test_elbo),
                    'current_train_elbo': float(train_elbo),
                    'current_train_vs_test_diff': float(train_elbo - test_elbo),
                    'test_donors': test_donors,
                    'test_shape': adata_test.shape,
                    'train_sample_size': train_sample_size,
                    'n_test_donors': adata_test.obs['donor_id'].nunique()
                }
                
                test_eval_file = f"{output_dir}/test_evaluation_{timestamp}.json"
                with open(test_eval_file, 'w') as f:
                    json.dump(test_evaluation, f, indent=2)
                print(f"Saved test evaluation to {test_eval_file}")
            
        except Exception as e:
            print(f"Error evaluating variational model on test set: {str(e)}")
    else:
        print("Variational model not found, skipping test evaluation")
    
    # Extract variational representations for ALL donors (training + test)
    if 'variational' in latent_results:
        print(f"\nExtracting variational representations for ALL donors...")
        try:
            # Use the full dataset (all donors)
            print(f"Full dataset shape: {adata.shape}")
            print(f"Total unique donors: {adata.obs['donor_id'].nunique()}")
            
            # Load the variational model with training data first
            variational_model_path = latent_results['variational']['model_path']
            variational_model_all = scvi.model.SCVI.load(variational_model_path, adata_train)
            
            # Transfer field mappings to handle all donors (including test donors)
            try:
                all_transferred_manager = variational_model_all.adata_manager.transfer_fields(adata, extend_categories=True)
            except TypeError:
                print("extend_categories parameter not available, trying alternative approach...")
                all_transferred_manager = variational_model_all.adata_manager.transfer_fields(adata)
            
            variational_model_all._adata_manager = all_transferred_manager
            variational_model_all._adata = adata
            
            # Get variational representation for all data
            print("Generating variational representation for all donors...")
            all_latent = variational_model_all.get_latent_representation()
            
            print(f"Variational representation shape for all donors: {all_latent.shape}")
            
            # Create AnnData object with variational representations
            adata_var_all = anndata.AnnData(
                X=all_latent,
                obs=adata.obs.copy(),
                var=pd.DataFrame(index=[f'variational_{i}' for i in range(all_latent.shape[1])])
            )
            
            # Add held-out status to obs
            adata_var_all.obs['is_held_out'] = adata_var_all.obs['donor_id'].isin(test_donors)
            adata_var_all.obs['donor_split'] = adata_var_all.obs['is_held_out'].map({True: 'test', False: 'train'})
            
            # Add summary statistics
            n_train_cells = (~adata_var_all.obs['is_held_out']).sum()
            n_test_cells = adata_var_all.obs['is_held_out'].sum()
            n_train_donors = adata_var_all.obs[~adata_var_all.obs['is_held_out']]['donor_id'].nunique()
            n_test_donors = adata_var_all.obs[adata_var_all.obs['is_held_out']]['donor_id'].nunique()
            
            print(f"Training: {n_train_cells} cells from {n_train_donors} donors")
            print(f"Test: {n_test_cells} cells from {n_test_donors} donors")
            
            # Add metadata to uns
            adata_var_all.uns['variational_all_metadata'] = {
                'batch_method': 'variational',
                'timestamp': timestamp,
                'model_path': variational_model_path,
                'latent_shape': list(all_latent.shape),
                'n_train_cells': int(n_train_cells),
                'n_test_cells': int(n_test_cells),
                'n_train_donors': int(n_train_donors),
                'n_test_donors': int(n_test_donors),
                'test_donors': test_donors,
                'input_file': input_file
            }
            
            # Save variational representations for all donors
            var_all_file = f"{output_dir}/variational_all_donors_{timestamp}.h5ad"
            adata_var_all.write(var_all_file)
            print(f"Saved variational representations for all donors to {var_all_file}")
            
            # Also create a donor-level summary
            donor_summary = adata_var_all.obs.groupby('donor_id').agg({
                'is_held_out': 'first',
                'donor_split': 'first'
            }).reset_index()
            
            donor_summary_file = f"{output_dir}/donor_split_summary_{timestamp}.csv"
            donor_summary.to_csv(donor_summary_file, index=False)
            print(f"Saved donor split summary to {donor_summary_file}")
            
        except Exception as e:
            print(f"Error extracting variational representations for all donors: {str(e)}")
    else:
        print("Variational model not found, skipping all-donor variational representation extraction")
    
    print(f"\nDownstream analysis completed!")
    print(f"Generated latent representations for {len(latent_results)} models:")
    for method, result in latent_results.items():
        print(f"  {method}: {result['shape']}")

if __name__ == "__main__":
    main()
"""
train_scvi.py
Train scVI model with negative binomial likelihood and epoch-level loss logging.
"""
import argparse
import json
import os
from datetime import datetime

import anndata as ad
import pandas as pd
import torch
from pytorch_lightning.callbacks import Callback
from scvi.model import SCVI

# Enable TF32 for faster training on Ampere+ GPUs (A100, L4, etc.)
torch.set_float32_matmul_precision("high")


class EpochLossLogger(Callback):
    """Custom callback to log losses at every epoch end."""

    def __init__(self, log_file: str = "epoch_losses.csv"):
        super().__init__()
        self.log_file = log_file
        self.epoch_data = []

    def on_train_epoch_end(self, trainer, pl_module):
        """Log training loss at epoch end."""
        epoch = trainer.current_epoch
        metrics = trainer.callback_metrics

        epoch_record = {
            "epoch": epoch,
            "train_loss": metrics.get("train_loss_epoch", None),
            "elbo_train": metrics.get("elbo_train", None),
            "reconstruction_loss_train": metrics.get("reconstruction_loss_train", None),
            "kl_local_train": metrics.get("kl_local_train", None),
        }

        # Convert tensors to floats
        for k, v in epoch_record.items():
            if hasattr(v, "item"):
                epoch_record[k] = v.item()

        self.epoch_data.append(epoch_record)
        pd.DataFrame(self.epoch_data).to_csv(self.log_file, index=False)

        train_loss = epoch_record.get("elbo_train")
        if train_loss:
            print(f"Epoch {epoch}: elbo_train={train_loss:.4f}")

    def on_validation_epoch_end(self, trainer, pl_module):
        """Log validation loss at epoch end."""
        epoch = trainer.current_epoch
        metrics = trainer.callback_metrics

        if self.epoch_data and self.epoch_data[-1]["epoch"] == epoch:
            self.epoch_data[-1].update(
                {
                    "val_loss": metrics.get("validation_loss", None),
                    "elbo_validation": metrics.get("elbo_validation", None),
                    "reconstruction_loss_validation": metrics.get(
                        "reconstruction_loss_validation", None
                    ),
                    "kl_local_validation": metrics.get("kl_local_validation", None),
                }
            )

            for k, v in self.epoch_data[-1].items():
                if hasattr(v, "item"):
                    self.epoch_data[-1][k] = v.item()

            pd.DataFrame(self.epoch_data).to_csv(self.log_file, index=False)

            val_loss = self.epoch_data[-1].get("elbo_validation")
            if val_loss:
                print(f"Epoch {epoch}: elbo_validation={val_loss:.4f}")


def train_scvi(
    adata_path: str,
    output_dir: str,
    layer: str = "counts_raw",
    size_factor_key: str = "total_counts",
    max_epochs: int = 400,
    n_latent: int = 10,
    n_hidden: int = 128,
    n_layers: int = 1,
    dropout_rate: float = 0.1,
    gene_likelihood: str = "nb",
    early_stopping: bool = True,
    early_stopping_patience: int = 45,
    batch_size: int = 128,
    train_size: float = 0.9,
):
    """Train scVI model with proper configuration and logging."""
    os.makedirs(output_dir, exist_ok=True)

    print(f"Loading data from {adata_path}...")
    adata = ad.read_h5ad(adata_path)
    print(f"Data shape: {adata.shape}")
    print(f"Using layer: {layer}")
    print(f"Using size factor key: {size_factor_key}")

    if layer not in adata.layers:
        raise ValueError(f"Layer '{layer}' not found. Available: {list(adata.layers.keys())}")

    if size_factor_key not in adata.obs.columns:
        raise ValueError(
            f"Size factor '{size_factor_key}' not found. Available: {list(adata.obs.columns)}"
        )

    print("Setting up AnnData for scVI...")
    SCVI.setup_anndata(
        adata,
        layer=layer,
        size_factor_key=size_factor_key,
    )

    print("Initializing scVI model...")
    model = SCVI(
        adata,
        n_latent=n_latent,
        n_hidden=n_hidden,
        n_layers=n_layers,
        dropout_rate=dropout_rate,
        gene_likelihood=gene_likelihood,
    )

    print(f"Model architecture:")
    print(f"  - n_latent: {n_latent}")
    print(f"  - n_hidden: {n_hidden}")
    print(f"  - n_layers: {n_layers}")
    print(f"  - gene_likelihood: {gene_likelihood}")
    print(f"  - n_genes: {adata.n_vars}")
    print(f"  - n_cells: {adata.n_obs}")

    epoch_logger = EpochLossLogger(log_file=os.path.join(output_dir, "epoch_losses.csv"))

    print(f"\nStarting training for {max_epochs} epochs...")
    print(f"  - batch_size: {batch_size}")
    print(f"  - train_size: {train_size}")
    print(f"  - early_stopping: {early_stopping}")

    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        train_size=train_size,
        early_stopping=early_stopping,
        early_stopping_patience=early_stopping_patience,
        early_stopping_monitor="elbo_validation",
        check_val_every_n_epoch=1,
        callbacks=[epoch_logger],
        enable_progress_bar=True,
    )

    model_path = os.path.join(output_dir, "scvi_model")
    print(f"\nSaving model to {model_path}...")
    model.save(model_path)

    history_path = os.path.join(output_dir, "training_history.json")
    history = {}
    for key in model.history.keys():
        val = model.history[key]
        history[key] = val.tolist() if hasattr(val, "tolist") else list(val)

    with open(history_path, "w") as f:
        json.dump(history, f, indent=2)
    print(f"Saved training history to {history_path}")

    config = {
        "adata_path": adata_path,
        "layer": layer,
        "size_factor_key": size_factor_key,
        "max_epochs": max_epochs,
        "n_latent": n_latent,
        "n_hidden": n_hidden,
        "n_layers": n_layers,
        "dropout_rate": dropout_rate,
        "gene_likelihood": gene_likelihood,
        "batch_size": batch_size,
        "train_size": train_size,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "actual_epochs": len(model.history.get("elbo_train", [])),
        "timestamp": datetime.now().isoformat(),
    }

    config_path = os.path.join(output_dir, "training_config.json")
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    print(f"Saved config to {config_path}")

    print("\n=== Training Complete ===")
    print(f"Output directory: {output_dir}")

    return model


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train scVI model")
    parser.add_argument("--adata", required=True, help="Path to h5ad file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--layer", default="counts_raw", help="Layer with raw counts")
    parser.add_argument(
        "--size-factor-key", default="total_counts", help="Size factor column in obs"
    )
    parser.add_argument("--max-epochs", type=int, default=400, help="Max training epochs")
    parser.add_argument("--n-latent", type=int, default=10, help="Latent dimensions")
    parser.add_argument("--n-hidden", type=int, default=128, help="Hidden layer size")
    parser.add_argument("--n-layers", type=int, default=1, help="Number of hidden layers")
    parser.add_argument("--batch-size", type=int, default=128, help="Batch size")
    parser.add_argument(
        "--gene-likelihood",
        default="nb",
        choices=["nb", "zinb", "poisson"],
        help="Gene likelihood distribution",
    )

    args = parser.parse_args()

    train_scvi(
        adata_path=args.adata,
        output_dir=args.output,
        layer=args.layer,
        size_factor_key=args.size_factor_key,
        max_epochs=args.max_epochs,
        n_latent=args.n_latent,
        n_hidden=args.n_hidden,
        n_layers=args.n_layers,
        batch_size=args.batch_size,
        gene_likelihood=args.gene_likelihood,
    )

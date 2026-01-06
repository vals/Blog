"""
prepare_data.py
Prepares AnnData with null gene injection for scVI training.
"""
import argparse

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, hstack


def add_null_gene(
    adata: ad.AnnData,
    layer: str = "counts_raw",
    null_gene_name: str = "NULL_CONTROL",
) -> ad.AnnData:
    """
    Add a synthetic null gene (all zeros) to the AnnData object.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    layer : str
        Name of the layer containing raw counts
    null_gene_name : str
        Name for the null control gene

    Returns
    -------
    AnnData with the null gene added to both .X and specified layer
    """
    n_cells = adata.n_obs

    # Create sparse zero column for the null gene
    null_column = csr_matrix((n_cells, 1), dtype=adata.layers[layer].dtype)

    # Add null column to the counts layer
    new_layer = hstack([adata.layers[layer], null_column], format="csr")

    # Add null column to .X if it exists
    if adata.X is not None:
        if hasattr(adata.X, "toarray"):
            new_X = hstack([adata.X, null_column], format="csr")
        else:
            new_X = np.hstack([adata.X, np.zeros((n_cells, 1))])
    else:
        new_X = None

    # Create new var DataFrame with null gene
    # Initialize with appropriate default values for each column type to avoid h5ad write issues
    null_gene_row = {}
    for col in adata.var.columns:
        dtype = adata.var[col].dtype
        if dtype == object or dtype.name == "category":
            null_gene_row[col] = ""  # Empty string for string/category columns
        elif np.issubdtype(dtype, np.bool_):
            null_gene_row[col] = False
        elif np.issubdtype(dtype, np.integer):
            null_gene_row[col] = 0
        else:
            null_gene_row[col] = 0.0  # Float for numeric columns

    new_gene_df = pd.DataFrame([null_gene_row], index=[null_gene_name])
    new_var = pd.concat([adata.var, new_gene_df], axis=0)

    # Create new AnnData
    adata_new = ad.AnnData(
        X=new_X,
        obs=adata.obs.copy(),
        var=new_var,
        uns=adata.uns.copy() if adata.uns else {},
        obsm=adata.obsm.copy() if adata.obsm else {},
        varm={},  # Reset varm as dimensions changed
        layers={layer: new_layer},
    )

    # Copy other layers if they exist
    for lname, ldata in adata.layers.items():
        if lname != layer:
            null_col_layer = csr_matrix((n_cells, 1), dtype=ldata.dtype)
            adata_new.layers[lname] = hstack([ldata, null_col_layer], format="csr")

    print(f"Added null gene '{null_gene_name}' at position {adata_new.n_vars - 1}")
    print(f"New shape: {adata_new.shape}")

    return adata_new


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add null control gene to AnnData")
    parser.add_argument("--input", required=True, help="Input h5ad file")
    parser.add_argument("--output", required=True, help="Output h5ad file")
    parser.add_argument("--layer", default="counts_raw", help="Layer with raw counts")
    args = parser.parse_args()

    print(f"Loading {args.input}...")
    adata = ad.read_h5ad(args.input)
    print(f"Original shape: {adata.shape}")

    adata = add_null_gene(adata, layer=args.layer)

    print(f"Saving to {args.output}...")
    adata.write_h5ad(args.output, compression="gzip")
    print("Done!")

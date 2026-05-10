"""
Per-cell-type differential expression: 10x 3' v3 vs 10x FRP
using edgePython's single-cell mixed model (NEBULA-LN).
"""

import anndata as ad
import edgepython as ep
import numpy as np
import pandas as pd

adata_full = ad.read_h5ad("PRJNA1106903_combined.h5ad")

# Use raw counts (X has normalized/scaled values)
adata = ad.AnnData(
    X=adata_full.raw.X,
    obs=adata_full.obs,
    var=adata_full.raw.var,
)

cell_types = adata.obs["predicted_celltype"].unique()
all_results = []

for ct in sorted(cell_types):
    mask = adata.obs["predicted_celltype"] == ct
    sub = adata[mask].copy()

    # Skip if any kit has < 20 cells
    kit_counts = sub.obs["kit"].value_counts()
    if kit_counts.min() < 20:
        print(f"Skipping {ct}: min kit count = {kit_counts.min()}")
        continue

    print(f"\n{'='*60}")
    print(f"{ct}: {sub.n_obs} cells across {sub.obs['kit'].nunique()} kits")
    print(f"  Kit sizes: {dict(kit_counts)}")

    # Design matrix: intercept + FRP indicator
    is_frp = (sub.obs["assay"] == "10x gene expression flex").astype(float).values
    design = pd.DataFrame(
        {"Intercept": np.ones(sub.n_obs), "FRP": is_frp},
        columns=["Intercept", "FRP"],
    )

    sample_ids = sub.obs["kit"].values

    # Fit NEBULA-LN mixed model
    fit = ep.glm_sc_fit(sub, design=design, sample=sample_ids, norm_method="TMM", verbose=False)

    # Empirical Bayes shrinkage of dispersion
    fit = ep.shrink_sc_disp(fit, robust=True)

    # Wald test on FRP coefficient
    result = ep.glm_sc_test(fit, coef=1)
    table = result["table"].copy()
    table["ave_log_abundance"] = fit["ave_log_abundance"]
    table["cell_type"] = ct
    table["gene"] = table.index

    n_sig = (table["FDR"] < 0.05).sum()
    n_total = len(table)
    print(f"  {n_sig}/{n_total} genes DE at FDR < 0.05")

    all_results.append(table)

# Combine and save
df = pd.concat(all_results, ignore_index=True)
df.to_csv("PRJNA1106903_de_results.csv", index=False)
print(f"\nSaved {len(df)} results to PRJNA1106903_de_results.csv")

# Summary
print("\n" + "=" * 60)
print("Summary: DE genes per cell type (FDR < 0.05)")
print("=" * 60)
summary = df.groupby("cell_type").apply(lambda x: pd.Series({
    "total_genes": len(x),
    "DE_genes": (x["FDR"] < 0.05).sum(),
    "up_in_FRP": ((x["FDR"] < 0.05) & (x["logFC"] > 0)).sum(),
    "down_in_FRP": ((x["FDR"] < 0.05) & (x["logFC"] < 0)).sum(),
})).astype(int)
print(summary.to_string())

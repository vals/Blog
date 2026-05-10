"""
MA plot: per-cell-type differential expression between 10x Flex and 10x 3'
Highlights cell type markers and immune response genes by process category.
"""

import pandas as pd
import numpy as np
from plotnine import *
from theme_nxn import theme_nxn, get_nxn_palette

# --- Load and filter DE results ---
df = pd.read_csv("PRJNA1106903_de_results.csv")
df = df[(df["converged"] == 1) & (df["ave_log_abundance"] > -0.8)].copy()

# --- Cell type markers (shown only in their own facet) ---
MARKER_GENES = {
    "B cell": ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5"],
    "CD4+ T cell": ["CD4", "IL7R", "CCR7", "LEF1", "TCF7"],
    "Cytotoxic T cell": ["CD8A", "GZMB", "NKG7", "CST7", "CD3E"],
    "CD14+ monocyte": ["CD14", "LYZ", "S100A9", "S100A8", "VCAN"],
    "CD16+ monocyte": ["FCGR3A", "MS4A7", "LST1", "CDKN1C", "LILRB2"],
    "Natural killer cell": ["GNLY", "NKG7", "KLRD1", "NCAM1", "KLRF1"],
    "Dendritic cell": ["FCER1A", "CD1C", "CLEC10A", "ITGAX", "FLT3"],
    "Plasmacytoid dendritic cell": ["CLEC4C", "IRF7", "IL3RA", "TCF4", "LILRA4"],
}
marker_pairs = set()
for ct, genes in MARKER_GENES.items():
    for g in genes:
        marker_pairs.add((ct, g))

# --- Immune response genes by process category ---
PROCESS_MAP = {
    # Interferon response (I + II merged)
    "IFITM1": "Interferon response", "IFITM3": "Interferon response",
    "ISG15": "Interferon response", "ISG20": "Interferon response",
    "MX1": "Interferon response", "MX2": "Interferon response",
    "OAS1": "Interferon response", "IFIT1": "Interferon response",
    "IFIT3": "Interferon response", "IRF1": "Interferon response",
    "STAT1": "Interferon response", "IRF8": "Interferon response",
    "GBP1": "Interferon response", "GBP5": "Interferon response",
    "IDO1": "Interferon response",
    # Antigen presentation
    "CD74": "Antigen presentation", "B2M": "Antigen presentation",
    "TAP1": "Antigen presentation", "TAP2": "Antigen presentation",
    "CIITA": "Antigen presentation",
    # Chemotaxis
    "CCL2": "Chemotaxis", "CCL5": "Chemotaxis", "CXCR4": "Chemotaxis",
    "CCR2": "Chemotaxis", "CX3CR1": "Chemotaxis",
    # T cell activation
    "CD69": "T cell activation", "CD28": "T cell activation",
    "ICOS": "T cell activation", "LAG3": "T cell activation",
    "PDCD1": "T cell activation", "TIGIT": "T cell activation",
    # Apoptosis
    "BCL2": "Apoptosis", "BAX": "Apoptosis", "FAS": "Apoptosis",
    "CASP3": "Apoptosis", "CASP8": "Apoptosis",
    # Other immune (pro-inflammatory, regulatory, cytotoxicity, innate sensing)
    "TNF": "Other immune", "IL1B": "Other immune", "IL6": "Other immune",
    "NFKBIA": "Other immune", "CXCL8": "Other immune",
    "IL10": "Other immune", "TGFB1": "Other immune", "FOXP3": "Other immune",
    "CTLA4": "Other immune", "IL2RA": "Other immune",
    "GZMA": "Other immune", "GZMB": "Other immune", "GZMK": "Other immune",
    "PRF1": "Other immune", "FASLG": "Other immune",
    "TLR2": "Other immune", "TLR4": "Other immune", "NLRP3": "Other immune",
    "CD163": "Other immune", "MARCO": "Other immune",
}

# --- Filter to cell types of interest ---
keep_ct = list(MARKER_GENES.keys())
df = df[df["cell_type"].isin(keep_ct)].copy()

# --- Classify each gene ---
def classify(row):
    if (row["cell_type"], row["gene_symbol"]) in marker_pairs:
        return "Cell type marker"
    elif row["gene_symbol"] in PROCESS_MAP:
        return PROCESS_MAP[row["gene_symbol"]]
    else:
        return "bg"

df["gene_set"] = df.apply(classify, axis=1)

# --- Split background into NS / DE ---
bg = df[df["gene_set"] == "bg"].copy()
bg["bg_class"] = "NS"
bg.loc[(bg["FDR"] < 0.05) & (bg["logFC"].abs() > 1), "bg_class"] = "DE"

# --- Highlighted genes and DE rings ---
highlighted = df[df["gene_set"] != "bg"].copy()
de_ring = highlighted[(highlighted["FDR"] < 0.05) & (highlighted["logFC"].abs() > 1)].copy()

cat_order = [
    "Cell type marker", "Antigen presentation", "Chemotaxis",
    "Interferon response", "T cell activation", "Apoptosis", "Other immune",
]
highlighted["gene_set"] = pd.Categorical(highlighted["gene_set"], categories=cat_order, ordered=True)
highlighted = highlighted.sort_values("gene_set")
de_ring["gene_set"] = pd.Categorical(de_ring["gene_set"], categories=cat_order, ordered=True)

# DE ring subset
de_ring = highlighted[(highlighted["FDR"] < 0.05) & (highlighted["logFC"].abs() > 1)].copy()

# DE cell-type markers — text labels
de_marker_labels = de_ring[de_ring["gene_set"] == "Cell type marker"].copy()

# Most extreme DE process genes — one per facet, biggest |logFC|, |logFC| > 1.7
EXTREME_DE = {
    ("B cell", "CCR2"),
    ("CD4+ T cell", "FOXP3"),
    ("Cytotoxic T cell", "CCR2"),
    ("CD14+ monocyte", "CXCR4"),
    ("Dendritic cell", "IRF8"),
    ("Dendritic cell", "IFITM3"),
    ("Plasmacytoid dendritic cell", "CCR2"),
}
extreme_de_labels = de_ring[
    de_ring.apply(lambda r: (r["cell_type"], r["gene_symbol"]) in EXTREME_DE, axis=1)
].copy()

# Combine cell-type markers + extreme process genes
all_labels = pd.concat([de_marker_labels, extreme_de_labels], ignore_index=True)
all_labels["label_x"] = all_labels["ave_log_abundance"]
all_labels["label_y"] = all_labels["logFC"] + np.where(
    all_labels["logFC"] > 0, 0.55, -0.55
)

# Manual overrides where points are close enough that labels collide
LABEL_OVERRIDES = {
    ("CD4+ T cell", "CD4"): (1.04, 2.05),
    ("CD4+ T cell", "TCF7"): (2.14, 1.55),
    ("CD4+ T cell", "FOXP3"): (-0.17, 3.1),
    ("Plasmacytoid dendritic cell", "CLEC4C"): (3.6, 2.05),
    ("Plasmacytoid dendritic cell", "IRF7"): (5.4, 1.85),
    ("Plasmacytoid dendritic cell", "CCR2"): (1.99, 3.1),
    ("Dendritic cell", "FCER1A"): (4.3, -2.4),
    ("Dendritic cell", "CLEC10A"): (1.8, -2.3),
    ("Dendritic cell", "IFITM3"): (-0.5, -3.4),
    ("B cell", "CCR2"): (-0.77, 2.6),
    ("Cytotoxic T cell", "CCR2"): (-0.41, 3.1),
}
for i, row in all_labels.iterrows():
    key = (row["cell_type"], row["gene_symbol"])
    if key in LABEL_OVERRIDES:
        all_labels.at[i, "label_x"] = LABEL_OVERRIDES[key][0]
        all_labels.at[i, "label_y"] = LABEL_OVERRIDES[key][1]

# --- Colors ---
palette = get_nxn_palette()
set_colors = {
    "Cell type marker": palette[0],
    "Antigen presentation": "#E64B35",
    "Chemotaxis": "#3C8DAD",
    "Interferon response": "#D95F02",
    "T cell activation": "#5E4FA2",
    "Apoptosis": "#1B9E77",
    "Other immune": "#B07AA1",
}
bg_colors = {"NS": "#D5D5D5", "DE": "#999999"}

# --- Direction labels and circle explanation (first facet only) ---
label_df = pd.DataFrame({
    "cell_type": ["B cell", "B cell", "B cell"],
    "ave_log_abundance": [0.0, 0.0, 4.8],
    "logFC": [4.2, -4.2, -4.2],
    "label": ["Higher in 10x Flex", "Higher in 10x 3'", "\u25C9 DE   \u25CF Not DE"],
})

# --- Build plot ---
p = (
    ggplot()
    + geom_hline(yintercept=0, linetype="dashed", color="grey", size=0.5)
    + geom_point(
        data=bg,
        mapping=aes(x="ave_log_abundance", y="logFC", color="bg_class"),
        size=0.5, raster=True,
    )
    + geom_point(
        data=highlighted,
        mapping=aes(x="ave_log_abundance", y="logFC", color="gene_set"),
        size=2.5,
    )
    + geom_point(
        data=de_ring,
        mapping=aes(x="ave_log_abundance", y="logFC", color="gene_set"),
        size=5.5, fill="none", stroke=1.0,
    )
    + geom_text(
        data=label_df,
        mapping=aes(x="ave_log_abundance", y="logFC", label="label"),
        ha="left", va="center", size=10, color="#666666", fontstyle="italic",
        family="DejaVu Sans",
    )
    + geom_text(
        data=all_labels,
        mapping=aes(x="label_x", y="label_y", label="gene_symbol", color="gene_set"),
        size=9, fontweight="bold", family="DejaVu Sans",
        show_legend=False,
    )
    + facet_wrap("~cell_type", ncol=4)
    + scale_color_manual(
        values={**bg_colors, **set_colors},
        breaks=cat_order,
    )
    + scale_x_continuous(limits=(-1, 9))
    + scale_y_continuous(limits=(-5, 5))
    + labs(
        x="Mean log abundance",
        y="log FC (10x Flex vs 10x 3')",
        color="",
    )
    + theme_nxn()
    + theme(
        figure_size=(12, 7),
        strip_text=element_text(size=12),
        axis_text=element_text(size=10),
        axis_title=element_text(size=12),
        legend_text=element_text(size=10),
        legend_position="top",
    )
)

p.save("PRJNA1106903_marker_ma.png", dpi=300)
print("Saved PRJNA1106903_marker_ma.png")

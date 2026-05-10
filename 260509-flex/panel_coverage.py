"""
Flex panel coverage analysis.

Quantifies how well the 10x Flex probe panel covers curated immune-relevant
gene sets, for a blog post comparing Flex to 10x 3'.
"""

import warnings
from pathlib import Path

import anndata as ad
import gseapy as gp
import mygene
import numpy as np
import pandas as pd
import requests
from plotnine import *
from theme_nxn import get_nxn_palette, theme_nxn

warnings.filterwarnings("ignore", category=FutureWarning)

OUTDIR = Path("results")
OUTDIR.mkdir(exist_ok=True)

# ── 1. Load probe panel from FRP h5ad files ─────────────────────────────────

print("Loading Flex panel gene lists from h5ad files...")
frp1 = ad.read_h5ad("10X_FRP-rep1.h5ad", backed="r")
frp2 = ad.read_h5ad("10X_FRP-rep2.h5ad", backed="r")
flex_panel = set(frp1.var["gene_name"]) | set(frp2.var["gene_name"])

three1 = ad.read_h5ad("10X_3-rep1.h5ad", backed="r")
three2 = ad.read_h5ad("10X_3-rep2.h5ad", backed="r")
three_prime_genes = set(three1.var["gene_name"]) | set(three2.var["gene_name"])

print(f"  Flex panel: {len(flex_panel):,} genes")
print(f"  3' reference: {len(three_prime_genes):,} genes")
print()

# ── 2. Fetch gene sets ──────────────────────────────────────────────────────

gene_sets: dict[str, dict[str, list[str]]] = {}

# --- MSigDB Hallmark ---
print("Fetching MSigDB Hallmark gene sets...")
hallmark_lib = gp.get_library("MSigDB_Hallmark_2020", organism="Human")

HALLMARK_MAP = {
    "Inflammatory Response": "HALLMARK_INFLAMMATORY_RESPONSE",
    "Interferon Alpha Response": "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "Interferon Gamma Response": "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "TNF-alpha Signaling via NF-kB": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "Complement": "HALLMARK_COMPLEMENT",
    "IL-6/JAK/STAT3 Signaling": "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "IL-2/STAT5 Signaling": "HALLMARK_IL2_STAT5_SIGNALING",
    "Allograft Rejection": "HALLMARK_ALLOGRAFT_REJECTION",
    "Apoptosis": "HALLMARK_APOPTOSIS",
}

gene_sets["Hallmark"] = {}
for friendly_name, canonical_name in HALLMARK_MAP.items():
    if friendly_name in hallmark_lib:
        gene_sets["Hallmark"][canonical_name] = hallmark_lib[friendly_name]
        print(f"  {canonical_name}: {len(hallmark_lib[friendly_name])} genes")
    else:
        print(f"  WARNING: {friendly_name} not found in library")

# --- PanglaoDB ---
print("\nFetching PanglaoDB markers...")
pangla_url = "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
pangla = pd.read_csv(pangla_url, sep="\t")

# Filter: immune/blood, human, canonical
pangla_immune = pangla[
    (pangla["organ"].isin(["Immune system", "Blood"]))
    & (pangla["species"].isin(["Hs", "Mm Hs"]))
].copy()

# Use canonical_marker if available
if "canonical marker" in pangla_immune.columns:
    pangla_immune = pangla_immune[pangla_immune["canonical marker"] == 1]

WANTED_CELL_TYPES = [
    "B cells",
    "T cells",
    "NK cells",
    "Monocytes",
    "Dendritic cells",
    "Plasmacytoid dendritic cells",
    "Neutrophils",
    "Basophils",
    "Eosinophils",
    "Mast cells",
    "Plasma cells",
    "Megakaryocytes",
    "Macrophages",
    "T regulatory cells",
    "Gamma delta T cells",
]

gene_sets["PanglaoDB"] = {}
for ct in WANTED_CELL_TYPES:
    subset = pangla_immune[pangla_immune["cell type"] == ct]
    if len(subset) == 0:
        print(f"  WARNING: {ct} not found in PanglaoDB")
        continue
    # If >50 markers, take top 25 by ubiquitousness_index
    if len(subset) > 50 and "ubiquitousness index" in subset.columns:
        subset = subset.nlargest(25, "ubiquitousness index")
    markers = sorted(subset["official gene symbol"].unique().tolist())
    gene_sets["PanglaoDB"][ct] = markers
    print(f"  {ct}: {len(markers)} markers")

# --- HGNC gene families ---
print("\nFetching HGNC gene families...")
HGNC_HEADERS = {"Accept": "application/json"}


def fetch_hgnc_group(group_id: int) -> list[str]:
    r = requests.get(
        f"https://rest.genenames.org/fetch/gene_group_id/{group_id}",
        headers=HGNC_HEADERS,
    )
    r.raise_for_status()
    docs = r.json()["response"]["docs"]
    return sorted(d["symbol"] for d in docs)


# Chemokine ligands (group 483)
chemokine_ligands = fetch_hgnc_group(483)

# Chemokine receptors: CC (1091), CXC (1094), CX3C (1093), XC (1092), atypical (1090)
chemokine_receptors = []
for gid in [1091, 1094, 1093, 1092, 1090]:
    chemokine_receptors.extend(fetch_hgnc_group(gid))
chemokine_receptors = sorted(set(chemokine_receptors))

# Cytokines: interleukins (601) + interferons (598) + TNF superfamily (781) + colony-stimulating factors
interleukins = fetch_hgnc_group(601)
interferons = fetch_hgnc_group(598)
tnf_ligands = fetch_hgnc_group(781)

# Core cytokines not in the above families
extra_cytokines = [
    "CSF1",
    "CSF2",
    "CSF3",
    "TGFB1",
    "TGFB2",
    "TGFB3",
    "LIF",
    "OSM",
    "KITLG",
    "FLT3LG",
    "TPO",
    "EPO",
    "TSLP",
    "IL14",
]
cytokines = sorted(set(interleukins + interferons + tnf_ligands + extra_cytokines))
# Filter out pseudogenes (ending in P + digits)
import re

cytokines = [g for g in cytokines if not re.match(r".*P\d*$", g) or g in ("TSLP",)]

gene_sets["HGNC"] = {
    "Cytokines": cytokines,
    "Chemokine ligands": chemokine_ligands,
    "Chemokine receptors": chemokine_receptors,
}
for name, genes in gene_sets["HGNC"].items():
    print(f"  {name}: {len(genes)} genes")

# --- Biologic drug targets ---
print("\nDefining biologic drug targets...")
drug_targets = [
    "CTLA4",
    "PDCD1",
    "CD274",
    "PDCD1LG2",
    "LAG3",
    "TIGIT",
    "HAVCR2",
    "TNFRSF9",
    "TNFRSF4",
    "CD40",
    "CD40LG",
    "IL17A",
    "IL17F",
    "IL23A",
    "IL12B",
    "IL4",
    "IL5",
    "IL13",
    "IL6",
    "IL6R",
    "IL1B",
    "IL1R1",
    "IL1RN",
    "TNF",
    "TNFRSF1A",
    "JAK1",
    "JAK2",
    "JAK3",
    "TYK2",
    "MS4A1",
    "CD19",
    "CD22",
    "CD38",
    "CD52",
    "TNFRSF17",
    "TNFSF13B",
    "ITGA4",
    "ITGB7",
    "CCR4",
    "CXCR4",
    "S1PR1",
    "ITK",
    "BTK",
    "SYK",
]
gene_sets["Drug targets"] = {"Biologic drug targets": drug_targets}
print(f"  Biologic drug targets: {len(drug_targets)} genes")

# ── 3. Match gene sets to panel & alias resolution ──────────────────────────

print("\nMatching gene sets to Flex panel...")

all_coverage_rows = []
all_gene_status_rows = []
all_unmatched_rows = []


def classify_gene(gene: str, panel: set[str]) -> str:
    if gene in panel:
        return "in_panel"
    return "not_in_panel"


for group, sets in gene_sets.items():
    for set_name, genes in sets.items():
        n_total = len(genes)
        statuses = {g: classify_gene(g, flex_panel) for g in genes}
        unmatched = [g for g, s in statuses.items() if s == "not_in_panel"]
        pct_unmatched = len(unmatched) / n_total if n_total > 0 else 0

        # Try alias resolution if >5% unmatched
        # Conservative: only accept if the query is a known *previous* symbol
        # for the target gene (not just any alias hit — avoids false matches
        # like VIP -> CPAMD8 where VIP is a real gene, just absent).
        resolved_aliases = {}
        if pct_unmatched > 0.05 and len(unmatched) > 0:
            mg = mygene.MyGeneInfo()
            results = mg.querymany(
                unmatched, scopes="alias,retired", species="human", fields="symbol"
            )
            # Also check which query symbols are valid current genes
            valid_check = mg.querymany(
                unmatched, scopes="symbol", species="human", fields="symbol"
            )
            valid_current_symbols = {
                r["query"]
                for r in valid_check
                if "symbol" in r and r["symbol"] == r["query"]
            }
            for r in results:
                if "symbol" in r and r["symbol"] != r["query"]:
                    # Skip if the query is itself a valid current gene symbol
                    if r["query"] in valid_current_symbols:
                        continue
                    canonical = r["symbol"]
                    if canonical in flex_panel:
                        resolved_aliases[r["query"]] = canonical
                        statuses[r["query"]] = "in_panel"

        n_in_panel = sum(1 for s in statuses.values() if s == "in_panel")

        # Also check 3' coverage
        n_in_3prime = sum(1 for g in genes if g in three_prime_genes)

        all_coverage_rows.append(
            {
                "group": group,
                "set_name": set_name,
                "n_total": n_total,
                "n_in_panel": n_in_panel,
                "pct_in_panel": round(100 * n_in_panel / n_total, 1) if n_total else 0,
                "n_in_3prime": n_in_3prime,
                "pct_in_3prime": (
                    round(100 * n_in_3prime / n_total, 1) if n_total else 0
                ),
            }
        )

        for gene in genes:
            in_3p = gene in three_prime_genes
            all_gene_status_rows.append(
                {
                    "group": group,
                    "set_name": set_name,
                    "gene": gene,
                    "flex_status": statuses[gene],
                    "in_3prime": in_3p,
                    "resolved_alias": resolved_aliases.get(gene, ""),
                }
            )

        for gene in unmatched:
            if gene not in resolved_aliases:
                all_unmatched_rows.append(
                    {
                        "group": group,
                        "set_name": set_name,
                        "gene": gene,
                        "resolved_via_alias": False,
                        "alias_target": "",
                    }
                )
            else:
                all_unmatched_rows.append(
                    {
                        "group": group,
                        "set_name": set_name,
                        "gene": gene,
                        "resolved_via_alias": True,
                        "alias_target": resolved_aliases[gene],
                    }
                )

        print(
            f"  {group} / {set_name}: {n_in_panel}/{n_total} "
            f"({all_coverage_rows[-1]['pct_in_panel']}%) in Flex, "
            f"{n_in_3prime}/{n_total} ({all_coverage_rows[-1]['pct_in_3prime']}%) in 3'"
        )

# ── 4. Write CSV outputs ────────────────────────────────────────────────────

print("\nWriting output files...")

coverage_df = pd.DataFrame(all_coverage_rows)
coverage_df.to_csv(OUTDIR / "coverage_table.csv", index=False)

gene_status_df = pd.DataFrame(all_gene_status_rows)
gene_status_df.to_csv(OUTDIR / "gene_status.csv", index=False)

unmatched_df = pd.DataFrame(all_unmatched_rows)
unmatched_df.to_csv(OUTDIR / "unmatched.csv", index=False)

# ── 5. V(D)J gene check ────────────────────────────────────────────────────

print("\nV(D)J segment check:")
vdj_prefixes = ["IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ", "TRAV", "TRAJ", "TRBV", "TRBJ"]
vdj_in_panel = [g for g in flex_panel if any(g.startswith(p) for p in vdj_prefixes)]
vdj_in_3prime = [g for g in three_prime_genes if any(g.startswith(p) for p in vdj_prefixes)]
print(f"  V(D)J genes in Flex panel: {len(vdj_in_panel)}")
print(f"  V(D)J genes in 3' reference: {len(vdj_in_3prime)}")
if vdj_in_panel:
    print(f"  Unexpected Flex V(D)J genes: {sorted(vdj_in_panel)}")

# Also check constant-region IG/TR genes
ig_const = [g for g in flex_panel if re.match(r"^(IGH[ADEGM]|IGK[C]|IGL[CV]\d|TRA[CV]|TRB[CV]|TRD[CV]|TRG[CV])", g)]
print(f"  IG/TR constant-region genes in Flex panel: {len(ig_const)}")
if ig_const:
    print(f"    {sorted(ig_const)}")

# ── 6. Plot ─────────────────────────────────────────────────────────────────

print("\nGenerating coverage plot...")

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

plot_df = coverage_df.copy()

# Clean up set names for display
LABEL_MAP = {
    "HALLMARK_INFLAMMATORY_RESPONSE": "Inflammatory response",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN-\u03b1 response",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN-\u03b3 response",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "TNF\u03b1 / NF-\u03baB",
    "HALLMARK_COMPLEMENT": "Complement",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": "IL-6 / JAK / STAT3",
    "HALLMARK_IL2_STAT5_SIGNALING": "IL-2 / STAT5",
    "HALLMARK_ALLOGRAFT_REJECTION": "Allograft rejection",
    "HALLMARK_APOPTOSIS": "Apoptosis",
    "Biologic drug targets": "Drug targets",
}
plot_df["label"] = plot_df["set_name"].map(lambda x: LABEL_MAP.get(x, x))

# Build ordered list: within each group, sort by Flex pct ascending
group_order = ["Hallmark", "PanglaoDB", "HGNC", "Drug targets"]
ordered_rows = []
for g in group_order:
    grp = plot_df[plot_df["group"] == g].sort_values("pct_in_panel", ascending=False)
    ordered_rows.append(grp)
plot_df = pd.concat(ordered_rows).reset_index(drop=True)

pal = get_nxn_palette()
color_flex = pal[0]
color_3p = pal[2] if len(pal) > 2 else "#999999"

n_sets = len(plot_df)
bar_height = 0.35
fig_height = max(6, n_sets * 0.38 + 1.5)

fig, ax = plt.subplots(figsize=(6.5, fig_height))

# Compute y positions with gaps between groups
y_positions = []
y = 0
prev_group = None
group_label_positions = {}  # group -> (y_start, y_end) for bracket
for i, row in plot_df.iterrows():
    if prev_group is not None and row["group"] != prev_group:
        y += 0.8  # gap between groups
    if row["group"] not in group_label_positions:
        group_label_positions[row["group"]] = [y, y]
    group_label_positions[row["group"]][1] = y
    y_positions.append(y)
    prev_group = row["group"]
    y += 1

y_positions = np.array(y_positions)

# Draw bars — 3' behind (wider position), Flex in front
ax.barh(
    y_positions + bar_height / 2,
    plot_df["pct_in_3prime"],
    height=bar_height,
    color=color_3p,
    alpha=0.8,
    label="3\u2032",
    zorder=2,
)
ax.barh(
    y_positions - bar_height / 2,
    plot_df["pct_in_panel"],
    height=bar_height,
    color=color_flex,
    alpha=0.85,
    label="Flex",
    zorder=2,
)

# Y-axis labels
ax.set_yticks(y_positions)
ax.set_yticklabels(plot_df["label"], fontsize=8)
ax.invert_yaxis()

# X-axis
ax.set_xlim(0, 105)
ax.xaxis.set_major_locator(mticker.MultipleLocator(25))
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v)}%"))
ax.set_xlabel("Genes in panel", fontsize=10)

# Group labels on right side
for g, (y_start, y_end) in group_label_positions.items():
    y_mid = (y_start + y_end) / 2
    ax.text(
        107, y_mid, g, fontsize=8, fontstyle="italic",
        va="center", ha="left", color="#555555",
    )

# Light grid
ax.set_axisbelow(True)
ax.grid(axis="x", color="#e0e0e0", linewidth=0.5, zorder=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Legend
ax.legend(loc="lower right", frameon=False, fontsize=9)

fig.tight_layout()
fig.savefig(OUTDIR / "coverage_plot.pdf", dpi=300, bbox_inches="tight")
fig.savefig(OUTDIR / "coverage_plot.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  Saved coverage_plot.pdf and coverage_plot.png")

# ── 6b. Dumbbell plot — compact blog version ────────────────────────────────

print("\nGenerating dumbbell plot...")

dumb_df = coverage_df.copy()
dumb_df["label"] = dumb_df["set_name"].map(lambda x: LABEL_MAP.get(x, x))
dumb_df["gap"] = dumb_df["pct_in_3prime"] - dumb_df["pct_in_panel"]

# Select the ~10 most interesting sets: biggest gap OR lowest Flex coverage
# Take union of bottom-8 by Flex coverage and top-8 by gap, deduplicate
bottom_flex = dumb_df.nsmallest(8, "pct_in_panel")
top_gap = dumb_df.nlargest(8, "gap")
selected = pd.concat([bottom_flex, top_gap]).drop_duplicates(subset="set_name")
# Sort by Flex coverage ascending
selected = selected.sort_values("pct_in_panel", ascending=True).reset_index(drop=True)

n = len(selected)
from matplotlib.path import Path as MplPath

# Half-circle markers for overlapping dots
_t_left = np.linspace(np.pi / 2, 3 * np.pi / 2, 40)
_vl = list(zip(np.cos(_t_left), np.sin(_t_left))) + [(0, -1), (0, 1)]
marker_left = MplPath(_vl)

_t_right = np.linspace(-np.pi / 2, np.pi / 2, 40)
_vr = list(zip(np.cos(_t_right), np.sin(_t_right))) + [(0, 1), (0, -1)]
marker_right = MplPath(_vr)

OVERLAP_THRESHOLD = 3  # pp — use half-circles when gap is smaller than this

fig, ax = plt.subplots(figsize=(8, 3.5))

y = np.arange(n)

for i, (_, row) in enumerate(selected.iterrows()):
    gap = abs(row["pct_in_3prime"] - row["pct_in_panel"])

    if gap >= OVERLAP_THRESHOLD:
        # Connecting line
        ax.plot(
            [row["pct_in_panel"], row["pct_in_3prime"]],
            [i, i],
            color="#cccccc", linewidth=2.5, zorder=1,
        )
        # Full dots
        ax.scatter(row["pct_in_3prime"], i, color=color_3p, s=55, zorder=2)
        ax.scatter(row["pct_in_panel"], i, color=color_flex, s=55, zorder=3)
    else:
        # Half-circle pair: Flex on left, 3' on right
        midx = (row["pct_in_panel"] + row["pct_in_3prime"]) / 2
        ax.scatter(
            midx, i, marker=marker_left, s=55,
            color=color_flex, zorder=3,
        )
        ax.scatter(
            midx, i, marker=marker_right, s=55,
            color=color_3p, zorder=3,
        )

# Direct labels: pick the row with the widest gap for clear placement
widest_idx = selected["gap"].idxmax()
widest_pos = selected.index.get_loc(widest_idx)
ref_row = selected.loc[widest_idx]
ax.annotate(
    "Flex", (ref_row["pct_in_panel"], widest_pos),
    textcoords="offset points", xytext=(0, 11),
    ha="center", va="bottom", fontsize=8.5, fontweight="bold", color=color_flex,
)
ax.annotate(
    "3\u2032", (ref_row["pct_in_3prime"], widest_pos),
    textcoords="offset points", xytext=(0, 11),
    ha="center", va="bottom", fontsize=8.5, fontweight="bold", color=color_3p,
)

# Y-axis: label with group name in grey
y_labels = []
for _, row in selected.iterrows():
    y_labels.append(row["label"])
ax.set_yticks(y)
ax.set_yticklabels(y_labels, fontsize=9)
ax.invert_yaxis()

ax.set_xlim(50, 105)
ax.xaxis.set_major_locator(mticker.MultipleLocator(10))
ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v)}%"))
ax.set_xlabel("Genes in data", fontsize=10)

ax.set_axisbelow(True)
ax.grid(axis="x", color="#e8e8e8", linewidth=0.5, zorder=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig(OUTDIR / "dumbbell_plot.pdf", dpi=300, bbox_inches="tight")
fig.savefig(OUTDIR / "dumbbell_plot.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  Saved dumbbell_plot.pdf and dumbbell_plot.png")

# ── 7. README ───────────────────────────────────────────────────────────────

print("\nGenerating README...")

median_flex = coverage_df["pct_in_panel"].median()
worst = coverage_df.loc[coverage_df["pct_in_panel"].idxmin()]
best = coverage_df.loc[coverage_df["pct_in_panel"].idxmax()]

# Drug targets missing
drug_status = gene_status_df[gene_status_df["group"] == "Drug targets"]
missing_drugs = drug_status[drug_status["flex_status"] == "not_in_panel"]["gene"].tolist()

# PanglaoDB markers missing from Flex
pangla_status = gene_status_df[gene_status_df["group"] == "PanglaoDB"]
missing_markers = pangla_status[pangla_status["flex_status"] == "not_in_panel"]
# Highlight specific cell types with notable gaps
notable_missing = {}
for ct in WANTED_CELL_TYPES:
    ct_missing = missing_markers[missing_markers["set_name"] == ct]["gene"].tolist()
    if ct_missing:
        notable_missing[ct] = ct_missing

# Big gaps between Flex and 3'
coverage_df["gap"] = coverage_df["pct_in_3prime"] - coverage_df["pct_in_panel"]
big_gaps = coverage_df[coverage_df["gap"] > 5].sort_values("gap", ascending=False)

readme_lines = [
    "# Flex Panel Coverage Analysis — Results",
    "",
    f"**Median coverage across all gene sets:** {median_flex:.1f}% of genes present in the Flex panel.",
    "",
    f"**Best-covered set:** {best['set_name']} ({best['pct_in_panel']:.0f}% Flex, {best['pct_in_3prime']:.0f}% 3')",
    f"**Worst-covered set:** {worst['set_name']} ({worst['pct_in_panel']:.0f}% Flex, {worst['pct_in_3prime']:.0f}% 3')",
    "",
]

if len(big_gaps) > 0:
    readme_lines.append("**Sets with large Flex-vs-3' gap (>5 pp):**")
    for _, row in big_gaps.iterrows():
        readme_lines.append(
            f"- {row['set_name']}: {row['pct_in_panel']:.0f}% Flex vs {row['pct_in_3prime']:.0f}% 3' (gap {row['gap']:.0f} pp)"
        )
    readme_lines.append("")

if missing_drugs:
    readme_lines.append(f"**Drug targets absent from Flex panel ({len(missing_drugs)}):** {', '.join(sorted(missing_drugs))}")
    readme_lines.append("")

if notable_missing:
    readme_lines.append("**Notable missing PanglaoDB cell-type markers:**")
    for ct, genes in sorted(notable_missing.items()):
        if len(genes) <= 8:
            readme_lines.append(f"- {ct}: {', '.join(sorted(genes))}")
        else:
            readme_lines.append(f"- {ct}: {len(genes)} markers missing (see gene_status.csv)")
    readme_lines.append("")

readme_lines.append(f"**V(D)J genes:** {len(vdj_in_panel)} in Flex panel, {len(vdj_in_3prime)} in 3' — expected absent from Flex as V(D)J segments are not targeted by the probe panel.")
readme_lines.append("")
readme_lines.append(f"*Analysis run on {pd.Timestamp.now().strftime('%Y-%m-%d')}. Flex panel: {len(flex_panel):,} genes; 3' reference: {len(three_prime_genes):,} genes.*")

readme_text = "\n".join(readme_lines) + "\n"
(OUTDIR / "README.md").write_text(readme_text)
print("  Saved README.md")

print("\nDone! All outputs in ./results/")

#!/usr/bin/env python
"""Additional markers of the Th+/Adrb2+ neuron population.

Two analyses of the same population:

1. Correlation. Genes are correlated with Adrb2 across Th+ neurons, restricted to cells
   in which both the gene and Adrb2 are detected. This is the analysis that identified
   Cdh12 (Spearman rho = 0.419) and the other reported markers.

2. Differential expression. Th+ neurons are split by Adrb2 detection and compared with a
   Wilcoxon rank-sum test. That test runs over the full ~128k-gene matrix, so its result
   is precomputed by ``08_export_precomputed.py`` and read back here.

    python scripts/06_th_adrb2_markers.py
    python scripts/06_th_adrb2_markers.py --strict-filter

Outputs:
    results/th_adrb2_correlation_table.csv   reported markers and their correlations
    results/th_adrb2_marker_table.csv        differential expression selection
    figures/volcano_th_adrb2.png             volcano plot
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

# --- Correlation analysis -----------------------------------------------------------

GENE1 = "Th"
GENE2 = "Adrb2"

#: Markers of the Th+/Adrb2+ population reported in the paper.
MARKER_GENES = [
    "Cdh12", "Auts2l1", "Cntnap5c", "Akap2", "Pcdhac2", "Samd5", "Chsy3", "Supt3h",
]

# --- Selection thresholds, as used for the reported table --------------------------

#: Minimum -log10 adjusted p-value.
PADJ_THR = 3.0
#: Minimum absolute log fold change.
FC_THR = 1.5
#: Number of genes reported in each direction.
N_LABEL = 12

#: Predicted-gene models carry no interpretable annotation, so they are excluded.
EXCLUDE_GENE_REGEX = r"^(?:Gm\d+)"

#: The wider exclusion set: predicted genes, BAC clone identifiers, pseudogenes,
#: ribosomal, mitochondrial and haemoglobin genes. Applied with --strict-filter.
STRICT_EXCLUDE_REGEX = (
    r"^(?:Gm\d+|RP\d+-|Rpl|Rps|mt-|Hb[abpdqz]|Rn\d)|(?:-ps\d*$)|(?:Rik$)"
)

GROUP_POS = "Th+/Adrb2+"
GROUP_NEG = "Th+/Adrb2-"


def _expression(adata: sc.AnnData, gene: str) -> np.ndarray:
    """Return one gene's expression as a dense 1-D array."""
    col = adata[:, gene].X
    if hasattr(col, "toarray"):
        col = col.toarray()
    return np.asarray(col).ravel()


def correlation_table(adata: sc.AnnData, genes: list[str]) -> pd.DataFrame:
    """Spearman correlation of each gene with Adrb2 across Th+ neurons.

    Following the original analysis, each correlation is computed on the cells in which
    both the gene and Adrb2 are detected, so the number of cells varies per gene.
    """
    th = _expression(adata, GENE1)
    adrb2 = _expression(adata, GENE2)
    th_pos = th > 0
    print(f"  {GENE1}+ neurons: {th_pos.sum():,}")

    rows = []
    for gene in genes:
        if gene not in adata.var_names:
            print(f"  {gene}: not in the expression panel, skipped")
            continue
        values = _expression(adata, gene)
        mask = th_pos & (adrb2 > 0) & (values > 0)
        n = int(mask.sum())
        if n < 3:
            rows.append({"gene": gene, "n_cells": n, "spearman_rho": np.nan, "p_value": np.nan})
            continue
        rho, p = spearmanr(adrb2[mask], values[mask])
        rows.append({"gene": gene, "n_cells": n, "spearman_rho": rho, "p_value": p})

    return pd.DataFrame(rows).sort_values("spearman_rho", ascending=False, na_position="last")


def load_de(path: Path, strict: bool = False) -> pd.DataFrame:
    """Load the Wilcoxon result and apply the filters used for the reported table."""
    df = pd.read_csv(path)

    df["logfoldchanges"] = pd.to_numeric(df["logfoldchanges"], errors="coerce")
    df["pvals_adj"] = pd.to_numeric(df["pvals_adj"], errors="coerce")

    n0 = len(df)
    df = df[np.isfinite(df["logfoldchanges"]) & np.isfinite(df["pvals_adj"])]
    df = df[df["pvals_adj"] < 1.0]

    pattern = STRICT_EXCLUDE_REGEX if strict else EXCLUDE_GENE_REGEX
    df = df[~df["names"].str.contains(pattern, regex=True, na=False)]
    df = df.copy()

    df["neglog10_padj"] = -np.log10(df["pvals_adj"].clip(lower=1e-300))
    print(f"  {n0:,} genes tested -> {len(df):,} after filtering"
          f"{' (strict gene filter)' if strict else ''}")
    return df


def select_markers(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Top enriched and depleted genes passing both thresholds."""
    passing = df[df["neglog10_padj"] >= PADJ_THR]
    top_up = passing[passing["logfoldchanges"] > FC_THR].nlargest(N_LABEL, "logfoldchanges")
    top_dn = passing[passing["logfoldchanges"] < -FC_THR].nsmallest(N_LABEL, "logfoldchanges")
    return top_up, top_dn


def _spread(vals: np.ndarray, low: float, high: float, min_dy: float) -> np.ndarray:
    """Push label positions apart so annotation text does not overlap."""
    vals = np.asarray(vals, dtype=float).copy()
    if vals.size == 0:
        return vals
    vals.sort()
    out = vals.copy()
    out[0] = max(out[0], low)
    for i in range(1, len(out)):
        out[i] = max(out[i], out[i - 1] + min_dy)
    overflow = out[-1] - high
    if overflow > 0:
        out -= overflow
        if out[0] < low:
            out = np.linspace(low, high, len(out))
    return out


def volcano(df: pd.DataFrame, labels: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 6.5))

    x = df["logfoldchanges"].to_numpy()
    y = df["neglog10_padj"].to_numpy()

    m_up = (x > FC_THR) & (y >= PADJ_THR)
    m_dn = (x < -FC_THR) & (y >= PADJ_THR)
    m_mid = ~(m_up | m_dn)

    ax.scatter(x[m_mid], y[m_mid], s=8, alpha=0.25, c="lightgray")
    ax.scatter(x[m_up], y[m_up], s=8, alpha=0.6, c="red")
    ax.scatter(x[m_dn], y[m_dn], s=8, alpha=0.6, c="blue")

    x_min, x_max = float(np.nanmin(x)), float(np.nanmax(x))
    y_max = float(np.nanmax(y))
    x_span = x_max - x_min
    y_span = y_max if y_max > 0 else 1.0

    ax.set_xlim(x_min - 0.35 * x_span, x_max + 0.35 * x_span)
    ax.set_ylim(0.0, y_max + 0.06 * y_span)

    min_dy = 0.045 * y_span
    low = max(PADJ_THR, 0.03 * y_span)

    for side, colour, x_text, ha in (
        (labels[labels["logfoldchanges"] < 0], "blue", x_min - 0.18 * x_span, "right"),
        (labels[labels["logfoldchanges"] > 0], "red", x_max + 0.18 * x_span, "left"),
    ):
        side = side.sort_values("neglog10_padj")
        y_positions = _spread(side["neglog10_padj"].to_numpy(), low, y_max, min_dy)
        for (_, row), y_txt in zip(side.iterrows(), y_positions):
            xr, yr = float(row["logfoldchanges"]), float(row["neglog10_padj"])
            ax.scatter([xr], [yr], s=18, c=colour, zorder=3)
            ax.annotate(
                str(row["names"]),
                xy=(xr, yr),
                xytext=(x_text, y_txt),
                ha=ha,
                va="center",
                fontsize=9,
                arrowprops=dict(arrowstyle="-", lw=0.6, color=colour, shrinkA=2, shrinkB=2),
            )

    ax.axvline(-FC_THR, color="k", lw=1, ls="--")
    ax.axvline(0, color="k", lw=2)
    ax.axvline(FC_THR, color="k", lw=1, ls="--")
    ax.axhline(PADJ_THR, color="k", lw=1, ls="--")
    ax.set_xlabel(f"logFC: {GROUP_POS} vs {GROUP_NEG}")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(f"Differential gene expression: {GROUP_POS} vs {GROUP_NEG}")
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indir", default="data/precomputed")
    parser.add_argument("--figdir", default="figures")
    parser.add_argument("--outdir", default="results")
    parser.add_argument(
        "--strict-filter",
        action="store_true",
        help="also exclude BAC clone identifiers, pseudogenes, ribosomal, "
             "mitochondrial and haemoglobin genes from the differential expression table",
    )
    args = parser.parse_args()

    indir = Path(args.indir)
    figdir, outdir = Path(args.figdir), Path(args.outdir)
    figdir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Correlation with Adrb2 -----------------------------------------------------

    print("correlation of reported markers with Adrb2 ...")
    adata = sc.read_h5ad(indir / "expression_panel.h5ad")
    corr = correlation_table(adata, MARKER_GENES)

    corr_out = outdir / "th_adrb2_correlation_table.csv"
    corr.to_csv(corr_out, index=False, float_format="%.4g")
    print(f"  wrote {corr_out}\n")
    print(corr.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print()

    # --- Differential expression ----------------------------------------------------

    print("loading differential expression result ...")
    df = load_de(indir / "th_adrb2_markers.csv.gz", strict=args.strict_filter)

    top_up, top_dn = select_markers(df)
    print(f"  {len(top_up)} enriched, {len(top_dn)} depleted at "
          f"|logFC| > {FC_THR}, -log10 padj >= {PADJ_THR}")

    table = pd.concat(
        [top_up.assign(direction="enriched"), top_dn.assign(direction="depleted")],
        ignore_index=True,
    ).drop_duplicates("names")

    columns = [
        "names", "direction", "logfoldchanges", "scores",
        "pvals", "pvals_adj", "neglog10_padj",
    ]
    columns += [c for c in ("pct_nz_group", "pct_nz_reference") if c in table]
    table = table[columns].rename(
        columns={
            "names": "gene",
            "logfoldchanges": "log2FC",
            "scores": "wilcoxon_score",
            "pvals": "p_value",
            "pvals_adj": "p_adjusted",
            "pct_nz_group": f"pct_expressing_{GROUP_POS}",
            "pct_nz_reference": f"pct_expressing_{GROUP_NEG}",
        }
    )

    out = outdir / "th_adrb2_marker_table.csv"
    table.to_csv(out, index=False, float_format="%.4g")
    print(f"  wrote {out}")

    print()
    print(table.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print()

    volcano(df, pd.concat([top_up, top_dn]).drop_duplicates("names"), figdir / "volcano_th_adrb2.png")
    print("done.")


if __name__ == "__main__":
    main()

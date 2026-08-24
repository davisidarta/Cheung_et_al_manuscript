#!/usr/bin/env python
"""Th and Adrb2 expression across the whole mouse brain.

Queries the Allen Institute whole mouse brain atlas for co-expression of Th and Adrb2,
to establish where in the brain the two are found together. The atlas is downloaded
through the ABC atlas cache; see scripts/00_fetch_source_data.sh.

    python scripts/07_wholebrain_abc_atlas.py --cache data/abc_atlas

Outputs:
    results/wholebrain_th_adrb2_by_region.csv
    figures/wholebrain_th_adrb2_by_region.png
    figures/wholebrain_coexpression.png
"""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# --- Analysis constants -------------------------------------------------------------

GENE1 = "Th"
GENE2 = "Adrb2"

#: Detection threshold on log-normalised expression.
THRESHOLD = 0.5

#: Column in the atlas metadata used to group double-positive cells.
REGION_KEY = "anatomical_division_label"

TARGET_SUM = 1e4


def concatenate_h5ad(directory: Path) -> sc.AnnData:
    """Read and concatenate the atlas expression files one at a time.

    The atlas ships as many per-region files; concatenating sequentially and collecting
    as we go keeps peak memory to roughly one file plus the running total.
    """
    files = sorted(p for p in directory.rglob("*.h5ad"))
    if not files:
        raise FileNotFoundError(f"no .h5ad files under {directory}")

    print(f"  {len(files)} expression files")
    combined = None
    for path in files:
        current = sc.read_h5ad(path)
        print(f"    {path.name}: {current.n_obs:,} cells")
        if combined is None:
            combined = current
        else:
            combined = sc.concat([combined, current], axis=0, join="outer")
            del current
            gc.collect()
    return combined


def to_gene_symbols(adata: sc.AnnData, cache: Path) -> sc.AnnData:
    """Relabel Ensembl gene IDs with symbols using the atlas gene table."""
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

    genes = AbcProjectCache.from_cache_dir(cache).get_metadata_dataframe(
        directory="WMB-10X", file_name="gene"
    )
    genes = genes.set_index("gene_identifier")

    mapping = genes["gene_symbol"].to_dict()
    adata.var["gene_symbol"] = [mapping.get(g, g) for g in adata.var_names]
    adata.var_names = pd.Index(adata.var["gene_symbol"].astype(str))
    adata.var_names_make_unique()
    return adata


def expression(adata: sc.AnnData, gene: str) -> np.ndarray:
    col = adata[:, gene].X
    if hasattr(col, "toarray"):
        col = col.toarray()
    return np.asarray(col).ravel()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", default="data/abc_atlas")
    parser.add_argument("--figdir", default="figures")
    parser.add_argument("--outdir", default="results")
    args = parser.parse_args()

    cache = Path(args.cache)
    figdir, outdir = Path(args.figdir), Path(args.outdir)
    figdir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    sc.settings.set_figure_params(dpi=150, facecolor="white")
    sc.settings.n_jobs = min(30, os.cpu_count() or 1)

    print("reading atlas ...")
    adata = concatenate_h5ad(cache)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    adata = to_gene_symbols(adata, cache)

    print("attaching cell metadata ...")
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

    meta = AbcProjectCache.from_cache_dir(cache).get_metadata_dataframe(
        directory="WMB-10X", file_name="cell_metadata_with_cluster_annotation"
    )
    meta = meta.set_index("cell_label")
    adata = adata[adata.obs_names.isin(meta.index)].copy()

    columns = [c for c in ("neurotransmitter", "class", "subclass", REGION_KEY) if c in meta]
    for column in columns:
        adata.obs[column] = meta.loc[adata.obs_names, column].values
    print(f"  {adata.n_obs:,} cells retained with annotation")

    print("normalising ...")
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)

    missing = [g for g in (GENE1, GENE2) if g not in adata.var_names]
    if missing:
        raise KeyError(f"not found in the atlas: {', '.join(missing)}")

    x1 = expression(adata, GENE1)
    x2 = expression(adata, GENE2)
    pos1, pos2 = x1 > THRESHOLD, x2 > THRESHOLD
    both = pos1 & pos2

    n = adata.n_obs
    print(
        f"\n  {GENE1}+          {pos1.sum():>9,}  ({100 * pos1.mean():.2f}% of all cells)\n"
        f"  {GENE2}+       {pos2.sum():>9,}  ({100 * pos2.mean():.2f}% of all cells)\n"
        f"  {GENE1}+/{GENE2}+ {both.sum():>9,}  ({100 * both.mean():.2f}% of all cells, "
        f"{100 * both.sum() / max(pos1.sum(), 1):.2f}% of {GENE1}+ cells)"
    )

    # --- Double-positive cells by anatomical region ---------------------------------

    if REGION_KEY in adata.obs:
        counts = adata.obs.loc[both, REGION_KEY].value_counts()
        counts.to_csv(outdir / "wholebrain_th_adrb2_by_region.csv", header=["n_cells"])
        print(f"\n  wrote {outdir / 'wholebrain_th_adrb2_by_region.csv'}")

        fig, ax = plt.subplots(figsize=(8, 5))
        counts.plot(kind="bar", color="skyblue", edgecolor="black", ax=ax)
        ax.set_title(f"{GENE1}+/{GENE2}+ cells by anatomical region")
        ax.set_ylabel("Number of cells")
        ax.set_xlabel("Anatomical region")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(figdir / "wholebrain_th_adrb2_by_region.png", dpi=200, facecolor="white")
        plt.close(fig)
        print(f"  wrote {figdir / 'wholebrain_th_adrb2_by_region.png'}")

    # --- Co-expression summary ------------------------------------------------------

    fig, ax = plt.subplots(figsize=(5, 5))
    fractions = [
        100 * pos1.mean(),
        100 * pos2.mean(),
        100 * both.mean(),
    ]
    ax.bar(
        [f"{GENE1}+", f"{GENE2}+", f"{GENE1}+/{GENE2}+"],
        fractions,
        color=["#c44e52", "#4c72b0", "#8172b2"],
        edgecolor="black",
    )
    for i, v in enumerate(fractions):
        ax.text(i, v, f"{v:.2f}%", ha="center", va="bottom", fontsize=10)
    ax.set_ylabel(f"% of all brain cells (n = {n:,})")
    ax.set_title(f"{GENE1} and {GENE2} in the whole mouse brain")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(figdir / "wholebrain_coexpression.png", dpi=200, facecolor="white")
    plt.close(fig)
    print(f"  wrote {figdir / 'wholebrain_coexpression.png'}")

    print("\ndone.")


if __name__ == "__main__":
    main()

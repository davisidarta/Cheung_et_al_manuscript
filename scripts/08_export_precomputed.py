#!/usr/bin/env python
"""Build the small, version-controllable inputs used by the figure scripts.

The processed neuron object is far too large to distribute through git, and the
differential expression test behind the Th+/Adrb2+ marker table runs over ~128k genes.
This script derives everything the figure scripts need from that object once, and writes
it to ``data/precomputed/`` so that scripts 05 and 06 run in seconds from a clean clone.

Run this only if you are regenerating the precomputed inputs; ordinary use of the
repository does not require it.

    python scripts/08_export_precomputed.py --input SympAtlas_Neurons.h5ad

Outputs (all under data/precomputed/):
    cell_metadata.csv.gz          per-cell annotations and cluster labels
    embeddings.npz                topoPaCMAP / topoMAP / UMAP / PCA coordinates
    expression_panel.h5ad         log1p expression for the genes the figures plot
    th_adrb2_markers.csv.gz       Wilcoxon table, Th+/Adrb2+ vs Th+/Adrb2-
    cluster_markers_lowres.csv.gz Wilcoxon table, per low-resolution cluster
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

# --- Analysis constants, kept identical to the original analysis -------------------

GENE1 = "Th"
GENE2 = "Adrb2"
THR1 = 0.0
THR2 = 0.0

CLUSTER_KEY = "Leiden on topological graph (low-res)"

#: Columns carried into the shipped metadata table.
OBS_COLUMNS = [
    "Study",
    "Ganglion",
    "Species",
    "internal_batch",
    "batch",
    "Leiden on topological graph (low-res)",
    "Leiden on topological graph (mid-res)",
    "Leiden on topological graph",
    "topoleiden_coarse",
    "phase",
    "S_score",
    "G2M_score",
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
    "doublet_score",
]

#: Embeddings carried into the shipped .npz.
OBSM_KEYS = [
    "X_topoPaCMAP",
    "X_topoPaCMAP3D",
    "X_topoMAP",
    "X_topoMAP3D",
    "X_umap",
    "X_UMAP_on_PCA",
    "X_pca",
]

#: Genes the figure scripts plot directly.
FIGURE_GENES = [
    # adrenergic receptors and the sympathetic / pan-neuronal markers
    "Th", "Adrb1", "Adrb2", "Adrb3", "Dbh", "Slc18a2", "Slc6a2", "Chat", "Slc5a7",
    "Tubb3", "Rbfox3", "Snap25", "Elavl4", "Phox2b", "Ret", "Grm7",
    # Th+/Adrb2+ marker candidates reported in the paper
    "Cdh12", "Auts2l1", "Cntnap5c", "Akap2", "Pcdhac2", "Samd5", "Chsy3", "Supt3h",
    # non-neuronal markers, for the cluster-identity panel
    "Mpz", "Apoe", "Cd74", "Dcn", "Sst", "Npy", "Calca", "Vip",
]

#: Top N markers per cluster added to the expression panel and shipped marker table.
N_PANEL_MARKERS = 15
N_TABLE_MARKERS = 200


def build_de_table(adata: sc.AnnData) -> pd.DataFrame:
    """Wilcoxon test of Th+/Adrb2+ against Th+/Adrb2- neurons.

    Mirrors the original notebook: cells are gated on raw (log1p) expression, the
    comparison is restricted to Th+ neurons, and Adrb2 detection splits them.
    """
    expr = adata.raw.to_adata() if adata.raw is not None else adata
    x1 = _dense_column(expr, GENE1)
    x2 = _dense_column(expr, GENE2)

    th_pos = x1 > THR1
    sub = adata[th_pos].copy()

    group_pos = f"{GENE1}+/{GENE2}+"
    group_neg = f"{GENE1}+/{GENE2}-"
    sub.obs["comparison"] = pd.Categorical(
        np.where(x2[th_pos] > THR2, group_pos, group_neg),
        categories=[group_pos, group_neg],
    )
    n_pos = int((sub.obs["comparison"] == group_pos).sum())
    print(
        f"  Th+ neurons: {th_pos.sum():,}  "
        f"({group_pos}: {n_pos:,} | {group_neg}: {th_pos.sum() - n_pos:,})"
    )

    key = f"{GENE1}_pos_{GENE2}_de"
    sc.tl.rank_genes_groups(
        sub,
        groupby="comparison",
        groups=[group_pos],
        reference=group_neg,
        method="wilcoxon",
        use_raw=True,
        pts=True,
        key_added=key,
    )
    df = sc.get.rank_genes_groups_df(sub, key=key, group=group_pos).copy()
    df["group"] = group_pos
    df["reference"] = group_neg
    return df


def _dense_column(adata: sc.AnnData, gene: str) -> np.ndarray:
    """Return one gene's expression as a dense 1-D array."""
    col = adata[:, gene].X
    if hasattr(col, "toarray"):
        col = col.toarray()
    return np.asarray(col).ravel()


def build_cluster_markers(adata: sc.AnnData) -> pd.DataFrame:
    """Per-cluster Wilcoxon markers, taken from the stored result if present."""
    stored = f"{CLUSTER_KEY}_wilcoxon"
    if stored in adata.uns:
        print(f"  using stored result: uns['{stored}']")
        key = stored
    else:
        print("  no stored result; recomputing")
        key = stored
        sc.tl.rank_genes_groups(
            adata, groupby=CLUSTER_KEY, method="wilcoxon", use_raw=True, key_added=key
        )

    frames = []
    for group in adata.uns[key]["names"].dtype.names:
        part = sc.get.rank_genes_groups_df(adata, key=key, group=group)
        part = part.head(N_TABLE_MARKERS).copy()
        part.insert(0, "cluster", group)
        frames.append(part)
    return pd.concat(frames, ignore_index=True)


def build_expression_panel(adata: sc.AnnData, markers: pd.DataFrame) -> sc.AnnData:
    """Subset raw log1p expression to the genes the figures actually plot.

    Built as a fresh AnnData rather than by slicing: the raw view carries the parent's
    uns and obsp, which hold the full-transcriptome test results and the neighbour
    graphs, and neither belongs in a distributed panel.
    """
    raw = adata.raw.to_adata()
    available = set(raw.var_names)

    panel = [g for g in FIGURE_GENES if g in available]
    missing = sorted(set(FIGURE_GENES) - available)
    if missing:
        print(f"  note: not present in the data, skipped -> {', '.join(missing)}")

    top = (
        markers.groupby("cluster", observed=True)
        .head(N_PANEL_MARKERS)["names"]
        .tolist()
    )
    panel += [g for g in top if g in available]

    panel = sorted(set(panel))
    print(f"  panel: {len(panel)} genes")

    idx = raw.var_names.get_indexer(panel)
    X = raw.X[:, idx]
    X = sp.csr_matrix(X, dtype="float32")

    out = sc.AnnData(
        X=X,
        obs=adata.obs[[c for c in OBS_COLUMNS if c in adata.obs]].copy(),
        var=pd.DataFrame(index=pd.Index(panel, name=None)),
    )
    for k in OBSM_KEYS:
        if k in adata.obsm:
            out.obsm[k] = np.asarray(adata.obsm[k], dtype="float32")
    print(f"  {X.nnz:,} non-zero values")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="SympAtlas_1-0_Neurons_Samson_Ji.h5ad",
        help="processed neuron object to derive the precomputed inputs from",
    )
    parser.add_argument("--outdir", default="data/precomputed")
    parser.add_argument(
        "--steps",
        default="all",
        help="comma-separated subset of: metadata,embeddings,markers,de,panel",
    )
    args = parser.parse_args()

    steps = {"metadata", "embeddings", "markers", "de", "panel"}
    if args.steps != "all":
        requested = {s.strip() for s in args.steps.split(",")}
        unknown = requested - steps
        if unknown:
            parser.error(f"unknown step(s): {', '.join(sorted(unknown))}")
        steps = requested
    # The panel needs the marker table to choose its genes.
    need_markers = "markers" in steps or "panel" in steps

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"reading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} integrated features")
    print(f"  raw: {adata.raw.shape[1]:,} genes")

    markers = None

    if "metadata" in steps:
        print("\n[metadata]")
        meta = adata.obs[[c for c in OBS_COLUMNS if c in adata.obs]].copy()
        meta.index.name = "cell_id"
        meta.to_csv(outdir / "cell_metadata.csv.gz", compression="gzip")
        print(f"  {meta.shape[0]:,} x {meta.shape[1]} -> cell_metadata.csv.gz")

    if "embeddings" in steps:
        print("\n[embeddings]")
        emb = {
            k: np.asarray(adata.obsm[k], dtype="float32")
            for k in OBSM_KEYS
            if k in adata.obsm
        }
        np.savez_compressed(outdir / "embeddings.npz", **emb)
        print(f"  {', '.join(emb)} -> embeddings.npz")

    if need_markers:
        print("\n[per-cluster markers]")
        markers = build_cluster_markers(adata)
        if "markers" in steps:
            markers.to_csv(
                outdir / "cluster_markers_lowres.csv.gz", index=False, compression="gzip"
            )
            print(f"  {markers.shape[0]:,} rows -> cluster_markers_lowres.csv.gz")

    if "de" in steps:
        print("\n[Th+/Adrb2+ differential expression]")
        de = build_de_table(adata)
        de.to_csv(outdir / "th_adrb2_markers.csv.gz", index=False, compression="gzip")
        print(f"  {de.shape[0]:,} rows -> th_adrb2_markers.csv.gz")

    if "panel" in steps:
        print("\n[expression panel]")
        panel = build_expression_panel(adata, markers)
        panel.write_h5ad(outdir / "expression_panel.h5ad", compression="gzip")
        print(f"  {panel.n_obs:,} x {panel.n_vars} -> expression_panel.h5ad")

    provenance = {
        "source_object": Path(args.input).name,
        "n_cells": int(adata.n_obs),
        "n_raw_genes": int(adata.raw.shape[1]),
        "cluster_key": CLUSTER_KEY,
        "de_gate": {"gene1": GENE1, "thr1": THR1, "gene2": GENE2, "thr2": THR2},
        "scanpy": sc.__version__,
    }
    (outdir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    print("\ndone.")
    for f in sorted(outdir.iterdir()):
        print(f"  {f.name:<32} {f.stat().st_size / 1e6:8.1f} MB")


if __name__ == "__main__":
    main()

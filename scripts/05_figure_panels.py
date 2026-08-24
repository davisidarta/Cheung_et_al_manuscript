#!/usr/bin/env python
"""Atlas panels: embeddings and adrenergic receptor expression across neuron populations.

Reads the precomputed inputs in ``data/precomputed/`` and writes the panels to
``figures/``. Nothing here is expensive; a clean clone runs it in seconds.

    python scripts/05_figure_panels.py

Panels produced:
    embedding_clusters.png      topoPaCMAP coloured by cluster
    embedding_study.png         topoPaCMAP coloured by dataset of origin
    embedding_ganglion.png      topoPaCMAP coloured by ganglion
    embedding_receptors.png     topoPaCMAP coloured by Th, Adrb2, Adrb3, Adrb1
    dotplot_adrb2_ganglion.png  Adrb2 across ganglia
    dotplot_receptors.png       Adrb1/2/3 across clusters
    violin_th_grm7.png          Th and Grm7 across ganglia
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

# --- Constants matching the original analysis --------------------------------------

CLUSTER_KEY = "Leiden on topological graph (low-res)"
BASIS = "topoPaCMAP"

#: Anatomical ordering, rostral to caudal.
GANGLION_ORDER = ["SCG", "Stellate", "ThoracicPool", "Lumbar", "Coeliac", "Pelvic"]

RECEPTOR_GENES = ["Adrb2", "Adrb3", "Adrb1"]
EMBEDDING_GENES = ["Th", "Adrb2", "Adrb3", "Adrb1"]


def load(indir: Path) -> sc.AnnData:
    """Load the expression panel and attach its embeddings."""
    adata = sc.read_h5ad(indir / "expression_panel.h5ad")

    if adata.obs[CLUSTER_KEY].dtype.name != "category":
        adata.obs[CLUSTER_KEY] = adata.obs[CLUSTER_KEY].astype("category")

    present = [g for g in GANGLION_ORDER if g in set(adata.obs["Ganglion"])]
    adata.obs["Ganglion"] = pd.Categorical(
        adata.obs["Ganglion"], categories=present, ordered=True
    )
    return adata


def save(path: Path) -> None:
    plt.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close("all")
    print(f"  wrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indir", default="data/precomputed")
    parser.add_argument("--outdir", default="figures")
    args = parser.parse_args()

    indir, outdir = Path(args.indir), Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Set through matplotlib rather than sc.settings.set_figure_params, which recurses
    # into itself on some scanpy builds.
    plt.rcParams.update({"figure.dpi": 150, "figure.facecolor": "white", "savefig.facecolor": "white"})

    adata = load(indir)
    print(f"{adata.n_obs:,} neurons x {adata.n_vars} genes")

    # --- Embeddings -----------------------------------------------------------------

    sc.pl.embedding(
        adata,
        basis=BASIS,
        color=[CLUSTER_KEY],
        legend_loc="on data",
        legend_fontsize=12,
        legend_fontoutline=2,
        title="Sympathetic neurons",
        frameon=False,
        show=False,
    )
    save(outdir / "embedding_clusters.png")

    sc.pl.embedding(
        adata,
        basis=BASIS,
        color=["Study"],
        title="Dataset of origin",
        frameon=False,
        show=False,
    )
    save(outdir / "embedding_study.png")

    sc.pl.embedding(
        adata,
        basis=BASIS,
        color=["Ganglion"],
        title="Ganglion",
        frameon=False,
        show=False,
    )
    save(outdir / "embedding_ganglion.png")

    sc.pl.embedding(
        adata,
        basis=BASIS,
        color=EMBEDDING_GENES,
        cmap="Reds",
        vmin=0,
        vmax=1.5,
        ncols=2,
        frameon=False,
        use_raw=False,
        show=False,
    )
    save(outdir / "embedding_receptors.png")

    # --- Receptor expression across populations -------------------------------------

    sc.pl.dotplot(
        adata,
        var_names=["Adrb2"],
        groupby="Ganglion",
        standard_scale="var",
        cmap="Reds",
        use_raw=False,
        show=False,
    )
    save(outdir / "dotplot_adrb2_ganglion.png")

    sc.pl.dotplot(
        adata,
        var_names=RECEPTOR_GENES,
        groupby=CLUSTER_KEY,
        standard_scale="var",
        swap_axes=True,
        vmax=0.3,
        mean_only_expressed=True,
        smallest_dot=0.1,
        use_raw=False,
        show=False,
    )
    save(outdir / "dotplot_receptors.png")

    sc.pl.stacked_violin(
        adata,
        var_names=["Th", "Grm7"],
        groupby="Ganglion",
        standard_scale="var",
        swap_axes=False,
        use_raw=False,
        density_norm="area",
        show=False,
    )
    save(outdir / "violin_th_grm7.png")

    print("done.")


if __name__ == "__main__":
    main()

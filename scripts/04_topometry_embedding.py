#!/usr/bin/env python
"""Dimensionality reduction and clustering of the integrated neuron atlas.

By default this attaches the embeddings distributed in ``data/precomputed/`` to an
integrated object, which is what the figure scripts consume. Passing ``--refit`` runs the
TopOMetry fit instead; that takes hours on 35k cells and needs a large machine, and
because the projection and the Leiden partition are both stochastic it produces an
equivalent but not coordinate-identical embedding.

    python scripts/04_topometry_embedding.py --input neurons.h5ad --output neurons_embedded.h5ad
    python scripts/04_topometry_embedding.py --input neurons.h5ad --output neurons_embedded.h5ad --refit

Parameters below are those used for the published atlas.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc

# --- TopOMetry configuration --------------------------------------------------------

#: Neighbours for the base kernel.
BASE_KNN = 20
#: Neighbours for the topological graph.
GRAPH_KNN = 15
#: Eigencomponents retained from the diffusion operator.
N_EIGS = 200
#: Kernel used to build the affinity matrix.
KERNEL = "bw_adaptive"

#: Leiden resolutions, and the obs key each one is stored under.
LEIDEN_RESOLUTIONS = {
    "Leiden on topological graph (low-res)": 0.2,
    "Leiden on topological graph (mid-res)": 0.5,
    "Leiden on topological graph": 3.0,
}

#: Principal components computed on the integrated matrix, for comparison plots.
N_PCS = 50


def attach_precomputed(adata: sc.AnnData, path: Path) -> sc.AnnData:
    """Attach distributed embeddings, matching them to the object by cell order."""
    npz = np.load(path)
    print(f"  loading {path}")
    for key in npz.files:
        coords = npz[key]
        if coords.shape[0] != adata.n_obs:
            raise ValueError(
                f"{key}: {coords.shape[0]} rows but the object has {adata.n_obs} cells. "
                "The precomputed embeddings correspond to the published neuron set; "
                "re-run with --refit if you have rebuilt the object from scratch."
            )
        adata.obsm[key] = coords
        print(f"    {key:<20} {coords.shape}")
    return adata


def refit(adata: sc.AnnData) -> sc.AnnData:
    """Run the TopOMetry fit from the integrated matrix."""
    import topo as tp

    print("  running PCA")
    sc.pp.pca(adata, n_comps=N_PCS, zero_center=True, use_highly_variable=False)

    print(f"  fitting TopOGraph (base_knn={BASE_KNN}, graph_knn={GRAPH_KNN}, n_eigs={N_EIGS})")
    tg = tp.TopOGraph(
        base_knn=BASE_KNN,
        graph_knn=GRAPH_KNN,
        n_eigs=N_EIGS,
        base_metric="cosine",
        graph_metric="euclidean",
        bases="diffusion",
        graphs="diff",
        projections=["PaCMAP", "MAP"],
        n_jobs=-1,
        verbosity=1,
    )
    tg.run_layouts(adata.X, bases=["diffusion"], graphs=["diff"], projections=["PaCMAP", "MAP"])
    tg.to_adata(adata)

    print("  clustering")
    for key, resolution in LEIDEN_RESOLUTIONS.items():
        sc.tl.leiden(adata, resolution=resolution, key_added=key, flavor="igraph", n_iterations=2)
        n = adata.obs[key].nunique()
        print(f"    {key:<45} resolution {resolution} -> {n} clusters")

    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="integrated neuron object (.h5ad)")
    parser.add_argument("--output", required=True, help="where to write the result")
    parser.add_argument(
        "--precomputed",
        default="data/precomputed/embeddings.npz",
        help="distributed embeddings to attach (default mode)",
    )
    parser.add_argument(
        "--refit",
        action="store_true",
        help="run the TopOMetry fit instead of attaching the distributed embeddings",
    )
    args = parser.parse_args()

    print(f"reading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} features")

    if args.refit:
        print("\nre-fitting embeddings")
        adata = refit(adata)
    else:
        print("\nattaching precomputed embeddings")
        adata = attach_precomputed(adata, Path(args.precomputed))

    adata.write_h5ad(args.output, compression="gzip")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()

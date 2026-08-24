# Single-cell analyses for Cheung _et al._, 2026

Code for the single-cell RNA sequencing analyses presented in "β₂AR Agonists Sustain
Thermogenesis and Leanness via Sympathofacilitation", by Samson W. Cheung and colleagues.

Two datasets are analysed:

- **Sympathetic ganglia** — a meta-analysis integrating ten public single-cell and
  single-nucleus RNA-seq datasets of mouse sympathetic ganglia, using the Seurat v3
  canonical correlation analysis anchor workflow followed by TopOMetry dimensionality
  reduction. 35,333 of 239,202 cells were classified as neurons by expression of the
  pan-neuronal markers `Tubb3` and `Rbfox3`.
- **Whole mouse brain** — the Allen Institute whole mouse brain atlas, queried
  programmatically for `Th` and `Adrb2` expression.

## Layout

```
scripts/
  00_fetch_source_data.sh        download the source datasets
  01_qc_and_normalise.R          per-dataset QC and log-normalisation
  02_integrate_seurat_cca.R      anchor-based CCA integration
  03_select_neurons_export.R     pan-neuronal marker gating, export to h5ad
  04_topometry_embedding.py      dimensionality reduction and clustering
  05_figure_panels.py            embeddings and receptor expression panels
  06_th_adrb2_markers.py         Th+/Adrb2+ marker table and volcano
  07_wholebrain_abc_atlas.py     whole-brain Th/Adrb2 query
  08_export_precomputed.py       regenerate the inputs under data/precomputed/
data/
  sample_manifest.csv            the ten source datasets and their accessions
  precomputed/                   small inputs consumed by the figure scripts
notebooks/                       the exploratory notebooks the analyses grew from
```

## Running it

```bash
conda env create -f environment.yml
conda activate cheung-symp
```

The figure scripts read the precomputed inputs in `data/precomputed/` and run in seconds
from a clean clone:

```bash
python scripts/05_figure_panels.py      # atlas panels -> figures/
python scripts/06_th_adrb2_markers.py   # marker table -> results/
```

Rebuilding the atlas from the source data is a much longer path. `00` downloads a few
hundred gigabytes, `02` needs a large-memory machine and several hours to find anchors
across ten datasets, and `04` fits the TopOMetry embedding:

```bash
bash scripts/00_fetch_source_data.sh ganglia
Rscript scripts/01_qc_and_normalise.R
Rscript scripts/02_integrate_seurat_cca.R
Rscript scripts/03_select_neurons_export.R
python scripts/04_topometry_embedding.py --input data/neurons.h5ad --output data/neurons_embedded.h5ad
```

The whole-brain analysis is independent of the ganglia pipeline:

```bash
bash scripts/00_fetch_source_data.sh brain
python scripts/07_wholebrain_abc_atlas.py
```

## Data

Accessions for the ten source datasets are in `data/sample_manifest.csv`. Two need
manual steps: the Lee _et al._ data is distributed through the SPARC portal, which asks
you to accept its licence before download, and the processed loom files for Zeisel
_et al._ are no longer served from the URLs given in that paper, so that dataset has to
be started from raw sequence data in SRA. `00_fetch_source_data.sh` leaves a note in the
relevant directory in both cases.

`data/precomputed/` holds the small inputs the figure scripts consume, derived from the
processed neuron object by `08_export_precomputed.py`:

| file | contents |
| --- | --- |
| `cell_metadata.csv.gz` | per-cell study, ganglion, cluster labels and QC statistics |
| `embeddings.npz` | topoPaCMAP, topoMAP, UMAP and PCA coordinates |
| `expression_panel.h5ad` | log-normalised expression for the genes the figures plot |
| `th_adrb2_markers.csv.gz` | Wilcoxon result, Th+/Adrb2+ against Th+/Adrb2− |
| `cluster_markers_lowres.csv.gz` | Wilcoxon markers per cluster |

The integration and the TopOMetry fit are expensive and depend on package versions, so
the embeddings and the full-transcriptome test results are distributed precomputed rather
than recomputed on every run; `04_topometry_embedding.py --refit` re-runs the fit, which
produces an equivalent but not coordinate-identical embedding, as any stochastic
projection does. The complete processed object is too large to distribute through git and
is deposited separately.

## Citation

```bibtex
@article{cheung2025b2ar,
  author  = {Cheung, Samson W. and Sidarta-Oliveira, David and Sarker, Gitalee and
             Qu, Ji and Yao, Lu and Zhu, Yitao and Lundh, Sofia and Griebel, Alina and
             Morgan, Donald A. and Martinez-Sanchez, Noelia and Liu, Kun and
             Goulding, Joelle and Wilcox, Sian and Hill, Stephen J. and
             Vyazovskiy, Vladyslav and Ziegler, Karin A. and Engelhardt, Stefan and
             Paterson, David J. and Li, Dan and Rahmouni, Kamal and Domingos, Ana I.},
  title   = {{\ss}2AR Agonists Sustain Thermogenesis and Leanness via Sympathofacilitation},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.06.13.659468},
  url     = {https://www.biorxiv.org/content/10.1101/2025.06.13.659468}
}
```

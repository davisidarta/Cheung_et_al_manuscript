#!/usr/bin/env bash
#
# Download the source data for both analyses.
#
#   bash scripts/00_fetch_source_data.sh ganglia    # the ten sympathetic ganglia datasets
#   bash scripts/00_fetch_source_data.sh brain      # the whole mouse brain atlas
#   bash scripts/00_fetch_source_data.sh all
#
# Accessions are listed in data/sample_manifest.csv. Downloads land in data/raw/<study>/.
# Expect this to take a long time and a few hundred gigabytes.

set -euo pipefail

RAW_DIR="${RAW_DIR:-data/raw}"
BRAIN_DIR="${BRAIN_DIR:-data/abc_atlas}"
TARGET="${1:-all}"

geo_supplementary() {
  # geo_supplementary <GSE accession> <destination directory>
  local acc="$1" dest="$2"
  local stub="${acc:0:$(( ${#acc} - 3 ))}nnn"
  mkdir -p "$dest"
  echo "  -> $acc"
  curl -fsSL --retry 3 \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${acc}/suppl/${acc}_RAW.tar" \
    -o "${dest}/${acc}_RAW.tar" \
    || echo "     no bundled RAW.tar; browse ${dest} manually from the accession page"
  if [ -f "${dest}/${acc}_RAW.tar" ]; then
    tar -xf "${dest}/${acc}_RAW.tar" -C "$dest" && rm "${dest}/${acc}_RAW.tar"
  fi
}

fetch_ganglia() {
  echo "== sympathetic ganglia datasets =="
  mkdir -p "$RAW_DIR"

  geo_supplementary GSE78845  "$RAW_DIR/Furlan_NatNeuro16"
  geo_supplementary GSE175421 "$RAW_DIR/MappsThomsen_CellRep22"
  geo_supplementary GSE232286 "$RAW_DIR/Ge_JVisExp22"
  geo_supplementary GSE231924 "$RAW_DIR/Sharma_eLife23"
  geo_supplementary GSE231767 "$RAW_DIR/Ziegler_Science23"
  geo_supplementary GSE233163 "$RAW_DIR/Sarker_bioRxiv24"
  geo_supplementary GSE232789 "$RAW_DIR/Sivori_eLife24"
  geo_supplementary GSE278457 "$RAW_DIR/Wang_Nature24"

  # Zeisel et al. 2018 is deposited as raw sequence data only. The processed loom files
  # that the Linnarsson lab previously served alongside it are no longer reachable, so
  # this study has to be started from SRA and quantified locally.
  echo "  -> SRP135960 (Zeisel_Cell18)"
  mkdir -p "$RAW_DIR/Zeisel_Cell18"
  cat > "$RAW_DIR/Zeisel_Cell18/README.txt" <<'EOF'
Zeisel et al. 2018, Cell 174:999.

Raw sequence data: SRA SRP135960 / BioProject PRJNA438862
  https://www.ncbi.nlm.nih.gov/sra/SRP135960

The processed per-taxon loom files (l6_r4_sympathetic_noradrenergic_neurons.loom and
related) were distributed from the linnarsson-lab-loom storage bucket, which no longer
resolves. Obtain the sympathetic subset either by quantifying SRP135960 with the
authors' pipeline, or by requesting the processed file from the original authors.
EOF
  echo "     see $RAW_DIR/Zeisel_Cell18/README.txt"

  # Lee, Thaker & Zeltser 2022 is distributed through SPARC rather than GEO. The portal
  # requires accepting the dataset licence, so fetch it interactively.
  echo "  -> SPARC dataset 263 (Lee_Protocols22)"
  mkdir -p "$RAW_DIR/Lee_Protocols22"
  cat > "$RAW_DIR/Lee_Protocols22/README.txt" <<'EOF'
Lee, Thaker & Zeltser 2022. SPARC dataset 263, DOI 10.26275/bvpu-cuz7
  https://sparc.science/datasets/263

Download the per-plate count matrices from the portal (CC-BY; the portal asks you to
accept the licence before download, which is why this is not automated).
EOF
  echo "     see $RAW_DIR/Lee_Protocols22/README.txt"

  echo "done. Source layout under $RAW_DIR/"
}

fetch_brain() {
  echo "== whole mouse brain atlas =="
  python - "$BRAIN_DIR" <<'EOF'
import sys
from pathlib import Path

from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

dest = Path(sys.argv[1])
dest.mkdir(parents=True, exist_ok=True)

cache = AbcProjectCache.from_cache_dir(dest)
print(f"  manifest: {cache.current_manifest}")

cache.get_metadata_dataframe(
    directory="WMB-10X", file_name="cell_metadata_with_cluster_annotation"
)
cache.get_metadata_dataframe(directory="WMB-10X", file_name="gene")
cache.get_directory_data("WMB-10Xv3")
print(f"  cached under {dest}")
EOF
}

case "$TARGET" in
  ganglia) fetch_ganglia ;;
  brain)   fetch_brain ;;
  all)     fetch_ganglia; fetch_brain ;;
  *)       echo "usage: $0 [ganglia|brain|all]" >&2; exit 1 ;;
esac

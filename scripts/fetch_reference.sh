#!/usr/bin/env bash
# Download the three reference inputs: host genome, 10x whitelists, viral
# FASTA+GTF. Each download is verified by a count, not by exit status --
# a truncated or redirected file downloads "successfully" all the time.
set -euo pipefail
cd "$(dirname "$0")/.."
source config/config.sh
mkdir -p "$REF_DIR/GRCh38" "$WHITELIST_DIR" "$REF_DIR/hsv1"

get() { # get <url> <dest>
  [ -s "$2" ] && { echo "  have $(basename "$2")"; return; }
  echo "  -> $(basename "$2")"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x"$ARIA_CONN" -s"$ARIA_CONN" -c -d "$(dirname "$2")" -o "$(basename "$2")" "$1"
  else
    curl -fL --retry 3 -o "$2" "$1"
  fi
}

echo "=== 1. host reference (CellRanger GRCh38 2020-A) ==="
# 2020-A specifically: the PDAC atlas was built on it, so gene symbols join
# without an alias map. A different build reintroduces symbol drift.
TAR="$REF_DIR/refdata-gex-GRCh38-2020-A.tar.gz"
get "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz" "$TAR"
if [ ! -s "$REF_DIR/GRCh38/fasta/genome.fa" ]; then
  echo "  extracting (~11 GB, a few minutes)..."
  tar -xzf "$TAR" -C "$REF_DIR" --strip-components=1 \
      refdata-gex-GRCh38-2020-A/fasta/genome.fa \
      refdata-gex-GRCh38-2020-A/genes/genes.gtf
  mkdir -p "$REF_DIR/GRCh38/fasta" "$REF_DIR/GRCh38/genes"
  mv "$REF_DIR/fasta/genome.fa" "$REF_DIR/GRCh38/fasta/" 2>/dev/null || true
  mv "$REF_DIR/genes/genes.gtf" "$REF_DIR/GRCh38/genes/" 2>/dev/null || true
fi
n=$(grep -c '^>' "$REF_DIR/GRCh38/fasta/genome.fa")
echo "  contigs: $n"
[ "$n" -ge 24 ] || { echo "  FATAL: expected >=24 contigs, got $n" >&2; exit 1; }

echo
echo "=== 2. 10x barcode whitelists ==="
# These ship inside CellRanger. If you have CellRanger installed, copy from
# <cellranger>/lib/python/cellranger/barcodes/ instead of downloading.
declare -A WL_LINES=( ["737K-august-2016.txt"]=737280 ["3M-february-2018.txt"]=6794880 )
BASE="https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes"
for wl in "${!WL_LINES[@]}"; do
  dest="$WHITELIST_DIR/$wl"
  if [ ! -s "$dest" ]; then
    get "$BASE/$wl" "$dest" || get "$BASE/${wl}.gz" "${dest}.gz" || true
    [ -s "${dest}.gz" ] && gunzip -f "${dest}.gz"
  fi
  if [ ! -s "$dest" ]; then
    echo "  COULD NOT FETCH $wl." >&2
    echo "  Copy it from a CellRanger install:" >&2
    echo "    <cellranger>/lib/python/cellranger/barcodes/$wl -> $WHITELIST_DIR/" >&2
    continue
  fi
  got=$(wc -l < "$dest"); want="${WL_LINES[$wl]}"
  echo "  $wl: $got lines (expected $want)"
  [ "$got" -eq "$want" ] || { echo "  FATAL: $wl is wrong size -- truncated or not the real list" >&2; exit 1; }
done

echo
echo "=== 3. viral FASTA + GTF (${GEO_SERIES}) ==="
PRE="${GEO_SERIES:0:6}nnn"
GEOBASE="https://ftp.ncbi.nlm.nih.gov/geo/series/${PRE}/${GEO_SERIES}/suppl"
for f in "${GEO_SERIES}_HSV1GFPtrunc_corr.fa.gz" "${GEO_SERIES}_HSV1GFPtrunc_corr.gtf.gz"; do
  get "$GEOBASE/$f" "$REF_DIR/hsv1/$f"
  [ -s "$REF_DIR/hsv1/$f" ] && gunzip -kf "$REF_DIR/hsv1/$f"
done
VF="$REF_DIR/hsv1/${GEO_SERIES}_HSV1GFPtrunc_corr.fa"
VG="$REF_DIR/hsv1/${GEO_SERIES}_HSV1GFPtrunc_corr.gtf"
if [ -s "$VF" ] && [ -s "$VG" ]; then
  echo "  viral contigs : $(grep '^>' "$VF" | sed 's/^>//' | cut -d' ' -f1 | tr '\n' ' ')"
  echo "  viral GTF rows: $(grep -vc '^#' "$VG")"
  echo "  feature types : $(awk -F'\t' '!/^#/{print $3}' "$VG" | sort | uniq -c | tr '\n' ' ')"
  echo
  echo "  Compare the contig name above against the GTF's first column."
  echo "  If they differ, set VIRAL_CONTIG_RENAME=\"<gtf>:<fasta>\" in config.sh."
  awk -F'\t' '!/^#/{print $1}' "$VG" | sort -u | sed 's/^/  GTF contig: /'
fi

echo
echo "reference inputs ready. next: scripts/build_reference.sh"

#!/usr/bin/env bash
# Align one single-cell library against the concatenated reference.
#
# Chemistry comes from config (CB_LEN/UMI_LEN/SOLO_WHITELIST). If CB_LEN is
# empty the script detects it and STOPS -- it will not align on a guess.
set -euo pipefail
ORIG_PWD="$PWD"
cd "$(dirname "$0")/.."
source config/config.sh

SAMPLE="${1:?usage: run_starsolo.sh <sample> <R1.fastq.gz> <R2.fastq.gz>}"
# Resolve input paths BEFORE the cd above takes effect for the caller's
# shell. A relative path typed at the prompt is relative to THEIR cwd, not
# to the repo root this script works from.
abspath() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *)  printf '%s\n' "$ORIG_PWD/$1" ;;
  esac
}
R1="$(abspath "${2:?}")"; R2="$(abspath "${3:?}")"
for f in "$R1" "$R2"; do
  [ -r "$f" ] || { echo "ERROR: cannot read $f" >&2; exit 1; }
done
OUT="$RUN_DIR/$SAMPLE"; mkdir -p "$OUT"
EXTRA_STAR=()

# STAR opens runThreadN x outBAMsortingBinsN temp files at once while
# sorting. With 32 threads and the default 50 bins that is 1,600 against a
# typical soft limit of 1,024, and STAR dies partway through sorting.
NEED=$(( THREADS * 50 + 256 ))
CUR=$(ulimit -n)
if [ "$CUR" != "unlimited" ] && [ "$CUR" -lt "$NEED" ]; then
  HARD=$(ulimit -Hn)
  if [ "$HARD" = "unlimited" ] || [ "$HARD" -ge "$NEED" ]; then
    ulimit -n "$NEED" && echo "raised open-file limit $CUR -> $NEED"
  else
    BINS=$(( (HARD - 256) / THREADS ))
    [ "$BINS" -lt 1 ] && BINS=1
    echo "open-file limit $CUR (hard $HARD) is below $NEED;"
    echo "using --outBAMsortingBinsN $BINS instead"
    EXTRA_STAR+=( --outBAMsortingBinsN "$BINS" )
  fi
fi

# A leftover _STARtmp from a crashed run makes the next run fail for an
# unrelated-looking reason.
rm -rf "$OUT/_STARtmp"

if [ -z "${CB_LEN:-}" ]; then
  echo "CB_LEN is not set in config/config.sh. Detecting..."
  echo
  python3 scripts/barcode_structure.py --r1 "$R1" | tee "$OUT/barcode_structure.log"
  echo
  echo "Set CB_LEN, UMI_LEN and SOLO_WHITELIST in config/config.sh from the"
  echo "verdict above, then re-run. Not aligning on an auto-detected value:"
  echo "chemistry belongs in config where it is visible and reviewable."
  exit 2
fi

# --soloCBmatchWLtype and the CellRanger-flavoured UMI options are only
# valid when a whitelist is in use. With --soloCBwhitelist None, passing
# 1MM_multi_* is an error, not a no-op.
WL_ARGS=( --soloCBwhitelist "$SOLO_WHITELIST" )
if [ "$SOLO_WHITELIST" = "None" ]; then
  WL_ARGS+=( --soloUMIdedup 1MM_All --soloUMIfiltering MultiGeneUMI )
else
  WL_ARGS+=( --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
             --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR )
fi

echo "sample      : $SAMPLE"
echo "chemistry   : CB ${CB_LEN} + UMI ${UMI_LEN}, whitelist ${SOLO_WHITELIST}"
echo "strand      : ${SOLO_STRAND}"
echo

# readFilesIn takes cDNA FIRST, then the barcode read. Reversing these is
# the most common single-cell alignment error and it does not error --
# it produces ~0 valid barcodes.
STAR --runThreadN "$THREADS" \
     --genomeDir "$REF_DIR/star_index" \
     --readFilesIn "$R2" "$R1" \
     --readFilesCommand zcat \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen "$CB_LEN" \
     --soloUMIstart $((CB_LEN+1)) --soloUMIlen "$UMI_LEN" \
     --soloBarcodeReadLength 0 \
     --soloStrand "$SOLO_STRAND" \
     "${WL_ARGS[@]}" \
     --soloFeatures Gene GeneFull \
     --soloCellFilter "$SOLO_CELL_FILTER" \
     --soloMultiMappers Unique EM \
     "${EXTRA_STAR[@]+"${EXTRA_STAR[@]}"}" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
     --outFileNamePrefix "$OUT/"

echo
echo "=== run QC ==="
# Barcodes.stats only populates when a whitelist is in use; with
# --soloCBwhitelist None every counter stays 0 and the block is noise.
# Summary.csv carries the numbers that matter on either path.
if [ -f "$OUT/Solo.out/Gene/Summary.csv" ]; then
  grep -E 'Number of Reads|Valid Barcodes|Saturation|Mapped to Genome: Unique,|Mapped to Gene: Unique Gene|Estimated Number of Cells|Median UMI|Median Gene' \
    "$OUT/Solo.out/Gene/Summary.csv" || cat "$OUT/Solo.out/Gene/Summary.csv"
elif [ -f "$OUT/Solo.out/Barcodes.stats" ]; then
  cat "$OUT/Solo.out/Barcodes.stats"
fi

echo
echo "=== convert to h5ad ==="
# Both feature sets: HSV-1's overlapping ORFs mean some viral reads are
# multi-gene and dropped under Gene. See validate_reference.py [5].
for FEAT in Gene GeneFull; do
  d="$OUT/Solo.out/$FEAT/filtered"
  [ -d "$d" ] || { echo "  no $FEAT/filtered -- skipping"; continue; }
  python3 scripts/solo_to_h5ad.py \
    --solo-dir "$d" \
    --viral-genes-file "$REF_DIR/viral_genes.txt" \
    --sample "$SAMPLE" \
    --out "$OUT/${SAMPLE}.${FEAT}.h5ad" | tee "$OUT/to_h5ad.${FEAT}.log"
done
echo "done: $OUT"

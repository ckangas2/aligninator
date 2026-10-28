#!/usr/bin/env bash
# Write a samplesheet for named samples in config/runs.tsv, using FASTQs
# already downloaded by fetch_runs.sh. Prints CSV to stdout.
#
#   bash scripts/make_samplesheet.sh exp2_iPSC-2_ctrl_rep1 exp2_iPSC-2_3dpi_rep1 > samplesheet.csv
#   bash scripts/make_samplesheet.sh --all > samplesheet.csv
#
# ENA's _1 is the barcode read and _2 the cDNA read for these runs;
# the pipeline's CHECK_READS step will stop if that is not true.
set -euo pipefail
cd "$(dirname "$0")/.."
source config/config.sh
MANIFEST="config/runs.tsv"
[ $# -gt 0 ] || { echo "usage: make_samplesheet.sh <sample>... | --all" >&2; exit 1; }
if [ "$1" = "--all" ]; then
  mapfile -t samples < <(awk -F'\t' '!/^#/ && NF && $1!="sample" {print $1}' "$MANIFEST")
else
  samples=("$@")
fi
# Check every sample before printing anything: a partial failure must not
# leave a header-only CSV behind in the redirect target.
rows=()
for s in "${samples[@]}"; do
  run=$(awk -F'\t' -v s="$s" '!/^#/ && $1==s {print $2}' "$MANIFEST")
  [ -n "$run" ] || { echo "no such sample in $MANIFEST: $s" >&2; exit 1; }
  r1="$FASTQ_DIR/${run}_1.fastq.gz"; r2="$FASTQ_DIR/${run}_2.fastq.gz"
  for f in "$r1" "$r2"; do
    [ -s "$f" ] || { echo "missing $f -- run scripts/fetch_runs.sh --sample $s, or set FASTQ_DIR" >&2; exit 1; }
  done
  rows+=("$s,$(cd "$(dirname "$r1")" && pwd)/$(basename "$r1"),$(cd "$(dirname "$r2")" && pwd)/$(basename "$r2")")
done
echo "sample,fastq_1,fastq_2"
printf '%s\n' "${rows[@]}"

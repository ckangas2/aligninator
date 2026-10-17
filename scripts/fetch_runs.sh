#!/usr/bin/env bash
# Download FASTQs for named runs from ENA.
#
# URLs are queried from ENA at run time rather than constructed. ENA's
# directory layout depends on accession length (an 11-character accession
# nests under a zero-padded last-two-digits directory, a 12-character one
# does not), and hand-rolling that rule is a silent 404 generator.
#
# HTTPS with many connections, not FTP: ENA's anonymous FTP is capped per
# stream at a few hundred KB/s. There is no account that lifts it.
set -euo pipefail
cd "$(dirname "$0")/.."
source config/config.sh

usage() { echo "usage: fetch_runs.sh <SRR...> | --sample <name>... | --all" >&2; exit 1; }
[ $# -gt 0 ] || usage

MANIFEST="config/runs.tsv"
runs=()
case "${1:-}" in
  --all)    mapfile -t runs < <(awk -F'\t' '!/^#/ && $2 && $1!="sample" {print $2}' "$MANIFEST") ;;
  --sample) shift; for s in "$@"; do
              r=$(awk -F'\t' -v s="$s" '!/^#/ && $1==s && $1!="sample" {print $2}' "$MANIFEST")
              [ -n "$r" ] || { echo "no such sample: $s" >&2; exit 1; }
              runs+=("$r"); done ;;
  *)        runs=("$@") ;;
esac

mkdir -p "$FASTQ_DIR"
echo "runs to fetch: ${runs[*]}"
echo

for run in "${runs[@]}"; do
  echo "=== $run ==="
  meta=$(curl -sS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp,fastq_md5,fastq_bytes&format=tsv")
  line=$(echo "$meta" | sed -n '2p')
  [ -n "$line" ] || { echo "  ENA returned nothing for $run -- skipping" >&2; continue; }
  IFS=$'\t' read -r ftp md5 bytes <<< "$line"
  IFS=';' read -ra urls  <<< "$ftp"
  IFS=';' read -ra md5s  <<< "$md5"

  if [ "${#urls[@]}" -lt 2 ]; then
    echo "  WARNING: only ${#urls[@]} file(s) for $run." >&2
    echo "  A 10x run needs both the barcode read and the cDNA read." >&2
    echo "  If the barcode read was not deposited this run is unusable." >&2
  fi

  for i in "${!urls[@]}"; do
    u="https://ftp.sra.ebi.ac.uk/${urls[$i]#ftp.sra.ebi.ac.uk/}"
    f="$FASTQ_DIR/$(basename "${urls[$i]}")"
    want="${md5s[$i]:-}"
    if [ -f "$f" ] && [ -n "$want" ]; then
      got=$(md5sum "$f" | cut -d' ' -f1)
      if [ "$got" = "$want" ]; then echo "  ok (verified): $(basename "$f")"; continue; fi
    fi
    echo "  -> $(basename "$f")"
    if command -v aria2c >/dev/null 2>&1; then
      aria2c -x"$ARIA_CONN" -s"$ARIA_CONN" -c --summary-interval=15 \
             -d "$FASTQ_DIR" -o "$(basename "$f")" "$u"
    else
      echo "  (aria2c not found -- falling back to wget, expect ~500 KB/s)" >&2
      wget -c -O "$f" "$u"
    fi
    if [ -n "$want" ]; then
      got=$(md5sum "$f" | cut -d' ' -f1)
      [ "$got" = "$want" ] || { echo "  MD5 MISMATCH on $(basename "$f")" >&2; exit 1; }
      echo "  md5 ok"
    fi
  done
done
echo
echo "fastq dir: $FASTQ_DIR"
echo "next: scripts/inspect_fastq.sh $FASTQ_DIR/<run>_1.fastq.gz $FASTQ_DIR/<run>_2.fastq.gz"

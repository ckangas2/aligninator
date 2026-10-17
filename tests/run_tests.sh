#!/usr/bin/env bash
# Synthetic end-to-end. No real data, no STAR, no network.
set -uo pipefail
cd "$(dirname "$0")/.."
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
pass=0; fail=0
chk() { # chk <name> <expected_exit> <cmd...>
  local n="$1" want="$2"; shift 2
  "$@" >"$T/out.txt" 2>&1; local got=$?
  if [ "$got" -eq "$want" ]; then echo "  PASS  $n"; pass=$((pass+1));
  else echo "  FAIL  $n (exit $got, wanted $want)"; sed 's/^/        /' "$T/out.txt" | tail -5; fail=$((fail+1)); fi
}

echo "building synthetic fixtures..."
python3 tests/make_synthetic_fastq.py "$T/fq"  >/dev/null
python3 tests/make_synthetic_ref.py   "$T/ref" >/dev/null
python3 tests/make_synthetic_solo.py  "$T/solo_prefixed" "HSV1-" >/dev/null
python3 tests/make_synthetic_solo.py  "$T/solo_bare" "" >/dev/null
printf 'HSV1-%s\n' RL1 RL2 LAT RS1 UL1 UL2 UL3 UL5 US1 US12 GFP > "$T/viral.txt"

echo "chemistry detection:"
chk "calls v2"                    0 python3 -I scripts/detect_chemistry.py --r1 "$T/fq/v2_R1.fastq.gz"        --whitelist-dir "$T/fq/whitelists"
chk "calls v3"                    0 python3 -I scripts/detect_chemistry.py --r1 "$T/fq/v3_R1.fastq.gz"        --whitelist-dir "$T/fq/whitelists"
chk "refuses trimmed R1"          2 python3 -I scripts/detect_chemistry.py --r1 "$T/fq/trimmed_R1.fastq.gz"   --whitelist-dir "$T/fq/whitelists"
chk "refuses R2 passed as R1"     2 python3 -I scripts/detect_chemistry.py --r1 "$T/fq/r2_by_mistake.fastq.gz" --whitelist-dir "$T/fq/whitelists"
chk "refuses wrong whitelist"     2 python3 -I scripts/detect_chemistry.py --r1 "$T/fq/crosslist_R1.fastq.gz" --whitelist-dir "$T/fq/whitelists"

echo "reference validation:"
chk "rejects as-deposited GTF"    1 python3 -I scripts/validate_reference.py --viral-fasta "$T/ref/viral.fa" --viral-gtf "$T/ref/viral_bad.gtf"     --host-gtf "$T/ref/host.gtf"
chk "rejects CDS-only GTF"        1 python3 -I scripts/validate_reference.py --viral-fasta "$T/ref/viral.fa" --viral-gtf "$T/ref/viral_cdsonly.gtf" --host-gtf "$T/ref/host.gtf"
chk "accepts clean GTF"           0 python3 -I scripts/validate_reference.py --viral-fasta "$T/ref/viral.fa" --viral-gtf "$T/ref/viral_good.gtf"    --host-gtf "$T/ref/host.gtf"

echo "gtf normalisation:"
python3 -I scripts/prepare_viral_gtf.py --in-gtf "$T/ref/viral_bad.gtf" --out-gtf "$T/ref/fixed.gtf" \
  --gene-list-out "$T/ref/genes.txt" --rename-contig HSV1:chrHSV1 --prefix HSV1- --add-gene-rows --copy-id-to-name >/dev/null
chk "bad GTF becomes valid"       0 python3 -I scripts/validate_reference.py --viral-fasta "$T/ref/viral.fa" --viral-gtf "$T/ref/fixed.gtf" --host-gtf "$T/ref/host.gtf"

echo "barcode structure (whitelist-free chemistry):"
python3 -I tests/make_synthetic_barcodes.py "$T/bc" >/dev/null
bcchk() { # bcchk <file> <want_cb> <want_umi>
  local got
  got=$(python3 -I scripts/barcode_structure.py --r1 "$T/bc/$1_R1.fastq.gz" --n-reads 200000 2>&1 \
        | grep -o 'cell barcode = [0-9]* bp, UMI = [0-9]* bp')
  if [ "$got" = "cell barcode = $2 bp, UMI = $3 bp" ]; then
    echo "  PASS  $1 -> CB$2/UMI$3"; pass=$((pass+1))
  else echo "  FAIL  $1 wanted CB$2/UMI$3, got: $got"; fail=$((fail+1)); fi
}
bcchk dropseq 12 8
bcchk tenx_v2 16 10
bcchk tenx_v3 16 12

echo "matrix conversion:"
chk "scores matching reference"   0 python3 -I scripts/solo_to_h5ad.py --solo-dir "$T/solo_prefixed" --viral-genes-file "$T/viral.txt" --out "$T/a.h5ad"
python3 -I scripts/solo_to_h5ad.py --solo-dir "$T/solo_bare" --viral-genes-file "$T/viral.txt" --out "$T/b.h5ad" >"$T/bare.txt" 2>&1
if grep -q "MISSING 11" "$T/bare.txt"; then echo "  PASS  reports mismatched reference"; pass=$((pass+1));
else echo "  FAIL  reports mismatched reference"; fail=$((fail+1)); fi

echo; echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]

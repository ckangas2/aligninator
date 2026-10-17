#!/usr/bin/env bash
# Report read lengths and the SJDB_OVERHANG to set. Run before building
# the reference: sjdbOverhang must be (cDNA read length - 1), and getting
# it wrong degrades splice-junction mapping without any error.
set -euo pipefail
R1="${1:?usage: inspect_fastq.sh <R1.fastq.gz> <R2.fastq.gz>}"
R2="${2:?}"
len() { zcat "$1" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -rn; }
echo "=== $(basename "$R1") (expect the SHORT read: barcode+UMI) ==="; len "$R1"
echo; echo "=== $(basename "$R2") (expect the LONG read: cDNA) ==="; len "$R2"
l1=$(zcat "$R1" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
l2=$(zcat "$R2" | head -4000 | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -rn | head -1 | awk '{print $2}')
echo
if [ "$l1" -gt "$l2" ]; then
  echo "WARNING: file 1 is LONGER than file 2. These are probably swapped."
  echo "STARsolo takes cDNA first, barcode second -- check before aligning."
fi
echo "R1 = ${l1} bp   R2 = ${l2} bp"
echo "set in config/config.sh:  SJDB_OVERHANG=$((l2-1))"
case "$l1" in
  26) echo "R1 26 bp -> 10x v2 (CB16+UMI10) or 5' v1/v2. detect_chemistry.py will disambiguate by whitelist." ;;
  28) echo "R1 28 bp -> 10x v3 (CB16+UMI12)." ;;
  *)  echo "R1 ${l1} bp is not a standard 10x barcode read length. Investigate before proceeding." ;;
esac

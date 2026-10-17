#!/usr/bin/env bash
# Build a concatenated host+virus STAR index.
#
# NOT YET RUN ON REAL DATA -- written against STAR 2.7.11 docs, validated
# only for argument structure. Run step 1 and 2 first; they are cheap and
# catch the errors that matter before genomeGenerate burns an hour.
set -euo pipefail
cd "$(dirname "$0")/.."
source config/config.sh

mkdir -p "$REF_DIR" "$LOG_DIR"

echo "=== 1. normalise the viral GTF ==="
python3 scripts/prepare_viral_gtf.py \
  --in-gtf   "$VIRAL_GTF_RAW" \
  --out-gtf  "$REF_DIR/viral.normalised.gtf" \
  --gene-list-out "$REF_DIR/viral_genes.txt" \
  ${VIRAL_CONTIG_RENAME:+--rename-contig "$VIRAL_CONTIG_RENAME"} \
  --prefix "$VIRAL_PREFIX" --add-gene-rows --copy-id-to-name \
  | tee "$LOG_DIR/prepare_viral_gtf.log"

echo
echo "=== 2. validate before building (fails fast) ==="
python3 scripts/validate_reference.py \
  --viral-fasta "$VIRAL_FASTA" \
  --viral-gtf   "$REF_DIR/viral.normalised.gtf" \
  --host-gtf    "$HOST_GTF" \
  --viral-prefix "$VIRAL_PREFIX" \
  | tee "$LOG_DIR/validate_reference.log"

echo
echo "=== 3. concatenate FASTA and GTF ==="
# Order matters only for reproducibility, not correctness -- but keep host
# first so contig order matches the host-only reference people compare to.
zcat -f "$HOST_FASTA" "$VIRAL_FASTA" > "$REF_DIR/genome.fa"
zcat -f "$HOST_GTF" > "$REF_DIR/genes.gtf"
cat "$REF_DIR/viral.normalised.gtf" >> "$REF_DIR/genes.gtf"
grep -c '^>' "$REF_DIR/genome.fa" | xargs echo "contigs in merged FASTA:"

echo
echo "=== 4. STAR genomeGenerate ==="
# sjdbOverhang should be (cDNA read length - 1). For 10x the cDNA read is R2.
STAR --runMode genomeGenerate \
     --runThreadN "$THREADS" \
     --genomeDir "$REF_DIR/star_index" \
     --genomeFastaFiles "$REF_DIR/genome.fa" \
     --sjdbGTFfile "$REF_DIR/genes.gtf" \
     --sjdbOverhang "$SJDB_OVERHANG" \
     --limitGenomeGenerateRAM "$GENOME_RAM" \
     --outFileNamePrefix "$LOG_DIR/genomeGenerate_"

echo
echo "=== 5. confirm viral genes survived into the index ==="
# If this count is 0 the reference is broken regardless of what STAR printed.
n=$(grep -c "^${VIRAL_PREFIX}" "$REF_DIR/star_index/geneInfo.tab" || true)
echo "viral genes in index: $n"
if [ "$n" -eq 0 ]; then
  echo "FATAL: no viral genes in the index. Do not align against this." >&2
  exit 1
fi
echo "reference ready: $REF_DIR/star_index"

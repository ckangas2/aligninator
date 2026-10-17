# All settings. Every value is overridable from the environment.
: "${THREADS:=32}"
: "${GENOME_RAM:=120000000000}"
: "${SJDB_OVERHANG:=62}"          # cDNA read length - 1 (GSE163952 R2 is 53-63bp, mode 63)

: "${PROJECT_DIR:=$PWD}"
: "${REF_DIR:=$PROJECT_DIR/reference}"
: "${RUN_DIR:=$PROJECT_DIR/runs}"
: "${LOG_DIR:=$PROJECT_DIR/logs}"
: "${WHITELIST_DIR:=$PROJECT_DIR/reference/whitelists}"

# --- host -------------------------------------------------------------
# Use the SAME build the PDAC atlas was made on (CellRanger GRCh38 2020-A)
# so gene symbols join without an alias map.
: "${HOST_FASTA:=$REF_DIR/GRCh38/fasta/genome.fa}"
: "${HOST_GTF:=$REF_DIR/GRCh38/genes/genes.gtf}"

# --- virus ------------------------------------------------------------
# GSE163952 ships its own FASTA + GTF. Verify the contig name in the FASTA
# header and set VIRAL_CONTIG_RENAME only if the GTF disagrees with it.
: "${VIRAL_FASTA:=$REF_DIR/hsv1/GSE163952_HSV1GFPtrunc_corr.fa}"
: "${VIRAL_GTF_RAW:=$REF_DIR/hsv1/GSE163952_HSV1GFPtrunc_corr.gtf}"
: "${VIRAL_CONTIG_RENAME:=}"      # e.g. "HSV1:chrHSV1" -- leave empty if they match
: "${VIRAL_PREFIX:=HSV1-}"        # prevents LAT/RS1 collision with host symbols

# --- downloads ---------------------------------------------------------
: "${FASTQ_DIR:=$PROJECT_DIR/fastq}"
: "${ARIA_CONN:=16}"              # parallel connections; ENA FTP single-stream is ~500 KB/s
: "${GEO_SERIES:=GSE163952}"

# --- chemistry ---------------------------------------------------------
# GSE163952 is Drop-seq, NOT 10x: R1 is 20 bp = 12 bp CB + 8 bp UMI, and
# the barcode-collapse sweep (scripts/barcode_structure.py) puts the
# plateau at 12. Drop-seq has no barcode whitelist; barcodes are called
# from the data, so SOLO_WHITELIST is the literal string None.
# Leave CB_LEN empty to make run_starsolo.sh detect and stop instead.
: "${CB_LEN:=12}"
: "${UMI_LEN:=8}"
: "${SOLO_WHITELIST:=None}"
: "${SOLO_STRAND:=Forward}"
: "${SOLO_CELL_FILTER:=EmptyDrops_CR}"

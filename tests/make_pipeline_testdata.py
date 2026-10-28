#!/usr/bin/env python3
"""Generate tests/data/: a tiny host + virus reference and two Drop-seq
libraries that the Nextflow `test` profile runs end to end.

Deterministic (fixed seed), so the committed files only change when this
script does. Re-generate with:

    python3 tests/make_pipeline_testdata.py tests/data

What is planted, and why:
  * The viral GTF is shipped broken the way real deposits are: its contig
    is `HSV1` while the FASTA says `chrHSV1`, rows carry gene_id but no
    gene_name, and the names LAT and RS1 also exist in the host GTF. The
    pipeline must normalise all three or the run fails / miscounts.
  * `mock` has no viral reads. `infected` has viral reads in ~70% of cells,
    at 3-50% of each cell's transcripts, plus a GFP reporter that must be
    excluded from viral_frac.
  * Each library has a few hundred near-empty barcodes, so cell calling
    has something to reject.
"""
import gzip, random, sys
from pathlib import Path

SEED      = 1
CB_LEN    = 12
UMI_LEN   = 8
READ_LEN  = 63           # cDNA read; sjdb_overhang = 62
N_CELLS   = 100
READS_PER_CELL = (180, 320)
N_EMPTY   = 400
INFECTED_SHARE = 0.70

rng = random.Random(SEED)
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]
def rand_seq(n): return "".join(rng.choice("ACGT") for _ in range(n))

def write_fasta(path, contigs):
    with gzip.open(path, "wt", compresslevel=9) as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

# ---------------------------------------------------------------- host
def make_host():
    contigs = {f"chr{i}": rand_seq(80_000) for i in (1, 2, 3)}
    names = [f"HOSTG{i}" for i in range(1, 27)] + ["LAT", "RS1"]   # collide with HSV-1
    genes, rows = [], []
    slot = 0
    for gi, name in enumerate(names):
        chrom = f"chr{slot % 3 + 1}"
        start = 2_000 + (slot // 3) * 7_000
        slot += 1
        strand = "+" if gi % 2 == 0 else "-"
        e1 = (start, start + 399)
        e2 = (start + 900, start + 1_499)            # 500 bp intron
        gid, tid = f"ENSGSYN{gi:05d}", f"ENSTSYN{gi:05d}"
        a = f'gene_id "{gid}"; gene_name "{name}"; gene_type "protein_coding";'
        t = a + f' transcript_id "{tid}";'
        rows.append((chrom, "gene", e1[0], e2[1], strand, a))
        rows.append((chrom, "transcript", e1[0], e2[1], strand, t))
        for n, (s, e) in enumerate((e1, e2), 1):
            rows.append((chrom, "exon", s, e, strand, t + f' exon_number "{n}";'))
        seq = contigs[chrom]
        mrna = seq[e1[0] - 1:e1[1]] + seq[e2[0] - 1:e2[1]]
        genes.append((name, mrna if strand == "+" else rc(mrna)))
    return contigs, rows, genes

# ---------------------------------------------------------------- virus
def make_virus():
    seq = rand_seq(25_000)
    # Single-exon, as HSV-1 genes essentially are. UL2 overlaps UL1 on the
    # same strand, like the real genome's dense ORF layout.
    layout = [("RL1", 500, 1_600, "+"), ("LAT", 2_000, 3_400, "+"),
              ("RS1", 3_800, 5_400, "-"), ("UL1", 6_000, 7_100, "-"),
              ("UL2", 6_900, 8_000, "-"), ("UL3", 8_500, 9_400, "+"),
              ("UL4", 9_800, 10_900, "+"), ("UL5", 11_300, 12_800, "-"),
              ("UL6", 13_200, 14_300, "+"), ("UL8", 14_800, 16_000, "-"),
              ("US1", 17_000, 18_200, "+"), ("US12", 19_000, 19_700, "-"),
              ("GFP", 21_000, 21_720, "+")]
    rows, genes = [], []
    for name, s, e, strand in layout:
        # As deposited: wrong contig name, gene_id only, no gene_name.
        rows.append(("HSV1", "exon", s, e, strand,
                     f'gene_id "{name}"; transcript_id "{name}.1";'))
        m = seq[s - 1:e]
        genes.append((name, m if strand == "+" else rc(m)))
    return {"chrHSV1": seq}, rows, genes

def write_gtf(path, rows, source):
    with gzip.open(path, "wt", compresslevel=9) as fh:
        for chrom, feat, s, e, strand, attrs in rows:
            fh.write(f"{chrom}\t{source}\t{feat}\t{s}\t{e}\t.\t{strand}\t.\t{attrs}\n")

# ---------------------------------------------------------------- reads
def sample_read(mrna):
    i = rng.randint(0, len(mrna) - READ_LEN)
    return mrna[i:i + READ_LEN]

def make_library(path_prefix, host_genes, viral_genes, infected):
    # Each cell has its own host expression profile; viral reads, if any,
    # replace a per-cell fraction of them.
    used, records = set(), []
    def new_cb():
        while True:
            cb = rand_seq(CB_LEN)
            if cb not in used:
                used.add(cb); return cb

    gfp = [g for g in viral_genes if g[0] == "GFP"]
    viral_no_gfp = [g for g in viral_genes if g[0] != "GFP"]

    for _ in range(N_CELLS):
        cb = new_cb()
        n = rng.randint(*READS_PER_CELL)
        weights = [rng.random() ** 2 for _ in host_genes]
        vfrac = rng.uniform(0.03, 0.50) if (infected and rng.random() < INFECTED_SHARE) else 0.0
        umis = []
        for _ in range(n):
            # ~15% PCR duplicates: reuse a UMI and the same molecule's sequence.
            if umis and rng.random() < 0.15:
                umi, seq = rng.choice(umis)
            else:
                if rng.random() < vfrac:
                    pool = gfp if rng.random() < 0.10 else viral_no_gfp
                    _, mrna = rng.choice(pool)
                else:
                    _, mrna = rng.choices(host_genes, weights)[0]
                umi, seq = rand_seq(UMI_LEN), sample_read(mrna)
                umis.append((umi, seq))
            records.append((cb + umi, seq))

    for _ in range(N_EMPTY):                       # ambient / empty droplets
        cb = new_cb()
        for _ in range(rng.randint(1, 8)):
            _, mrna = rng.choice(host_genes)
            records.append((cb + rand_seq(UMI_LEN), sample_read(mrna)))

    rng.shuffle(records)
    with gzip.open(f"{path_prefix}_R1.fastq.gz", "wt", compresslevel=9) as f1, \
         gzip.open(f"{path_prefix}_R2.fastq.gz", "wt", compresslevel=9) as f2:
        for i, (bc, cdna) in enumerate(records):
            f1.write(f"@r{i}\n{bc}\n+\n{'I' * len(bc)}\n")
            f2.write(f"@r{i}\n{cdna}\n+\n{'I' * len(cdna)}\n")
    return len(records)

def main(out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    host, host_rows, host_genes = make_host()
    virus, viral_rows, viral_genes = make_virus()
    write_fasta(out / "host.fa.gz", host)
    write_gtf(out / "host.gtf.gz", host_rows, "synthetic")
    write_fasta(out / "virus.fa.gz", virus)
    write_gtf(out / "virus.gtf.gz", viral_rows, "deposit")
    n_mock = make_library(out / "mock", host_genes, viral_genes, infected=False)
    n_inf  = make_library(out / "infected", host_genes, viral_genes, infected=True)
    (out / "samplesheet.csv").write_text(
        "sample,fastq_1,fastq_2\n"
        "mock,mock_R1.fastq.gz,mock_R2.fastq.gz\n"
        "infected,infected_R1.fastq.gz,infected_R2.fastq.gz\n")
    print(f"host: {len(host_genes)} genes, virus: {len(viral_genes)} genes")
    print(f"reads: mock {n_mock:,}, infected {n_inf:,}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/data")

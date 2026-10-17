# aligninator

STARsolo pipeline for aligning single-cell RNA-seq against a concatenated
host + virus reference, so viral transcripts are counted per cell instead
of discarded. Built for HSV-1; the host/virus logic is not HSV-specific.

Most published oncolytic-HSV single-cell data was aligned to a human-only
reference. The viral reads were sequenced; they were simply never counted,
and the deposited matrices cannot give a per-cell viral label. This
pipeline rebuilds that label from FASTQ.

**Validated end to end** on GSE163952: mock organoids give 0.0% of cells
above 1% viral load, 3 dpi gives 83.3%, through an identical reference and
code path. Full numbers in [RESULTS.md](RESULTS.md).

## Prior art, and what this is not

Counting viral reads in single-cell data is not an open problem.
[scPathoQuant](https://academic.oup.com/bioinformatics/article/40/4/btae145/7628127)
(Bioinformatics 2024) takes reads that failed to align to the host genome
and maps them to a pathogen genome, adding the counts to standard
single-cell matrices.
[Viral-Track](https://www.sciencedirect.com/science/article/pii/S0092867420305687)
(Bost et al., Cell 2020) does a related job. If you want a maintained,
published tool, use one of those.

This repository takes a different route — one concatenated reference, so
host-vs-virus assignment happens in a single alignment pass rather than
downstream of a host-only alignment that already rejected the reads. It
also handles whitelist-free chemistries (Drop-seq, Seq-Well), which the
10x-oriented tools above do not advertise.

The part worth borrowing is not the alignment. It is the validation layer:
every check in the table below exists because the corresponding failure
produces a successful run and a plausible, wrong matrix.

## Why this is not a thin wrapper around STAR

Every failure below produces a **successful run with a plausible matrix**.
None of them raises an error. Each is checked explicitly.

| failure | what you see | what it actually is |
|---|---|---|
| GTF contig name disagrees with FASTA header | `viral_frac == 0` in every cell | STAR dropped every viral annotation; reads became intergenic |
| viral GTF has CDS rows but no `exon` rows | zero viral features in the matrix | STARsolo counts to exons |
| viral `gene_name` collides with a host symbol | one gene silently renamed `LAT-1` | HSV-1 carries **LAT** and **RS1**; so does the human genome. `var_names_make_unique()` renames by *position*, so which one keeps the bare symbol depends on file order |
| viral genes matched by string prefix | no infection detected, or host genes counted as viral | on an unprefixed reference `startswith("HSV1-")` finds nothing, and `startswith("LAT")` finds the host gene |
| wrong 10x chemistry assumed | a matrix with a believable cell count | read length alone cannot separate 3' v2 from 5' v1 |
| chemistry assumed to be 10x at all | every barcode invalid, or an empty matrix | GSE163952 is **Drop-seq** (R1 = 20 bp = CB12+UMI8), not 10x, despite being described as 10x Chromium in several places |
| `--soloCBmatchWLtype` passed without a whitelist | STAR errors, or silently misbehaves | those options are only valid alongside `--soloCBwhitelist <file>` |
| `readFilesIn` given R1 before R2 | ~0 valid barcodes | STARsolo takes **cDNA first**, barcode read second |
| overlapping viral ORFs | viral counts below true depth | STARsolo assigns a read only when exactly one gene matches; HSV-1 is gene-dense and overlapping |

## Layout

```
scripts/
  detect_chemistry.py     infer 10x chemistry by whitelist match; refuses to guess
  barcode_structure.py    infer CB/UMI length with NO whitelist, from read collapse
  prepare_viral_gtf.py    rename contig, prefix gene names, synthesise gene rows
  validate_reference.py   pre-flight checks; run BEFORE genomeGenerate
  build_reference.sh      normalise -> validate -> concatenate -> STAR index
  run_starsolo.sh         detect chemistry -> align -> convert
  solo_to_h5ad.py         matrix -> .h5ad with per-cell viral_frac
config/config.sh          all settings, environment-overridable
tests/                    synthetic data generators, one per failure mode
```

## Quickstart — minimal working example

Two libraries from GSE163952, mock and 3 dpi from the same iPSC line.
~2.5 GB of FASTQ. The mock arm is the negative control: if it does not
come back at `viral_frac` near zero while 3 dpi is clearly positive,
something is wrong and you know it before committing to the full 83 GB.

```bash
bash tests/run_tests.sh                       # synthetic, no data needed
bash scripts/fetch_reference.sh               # GRCh38 2020-A + whitelists + viral FASTA/GTF
bash scripts/fetch_runs.sh --sample exp2_iPSC-2_ctrl_rep1 exp2_iPSC-2_3dpi_rep1
bash scripts/inspect_fastq.sh fastq/SRR21855376_1.fastq.gz fastq/SRR21855376_2.fastq.gz
#   -> set SJDB_OVERHANG in config/config.sh from the reported R2 length
bash scripts/build_reference.sh
bash scripts/run_starsolo.sh exp2_iPSC-2_ctrl_rep1 fastq/SRR21855376_1.fastq.gz fastq/SRR21855376_2.fastq.gz
bash scripts/run_starsolo.sh exp2_iPSC-2_3dpi_rep1 fastq/SRR21855380_1.fastq.gz fastq/SRR21855380_2.fastq.gz
```

`config/runs.tsv` lists all 22 scRNA-seq runs with condition, replicate and
download size. `fetch_runs.sh --all` takes the lot (83.2 GB).

### Downloading

ENA's anonymous FTP is capped per stream at roughly 500 KB/s and no account
changes that. `fetch_runs.sh` uses HTTPS with 16 parallel connections via
`aria2c` (`conda install -c conda-forge aria2`), which is typically 10–30×
faster, and verifies every file against ENA's MD5.

## Finding the barcode without a whitelist

`detect_chemistry.py` identifies 10x libraries against a published
whitelist. Drop-seq, Seq-Well and inDrop have no whitelist, so that
method correctly reports no call and stops — useful as a guard, useless
as an answer.

`barcode_structure.py` uses the signal that does not need a list: a real
barcode **collapses reads onto itself**. Sweep the prefix length and
measure the share of reads held by the top K prefixes. Inside the
barcode the share sits on a plateau — adding a base changes nothing,
the same cells still hold the same reads. The first base of UMI shatters
every barcode and the share falls. The barcode is the last length on the
plateau.

On GSE163952, 100k reads:

```
 len    distinct   reads/prefix   reads in topK     drop
  11      71,399           2.80          55.6%    1.000
  12      71,858           2.78          55.6%    1.000   <- barcode ends here
  13      77,923           2.57          39.6%    0.713
  16     150,966           1.32          11.1%    0.645
```

Note the trap: the *steepest* fall is not at the barcode boundary. The
share keeps decaying toward K/n, so the global minimum lands near the
read end and would call CB = 17. The plateau edge is the answer, not the
biggest drop.

## Design decisions worth stating

**Viral features come from a file, not a prefix.** `prepare_viral_gtf.py`
writes `viral_genes.txt` from the GTF itself, and `solo_to_h5ad.py` reads
it. If a listed feature is missing from the matrix the script says so
loudly rather than scoring on a subset — a mismatch means the reference
used for alignment is not the one the list came from.

**GFP is excluded from `viral_frac`.** It is driven by an artificial
promoter and does not track viral transcription. It stays in the matrix
and is reported separately as a correlation check.

**Both `Gene` and `GeneFull` are produced.** HSV-1's overlapping ORFs mean
a real share of viral reads are multi-gene and dropped under `Gene`.
`validate_reference.py` reports the overlap fraction for the specific
reference being built, so the choice is made against a number.

## Status

Python components are tested against synthetic data that plants each
failure mode (14 tests, `bash tests/run_tests.sh`). The full pipeline has
been run end to end on GSE163952 with STAR 2.7.11b — see
[RESULTS.md](RESULTS.md).

Known limitations:

- `fetch_reference.sh` cannot download the 10x barcode whitelists; the URL
  it tries is wrong. They ship inside CellRanger. This does not affect
  whitelist-free chemistries, which need no list.
- Reconciliation against the deposited Seurat object is incomplete: that
  object holds 7,226 cells across 22 sequencing runs, so it is a subset,
  and which runs it covers has not been established.
- Only Drop-seq and 10x v2/v3 are in the chemistry table.

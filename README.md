# aligninator

[![CI](https://github.com/ckangas2/aligninator/actions/workflows/ci.yml/badge.svg)](https://github.com/ckangas2/aligninator/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-23aa62.svg)](https://www.nextflow.io/)
[![Container](https://img.shields.io/badge/container-ghcr.io%2Fckangas2%2Faligninator-blue)](https://github.com/ckangas2/aligninator/pkgs/container/aligninator)

A Nextflow pipeline that aligns single-cell RNA-seq against a concatenated
host + virus reference with STARsolo, so viral transcripts are counted per
cell instead of discarded. Built for HSV-1; nothing in the host/virus logic
is HSV-specific.

Most published oncolytic-HSV single-cell data was aligned to a human-only
reference. The viral reads were sequenced; they were simply never counted,
and the deposited matrices cannot give a per-cell viral label. This
pipeline rebuilds that label from FASTQ.

**Validated end to end** on GSE163952: mock organoids give 0.0% of cells
above 1% viral load, 3 dpi gives 83.3%, through an identical reference and
code path. Full numbers in [RESULTS.md](RESULTS.md).

## Try it in two minutes

You need [Nextflow](https://www.nextflow.io/docs/latest/install.html)
(≥ 24.04, which needs Java 17+) and Docker. Nothing else; STAR and Python
come in the container.

```bash
nextflow run ckangas2/aligninator -r v0.2.0 -profile test,docker
```

This builds a small host + virus reference, aligns a mock library and an
infected one (synthetic Drop-seq data in `tests/data/`), and writes
`results/summary.tsv`:

```
  sample feature_set  cells  median_viral_frac  pct_cells_gt_1pct  pct_cells_gt_5pct
infected        Gene    100             0.1330               66.0               60.0
    mock        Gene    100             0.0000                0.0                0.0
```

The test data is deliberately broken the way real deposits are — the viral
GTF names its contig differently from the FASTA, has no `gene_name`, and
uses symbols (LAT, RS1) that also exist in the host — so the run exercises
the validation layer rather than a happy path. CI runs it on every push and
checks the outcome with `tests/check_test_results.py`.

No Docker? Use `-profile test,singularity`, `-profile test,apptainer` or
`-profile test,conda` instead.

## Run it on your data

**1. A samplesheet.** CSV, one row per library. `fastq_1` is the barcode
read (cell barcode + UMI), `fastq_2` is the cDNA read. Relative paths are
resolved against the samplesheet's own folder.

```csv
sample,fastq_1,fastq_2
mock_rep1,fastq/SRR21855376_1.fastq.gz,fastq/SRR21855376_2.fastq.gz
hsv_3dpi_rep1,fastq/SRR21855380_1.fastq.gz,fastq/SRR21855380_2.fastq.gz
```

**2. A host and a viral reference**, each as FASTA + GTF. Use the host
build your downstream atlas was made on, so gene symbols join without an
alias map.

**3. The chemistry, stated explicitly.** The pipeline never aligns on a
guessed barcode layout.

```bash
nextflow run ckangas2/aligninator -r v0.2.0 -profile docker \
    --input samplesheet.csv \
    --host_fasta GRCh38/genome.fa   --host_gtf GRCh38/genes.gtf \
    --viral_fasta hsv1.fa           --viral_gtf hsv1.gtf \
    --cb_len 12 --umi_len 8 --whitelist None \
    --sjdb_overhang 62 \
    --outdir results
```

Don't know the chemistry? Leave `--cb_len` off. The pipeline runs
`barcode_structure.py` on every library, writes the verdict to
`results/barcode_structure/`, and stops. It works without a whitelist, so
it covers Drop-seq and Seq-Well as well as 10x
([how](#finding-the-barcode-without-a-whitelist)). For 10x, pass the
whitelist file that ships with CellRanger, e.g.
`--cb_len 16 --umi_len 12 --whitelist 3M-february-2018.txt`.

`--sjdb_overhang` should be the cDNA read length minus one. Each library's
`read_check.log` reports the ideal value for its reads.

### Key parameters

| parameter | default | |
|---|---|---|
| `--input` | — | samplesheet CSV |
| `--host_fasta`, `--host_gtf`, `--viral_fasta`, `--viral_gtf` | — | reference inputs (plain or gzipped) |
| `--viral_contig_rename` | — | `OLD:NEW` when the viral GTF's contig name differs from its FASTA |
| `--viral_prefix` | `HSV1-` | prefixed onto viral gene names to prevent host symbol collisions |
| `--star_index`, `--viral_genes` | — | reuse a built index; both together, from the same build |
| `--cb_len`, `--umi_len`, `--whitelist` | — | chemistry; `--whitelist None` for whitelist-free |
| `--solo_strand` | `Forward` | |
| `--solo_cell_filter` | `EmptyDrops_CR` | passed to `--soloCellFilter` |
| `--solo_features` | `Gene GeneFull` | see [Gene vs GeneFull](RESULTS.md#gene-vs-genefull) |
| `--exclude_from_viral` | `GFP` | reporter genes kept in the matrix but left out of `viral_frac` |
| `--save_bam` / `--save_reference` | `false` | publish BAMs / the STAR index |
| `--max_cpus`, `--max_memory`, `--max_time` | 16, 64 GB, 24 h | ceilings for any single task |

`nextflow run ckangas2/aligninator --help` prints the full list.

**Resources.** Indexing GRCh38 + virus needs about 32 GB RAM and takes
~16 min on 32 cores; each library then aligns in minutes. To run on a
cluster, add an executor profile to a config file and pass it with `-c`
(see the [Nextflow executor docs](https://www.nextflow.io/docs/latest/executor.html)).
Pass `-resume` after any failure to pick up where it stopped.

To skip re-indexing on later runs, keep the index with
`--save_reference` and pass it back as
`--star_index results/reference/star_index --viral_genes results/reference/viral_genes.txt`.

### Outputs

```
results/
  summary.tsv                       one row per sample x feature set: cells,
                                    viral_frac distribution, STARsolo QC
  reference/
    viral.normalised.gtf            viral GTF after rename / prefix / gene rows
    viral_genes.txt                 authoritative list of viral feature names
    validate_reference.log          pre-flight checks
    star_index_check.log            viral genes confirmed inside the index
  samples/<sample>/
    <sample>.Gene.h5ad              cells x genes, obs: viral_frac, viral_counts,
    <sample>.GeneFull.h5ad            total_counts, gfp_counts, sample
    Solo.out/                       raw STARsolo output (raw + filtered matrices)
    <sample>.star_summary.csv       STARsolo Summary.csv
    <sample>.read_check.log         read lengths vs configured chemistry
  pipeline_info/                    Nextflow report, timeline, trace, DAG
```

## Reproduce the validation run (GSE163952)

Two libraries, mock and 3 dpi from the same iPSC line, ~2.5 GB of FASTQ.
The mock arm is the negative control: if it does not come back at
`viral_frac` near zero while 3 dpi is clearly positive, something is wrong
and you know it before committing to the full 83 GB.

```bash
git clone https://github.com/ckangas2/aligninator && cd aligninator
bash scripts/fetch_reference.sh               # GRCh38 2020-A + HSV-1 FASTA/GTF from GEO
bash scripts/fetch_runs.sh --sample exp2_iPSC-2_ctrl_rep1 exp2_iPSC-2_3dpi_rep1
bash scripts/make_samplesheet.sh exp2_iPSC-2_ctrl_rep1 exp2_iPSC-2_3dpi_rep1 > samplesheet.csv
nextflow run . -profile gse163952,docker --input samplesheet.csv
```

The `gse163952` profile carries the reference paths and the chemistry
(Drop-seq, CB 12 + UMI 8, no whitelist —
[evidence](docs/CORRECTIONS.md)). `config/runs.tsv` lists all 22 runs with
condition, replicate and download size; `fetch_runs.sh --all` takes the lot.

ENA's anonymous FTP is capped at roughly 500 KB/s per stream.
`fetch_runs.sh` uses HTTPS with 16 parallel connections via `aria2c`
(`conda install -c conda-forge aria2`), typically 10–30× faster, and
verifies every file against ENA's MD5.

## Prior art, and what this is not

Counting viral reads in single-cell data is not an open problem.
[scPathoQuant](https://academic.oup.com/bioinformatics/article/40/4/btae145/7628127)
(Bioinformatics 2024) takes reads that failed to align to the host genome
and maps them to a pathogen genome, adding the counts to standard
single-cell matrices.
[Viral-Track](https://www.sciencedirect.com/science/article/pii/S0092867420305687)
(Bost et al., Cell 2020) does a related job. If you want a maintained,
published tool, use one of those.

This pipeline takes a different route — one concatenated reference, so
host-vs-virus assignment happens in a single alignment pass rather than
downstream of a host-only alignment that already rejected the reads. It
also handles whitelist-free chemistries (Drop-seq, Seq-Well), which the
10x-oriented tools above do not advertise.

The part worth borrowing is not the alignment. It is the validation layer.

## Why this is not a thin wrapper around STAR

Every failure below produces a **successful run with a plausible matrix**.
None of them raises an error. Each is checked explicitly, and each check
stops the pipeline rather than printing a warning nobody reads.

| failure | what you see | what it actually is | caught by |
|---|---|---|---|
| GTF contig name disagrees with FASTA header | `viral_frac == 0` in every cell | STAR dropped every viral annotation; reads became intergenic | `VALIDATE_REFERENCE` |
| viral GTF has CDS rows but no `exon` rows | zero viral features in the matrix | STARsolo counts to exons | `VALIDATE_REFERENCE` |
| viral `gene_name` collides with a host symbol | one gene silently renamed `LAT-1` | HSV-1 carries **LAT** and **RS1**; so does the human genome. `var_names_make_unique()` renames by *position*, so which one keeps the bare symbol depends on file order | `PREPARE_VIRAL_GTF` prefixes, `VALIDATE_REFERENCE` checks |
| viral annotations lost during indexing | `viral_frac == 0` | any of the above, after concatenation | `STAR_INDEX` counts viral genes in `geneInfo.tab` |
| viral genes matched by string prefix | no infection detected, or host genes counted as viral | on an unprefixed reference `startswith("HSV1-")` finds nothing, and `startswith("LAT")` finds the host gene | `viral_genes.txt` + `--strict` in `SOLO_TO_H5AD` |
| wrong 10x chemistry assumed | a matrix with a believable cell count | read length alone cannot separate 3' v2 from 5' v1 | chemistry is a required parameter; `detect_chemistry.py` for 10x |
| chemistry assumed to be 10x at all | every barcode invalid, or an empty matrix | GSE163952 is **Drop-seq** (R1 = 20 bp = CB12+UMI8), not 10x, despite being described as 10x Chromium in several places | `BARCODE_STRUCTURE` |
| `--soloCBmatchWLtype` passed without a whitelist | STAR errors, or silently misbehaves | those options are only valid alongside `--soloCBwhitelist <file>` | `STARSOLO` builds args from `--whitelist` |
| R1 and R2 swapped | ~0 valid barcodes | STARsolo takes **cDNA first**, barcode read second | `CHECK_READS` |
| barcode read shorter than CB + UMI | garbage barcodes | trimmed R1, or the wrong chemistry | `CHECK_READS` |
| overlapping viral ORFs | viral counts below true depth | STARsolo assigns a read only when exactly one gene matches; HSV-1 is gene-dense and overlapping | `VALIDATE_REFERENCE` reports overlap; `Gene` and `GeneFull` both produced |

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
biggest drop. K should be near the expected number of cells
(`--barcode_top_k`); a K larger than the real cell count lets empty
barcodes into the top set and blurs the plateau.

## Design decisions worth stating

**Viral features come from a file, not a prefix.** `prepare_viral_gtf.py`
writes `viral_genes.txt` from the GTF itself, and `solo_to_h5ad.py` reads
it. If a listed feature is missing from the matrix the pipeline fails
rather than scoring on a subset — a mismatch means the reference used for
alignment is not the one the list came from. This is also why
`--star_index` refuses to run without `--viral_genes`.

**GFP is excluded from `viral_frac`.** It is driven by an artificial
promoter and does not track viral transcription. It stays in the matrix
and is reported separately as `gfp_counts`.

**Both `Gene` and `GeneFull` are produced.** HSV-1's overlapping ORFs mean
a real share of viral reads are multi-gene and dropped under `Gene`.
`validate_reference.py` reports the overlap fraction for the specific
reference being built, so the choice is made against a number. Compare
`viral_frac` only within one feature set: GeneFull adds intronic reads to
the denominator ([why it matters](RESULTS.md#gene-vs-genefull)).

**Chemistry is a parameter, not an inference.** Auto-detection exists, but
the pipeline stops after it. A detected value used silently is a value
nobody reviewed.

**No BAM by default.** STARsolo does not need one, and a sorted BAM per
library is the largest and slowest output. `--save_bam` turns it on, with
the open-file-limit handling that coordinate sorting at high thread counts
needs.

## Layout

```
main.nf                   workflow: reference -> per-sample alignment -> summary
nextflow.config           parameters, resources, profiles
conf/
  test.config             synthetic end-to-end run (CI)
  gse163952.config        the validation dataset
modules/local/            one process per file
bin/                      Python on the pipeline's PATH (and in the container)
  prepare_viral_gtf.py    rename contig, prefix gene names, synthesise gene rows
  validate_reference.py   pre-flight checks before genomeGenerate
  barcode_structure.py    infer CB/UMI length with no whitelist, from read collapse
  detect_chemistry.py     infer 10x chemistry by whitelist match; refuses to guess
  solo_to_h5ad.py         matrix -> .h5ad with per-cell viral_frac
  collect_summary.py      per-run summary table
scripts/                  data download, and the original shell pipeline
Dockerfile                STAR 2.7.11b + Python, from environment.yml
environment.yml           pinned software, shared by Docker and -profile conda
tests/
  run_tests.sh            unit tests: one synthetic fixture per failure mode
  make_pipeline_testdata.py  regenerates tests/data/ (deterministic)
  check_test_results.py   asserts the test profile's outcome
```

### Without Nextflow

The original shell entry points still work if you have STAR and the Python
dependencies installed (`micromamba env create -f environment.yml`):
`scripts/build_reference.sh` then `scripts/run_starsolo.sh <sample> <R1> <R2>`,
both configured through `config/config.sh`. The Nextflow pipeline runs the
same scripts with the same checks, and adds the read-structure check and the
`--strict` viral-feature check.

## Development

```bash
bash tests/run_tests.sh                               # unit tests, seconds, no STAR
nextflow run . -profile test,docker                   # end to end, ~2 min
python3 tests/check_test_results.py results
```

CI runs both on every push and pull request, then on `main` and on `v*`
tags publishes the container to `ghcr.io/ckangas2/aligninator`. The
pipeline pins the image by version in `nextflow.config`, so a release is:
bump `manifest.version` and `process.container` together, then tag
`v<version>`.

## Status

Tested at three levels: 17 unit tests that plant each failure mode, an
end-to-end synthetic run with asserted outcomes and negative cases in CI,
and the GSE163952 validation run with STAR 2.7.11b
([RESULTS.md](RESULTS.md)).

Known limitations:

- `fetch_reference.sh` cannot download the 10x barcode whitelists; the URL
  it tries is wrong. They ship inside CellRanger. This does not affect
  whitelist-free chemistries, which need no list.
- Reconciliation against the deposited Seurat object is incomplete: that
  object holds 7,226 cells across 22 sequencing runs, so it is a subset,
  and which runs it covers has not been established.
- Only Drop-seq and 10x v2/v3 are in the chemistry table.
- One host and one virus per run. Multiple pathogens would need
  `viral_genes.txt` to carry a source column.

# Corrections to project notes

## GSE163952 is Drop-seq, not 10x Chromium  (2026-10-17)

`DATASETS.md` records GSE163952 as "10x Chromium". It is not.
`HANDOFF.md` §3 already said "Drop-seq despite what the search said" —
the two files disagreed and the reads settle it.

Evidence, from `SRR21855376_1.fastq.gz`:

- R1 is **20 bp**, uniform across 100,000 reads. No 10x chemistry has a
  20 bp barcode read: v2 is 26 (CB16+UMI10), v3 is 28 (CB16+UMI12).
- Barcode collapse: 40,653 distinct 12-mers (2.46 reads each) against
  92,639 distinct 16-mers (1.08 reads each) in the same reads. Bases
  13–16 add essentially no collapsing, so they are UMI.
- The plateau sweep in `barcode_structure.py` puts the barcode boundary
  at 12 bp.

20 bp = 12 bp cell barcode + 8 bp UMI, which is Drop-seq / Seq-Well.
Consistent with the lab: the Rajewsky/Landthaler group uses Drop-seq,
and Wyler (GSE123782, also Drop-seq) is a co-author.

R2 is **pre-trimmed**, 53–63 bp with mode 63, so `SJDB_OVERHANG=62`.

### What this changes

- No barcode whitelist exists. `--soloCBwhitelist None`, and
  `--soloCBmatchWLtype` / the `*_CR` UMI options must NOT be passed.
- `DATASETS.md` calls GSE163952 primary partly because it is "10x
  Chromium" and therefore deeper than the Drop-seq baseline. That
  rationale does not hold. The *observed* numbers it reports from the
  deposited object (NECTIN1 11.9%, median 2,007 UMI/cell) stand on their
  own and are still the reason to prefer it — but the mechanism is not
  platform depth.

### Still open

The deposited object has 7,226 cells across 22 sequencing runs — about
330 cells per library, implausibly low. `GSE163952_hg38HSV1GFPtc.Seurat.500.RDS.gz`
is therefore a subset, not the whole deposit. Which runs it covers must
be resolved (`table(obj$orig.ident)`) before any reconciliation claim.

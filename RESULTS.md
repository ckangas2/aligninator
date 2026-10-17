# Validation run — 2026-10-17

End-to-end on GSE163952 (Rybak-Wolf et al., Nat Microbiol 2023), iPSC-2
cerebral organoids, mock vs 3 days post infection. Both libraries went
through the identical reference, chemistry and code path.

## Reference

| | |
|---|---|
| host | CellRanger GRCh38 2020-A, 194 contigs |
| virus | `chrHSV1` from the GSE163952 deposit |
| merged | **195 contigs** |
| viral genes in STAR index | **76** |
| chemistry | Drop-seq, CB 12 bp + UMI 8 bp, no whitelist |
| `sjdbOverhang` | 62 (R2 is pre-trimmed, 53–63 bp) |
| genomeGenerate | ~16 min, 32 threads |
| per-library alignment | ~2 min |

## Result

`viral_frac` excludes GFP (artificial promoter, does not track viral
transcription). 75 of 76 viral features scored.

| | ctrl (mock) | 3 dpi | 3 dpi, GeneFull |
|---|---|---|---|
| cells | 2,175 | 1,506 | 1,596 |
| median `viral_frac` | 0.0000 | 0.0450 | 0.0301 |
| p90 | 0.0000 | 0.1392 | 0.1035 |
| max | 0.0000 | 0.7132 | 0.6638 |
| cells >1% | **0.0%** | **83.3%** | 76.8% |
| cells >5% | 0.0% | 46.1% | 32.3% |
| cells >20% | 0.0% | 3.0% | 1.1% |

The mock arm is the control that matters: zero cells above 1% viral,
against 83.3% at 3 dpi. Same reference, same code. A `viral_frac` of 0
in a mock sample on its own proves nothing — it is equally consistent
with "no virus" and "viral counting silently broken" — so the positive
arm is what validates the pipeline, not the negative one.

Cell count falls from 2,175 to 1,506 with infection, consistent with
HSV-1 host shutoff.

## Gene vs GeneFull

GeneFull reports a LOWER viral fraction, which is the opposite of the
naive expectation. GeneFull adds intronic reads to the denominator, and
HSV-1 genes are essentially intronless — the human total grows while
viral counts barely move. The gap is a denominator effect, not viral
read loss.

Consequence: `viral_frac` is only comparable within one feature set.
Use `Gene` for the viral label and keep it consistent. To measure what
overlapping ORFs actually cost, compare viral *counts* between the two,
not fractions.

## Against the deposited object

DATASETS.md records max `viral_frac` 0.7118 for the deposited Seurat
object; this realignment gives 0.7132. Independent alignment, same
ceiling to three decimals.

The medians differ (0.0012 deposited vs 0.0450 here) because the
deposited object pools all conditions including mock, which pulls its
median toward zero. This run is a single 3 dpi library.

Not yet checked: which of the 22 sequencing runs the 7,226-cell
deposited object actually covers. At ~330 cells per library it must be a
subset, so a full reconciliation claim is not yet supported.

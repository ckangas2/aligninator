#!/usr/bin/env python3
"""Assert the outcome of `nextflow run . -profile test`.

A run that finishes is not a run that worked. These are the properties that
fail silently when the pipeline is broken:

  * the mock library must come back with no cell above 1% viral
  * the infected library must come back clearly positive (~70% planted)
  * every viral feature in viral_genes.txt must be present in the matrix
  * the host's LAT and the viral HSV1-LAT must both survive as separate features
"""
import gzip, sys
from pathlib import Path
import pandas as pd

def main(outdir="results"):
    out = Path(outdir)
    fails = []
    def check(cond, msg):
        print(f"  [{'ok' if cond else 'FAIL':>4}] {msg}")
        if not cond:
            fails.append(msg)

    df = pd.read_csv(out / "summary.tsv", sep="\t")
    print(df.to_string(index=False), "\n")
    gene = df[df.feature_set == "Gene"].set_index("sample")

    check(set(gene.index) == {"mock", "infected"}, "both samples summarised")
    check((df.viral_features_missing == 0).all(), "no viral feature missing from any matrix")
    check((df.viral_features_found == 13).all(), "all 13 viral features found")
    check((df.viral_features_scored == 12).all(), "GFP excluded from viral_frac (12 scored)")
    for s in ("mock", "infected"):
        check(80 <= gene.loc[s, "cells"] <= 110, f"{s}: cell count {gene.loc[s, 'cells']} in [80, 110] (100 planted)")
    check(gene.loc["mock", "pct_cells_gt_1pct"] == 0.0, "mock: 0% of cells above 1% viral")
    check(gene.loc["mock", "max_viral_frac"] == 0.0, "mock: max viral_frac is 0")
    inf = gene.loc["infected", "pct_cells_gt_1pct"]
    check(55.0 <= inf <= 85.0, f"infected: {inf}% of cells above 1% viral (70% planted)")

    feats = out / "samples" / "infected" / "Solo.out" / "Gene" / "filtered" / "features.tsv"
    opener = gzip.open if not feats.exists() else open
    if not feats.exists():
        feats = feats.with_suffix(".tsv.gz")
    names = {l.split("\t")[1] for l in opener(feats, "rt").read().splitlines() if l}
    check({"LAT", "HSV1-LAT", "RS1", "HSV1-RS1"} <= names,
          "host LAT/RS1 and viral HSV1-LAT/HSV1-RS1 are distinct features")

    idx_log = (out / "reference" / "star_index_check.log").read_text()
    check("viral genes in index    : 13" in idx_log, "13 viral genes in the STAR index")

    print(f"\n{'PASSED' if not fails else f'FAILED ({len(fails)})'}")
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))

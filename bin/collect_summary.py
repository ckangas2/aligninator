#!/usr/bin/env python3
"""Join per-sample viral metrics with the STARsolo QC numbers that matter.

One row per (sample, feature set). The STAR columns come from the first
feature set's Summary.csv and repeat across that sample's rows.
"""
import argparse, sys
from pathlib import Path
import pandas as pd

STAR_KEYS = {
    "Number of Reads": "reads",
    "Reads With Valid Barcodes": "valid_barcodes",
    "Sequencing Saturation": "saturation",
    "Reads Mapped to Genome: Unique": "mapped_genome_unique",
    "Reads Mapped to Gene: Unique Gene": "mapped_gene_unique",
    "Estimated Number of Cells": "star_estimated_cells",
    "Median UMI per Cell": "median_umi_per_cell",
}

def read_star(path):
    out = {}
    for line in Path(path).read_text().splitlines():
        k, _, v = line.partition(",")
        if k in STAR_KEYS:
            out[STAR_KEYS[k]] = v
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--viral", nargs="+", required=True)
    ap.add_argument("--star", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    viral = pd.concat([pd.read_csv(p, sep="\t") for p in a.viral], ignore_index=True)
    star = pd.DataFrame([
        {"sample": Path(p).name.removesuffix(".star_summary.csv"), **read_star(p)}
        for p in a.star])
    df = viral.merge(star, on="sample", how="left")
    df = df.sort_values(["sample", "feature_set"]).reset_index(drop=True)
    df.to_csv(a.out, sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(df[["sample", "feature_set", "cells", "valid_barcodes",
                  "median_viral_frac", "pct_cells_gt_1pct", "pct_cells_gt_5pct"]]
              .to_string(index=False))
    return 0

if __name__ == "__main__":
    sys.exit(main())

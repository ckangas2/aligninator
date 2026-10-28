#!/usr/bin/env python3
"""
Convert a STARsolo raw/filtered matrix into .h5ad with a per-cell viral load.

Design note -- why viral genes come from a FILE, not a prefix
-------------------------------------------------------------
The obvious implementation is `adata.var_names.str.startswith("HSV1-")`.
That works only if the reference was built with prefixing. On a reference
built from the deposit as-is, viral features are bare UL1/RS1/LAT, and
prefix matching finds nothing -- viral_frac comes out 0.0 for every cell
and the run looks like a failed infection. Worse, a prefix like "LAT"
would capture host genes.

So the viral feature list is an explicit input, produced by
prepare_viral_gtf.py from the GTF itself. If a listed feature is absent
from the matrix we say so loudly rather than quietly scoring on a subset.

GFP is excluded from viral_frac by default: it is driven by an artificial
promoter and does not track viral transcription. It is retained in the
matrix and reported separately as a correlation check.
"""
import argparse, sys
from pathlib import Path
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
import pandas as pd

def read_solo_dir(d):
    d = Path(d)
    mtx = d / "matrix.mtx"
    if not mtx.exists(): mtx = d / "matrix.mtx.gz"
    X = sio.mmread(str(mtx)).tocsr()          # STARsolo: genes x cells
    bc = pd.read_csv(d / ("barcodes.tsv" if (d/"barcodes.tsv").exists()
                          else "barcodes.tsv.gz"), header=None, sep="\t")
    ft = pd.read_csv(d / ("features.tsv" if (d/"features.tsv").exists()
                          else "features.tsv.gz"), header=None, sep="\t")
    return X, bc, ft

def compute_viral(X_cells_by_genes, var_names, viral_names, exclude=("GFP",)):
    """
    X: CSR, cells x genes. Returns dict of per-cell vectors + a report.
    Pure scipy so it can be unit-tested without anndata.
    """
    var_index = {g: i for i, g in enumerate(var_names)}
    requested = list(viral_names)
    found  = [g for g in requested if g in var_index]
    missing = [g for g in requested if g not in var_index]

    excl = {g for g in found if any(g == e or g.endswith("-" + e) or
                                    g.endswith(e) and e in ("GFP",) for e in exclude)}
    scored = [g for g in found if g not in excl]

    total = np.asarray(X_cells_by_genes.sum(axis=1)).ravel()
    def colsum(names):
        if not names:
            return np.zeros(X_cells_by_genes.shape[0])
        idx = [var_index[g] for g in names]
        return np.asarray(X_cells_by_genes[:, idx].sum(axis=1)).ravel()

    viral = colsum(scored)
    gfp   = colsum(sorted(excl))
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(total > 0, viral / np.maximum(total, 1), 0.0)
    report = dict(requested=len(requested), found=len(found),
                  missing=missing, excluded=sorted(excl), scored=len(scored))
    return dict(total_counts=total, viral_counts=viral, gfp_counts=gfp,
                viral_frac=frac), report

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solo-dir", required=True)
    ap.add_argument("--viral-genes-file", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default="")
    ap.add_argument("--exclude-from-viral", default="GFP")
    ap.add_argument("--feature-name-col", type=int, default=1,
                    help="column in features.tsv to use as var_names (0=id,1=name)")
    ap.add_argument("--feature-set", default="",
                    help="label for the summary row, e.g. Gene or GeneFull")
    ap.add_argument("--summary-out", default="",
                    help="write a one-row TSV of per-sample viral metrics here")
    ap.add_argument("--strict", action="store_true",
                    help="exit 3 if any listed viral feature is absent from the matrix")
    args = ap.parse_args()

    X, bc, ft = read_solo_dir(args.solo_dir)
    X = X.T.tocsr()                                   # -> cells x genes
    col = min(args.feature_name_col, ft.shape[1]-1)
    var_names = ft[col].astype(str).tolist()
    viral = [l.strip() for l in Path(args.viral_genes_file).read_text().splitlines() if l.strip()]

    vals, rep = compute_viral(X, var_names, viral,
                              exclude=tuple(s for s in args.exclude_from_viral.split(",") if s))

    print(f"matrix            : {X.shape[0]:,} cells x {X.shape[1]:,} genes")
    print(f"viral requested   : {rep['requested']}")
    print(f"viral found       : {rep['found']}")
    print(f"viral scored      : {rep['scored']}  (excluded: {rep['excluded']})")
    if rep["missing"]:
        print(f"  MISSING {len(rep['missing'])} requested viral features from the matrix:")
        print(f"  {rep['missing'][:20]}")
        print("  A viral gene list that does not match the matrix means the")
        print("  reference used for alignment is not the one this list came from.")
    frac = vals["viral_frac"]
    print(f"\nviral_frac        : median {np.median(frac):.4f}  "
          f"p90 {np.percentile(frac,90):.4f}  max {frac.max():.4f}")
    for t in (0.01, 0.05, 0.20):
        print(f"  cells > {t:>4.0%}      : {(frac>t).sum():,} ({(frac>t).mean():.1%})")

    if args.summary_out:
        n = len(frac)
        row = dict(
            sample=args.sample or Path(args.solo_dir).name,
            feature_set=args.feature_set,
            cells=n,
            viral_features_listed=rep["requested"],
            viral_features_found=rep["found"],
            viral_features_scored=rep["scored"],
            viral_features_missing=len(rep["missing"]),
            median_viral_frac=round(float(np.median(frac)), 4) if n else 0.0,
            p90_viral_frac=round(float(np.percentile(frac, 90)), 4) if n else 0.0,
            max_viral_frac=round(float(frac.max()), 4) if n else 0.0,
            pct_cells_gt_1pct=round(100 * float((frac > 0.01).mean()), 1) if n else 0.0,
            pct_cells_gt_5pct=round(100 * float((frac > 0.05).mean()), 1) if n else 0.0,
            pct_cells_gt_20pct=round(100 * float((frac > 0.20).mean()), 1) if n else 0.0,
        )
        pd.DataFrame([row]).to_csv(args.summary_out, sep="\t", index=False)

    if args.strict and rep["missing"]:
        print("\nFATAL (--strict): viral gene list does not match the matrix.",
              file=sys.stderr)
        return 3

    try:
        import anndata as ad
    except ImportError:
        out = Path(args.out).with_suffix(".csv")
        pd.DataFrame({k: v for k, v in vals.items()},
                     index=bc[0].astype(str)).to_csv(out)
        print(f"\nanndata unavailable; wrote per-cell metrics to {out}")
        return 0

    adata = ad.AnnData(X=X)
    adata.obs_names = bc[0].astype(str).values
    adata.var_names = pd.Index(var_names)
    adata.var_names_make_unique()
    for k, v in vals.items():
        adata.obs[k] = v
    if args.sample:
        adata.obs["sample"] = args.sample
    adata.uns["viral_features"] = rep
    adata.write_h5ad(args.out)
    print(f"\nwrote {args.out}")
    return 0

if __name__ == "__main__":
    sys.exit(main())

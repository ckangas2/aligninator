#!/usr/bin/env python3
"""
Infer cell-barcode length from read structure alone -- no whitelist needed.

Why this is necessary
---------------------
detect_chemistry.py identifies 10x libraries by matching the barcode
prefix against a published whitelist. Whitelist-free chemistries
(Drop-seq, Seq-Well, inDrop) have no such list, so that method reports
"no confident call" and stops -- correct, but not useful.

The signal that does work: a real cell barcode COLLAPSES reads onto
itself. A handful of beads capture most of the library, so the number
of distinct prefixes stays low and the top prefixes hold most of the
reads. Extend the prefix one base past the end of the barcode and you
start reading UMI, which is random -- every barcode shatters into ~4
pieces and the read mass in the top prefixes drops hard.

So: sweep prefix length, measure the share of reads held by the top K
prefixes, and find the length after which that share falls off a cliff.
The cliff is the end of the barcode.

Worked example (GSE163952, 100k reads):
    12 bp -> 40,653 distinct   (2.46 reads each)
    16 bp -> 92,639 distinct   (1.08 reads each)
Bases 13-16 add almost no collapsing, so they are UMI. CB is 12.
"""
import argparse, gzip, sys
from collections import Counter

def open_maybe_gz(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

def read_seqs(path, n):
    out = []
    with open_maybe_gz(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                out.append(line.strip())
                if len(out) >= n:
                    break
    return out

# (name, cb_len, umi_len, r1_len, whitelist?)
KNOWN = [
    ("Drop-seq",      12,  8, 20, False),
    ("Seq-Well",      12,  8, 20, False),
    ("10x v2",        16, 10, 26, True),
    ("10x v3",        16, 12, 28, True),
    ("inDrop v3",     16,  6, 22, False),
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--n-reads", type=int, default=200000)
    ap.add_argument("--top-k", type=int, default=1000,
                    help="prefixes counted as 'top'; ~expected cell count")
    ap.add_argument("--min-len", type=int, default=8)
    ap.add_argument("--max-len", type=int, default=0)
    args = ap.parse_args()

    seqs = read_seqs(args.r1, args.n_reads)
    if not seqs:
        sys.exit(f"ERROR: no reads in {args.r1}")
    n = len(seqs)
    rlen = Counter(len(s) for s in seqs).most_common(1)[0][0]
    hi = args.max_len or rlen
    print(f"reads sampled : {n:,}")
    print(f"R1 length     : {rlen} bp")
    print(f"top-K         : {args.top_k}\n")

    print(f"{'len':>4}  {'distinct':>10}  {'reads/prefix':>13}  {'reads in topK':>14}  {'drop':>7}")
    print("-" * 60)
    stats = {}
    prev_f = None
    for L in range(args.min_len, hi + 1):
        c = Counter(s[:L] for s in seqs)
        top = sum(v for _, v in c.most_common(args.top_k))
        f = top / n
        drop = (f / prev_f) if prev_f else None
        stats[L] = dict(distinct=len(c), frac_top=f, ratio=drop)
        ds = f"{drop:6.3f}" if drop is not None else "     -"
        print(f"{L:>4}  {len(c):>10,}  {n/len(c):>13.2f}  {f:>13.1%}  {ds:>7}")
        prev_f = f

    # Within the barcode, adding a base changes nothing -- the same cells
    # still hold the same reads, so the top-K share sits on a plateau
    # (ratio ~1.0). The first base of UMI shatters each barcode and the
    # share falls. The barcode is the LAST length on the plateau.
    #
    # Taking the steepest fall anywhere is wrong: the share keeps
    # decaying past the barcode and the global minimum lands near the
    # read end.
    PLATEAU = 0.90
    cb = None
    for L in range(args.min_len, hi):
        r = stats.get(L + 1, {}).get("ratio")
        if r is not None and r < PLATEAU:
            cb = L
            break
    if cb is None:
        print("\nVERDICT: NO CLEAR BARCODE BOUNDARY.")
        print("  The top-K read share never falls -- no prefix length shows")
        print("  the collapse a real barcode produces. Possible causes: reads")
        print("  are not barcoded, R1/R2 swapped, a variable-length barcode")
        print("  (inDrop), or too few reads sampled.")
        return 2
    if cb == args.min_len:
        print(f"\n  WARNING: the plateau breaks at the first length tested "
              f"({args.min_len}).")
        print(f"  The true barcode may be shorter. Re-run with --min-len 4.")
    print(f"\nplateau ends at {cb} (top-K share falls to "
          f"{stats[cb+1]['ratio']:.1%} at {cb+1})")

    # cross-check: distinct prefixes should also stop growing at the barcode
    g_in  = stats[cb]["distinct"] / stats[cb-1]["distinct"] if cb-1 in stats else None
    g_out = stats[cb+1]["distinct"] / stats[cb]["distinct"]
    if g_in is not None:
        print(f"  cross-check: distinct prefixes grow {g_in:.3f}x into {cb}, "
              f"{g_out:.3f}x past it")
        if g_in > 1.25:
            print("  WARNING: distinct count was still climbing at the called "
                  "length; the two signals disagree. Inspect the table.")

    umi = rlen - cb
    print(f"\nVERDICT: cell barcode = {cb} bp, UMI = {umi} bp  (R1 {rlen} bp)")
    match = [k for k in KNOWN if k[1] == cb and k[2] == umi and k[3] == rlen]
    if match:
        names = " / ".join(m[0] for m in match)
        needs_wl = any(m[4] for m in match)
        print(f"  consistent with: {names}")
        print(f"  whitelist: {'required' if needs_wl else 'NONE -- use --soloCBwhitelist None'}")
    else:
        print("  does not match a chemistry in the table; verify before aligning.")

    print(f"\nSTARsolo arguments:")
    print(f"  --soloType CB_UMI_Simple \\")
    print(f"  --soloCBstart 1 --soloCBlen {cb} \\")
    print(f"  --soloUMIstart {cb+1} --soloUMIlen {umi} \\")
    if not (match and any(m[4] for m in match)):
        print(f"  --soloCBwhitelist None")
        print(f"  # NOTE: --soloCBmatchWLtype must be left at default (Exact)")
        print(f"  # when there is no whitelist. Passing 1MM_multi_* is an error.")
    return 0

if __name__ == "__main__":
    sys.exit(main())

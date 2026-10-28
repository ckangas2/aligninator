#!/usr/bin/env python3
"""
Infer 10x chemistry from R1 reads by testing barcode prefixes against
candidate whitelists.

Why this exists
---------------
STARsolo needs --soloCBlen/--soloUMIlen/--soloCBwhitelist to be *right*.
Guess wrong and you do not get an error: you get a run that completes and
produces a matrix with a plausible-looking number of cells and badly wrong
counts. Read length alone cannot separate 3' v2 from 5' v1 (both 26 bp), so
we also measure the whitelist hit rate, which does separate v2 from v3.

The output is a decision plus the evidence for it. Never just the decision.
"""
import argparse, gzip, json, sys
from pathlib import Path
from collections import Counter

# (name, cb_len, umi_len, whitelist basename, soloStrand, soloType)
CHEMISTRIES = [
    ("10xv2",        16, 10, "737K-august-2016.txt", "Forward", "CB_UMI_Simple"),
    ("10xv3",        16, 12, "3M-february-2018.txt", "Forward", "CB_UMI_Simple"),
    ("10x5prime_v1", 16, 10, "737K-august-2016.txt", "Reverse", "CB_UMI_Simple"),
    ("10x5prime_v2", 16, 10, "737K-august-2016.txt", "Reverse", "CB_UMI_Simple"),
    ("10x_multiome", 16, 12, "737K-arc-v1.txt",      "Forward", "CB_UMI_Simple"),
]

def open_maybe_gz(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

def read_r1(path, n):
    """Yield up to n sequence lines from a FASTQ."""
    seqs = []
    with open_maybe_gz(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
                if len(seqs) >= n:
                    break
    return seqs

def load_whitelist(path):
    wl = set()
    with open_maybe_gz(path) as fh:
        for line in fh:
            s = line.strip()
            if s:
                wl.add(s.split("\t")[0])
    return wl

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True, help="R1 fastq(.gz)")
    ap.add_argument("--whitelist-dir", required=True)
    ap.add_argument("--n-reads", type=int, default=100000)
    ap.add_argument("--min-match", type=float, default=0.50,
                    help="minimum whitelist hit rate to accept a call")
    ap.add_argument("--json-out", default="")
    args = ap.parse_args()

    seqs = read_r1(args.r1, args.n_reads)
    if not seqs:
        sys.exit(f"ERROR: no reads parsed from {args.r1}")

    lens = Counter(len(s) for s in seqs)
    modal_len, modal_n = lens.most_common(1)[0]
    frac_modal = modal_n / len(seqs)

    print(f"reads sampled      : {len(seqs):,}")
    print(f"R1 length (modal)  : {modal_len} bp  ({frac_modal:.1%} of reads)")
    if len(lens) > 1:
        print(f"  other lengths    : {dict(lens.most_common()[1:6])}")
    if frac_modal < 0.90:
        print("  WARNING: R1 lengths are heterogeneous. Reads may be trimmed.")
    print()

    wl_dir = Path(args.whitelist_dir)
    wl_cache = {}
    rows = []
    for name, cb, umi, wl_file, strand, solotype in CHEMISTRIES:
        need = cb + umi
        if modal_len < need:
            rows.append(dict(chemistry=name, cb_len=cb, umi_len=umi,
                             whitelist=wl_file, strand=strand, solo_type=solotype,
                             match_rate=None,
                             note=f"R1 too short ({modal_len} < {need})"))
            continue
        wl_path = wl_dir / wl_file
        if not wl_path.exists() and not (wl_dir / (wl_file + ".gz")).exists():
            rows.append(dict(chemistry=name, cb_len=cb, umi_len=umi,
                             whitelist=wl_file, strand=strand, solo_type=solotype,
                             match_rate=None, note="whitelist file not found"))
            continue
        if wl_file not in wl_cache:
            real = wl_path if wl_path.exists() else wl_dir / (wl_file + ".gz")
            wl_cache[wl_file] = load_whitelist(real)
        wl = wl_cache[wl_file]
        hits = sum(1 for s in seqs if s[:cb] in wl)
        rows.append(dict(chemistry=name, cb_len=cb, umi_len=umi,
                         whitelist=wl_file, strand=strand, solo_type=solotype,
                         match_rate=hits / len(seqs), note=""))

    print(f"{'chemistry':<16}{'CB':>4}{'UMI':>5}  {'match':>8}  note")
    print("-" * 70)
    for r in sorted(rows, key=lambda r: (r["match_rate"] is None, -(r["match_rate"] or 0))):
        mr = "    n/a" if r["match_rate"] is None else f"{r['match_rate']:7.2%}"
        print(f"{r['chemistry']:<16}{r['cb_len']:>4}{r['umi_len']:>5}  {mr}  {r['note']}")
    print()

    scored = [r for r in rows if r["match_rate"] is not None]
    best = max(scored, key=lambda r: r["match_rate"]) if scored else None

    if best is None or best["match_rate"] < args.min_match:
        got = f"{best['match_rate']:.2%}" if best else "none"
        print(f"VERDICT: NO CONFIDENT CALL (best {got} < {args.min_match:.0%})")
        print("  Do not proceed on a guess. Likely causes: non-10x library,")
        print("  R1 already trimmed, wrong file passed as R1, or a chemistry")
        print("  not in the table above.")
        result = dict(call=None, reason="below_threshold", modal_r1_len=modal_len,
                      candidates=rows)
        rc = 2
    else:
        # ambiguity check: anything within 5 points of the winner on a
        # *different* whitelist is a real tie, not a rounding artifact.
        rivals = [r for r in scored
                  if r["chemistry"] != best["chemistry"]
                  and r["match_rate"] > best["match_rate"] - 0.05
                  and (r["umi_len"] != best["umi_len"] or r["strand"] != best["strand"])]
        print(f"VERDICT: {best['chemistry']}  "
              f"(CB {best['cb_len']}, UMI {best['umi_len']}, "
              f"{best['match_rate']:.2%} on {best['whitelist']})")
        if rivals:
            print("  AMBIGUOUS — these differ in UMI length or strand and score close:")
            for r in rivals:
                print(f"    {r['chemistry']}  {r['match_rate']:.2%}  "
                      f"UMI {r['umi_len']}  strand {r['strand']}")
            print("  Read length cannot separate 3' from 5' chemistry. Confirm")
            print("  against the GEO/SRA library description before running.")
        result = dict(call=best["chemistry"], cb_len=best["cb_len"],
                      umi_len=best["umi_len"], whitelist=best["whitelist"],
                      strand=best["strand"], solo_type=best["solo_type"],
                      match_rate=best["match_rate"], modal_r1_len=modal_len,
                      ambiguous=[r["chemistry"] for r in rivals], candidates=rows)
        rc = 0

    if args.json_out:
        Path(args.json_out).write_text(json.dumps(result, indent=2))
        print(f"\nwrote {args.json_out}")
    return rc

if __name__ == "__main__":
    sys.exit(main())

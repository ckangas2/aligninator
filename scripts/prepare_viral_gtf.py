#!/usr/bin/env python3
"""
Normalise a deposited viral GTF so it can be safely concatenated with a
host reference.

Three transformations, each addressing a failure in validate_reference.py:

  --rename-contig OLD:NEW   fix contig-name disagreement with the FASTA
  --prefix HSV1-            prefix gene_id and gene_name so viral features
                            cannot collide with host symbols (LAT, RS1)
  --add-gene-rows           synthesise `gene` rows from exon extents, for
                            tools stricter than STAR (cellranger mkref)

Also emits viral_genes.txt -- the authoritative list of viral feature names
in the finished reference. Downstream code must read that file rather than
pattern-matching on a prefix: prefix matching on an unprefixed deposit
silently captures the host's LAT and RS1.
"""
import argparse, gzip, re, sys
from pathlib import Path
from collections import defaultdict

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')

def open_maybe_gz(p, mode="rt"):
    p = str(p)
    return gzip.open(p, mode) if p.endswith(".gz") else open(p, mode)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-gtf", required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--gene-list-out", default="")
    ap.add_argument("--rename-contig", default="", metavar="OLD:NEW")
    ap.add_argument("--prefix", default="HSV1-")
    ap.add_argument("--add-gene-rows", action="store_true")
    ap.add_argument("--copy-id-to-name", action="store_true",
                    help="if a row lacks gene_name, set it from gene_id")
    args = ap.parse_args()

    old = new = None
    if args.rename_contig:
        if ":" not in args.rename_contig:
            sys.exit("--rename-contig must be OLD:NEW")
        old, new = args.rename_contig.split(":", 1)

    n_in = n_renamed = n_prefixed = n_named = 0
    rows, extents = [], defaultdict(lambda: [None, None, None, None])  # gid -> [seq,start,end,strand]

    with open_maybe_gz(args.in_gtf) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            n_in += 1
            if old is not None and f[0] == old:
                f[0] = new; n_renamed += 1
            attrs = dict(ATTR_RE.findall(f[8]))

            if "gene_name" not in attrs and args.copy_id_to_name and "gene_id" in attrs:
                attrs["gene_name"] = attrs["gene_id"]; n_named += 1

            for key in ("gene_id", "gene_name"):
                if key in attrs and args.prefix and not attrs[key].startswith(args.prefix):
                    attrs[key] = args.prefix + attrs[key]
                    if key == "gene_id":
                        n_prefixed += 1
            if "transcript_id" in attrs and args.prefix \
               and not attrs["transcript_id"].startswith(args.prefix):
                attrs["transcript_id"] = args.prefix + attrs["transcript_id"]

            f[8] = " ".join(f'{k} "{v}";' for k, v in attrs.items())
            rows.append((f, attrs))

            gid = attrs.get("gene_id")
            if gid and f[2] in ("exon", "CDS"):
                e = extents[gid]
                s_, e_ = int(f[3]), int(f[4])
                e[0] = f[0]
                e[1] = s_ if e[1] is None else min(e[1], s_)
                e[2] = e_ if e[2] is None else max(e[2], e_)
                e[3] = f[6]

    have_gene = {a.get("gene_id") for f, a in rows if f[2] == "gene"}
    out_lines = []
    n_gene_rows = 0
    if args.add_gene_rows:
        for gid, (seq, s_, e_, st) in extents.items():
            if gid in have_gene or seq is None:
                continue
            name = next((a.get("gene_name", gid) for f, a in rows
                         if a.get("gene_id") == gid), gid)
            out_lines.append(f'{seq}\tprepare_viral_gtf\tgene\t{s_}\t{e_}\t.\t{st}\t.\t'
                             f'gene_id "{gid}"; gene_name "{name}";\n')
            n_gene_rows += 1

    for f, _ in rows:
        out_lines.append("\t".join(f) + "\n")

    with open_maybe_gz(args.out_gtf, "wt") as fh:
        fh.writelines(out_lines)

    names = sorted({a.get("gene_name") or a.get("gene_id")
                    for f, a in rows if f[2] in ("exon", "CDS", "gene")} - {None})
    if args.gene_list_out:
        Path(args.gene_list_out).write_text("\n".join(names) + "\n")

    print(f"input rows        : {n_in:,}")
    print(f"contig renamed    : {n_renamed:,} rows  ({old} -> {new})" if old else
          f"contig renamed    : (not requested)")
    print(f"gene_name set     : {n_named:,} rows from gene_id")
    print(f"prefixed          : {n_prefixed:,} features with '{args.prefix}'")
    print(f"gene rows added   : {n_gene_rows:,}")
    print(f"viral features    : {len(names):,}")
    print(f"wrote             : {args.out_gtf}")
    if args.gene_list_out:
        print(f"wrote             : {args.gene_list_out}")
    print(f"\nfirst 8 features  : {names[:8]}")

if __name__ == "__main__":
    sys.exit(main())

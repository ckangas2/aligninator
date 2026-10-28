#!/usr/bin/env python3
"""
Pre-flight checks on a concatenated host+virus reference, run BEFORE
STAR genomeGenerate (which takes hours and will not warn you about any
of this).

The failures below all share one property: STAR completes successfully
and hands you a matrix. Nothing errors. The counts are just wrong.

  1. Contig name mismatch between FASTA and GTF.
     STAR silently drops GTF rows whose sequence name is absent from the
     FASTA. A viral GTF on "HSV1" against a FASTA headed ">chrHSV1"
     yields a reference with the viral genome present but ZERO viral
     genes annotated. Every viral read becomes intergenic and vanishes
     from the gene matrix. Result: viral_frac == 0 everywhere, which
     reads as "no infection" rather than "broken reference".

  2. Missing exon rows.
     STARsolo counts to exons. A GTF with only CDS rows annotates
     nothing. GSE163952's deposited GTF has no `gene` rows at all --
     that is fine (STAR builds gene records from exon gene_id), but the
     absence of exon rows is fatal.

  3. gene_name collisions between host and virus.
     HSV-1 carries LAT and RS1. Human carries LAT (linker for activation
     of T cells) and RS1 (retinoschisin 1). Concatenate naively and the
     merged feature list has two genes named LAT. scanpy's
     var_names_make_unique() then renames one to "LAT-1" by position --
     silently, and the assignment depends on file order. Downstream,
     "LAT" may be the human gene or the viral one and nothing tells you
     which. Prefixing viral gene_names at build time is the only fix
     that survives every downstream tool.

  4. Same-strand overlapping features.
     HSV-1 is gene-dense with extensive overlap. STARsolo's default
     behaviour assigns a read only when it maps to exactly ONE gene;
     reads landing in an overlap are discarded. On a compact viral
     genome this can silently remove a large share of viral signal. We
     quantify it here so the loss is a known number, not a surprise.
"""
import argparse, gzip, re, sys
from pathlib import Path
from collections import defaultdict, Counter

def open_maybe_gz(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "rt")

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')

def parse_attrs(field):
    return dict(ATTR_RE.findall(field))

def fasta_contigs(path):
    names = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].strip().split()[0])
    return names

def parse_gtf(path):
    rows = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            rows.append(dict(seqname=f[0], feature=f[2], start=int(f[3]),
                             end=int(f[4]), strand=f[6], attrs=parse_attrs(f[8])))
    return rows

GENE_NAME_RE = re.compile(r'gene_name\s+"([^"]*)"')

def host_exon_gene_names(path):
    """gene_name values on exon rows, streamed.

    Only the set of names is needed for the collision check. Holding every
    row of a human GTF (~2.5M rows) as dicts takes many GB and gets the task
    OOM-killed; streaming keeps this to a few MB.
    """
    names = set()
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "exon":
                continue
            m = GENE_NAME_RE.search(f[8])
            if m:
                names.add(m.group(1))
    return names

def overlap_report(rows, contig, gene_key="gene_id"):
    """Fraction of a contig's exonic bases covered by >1 gene on the same strand."""
    exons = [r for r in rows if r["seqname"] == contig and r["feature"] == "exon"]
    out = {}
    for strand in ("+", "-"):
        iv = [(r["start"], r["end"], r["attrs"].get(gene_key, "?"))
              for r in exons if r["strand"] == strand]
        if not iv:
            out[strand] = dict(bases=0, multi=0, frac=0.0, genes=0)
            continue
        cov = defaultdict(set)
        for s, e, g in iv:
            for pos in range(s, e + 1):
                cov[pos].add(g)
        multi = sum(1 for p, gs in cov.items() if len(gs) > 1)
        out[strand] = dict(bases=len(cov), multi=multi,
                           frac=multi / len(cov) if cov else 0.0,
                           genes=len({g for _, _, g in iv}))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host-gtf"); ap.add_argument("--host-fasta")
    ap.add_argument("--viral-gtf", required=True)
    ap.add_argument("--viral-fasta", required=True)
    ap.add_argument("--viral-prefix", default="HSV1-",
                    help="prefix applied to viral gene_name (set '' to disable)")
    ap.add_argument("--overlap-contig", default="",
                    help="contig to run the overlap analysis on (default: first viral contig)")
    args = ap.parse_args()

    fail, warn = [], []
    def ok(m):   print(f"  [ ok ] {m}")
    def bad(m):  print(f"  [FAIL] {m}"); fail.append(m)
    def wrn(m):  print(f"  [warn] {m}"); warn.append(m)

    # ---------------------------------------------------------------- 1
    print("\n[1] contig names: viral FASTA vs viral GTF")
    vf = fasta_contigs(args.viral_fasta)
    vrows = parse_gtf(args.viral_gtf)
    vgtf_contigs = sorted({r["seqname"] for r in vrows})
    print(f"      FASTA contigs : {vf}")
    print(f"      GTF   contigs : {vgtf_contigs}")
    missing = [c for c in vgtf_contigs if c not in set(vf)]
    if missing:
        bad(f"GTF references contigs absent from FASTA: {missing} "
            f"-> these annotations are silently dropped by STAR")
    else:
        ok("every GTF contig is present in the FASTA")
    unann = [c for c in vf if c not in set(vgtf_contigs)]
    if unann:
        wrn(f"FASTA contigs with no annotation: {unann}")

    # ---------------------------------------------------------------- 2
    print("\n[2] feature types in viral GTF")
    feats = Counter(r["feature"] for r in vrows)
    print(f"      {dict(feats)}")
    if feats.get("exon", 0) == 0:
        bad("no `exon` rows -- STARsolo counts to exons and will annotate nothing")
    else:
        ok(f"{feats['exon']} exon rows present")
    if feats.get("gene", 0) == 0:
        wrn("no `gene` rows (OK for STAR: gene records are built from exon gene_id)")

    print("\n[3] required attributes on viral exon rows")
    exons = [r for r in vrows if r["feature"] == "exon"]
    for key in ("gene_id", "gene_name"):
        n = sum(1 for r in exons if key in r["attrs"])
        if n == 0:
            (bad if key == "gene_id" else wrn)(
                f"no exon row carries {key}" +
                ("" if key == "gene_id" else " -- features.tsv will fall back to gene_id"))
        elif n < len(exons):
            wrn(f"{len(exons)-n}/{len(exons)} exon rows lack {key}")
        else:
            ok(f"all {len(exons)} exon rows carry {key}")

    # ---------------------------------------------------------------- 3
    print("\n[4] host/virus gene_name collisions")
    vnames = {r["attrs"].get("gene_name") or r["attrs"].get("gene_id")
              for r in exons}
    vnames.discard(None)
    print(f"      viral gene names: {len(vnames)}")
    if args.host_gtf:
        hnames = host_exon_gene_names(args.host_gtf)
        hnames.discard(None)
        clash = sorted(vnames & hnames)
        if clash:
            bad(f"{len(clash)} viral gene_name(s) also exist in the host GTF: "
                f"{clash}")
            print("         These become duplicate rows in features.tsv.")
            print("         var_names_make_unique() will rename one by POSITION,")
            print("         so which gene ends up as the bare symbol depends on")
            print("         file concatenation order. Prefix the viral names.")
        else:
            ok("no host/virus gene_name collisions")
    else:
        known = {"LAT", "RS1", "UL16", "US1", "GFP"}
        risky = sorted(vnames & known)
        wrn("no --host-gtf given; cannot check collisions directly")
        if risky:
            print(f"         viral names that are known human symbols: {risky}")

    if args.viral_prefix:
        would = sorted(n for n in vnames if not n.startswith(args.viral_prefix))
        print(f"      with --viral-prefix '{args.viral_prefix}': "
              f"{len(would)} names would be renamed, collisions resolved")

    # ---------------------------------------------------------------- 4
    contig = args.overlap_contig or (vgtf_contigs[0] if vgtf_contigs else "")
    print(f"\n[5] same-strand exon overlap on {contig}")
    if contig:
        rep = overlap_report(vrows, contig)
        for strand, d in rep.items():
            if d["bases"]:
                print(f"      strand {strand}: {d['genes']:>3} genes, "
                      f"{d['bases']:>6,} exonic bp, "
                      f"{d['multi']:>6,} bp multi-gene ({d['frac']:.1%})")
        worst = max((d["frac"] for d in rep.values()), default=0.0)
        if worst > 0.20:
            wrn(f"{worst:.1%} of exonic bases are multi-gene on one strand. "
                f"STARsolo discards reads that cannot be assigned to exactly "
                f"one gene -- expect viral counts well below true depth.")
            print("         Mitigations: --soloFeatures Gene GeneFull (GeneFull")
            print("         counts intronic/ambiguous reads too), or score viral")
            print("         load from the BAM by contig rather than from the")
            print("         gene matrix. Decide deliberately; do not default.")
        else:
            ok(f"max same-strand overlap {worst:.1%}")

    # ---------------------------------------------------------------- done
    print("\n" + "=" * 62)
    print(f"RESULT: {len(fail)} failure(s), {len(warn)} warning(s)")
    if fail:
        print("\nDo not build the reference until the failures are fixed:")
        for m in fail:
            print(f"  - {m}")
    print("=" * 62)
    return 1 if fail else 0

if __name__ == "__main__":
    sys.exit(main())

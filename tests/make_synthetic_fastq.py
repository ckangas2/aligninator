#!/usr/bin/env python3
"""Build synthetic 10x FASTQ + whitelists to exercise detect_chemistry.py."""
import gzip, random, sys
from pathlib import Path

random.seed(0)
BASES = "ACGT"

def rand_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))

def make_whitelist(path, n, seed):
    rng = random.Random(seed)
    bcs = set()
    while len(bcs) < n:
        bcs.add("".join(rng.choice(BASES) for _ in range(16)))
    bcs = sorted(bcs)
    with open(path, "w") as fh:
        fh.write("\n".join(bcs) + "\n")
    return bcs

def write_fastq(path, reads):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt") as fh:
        for i, s in enumerate(reads):
            fh.write(f"@read{i}\n{s}\n+\n{'I'*len(s)}\n")

def main(outdir):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    wl = out / "whitelists"; wl.mkdir(exist_ok=True)

    # Two disjoint whitelists, mimicking the real v2/v3 lists.
    v2 = make_whitelist(wl / "737K-august-2016.txt", 5000, seed=1)
    v3 = make_whitelist(wl / "3M-february-2018.txt", 20000, seed=2)
    make_whitelist(wl / "737K-arc-v1.txt", 5000, seed=3)

    n = 20000
    # --- a true v2 library: 16 CB (from v2 list) + 10 UMI = 26 bp
    reads = [random.choice(v2) + rand_seq(10) for _ in range(int(n*0.92))]
    reads += [rand_seq(26) for _ in range(n - len(reads))]   # 8% junk barcodes
    random.shuffle(reads); write_fastq(out / "v2_R1.fastq.gz", reads)

    # --- a true v3 library: 16 CB (from v3 list) + 12 UMI = 28 bp
    reads = [random.choice(v3) + rand_seq(12) for _ in range(int(n*0.90))]
    reads += [rand_seq(28) for _ in range(n - len(reads))]
    random.shuffle(reads); write_fastq(out / "v3_R1.fastq.gz", reads)

    # --- FAILURE MODE 1: R1 pre-trimmed to 16 bp (barcode only, UMI gone).
    # A naive length check would call this "too short" and stop; we want to
    # confirm the script refuses rather than silently picking something.
    reads = [random.choice(v3) for _ in range(n)]
    write_fastq(out / "trimmed_R1.fastq.gz", reads)

    # --- FAILURE MODE 2: R2 handed in by mistake (cDNA, 91 bp, no barcodes).
    # Length passes every check. Only the whitelist rate catches it.
    write_fastq(out / "r2_by_mistake.fastq.gz", [rand_seq(91) for _ in range(n)])

    # --- FAILURE MODE 3: right length, wrong whitelist (v2-length reads whose
    # barcodes come from the v3 list). Mimics a non-10x or mismatched kit.
    reads = [random.choice(v3) + rand_seq(10) for _ in range(n)]
    write_fastq(out / "crosslist_R1.fastq.gz", reads)

    print(f"synthetic data in {out}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/synth")

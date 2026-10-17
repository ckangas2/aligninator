#!/usr/bin/env python3
"""Synthetic R1 with a KNOWN barcode length, for barcode_structure.py."""
import gzip, sys, numpy as np
from pathlib import Path

def make(path, cb_len, umi_len, n_cells, n_reads, ambient_frac, seed):
    rng = np.random.default_rng(seed)
    B = np.array(list("ACGT"))
    def rand(n, L): return ["".join(x) for x in B[rng.integers(0,4,(n,L))]]
    cells = rand(n_cells, cb_len)
    # heavy-tailed: a few beads dominate, as in a real droplet library
    w = rng.pareto(1.2, n_cells) + 1.0; w /= w.sum()
    n_amb = int(n_reads * ambient_frac)
    n_cell = n_reads - n_amb
    idx = rng.choice(n_cells, size=n_cell, p=w)
    bcs = [cells[i] for i in idx] + rand(n_amb, cb_len)   # ambient = junk beads
    umis = rand(n_reads, umi_len)
    rng.shuffle(bcs)
    with gzip.open(path, "wt") as fh:
        for i,(b,u) in enumerate(zip(bcs, umis)):
            s = b+u
            fh.write(f"@r{i}\n{s}\n+\n{'I'*len(s)}\n")

if __name__ == "__main__":
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    make(out/"dropseq_R1.fastq.gz", 12, 8, 2000, 200000, 0.35, 1)
    make(out/"tenx_v2_R1.fastq.gz", 16,10, 3000, 200000, 0.25, 2)
    make(out/"tenx_v3_R1.fastq.gz", 16,12, 3000, 200000, 0.25, 3)
    print(f"wrote {out}")

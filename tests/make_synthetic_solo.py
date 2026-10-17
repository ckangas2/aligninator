#!/usr/bin/env python3
"""Synthetic STARsolo output: genes x cells matrix with a planted infection."""
import sys, numpy as np, scipy.io as sio, scipy.sparse as sp
from pathlib import Path
rng = np.random.default_rng(11)

def main(out, prefix="HSV1-"):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    host = ["GAPDH","ACTB","NECTIN1","LAT","RS1","CSF2","IFIT2","MKI67"] + \
           [f"HOSTG{i}" for i in range(200)]
    viral = [f"{prefix}{g}" for g in
             ["RL1","RL2","LAT","RS1","UL1","UL2","UL3","UL5","US1","US12","GFP"]]
    genes = host + viral
    n_cells = 500

    H = sp.random(len(host), n_cells, density=0.25, random_state=3,
                  data_rvs=lambda n: rng.integers(1, 40, n)).tocsr()
    # 30% of cells infected, with a wide viral load range
    V = np.zeros((len(viral), n_cells))
    inf = rng.choice(n_cells, size=int(n_cells*0.30), replace=False)
    for c in inf:
        load = rng.lognormal(3.0, 1.2)
        V[:, c] = rng.poisson(load, size=len(viral))
    V[-1, :] = 0; V[-1, inf] = rng.poisson(20, size=len(inf))   # GFP
    X = sp.vstack([H, sp.csr_matrix(V)]).tocsr()

    sio.mmwrite(str(out/"matrix.mtx"), X)
    (out/"barcodes.tsv").write_text("\n".join(f"CELL{i:05d}-1" for i in range(n_cells))+"\n")
    (out/"features.tsv").write_text("\n".join(f"ID_{g}\t{g}\tGene Expression" for g in genes)+"\n")
    print(f"solo dir: {out}  ({X.shape[0]} genes x {X.shape[1]} cells, {len(inf)} infected)")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv)>2 else "HSV1-")

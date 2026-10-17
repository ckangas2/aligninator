#!/usr/bin/env python3
"""Synthetic references that plant each failure validate_reference.py must catch."""
import random, sys
from pathlib import Path
random.seed(7)

def fasta(path, name, length):
    seq = "".join(random.choice("ACGT") for _ in range(length))
    with open(path, "w") as fh:
        fh.write(f">{name} synthetic\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + "\n")

def gtf_row(seq, feat, s, e, strand, gid, gname=None):
    a = f'gene_id "{gid}";'
    if gname: a += f' gene_name "{gname}";'
    a += f' transcript_id "{gid}-1";'
    return f"{seq}\tsynth\t{feat}\t{s}\t{e}\t.\t{strand}\t.\t{a}\n"

def main(out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    fasta(out/"viral.fa", "chrHSV1", 152000)

    # Realistic-ish HSV-1 layout: dense, overlapping ORFs on both strands,
    # including the two names that collide with human symbols.
    genes = [("RL1",1000,2200,"+"),("RL2",2100,4500,"-"),("LAT",3000,5500,"+"),
             ("UL1",6000,7200,"-"),("UL2",7100,8400,"-"),("UL3",8300,9100,"+"),
             ("UL4",9000,9900,"+"),("UL5",9800,12500,"-"),("UL6",12400,14300,"-"),
             ("UL7",14200,15100,"+"),("UL8",15000,17300,"+"),("UL9",17200,19800,"-"),
             ("US1",120000,121500,"+"),("US12",128000,129000,"-"),
             ("RS1",130000,134000,"-"),("GFP",140000,140800,"+")]

    # --- BAD deposit: contig "HSV1" in GTF vs ">chrHSV1" in FASTA,
    #     no gene rows, bare gene_names, heavy overlap.
    with open(out/"viral_bad.gtf","w") as fh:
        for g,s,e,st in genes:
            fh.write(gtf_row("HSV1","exon",s,e,st,g,g))
            fh.write(gtf_row("HSV1","CDS",s,e,st,g,g))

    # --- GOOD: contig matches, exons present, names prefixed.
    with open(out/"viral_good.gtf","w") as fh:
        for g,s,e,st in genes:
            fh.write(gtf_row("chrHSV1","gene",s,e,st,f"HSV1-{g}",f"HSV1-{g}"))
            fh.write(gtf_row("chrHSV1","exon",s,e,st,f"HSV1-{g}",f"HSV1-{g}"))

    # --- CDS-only: the fatal "no exon rows" case.
    with open(out/"viral_cdsonly.gtf","w") as fh:
        for g,s,e,st in genes:
            fh.write(gtf_row("chrHSV1","CDS",s,e,st,g,g))

    # --- minimal host GTF carrying the real human LAT and RS1.
    host = [("LAT","chr16",1000,5000,"+"),("RS1","chrX",2000,6000,"-"),
            ("NECTIN1","chr11",7000,9000,"+"),("CSF2","chr5",10000,11000,"+")]
    with open(out/"host.gtf","w") as fh:
        for g,c,s,e,st in host:
            fh.write(gtf_row(c,"exon",s,e,st,f"ENSG_{g}",g))
    print(f"synthetic refs in {out}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv)>1 else "tests/ref")

#!/usr/bin/env python3
"""
iSTART: does HM13 end inside the intron that holds MCTS2? (intronic polyadenylation, IPA)

Evidence, per cell line (total RNA-seq BAMs):
 1. Poly(A)-tailed reads: reads whose alignment ends in a soft-clipped run of >= --min-a nucleotides that is
    >= 85% A (for genes on the + strand; T at the start for - strand). Ends are dropped when the genome right after
    the end is A-rich (>= 7 A in 10 nt: internal priming). Ends within 30 nt are clustered.
 2. PolyASite 2.0 clusters and annotated 3' ends of HM13 transcripts (GENCODE) in the same region.
 3. RNA-seq coverage (RNA strand) across the MCTS2 intron in 500-nt bins.
Poly(A) read ends are counted per HM13 intron and normalised by intron length, so the MCTS2 intron can be
compared with the host's other introns.
"""
import argparse
import gzip
import importlib.util
import os
import re
import subprocess
from collections import defaultdict

import numpy as np
import pyBigWig

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, fname):
    s = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


kg = load("kg", "30_known_genes.py")
u1 = load("u1", "25_u1pas.py")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def polya_ends(sam_lines, strand, min_a=8, frac=0.85):
    """Positions (0-based, last transcribed base) of reads ending in a non-genomic poly(A) tail."""
    out = []
    for line in sam_lines:
        f = line.split("\t")
        if len(f) < 10:
            continue
        pos, cig, seq = int(f[3]) - 1, f[5], f[9]
        ops = [(int(n), o) for n, o in CIG.findall(cig)]
        if not ops:
            continue
        ref_len = sum(n for n, o in ops if o in "MDN=X")
        if strand == "+" and ops[-1][1] == "S" and ops[-1][0] >= min_a:
            clip = seq[-ops[-1][0]:]
            if clip.count("A") >= frac * len(clip):
                out.append(pos + ref_len - 1)
        if strand == "-" and ops[0][1] == "S" and ops[0][0] >= min_a:
            clip = seq[:ops[0][0]]
            if clip.count("T") >= frac * len(clip):
                out.append(pos)
    return out


def cluster(positions, w=30):
    out = []
    for p in sorted(positions):
        if out and p - out[-1][-1] <= w:
            out[-1].append(p)
        else:
            out.append([p])
    return [(int(np.median(c)), len(c)) for c in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", default="ref/gencode.v47.annotation.gtf.gz")
    ap.add_argument("--fasta", default="ref/GRCh38.primary_assembly.genome.fa")
    ap.add_argument("--pas", default="ref/atlas.clusters.2.0.GRCh38.96.bed.gz")
    ap.add_argument("--host", default="HM13")
    ap.add_argument("--nested", default="MCTS2")
    ap.add_argument("--line", action="append", required=True,
                    help="NAME|bam1,bam2|rna_plus.bw,rna_minus.bw")
    ap.add_argument("--min-a", type=int, default=8)
    a = ap.parse_args()

    genes, tx, exons = kg.parse_gtf(a.gtf)
    by_name = {g["name"]: g for g in genes.values()}
    H, N = by_name[a.host], by_name[a.nested]
    chrom, st = H["chrom"], H["strand"]
    htx = [t for tid, t in tx.items() if t["gene_id"] == H["gene_id"]]
    hids = [tid for tid, t in tx.items() if t["gene_id"] == H["gene_id"]]
    ntx = [tid for tid, t in tx.items() if t["gene_id"] == N["gene_id"]]
    n_tss = sorted({(tx[t]["start"] if st == "+" else tx[t]["end"] - 1) for t in ntx})
    can = [t for t in hids if "Ensembl_canonical" in tx[t]["tags"]][0]
    cex = sorted(exons[can])
    introns = [(e1, s2) for (s1, e1), (s2, e2) in zip(cex[:-1], cex[1:])]
    mi = [i for i, (s, e) in enumerate(introns) if s <= min(n_tss) < e]
    mi = mi[0] if mi else None
    print(f"{a.host}: {chrom}:{H['start']}-{H['end']} ({st}), {len(hids)} transcripts; {a.nested} starts at "
          f"{', '.join(map(str, n_tss))}; canonical intron holding it: "
          f"{'#' + str(mi + 1) + ' ' + str(introns[mi]) if mi is not None else 'none'}")
    all_ex = sorted({e for t in hids for e in exons[t]})

    def where(p):
        if any(s <= p < e for s, e in all_ex):
            return f"{a.host} exon"
        for i, (s, e) in enumerate(introns):
            if s <= p < e:
                tag = f"intron {i + 1}"
                if i == mi:
                    first = min(n_tss) if st == "+" else max(n_tss)
                    tag += (" (before " if (p < first) == (st == "+") else " (after ") + f"{a.nested} start)"
                return tag
        return "outside canonical introns"

    # annotated host 3' ends and PolyASite clusters
    ends = sorted({(t["end"] - 1 if st == "+" else t["start"]) for t in htx})
    print(f"\nAnnotated 3' ends of {a.host} transcripts: " + "; ".join(f"{p} [{where(p)}]" for p in ends))
    print("PolyASite 2.0 clusters on the host strand:")
    with gzip.open(a.pas, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            c = f[0] if f[0].startswith("chr") else "chr" + f[0]
            if c == chrom and f[5].strip() == st and H["start"] <= int(f[1]) < H["end"]:
                p = int(f[2]) - 1 if st == "+" else int(f[1])
                print(f"  {p}  [{where(p)}]  {f[3]}")

    fa = u1.Fasta(a.fasta)
    for spec in a.line:
        name, bams, rna = spec.split("|")
        ends_all = []
        for b in bams.split(","):
            r = subprocess.run(["samtools", "view", "-F", "0x904", "-q", "10", b, f"{chrom}:{H['start'] + 1}-{H['end']}"],
                               capture_output=True, text=True, check=True)
            ends_all += polya_ends(r.stdout.splitlines(), st, a.min_a)
        kept, primed = [], 0
        for p in ends_all:
            down = fa.get(chrom, p + 1, p + 11) if st == "+" else fa.get(chrom, p - 10, p).translate(u1.COMP)[::-1]
            if down.count("A") >= 7:
                primed += 1
            else:
                kept.append(p)
        cl = cluster(kept)
        print(f"\n===== {name} =====")
        print(f"reads ending in a poly(A) tail: {len(ends_all)}; dropped as internal priming: {primed}; kept {len(kept)} in {len(cl)} clusters")
        for p, n in sorted(cl, key=lambda x: -x[1])[:15]:
            print(f"  {chrom}:{p}  reads {n:>4}  [{where(p)}]")
        per = defaultdict(int)
        for p in kept:
            for i, (s, e) in enumerate(introns):
                if s <= p < e:
                    per[i] += 1
        print("poly(A) read ends per kb of each canonical intron:")
        print("  " + "  ".join(f"#{i + 1}{'*' if i == mi else ''}: {per[i] / ((e - s) / 1000):.2f}" for i, (s, e) in enumerate(introns)))
        if mi is not None:
            s, e = introns[mi]
            p_, m_ = rna.split(",")
            bw = pyBigWig.open(p_ if st == "+" else m_)
            v = np.abs(np.nan_to_num(bw.values(chrom, s, e, numpy=True)))
            bins = [(s + k, v[k:k + 500].mean()) for k in range(0, len(v), 500)]
            print(f"RNA-seq coverage across intron #{mi + 1} (* = bin containing the {a.nested} start), 500-nt bins:")
            print("  " + "  ".join(f"{b0}{'*' if any(b0 <= t < b0 + 500 for t in n_tss) else ''}:{m:.2f}" for b0, m in bins))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
iSTART: does host transcription suppress initiation on its own strand? (sense–antisense asymmetry)

Almost all PRO-cap start sites are divergent pairs: one arm on the plus strand, the other on minus, about 100–300 nt
apart, in the same surroundings. For each pair lying in an intron of a protein-coding gene (both arms
at least --min-exon-dist from its exons) we compare the strength of the arm on the host strand (sense)
with the opposite arm (anti):
    DI = log2[(sense signal + 1) / (anti signal + 1)]     (capped signal within ±50 nt of the arm)
DI < 0: initiation on the host strand weaker than on the opposite strand. Pairs between genes (arbitrary orientation)
give the reference distribution. The suppression hypothesis predicts DI < 0 in active genes, the more so
the stronger the local host transcription (Bru-seq on the host strand, 5–1 kb upstream of the pair)
and when H3K36me3 is nearby; in silent genes DI ≈ 0.
"""
import argparse
import gzip
import importlib.util
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, fname):
    s = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


c10 = load("c10", "10_candidates.py")
c11 = load("c11", "11_evidence.py")
c12 = load("c12", "12_itss.py")
c21 = load("c21", "21_describe.py")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--bi", required=True, help="PINTS bidirectional peaks (col. 5 = plus strand)")
    ap.add_argument("--cap-plus", required=True)
    ap.add_argument("--cap-minus", required=True)
    ap.add_argument("--host-tss", required=True, help="host_tss.tsv (step 10)")
    ap.add_argument("--bru-plus", required=True, help="Bru10, bigWig of the plus READ strand")
    ap.add_argument("--bru-minus", required=True)
    ap.add_argument("--bru-flip", choices=["yes", "no"], default="yes", help="reads are on the strand opposite to the RNA")
    ap.add_argument("--h3k36me3", help="H3K36me3 peaks (optional)")
    ap.add_argument("--min-exon-dist", type=int, default=1000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    tx, exons = c10.parse_gtf(a.gtf)
    hosts = c10.build_hosts(tx, exons, {"protein_coding"}, set(c10.HOST_EXON_TX_TYPES_DEFAULT.split(",")))
    hidx = c10.index_hosts(hosts)
    gene_any = defaultdict(list)
    for t in tx.values():
        if t["gene_type"] == "protein_coding":
            gene_any[t["chrom"]].append((t["start"], t["end"]))
    GA = {c: c10.IntervalSet(v) for c, v in gene_any.items()}
    E = c10.IntervalSet([])
    cap = c11.Track(a.cap_plus, a.cap_minus, "capped")
    bru = c11.Track(a.bru_plus, a.bru_minus, "bru")
    bru.flip = a.bru_flip == "yes"
    k36 = c21.load_peaks(a.h3k36me3) if a.h3k36me3 else None

    # host activity: the strongest of its annotated start sites (±100 nt)
    ht = pd.read_csv(a.host_tss, sep="\t")
    act = {}
    for g, df in ht.drop_duplicates(["host_gene_id", "pos"]).groupby("host_gene_id"):
        act[g] = max(cap.values(r.chrom, r.pos - 100, r.pos + 101, r.strand).sum() for r in df.itertuples())

    def sig(chrom, p, strand):
        return cap.values(chrom, p - 50, p + 51, strand).sum()

    rows = []
    for c, p5, p6, conf in c12.read_bi(a.bi):
        pp, pm = p5, p6  # col. 5 = plus strand (checked in step 12)
        lo, hi = min(pp, pm), max(pp, pm)
        mid = (pp + pm) // 2
        cands = []
        for s in ("+", "-"):
            for h in c10.containing_hosts(hidx, c, s, lo, hi + 1):
                ri = c10.ref_intron_of(h, lo, hi + 1)
                if ri and h["exon_set"].dist(lo) >= a.min_exon_dist and h["exon_set"].dist(hi) >= a.min_exon_dist:
                    cands.append(h)
        strands = {h["strand"] for h in cands}
        if len(strands) == 2:
            continue  # genes on both strands: ambiguous
        if cands:
            h = min(cands, key=lambda x: x["end"] - x["start"])
            s = h["strand"]
            ps, pa = (pp, pm) if s == "+" else (pm, pp)
            a_s, a_a = sig(c, ps, s), sig(c, pa, "-" if s == "+" else "+")
            # local host transcription: Bru10 on the host strand, 5–1 kb upstream of the sense arm
            ws, we = (ps - 5000, ps - 1000) if s == "+" else (ps + 1001, ps + 5001)
            host_tx = bru.values(c, ws, we, s).mean()
            # klasa values written to the table: gene intron / intergenic (below)
            rows.append(dict(chrom=c, mid=mid, klasa="intron genu", host=h["gene_id"], host_name=h["gene_name"],
                             host_signal=act.get(h["gene_id"], 0.0), sens=a_s, anty=a_a,
                             DI=np.log2((a_s + 1) / (a_a + 1)), host_tx=host_tx,
                             k36=c21.has_peak(k36, c, mid, 1000) if k36 else np.nan))
        elif not GA.get(c, E).overlaps(lo - 5000, hi + 5000):
            a_s, a_a = sig(c, pp, "+"), sig(c, pm, "-")
            rows.append(dict(chrom=c, mid=mid, klasa="międzygenowa", host="NA", host_name="NA", host_signal=np.nan,
                             sens=a_s, anty=a_a, DI=np.log2((a_s + 1) / (a_a + 1)), host_tx=np.nan,
                             k36=c21.has_peak(k36, c, mid, 1000) if k36 else np.nan))
    d = pd.DataFrame(rows)
    d.to_csv(a.out, sep="\t", index=False)

    def line(name, g):
        if len(g) < 5:
            return f"  {name:<46}{len(g):>6}"
        p = wilcoxon(g.DI).pvalue if (g.DI != 0).any() else np.nan
        return (f"  {name:<46}{len(g):>6}{g.DI.median():>12.2f}{(g.DI < 0).mean() * 100:>9.0f}%"
                f"{2 ** g.DI.median():>12.2f}{p:>11.1e}")

    gi = d[d.klasa == "intron genu"]
    gi_act = gi[gi.host_signal >= 5]
    gi_sil = gi[gi.host_signal < 5]
    print(f"Divergent pairs in protein-coding gene introns (≥{a.min_exon_dist} nt from exons) and between genes (≥5 kb):")
    print(f"  {'group':<46}{'n':>6}{'median DI':>12}{'DI<0':>9}{'sense/anti':>12}{'p (DI=0)':>11}")
    print(line("intergenic pairs (arbitrary orientation)", d[d.klasa == "międzygenowa"]))
    print(line("intron of silent gene", gi_sil))
    print(line("intron of active gene", gi_act))
    if len(gi_act) >= 15:
        # thirds of local host transcription
        q = pd.qcut(gi_act.host_tx.rank(method="first"), 3, labels=["weak", "medium", "strong"])
        for lev in ["weak", "medium", "strong"]:
            print(line(f"  local host transcription: {lev}", gi_act[q == lev]))
        r, pr = spearmanr(gi_act.host_tx, gi_act.DI)
        print(f"  correlation of DI with local host transcription (Spearman): r = {r:.2f}, p = {pr:.1e}")
    if k36 is not None:
        print(line("intron of active gene, with H3K36me3 (±1 kb)", gi_act[gi_act.k36 == True]))
        print(line("intron of active gene, no H3K36me3", gi_act[gi_act.k36 == False]))
    print("  DI = log2[(signal of the arm on the host strand + 1) / (signal of the opposite arm + 1)]; "
          "sense/anti = 2^median DI")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
iSTART: does an intronic location change the fate of a start site? Intronic starts vs starts between genes.

Groups (column "grupa"; only UNANNOTATED starts: farther than ±500 nt from any GENCODE start):
  intron_akt  – capped iTSS in an intron of an active host (from the describe table, step 21)
  intron_cichy– capped iTSS in an intron of a silent host
  miedzygen   – starts from PRO-cap peaks at least --gap-pc nt from protein-coding genes (both strands),
                with the same cap test (≥5 reads, Poisson vs flanks, ratio to uncapped RNA)
  wewn_prom   – capped internal promoters of the host (reference)
  kontrole    – random positions in introns of active hosts (background)
Measures (on the strand transcribed from the start):
  RNA_nadw  – stable RNA excess (RNA-seq): mean coverage 0..3 kb minus mean -3.2..-0.2 kb, per PRO-cap read
  RNA_skok  – stable RNA rise: log2 (0..3 kb / -3.2..-0.2 kb); > 95th percentile of controls
  DRB_dal   – (optional) excess 10-min signal in the far 2..8 kb window over the level before the start, per
              PRO-cap read; plus fractions with a near and a far burst (ratio to 40 min, as in step 27)
  PAS sens/anty – U1–PAS axis, sense/antisense (1 kb)
  MTREX_d   – (optional) change of the rise after MTREX depletion; > 95th percentile of controls
Comparisons within tertiles of start strength (PRO-cap, shared by the three start-site groups).
"""
import argparse
import bisect
import importlib.util
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import mannwhitneyu, poisson

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name, fname):
    s = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


c10 = load("c10", "10_candidates.py")
c12 = load("c12", "12_itss.py")
c25 = load("c25", "25_u1pas.py")


class Cov:
    def __init__(self, plus, minus, flip=False):
        self.bw = {"+": pyBigWig.open(plus), "-": pyBigWig.open(minus)}
        self.flip = flip

    def mean(self, chrom, s, e, strand):
        if self.flip:
            strand = "-" if strand == "+" else "+"
        bw = self.bw[strand]
        L = bw.chroms(chrom)
        if L is None:
            return 0.0
        s, e = max(0, s), min(L, e)
        if e <= s:
            return 0.0
        v = bw.values(chrom, s, e, numpy=True)
        return float(np.abs(np.nan_to_num(v)).mean()) if len(v) else 0.0


def win(p, st, a, b):
    return (p + a, p + b) if st == "+" else (p - b + 1, p - a + 1)


def pair(spec, flip=False):
    if not spec:
        return None
    p, m = spec.split(",")
    return Cov(p, m, flip)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--describe", required=True)
    ap.add_argument("--uni", required=True)
    ap.add_argument("--bi", required=True)
    ap.add_argument("--cap", required=True, help="plus,minus")
    ap.add_argument("--uncap", help="plus,minus (optional)")
    ap.add_argument("--rna", required=True, help="plus,minus (RNA-seq, RNA strand)")
    ap.add_argument("--drb", help="plus,minus of the 10-min sample (reads on the opposite strand)")
    ap.add_argument("--drb-ref", help="plus,minus of the 40-min sample")
    ap.add_argument("--drb-front", type=int, default=16000)
    ap.add_argument("--ko", help="plus,minus RNA-seq after MTREX depletion")
    ap.add_argument("--ko-ctrl", help="plus,minus RNA-seq of the control")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gap-pc", type=int, default=5000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # --- annotation: starts of all transcripts and spans of protein-coding genes ---
    tx, exons = c10.parse_gtf(a.gtf)
    tss = defaultdict(list)
    pcs = defaultdict(list)
    for t in tx.values():
        tss[(t["chrom"], t["strand"])].append(c10.tss_of(t))
        if t["gene_type"] == "protein_coding":
            pcs[t["chrom"]].append((t["start"], t["end"]))
    tss = {k: np.sort(np.array(v)) for k, v in tss.items()}
    PC = {c: c10.IntervalSet(v) for c, v in pcs.items()}
    E = c10.IntervalSet([])

    def annotated(c, st, p):
        v = tss.get((c, st))
        if v is None or not len(v):
            return False
        j = np.searchsorted(v, p)
        return any(0 <= jj < len(v) and abs(v[jj] - p) <= 500 for jj in (j - 1, j))

    cap, uncap = pair(a.cap), pair(a.uncap)
    rna = pair(a.rna)
    drb, drbr = pair(a.drb, True), pair(a.drb_ref, True)
    ko, koc = pair(a.ko), pair(a.ko_ctrl)
    fa = c25.Fasta(a.fasta)
    cap_tot = sum(abs(cap.bw[s].header()["sumData"]) for s in "+-")
    unc_tot = sum(abs(uncap.bw[s].header()["sumData"]) for s in "+-") if uncap else None

    def capsum(c, p, st):
        return cap.mean(c, p - 100, p + 101, st) * 201

    def cap_test(c, p, st):
        s = capsum(c, p, st)
        fl = (cap.mean(c, p - 2000, p - 500, st) + cap.mean(c, p + 500, p + 2000, st)) / 2
        lam = max(fl, 1e-3) * 201
        pv = poisson.sf(s - 1, lam) if s > 0 else 1.0
        ok = s >= 5 and pv < 1e-3
        if ok and uncap:
            u = uncap.mean(c, p - 100, p + 101, st) * 201
            ok = ((s + 1) / cap_tot) / ((u + 1) / unc_tot) >= 2
        return s, ok

    # --- start sites from the describe table (intron) ---
    d = pd.read_csv(a.describe, sep="\t", low_memory=False)
    for col in ("host_active", "tss_supported"):
        d[col] = d[col].astype(str) == "True"
    up = d.get("upstream_of_used_host_tss", pd.Series(False, index=d.index)).astype(str) == "True"
    sites = []
    for r in d[~up].itertuples():
        if r.site_type == "candidate" and r.tss_supported:
            g = "intron_akt" if r.host_active else "intron_cichy"
        elif r.site_type == "internal_promoter" and r.tss_supported and r.host_active:
            g = "wewn_prom"
        elif r.site_type == "random_control" and r.host_active:
            g = "kontrole"
        else:
            continue
        if g.startswith("intron") and annotated(r.chrom, r.strand, int(r.pos)):
            continue
        dist = r.dist_used_tss if r.host_active else np.nan
        sites.append(dict(site_id=r.site_id, grupa=g, chrom=r.chrom, strand=r.strand, pos=int(r.pos), dist=dist,
                          cap_sum=capsum(r.chrom, int(r.pos), r.strand)))

    # --- start sites between genes, from PRO-cap peaks ---
    starts = [(c, s, p) for c, s, p, _, _ in c12.read_uni(a.uni)]
    for c, p5, p6, _ in c12.read_bi(a.bi):
        starts += [(c, "+", p5), (c, "-", p6)]
    starts.sort()
    seen = defaultdict(list)
    n_ig = 0
    for c, s, p in starts:
        if PC.get(c, E).overlaps(p - a.gap_pc, p + a.gap_pc) or annotated(c, s, p):
            continue
        if any(abs(p - q) <= 50 for q in seen[(c, s)][-3:]):
            continue
        seen[(c, s)].append(p)
        cs, ok = cap_test(c, p, s)
        if ok:
            sites.append(dict(site_id=f"IG:{c}:{p}:{s}", grupa="miedzygen", chrom=c, strand=s, pos=p, dist=np.nan, cap_sum=cs))
            n_ig += 1
    print(f"[{a.name}] start sites between genes (≥{a.gap_pc} nt from protein-coding genes, unannotated, capped): {n_ig}",
          file=sys.stderr)
    S = pd.DataFrame(sites)

    # --- measures ---
    rows = []
    for r in S.itertuples():
        c, p, st = r.chrom, r.pos, r.strand
        rec = {}
        dn, upv = rna.mean(c, *win(p, st, 0, 3000), st), rna.mean(c, *win(p, st, -3200, -200), st)
        rec["RNA_skok"] = np.log2((dn + 1e-3) / (upv + 1e-3))
        rec["RNA_nadw"] = (dn - upv) / max(r.cap_sum, 1)
        if drb is not None:
            T = {k: drb.mean(c, *win(p, st, *w), st) for k, w in (("up", (-5200, -200)), ("near", (0, 2000)), ("far", (2000, 8000)))}
            R = {k: drbr.mean(c, *win(p, st, *w), st) for k, w in (("up", (-5200, -200)), ("near", (0, 2000)), ("far", (2000, 8000)))}
            rec["DRB_dal"] = (T["far"] - T["up"]) / max(r.cap_sum, 1)
            rec["DRB_far_ratio"] = np.log2((T["far"] + 1e-3) / (T["up"] + 1e-3)) - np.log2((R["far"] + 1e-3) / (R["up"] + 1e-3))
            rec["DRB_near_ratio"] = np.log2((T["near"] + 1e-3) / (T["up"] + 1e-3)) - np.log2((R["near"] + 1e-3) / (R["up"] + 1e-3))
        if ko is not None:
            kd, ku = ko.mean(c, *win(p, st, 0, 3000), st), ko.mean(c, *win(p, st, -3200, -200), st)
            cd, cu = koc.mean(c, *win(p, st, 0, 3000), st), koc.mean(c, *win(p, st, -3200, -200), st)
            rec["MTREX_d"] = np.log2((kd + .01) / (ku + .01)) - np.log2((cd + .01) / (cu + .01))
        seq_s = c25.downstream(fa, c, p, st, 1000)
        seq_a = c25.downstream(fa, c, p - 300 if st == "+" else p + 300, "-" if st == "+" else "+", 1000)
        rec["pas_sens"] = c25.scan(seq_s)[0] if len(seq_s) >= 900 else np.nan
        rec["pas_anty"] = c25.scan(seq_a)[0] if len(seq_a) >= 900 else np.nan
        rows.append(rec)
    S = pd.concat([S.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    # DRB: only sites ahead of the host's wave front
    if drb is not None:
        bad = (S.grupa.isin(["intron_akt", "wewn_prom", "kontrole"])) & ~(S.dist >= a.drb_front + 6200)
        for col in ("DRB_dal", "DRB_far_ratio", "DRB_near_ratio"):
            S.loc[bad, col] = np.nan
    S.to_csv(a.out, sep="\t", index=False)

    ctrl = S[S.grupa == "kontrole"]
    thr = {col: np.nanpercentile(ctrl[col], 95) for col in ("RNA_skok", "DRB_far_ratio", "DRB_near_ratio", "MTREX_d") if col in S}
    st3 = S[S.grupa.isin(["intron_akt", "intron_cichy", "miedzygen"])].copy()
    st3["sila"] = pd.qcut(st3.cap_sum.rank(method="first"), 3, labels=["weak", "medium", "strong"])
    S = S.merge(st3[["site_id", "sila"]], on="site_id", how="left")

    names = {"intron_akt": "intronic, active host", "intron_cichy": "intronic, silent host",
             "miedzygen": "between genes", "wewn_prom": "internal promoters", "kontrole": "random intronic positions"}
    print(f"\n===== {a.name} =====")
    print("Counts and strength (median PRO-cap): " + "; ".join(
        f"{names[g]} {int((S.grupa == g).sum())} ({S[S.grupa == g].cap_sum.median():.0f})" for g in names))

    def block(title, col, kind):
        if col not in S:
            return
        print(f"\n{title}")
        print(f"  {'tertile':<14}" + "".join(f"{names[g][:24]:>27}" for g in ("intron_akt", "intron_cichy", "miedzygen"))
              + f"{'p (act. vs intergenic)':>24}")
        for s in ("weak", "medium", "strong"):
            cells, vals = [], {}
            for g in ("intron_akt", "intron_cichy", "miedzygen"):
                x = S[(S.grupa == g) & (S.sila == s)][col].dropna()
                vals[g] = x
                if kind == "frac":
                    cells.append(f"{(x > thr[col]).mean() * 100:>22.0f}% n={len(x):<4}" if len(x) else f"{'–':>27}")
                else:
                    cells.append(f"{x.median():>22.4f} n={len(x):<4}" if len(x) else f"{'–':>27}")
            a1, b1 = vals["intron_akt"], vals["miedzygen"]
            p = mannwhitneyu(a1, b1).pvalue if len(a1) > 5 and len(b1) > 5 else np.nan
            print(f"  {s:<14}" + "".join(cells) + f"{p:>24.1e}")
        ref = S[S.grupa == "wewn_prom"][col].dropna()
        if len(ref):
            v = f"{(ref > thr[col]).mean() * 100:.0f}%" if kind == "frac" else f"{ref.median():.4f}"
            print(f"  reference: internal promoters {v} (n={len(ref)}); random controls " +
                  (f"{(ctrl[col].dropna() > thr[col]).mean() * 100:.0f}%" if kind == "frac" else f"{ctrl[col].median():.4f}"))

    block("Stable RNA: fraction with a rise above the 95th percentile of controls", "RNA_skok", "frac")
    block("Stable RNA: excess after the start per PRO-cap read (median; independent of upstream background)", "RNA_nadw", "med")
    block("DRB 10 min: fraction with a far burst (2–8 kb) above the control threshold", "DRB_far_ratio", "frac")
    block("DRB 10 min: fraction with a near burst (0–2 kb) above the control threshold", "DRB_near_ratio", "frac")
    block("DRB 10 min: excess in the 2–8 kb window per PRO-cap read (median)", "DRB_dal", "med")
    block("MTREX depletion: fraction with an increase of the rise above the control threshold", "MTREX_d", "frac")
    print("\nU1–PAS axis (PAS per kb in the first kb; sense / antisense / ratio):")
    for g in names:
        x = S[S.grupa == g]
        ps, pa = x.pas_sens.mean(), x.pas_anty.mean()
        print(f"  {names[g]:<34} {ps:.2f} / {pa:.2f} / {ps / pa if pa else float('nan'):.2f}")


if __name__ == "__main__":
    main()

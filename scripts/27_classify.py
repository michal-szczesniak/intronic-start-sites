#!/usr/bin/env python3
"""
iSTART: classification of single start sites from DRB time courses with replicates (HCT116, DRB/TTchem-seq).

For each site and each replicate, in the test sample T (e.g. 10 min) and the reference R (e.g. 40 min):
  up   = -5.2..-0.2 kb, near = 0..2 kb, far = 2..8 kb (on the transcribed strand)
  near_score = log2[(T_near+0.5)/(T_up+0.5)] - log2[(R_near+0.5)/(R_up+0.5)]   (read counts)
  far_score  = the same for the far window
Thresholds: 95th percentile of random controls (active host), separately per replicate and window.
Classes (column "klasa"; must agree in BOTH replicates):
  produktywny (productive)        = far score above threshold
  krótki (short)                  = near score above threshold, far score not
  bez odpowiedzi (no response)    = neither above threshold
  nierozstrzygnięty (unresolved)  = too few reads in either replicate
Eligibility: silent host, or start farther from the host start than the wave front + 6.2 kb.
Validation (optional): U1–PAS protection, stable RNA, phyloP, start strength, chromatin type.
"""
import argparse

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import fisher_exact


class Cov:
    def __init__(self, plus, minus):
        self.bw = {"plus": pyBigWig.open(plus), "minus": pyBigWig.open(minus)}
        self.flip = False

    def sum(self, chrom, s, e, strand):
        k = "plus" if strand == "+" else "minus"
        if self.flip:
            k = "minus" if k == "plus" else "plus"
        bw = self.bw[k]
        L = bw.chroms(chrom)
        if L is None:
            return 0.0
        s, e = max(0, s), min(L, e)
        if e <= s:
            return 0.0
        v = bw.values(chrom, s, e, numpy=True)
        return float(np.abs(np.nan_to_num(v)).sum()) if len(v) else 0.0


def win(p, st, a, b):
    return (p + a, p + b) if st == "+" else (p - b + 1, p - a + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--test", nargs="+", required=True, help="plus,minus bigWig pairs of the test sample, one per replicate")
    ap.add_argument("--ref", nargs="+", required=True, help="plus,minus bigWig pairs of the reference, in the same order")
    ap.add_argument("--front", type=int, default=16000, help="wave front in the test sample (nt)")
    ap.add_argument("--read-len", type=float, default=75)
    ap.add_argument("--min-reads", type=float, default=20, help="min. reads in the windows (T+R) per replicate")
    ap.add_argument("--u1pas")
    ap.add_argument("--expr")
    ap.add_argument("--describe")
    ap.add_argument("--phylop")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    assert len(a.test) == len(a.ref)

    ev = pd.read_csv(a.evidence, sep="\t", low_memory=False)
    for c in ("host_active", "tss_supported"):
        ev[c] = ev[c].astype(str) == "True"
    up_flag = ev.get("upstream_of_used_host_tss", pd.Series(False, index=ev.index)).astype(str) == "True"
    min_d = a.front + 6200
    ev["eligible"] = ~up_flag & (~ev.host_active | (ev.dist_used_tss >= min_d))
    ev = ev[ev.eligible & ev.site_type.isin(["candidate", "internal_promoter", "random_control"])].copy()

    reps = []
    for t, r in zip(a.test, a.ref):
        T, R = Cov(*t.split(",")), Cov(*r.split(","))
        hp = ev[ev.host_active].drop_duplicates("host_gene_id").dropna(subset=["host_tss_used"]).head(1500)
        same = sum(R.sum(x.chrom, *win(int(x.host_tss_used), x.strand, 0, 5000), x.strand) for x in hp.itertuples())
        T.flip = R.flip = True
        opp = sum(R.sum(x.chrom, *win(int(x.host_tss_used), x.strand, 0, 5000), x.strand) for x in hp.itertuples())
        T.flip = R.flip = opp > same  # reads on the strand opposite to the RNA -> swap
        print(f"replicate {len(reps) + 1}: orientation {'swapped' if opp > same else 'not swapped'} (RNA strand {same:.0f}, opposite {opp:.0f})")
        reps.append((T, R))

    rows = []
    for x in ev.itertuples():
        p, st = int(x.pos), x.strand
        rec = dict(site_id=x.site_id)
        for i, (T, R) in enumerate(reps, 1):
            c = {}
            for nm, (s0, e0) in (("up", win(p, st, -5200, -200)), ("near", win(p, st, 0, 2000)), ("far", win(p, st, 2000, 8000))):
                c["T" + nm] = T.sum(x.chrom, s0, e0, st) / a.read_len
                c["R" + nm] = R.sum(x.chrom, s0, e0, st) / a.read_len
            lt = lambda k: np.log2((c["T" + k] + .5) / (c["Tup"] + .5)) - np.log2((c["R" + k] + .5) / (c["Rup"] + .5))
            rec[f"near{i}"], rec[f"far{i}"] = lt("near"), lt("far")
            rec[f"reads{i}"] = sum(c.values())
        rows.append(rec)
    d = ev.merge(pd.DataFrame(rows), on="site_id")
    n = len(reps)
    d["covered"] = np.all([d[f"reads{i}"] >= a.min_reads for i in range(1, n + 1)], axis=0)
    ctrl = d[(d.site_type == "random_control") & d.host_active & d.covered]
    thr = {f"{w}{i}": np.percentile(ctrl[f"{w}{i}"], 95) for w in ("near", "far") for i in range(1, n + 1)}
    near_pos = np.all([d[f"near{i}"] > thr[f"near{i}"] for i in range(1, n + 1)], axis=0)
    far_pos = np.all([d[f"far{i}"] > thr[f"far{i}"] for i in range(1, n + 1)], axis=0)
    d["klasa"] = np.where(~d.covered, "nierozstrzygnięty",
                          np.where(far_pos, "produktywny", np.where(near_pos, "krótki", "bez odpowiedzi")))

    # validation: attach features
    if a.u1pas:
        u = pd.read_csv(a.u1pas, sep="\t", low_memory=False)[["site_id", "prot_sens", "pas_sens"]]
        d = d.merge(u, on="site_id", how="left")
        d["chronione"] = d.prot_sens.astype(str) == "True"
    if a.expr:
        e = pd.read_csv(a.expr, sep="\t", low_memory=False)
        e["host_active"] = e.host_active.astype(str) == "True"
        thr_r = np.nanpercentile(e[(e.site_type == "random_control") & e.host_active].rna_step, 95)
        d = d.merge(e[["site_id", "rna_step"]], on="site_id", how="left")
        d["stabilne_RNA"] = d.rna_step > thr_r
    if a.describe and a.phylop:
        ds = pd.read_csv(a.describe, sep="\t", low_memory=False).reset_index()
        t = pd.read_csv(a.phylop, sep="\t", header=None, names=["n", "size", "covered", "sum", "mean0", "mean"])
        t["index"] = t.n.str[1:].astype(int)
        ds = ds.merge(t[["index", "mean"]], on="index", how="left")[["site_id", "typ", "mean"]].rename(columns={"mean": "phylop"})
        d = d.merge(ds, on="site_id", how="left")
        d["konserwowane"] = d.phylop > 1
    d.to_csv(a.out, sep="\t", index=False)

    order = ["produktywny", "krótki", "bez odpowiedzi", "nierozstrzygnięty"]
    print(f"Eligibility: silent host or start ≥ {min_d / 1000:.1f} kb from the host start; replicates: {n}")
    print("Thresholds (95th percentile of controls): " + ", ".join(f"{k} {v:.2f}" for k, v in thr.items()))
    groups = [("random controls", (d.site_type == "random_control") & d.host_active),
              ("capped internal promoters", (d.site_type == "internal_promoter") & d.host_active & d.tss_supported),
              ("capped iTSS, active host", (d.site_type == "candidate") & d.host_active & d.tss_supported),
              ("capped iTSS, inactive host", (d.site_type == "candidate") & ~d.host_active & d.tss_supported)]
    print(f"\n  {'group':<40}{'n':>6}" + "".join(f"{o[:14]:>16}" for o in order))
    for name, m in groups:
        g = d[m]
        if not len(g):
            continue
        cov = g[g.klasa != "nierozstrzygnięty"]
        cells = "".join(f"{(g.klasa == o).sum():>8} ({(g.klasa == o).sum() / max(len(cov), 1) * 100 if o != 'nierozstrzygnięty' else (g.klasa == o).mean() * 100:4.1f}%)" for o in order)
        print(f"  {name:<40}{len(g):>6}{cells}")
    print("  (class percentages among resolved sites; for 'nierozstrzygnięty' (unresolved) – among all sites)")

    it = d[(d.site_type == "candidate") & d.tss_supported & (d.klasa.isin(["produktywny", "krótki", "bez odpowiedzi"]))].copy()
    if len(it):
        print("\nFeatures of capped intronic start sites by class:")
        feats = [(c, lab) for c, lab in (("chronione", "U1-protected"), ("stabilne_RNA", "stable RNA"),
                                         ("konserwowane", "phyloP > 1")) if c in it]
        print(f"  {'class':<16}{'n':>6}{'median PRO-cap':>17}" + "".join(f"{lab:>15}" for _, lab in feats)
              + (f"{'promoter type':>15}" if "typ" in it else ""))
        for o in order[:3]:
            g = it[it.klasa == o]
            if not len(g):
                continue
            row = f"  {o:<16}{len(g):>6}{g.cap_sum.median():>17.0f}" + "".join(f"{g[c].mean() * 100:>14.0f}%" for c, _ in feats)
            if "typ" in it:
                row += f"{(g.typ == 'promotor-podobny').mean() * 100:>14.0f}%"
            print(row)
        for c, lab in feats:
            a_ = it[it.klasa == "produktywny"][c]; b_ = it[it.klasa == "krótki"][c]
            if len(a_) and len(b_):
                p = fisher_exact([[int(a_.sum()), len(a_) - int(a_.sum())], [int(b_.sum()), len(b_) - int(b_.sum())]])[1]
                print(f"  {lab}: productive {a_.mean() * 100:.0f}% vs short {b_.mean() * 100:.0f}%  (Fisher p = {p:.1e})")
        # the same within tertiles of start strength
        it["sila"] = pd.qcut(it.cap_sum.rank(method="first"), 3, labels=["weak", "medium", "strong"])
        if "chronione" in it:
            print("  U1 protection by start strength (productive vs short): " + "; ".join(
                f"{s}: {g[g.klasa == 'produktywny'].chronione.mean() * 100:.0f}% vs {g[g.klasa == 'krótki'].chronione.mean() * 100:.0f}%"
                for s, g in it.groupby("sila", observed=True)))


if __name__ == "__main__":
    main()

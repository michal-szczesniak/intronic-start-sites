#!/usr/bin/env python3
"""
iSTART: numbers for the supplementary figures, computed from the existing result files.

S1  Cap test: share of sites passing, by site type (random intronic positions, internal promoters, intronic
    start sites from PRO-cap peaks) and host activity, for each cell line; annotated intronic RNAs if given.
S2  Reproducibility in HCT116: intronic start sites confirmed in the CTCF degron line (active host) tested with
    the same cap test in the RAD21 degron line (no uncapped control), by thirds of strength in the CTCF line.
S3  The DRB wave along active host genes: nascent RNA at each time point divided by the reference, summed over
    genes, in 2-kb bins from the host start in use (genes extending at least --wave-len beyond it).
S6  Stable RNA by the rise measure: share of sites whose log2(after/before) exceeds the 95th percentile of random
    intronic positions, by thirds of strength, for intronic sites in active and silent hosts and sites between genes.

Library strand is detected per data set from the signal over the first 5 kb of active host genes.
"""
import argparse
import importlib.util
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import poisson, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))


class Cov:
    def __init__(self, spec):
        p, m = spec.split(",")
        self.bw = {"+": pyBigWig.open(p), "-": pyBigWig.open(m)}
        self.flip = False

    def vals(self, chrom, s, e, strand):
        if self.flip:
            strand = "-" if strand == "+" else "+"
        bw = self.bw[strand]
        L = bw.chroms(chrom)
        out = np.zeros(max(e - s, 0))
        if L is None or e <= 0:
            return out
        s2, e2 = max(0, s), min(L, e)
        if e2 > s2:
            v = bw.values(chrom, s2, e2, numpy=True)
            out[s2 - s:e2 - s] = np.abs(np.nan_to_num(v))
        return out

    def total(self):
        return sum(abs(self.bw[s].header()["sumData"]) for s in "+-")


def along(p, st, a, b):
    """genomic interval [s, e) covering offsets a..b from p in the direction of strand st"""
    return (p + a, p + b) if st == "+" else (p - b + 1, p - a + 1)


def detect_flip(cov, hosts):
    same = opp = 0.0
    for chrom, pos, st in hosts[:300]:
        s, e = along(pos, st, 0, 5000)
        same += cov.vals(chrom, s, e, st).sum()
        opp += cov.vals(chrom, s, e, "-" if st == "+" else "+").sum()
    cov.flip = opp > same
    return cov.flip


def tobool(d, cols):
    for c in cols:
        if c in d:
            d[c] = d[c].astype(str) == "True"
    return d


def kv(items):
    out = {}
    for it in items or []:
        k, v = it.split("=", 1)
        out[k] = v
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence", nargs="+", required=True, help="NAME=itss.evidence.tsv")
    ap.add_argument("--annot", nargs="*", help="NAME=annot.evidence.tsv (annotated intronic RNAs)")
    ap.add_argument("--rad21", help="plus,minus bigWig of the RAD21-line PRO-cap (HCT116)")
    ap.add_argument("--gtf", default="ref/gencode.v47.annotation.gtf.gz")
    ap.add_argument("--wave", nargs="*", help="NAME:LABEL=plus,minus ; LABEL 'ref' is the reference")
    ap.add_argument("--wave-len", type=int, default=100000)
    ap.add_argument("--wave-genes", type=int, default=1000)
    ap.add_argument("--context", nargs="*", help="NAME=context_NAME.tsv")
    ap.add_argument("--out-prefix", default="results/supp/supp")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out_prefix) or ".", exist_ok=True)
    ev = {k: tobool(pd.read_csv(v, sep="\t", low_memory=False), ["host_active", "tss_supported"])
          for k, v in kv(a.evidence).items()}

    # ---------------- S1 ----------------
    rows = []
    print("S1  cap test: % of sites passing (n)")
    lab = {"random_control": "random intronic positions", "internal_promoter": "internal promoters",
           "candidate": "intronic start sites (PRO-cap peaks)"}
    for name, d in ev.items():
        for st in ("random_control", "internal_promoter", "candidate"):
            for act in (True, False):
                x = d[(d.site_type == st) & (d.host_active == act)]
                if len(x):
                    rows.append(dict(cell_line=name, group=lab[st], host="active" if act else "silent",
                                     n=len(x), passed=int(x.tss_supported.sum()), pct=100 * x.tss_supported.mean()))
    for name, path in kv(a.annot).items():
        d = tobool(pd.read_csv(path, sep="\t", low_memory=False), ["host_active", "tss_supported"])
        d = d[d.site_type == "candidate"]
        for act in (True, False):
            x = d[d.host_active == act]
            if len(x):
                rows.append(dict(cell_line=name, group="annotated intronic RNAs", host="active" if act else "silent",
                                 n=len(x), passed=int(x.tss_supported.sum()), pct=100 * x.tss_supported.mean()))
    s1 = pd.DataFrame(rows)
    for r in s1.itertuples():
        print(f"  {r.cell_line:<7} {r.group:<40} {r.host:<7} {r.pct:6.2f}%  ({r.passed}/{r.n})")
    s1.to_csv(f"{a.out_prefix}.S1_captest.tsv", sep="\t", index=False)

    # active host start positions (for strand detection and the wave)
    def hosts_of(d):
        h = d[d.host_active].dropna(subset=["host_tss_used"]).drop_duplicates("host_gene_id")
        return h

    # ---------------- S2 ----------------
    if a.rad21 and "HCT116" in ev:
        d = ev["HCT116"]
        R = Cov(a.rad21)
        hh = hosts_of(d)
        detect_flip(R, [(r.chrom, int(r.host_tss_used), r.strand) for r in hh.itertuples()])
        x = d[(d.site_type == "candidate") & d.tss_supported & d.host_active].copy()
        res = []
        for r in x.itertuples():
            p = int(r.pos)
            v = R.vals(r.chrom, p - 2000, p + 2001, r.strand)
            s = v[1900:2101].sum()
            fl = np.concatenate([v[0:1500], v[2501:4001]]).mean()
            pv = poisson.sf(s - 1, max(fl, 1e-3) * 201) if s > 0 else 1.0
            res.append((s, s >= 5 and pv < 1e-3))
        x["rad21_reads"] = [t[0] for t in res]
        x["rad21_pass"] = [t[1] for t in res]
        x["third"] = pd.qcut(x.cap_sum.rank(method="first"), 3, labels=["weak", "medium", "strong"])
        print("\nS2  HCT116 reproducibility: CTCF-line intronic start sites (confirmed, active host) passing in the RAD21 line")
        out = []
        for th in ("weak", "medium", "strong"):
            g = x[x.third == th]
            out.append(dict(third=th, n=len(g), passed=int(g.rad21_pass.sum()), pct=100 * g.rad21_pass.mean(),
                            ctcf_median=g.cap_sum.median()))
            print(f"  {th:<7} n={len(g):>5}  passed {g.rad21_pass.mean() * 100:5.1f}%  (median CTCF-line reads {g.cap_sum.median():.0f})")
        print(f"  all     n={len(x):>5}  passed {x.rad21_pass.mean() * 100:5.1f}%")
        sh = x[x.rad21_pass]
        if len(sh) > 10:
            r_, p_ = spearmanr(sh.cap_sum, sh.rad21_reads)
            print(f"  PRO-cap reads at sites passing in both lines: Spearman r = {r_:.2f} (p = {p_:.1e}, n = {len(sh)})")
        pd.DataFrame(out).to_csv(f"{a.out_prefix}.S2_reproducibility.tsv", sep="\t", index=False)
        x[["site_id", "chrom", "strand", "pos", "cap_sum", "rad21_reads", "rad21_pass", "third"]].to_csv(
            f"{a.out_prefix}.S2_sites.tsv", sep="\t", index=False)

    # ---------------- S3 ----------------
    if a.wave:
        spec = defaultdict(dict)
        for it in a.wave:
            k, v = it.split("=", 1)
            name, label = k.split(":", 1)
            spec[name][label] = v
        s = importlib.util.spec_from_file_location("kg", os.path.join(HERE, "30_known_genes.py"))
        kg = importlib.util.module_from_spec(s)
        s.loader.exec_module(kg)
        genes, _, _ = kg.parse_gtf(a.gtf)
        nb = a.wave_len // 2000
        rows = []
        for name, files in spec.items():
            d = ev[name]
            hh = hosts_of(d)
            sel = []
            for r in hh.itertuples():
                g = genes.get(str(r.host_gene_id).split(".")[0])
                if g is None:
                    continue
                p = int(r.host_tss_used)
                span = (g["end"] - p) if r.strand == "+" else (p - g["start"])
                if span >= a.wave_len:
                    sel.append((r.chrom, p, r.strand))
            sel = sel[:a.wave_genes]
            covs = {lab: Cov(f) for lab, f in files.items()}
            for lab, c in covs.items():
                detect_flip(c, sel)
            sums = {lab: np.zeros(nb) for lab in covs}
            for chrom, p, st in sel:
                s0, e0 = along(p, st, 0, nb * 2000)
                for lab, c in covs.items():
                    v = c.vals(chrom, s0, e0, st)
                    if st == "-":
                        v = v[::-1]
                    sums[lab] += v.reshape(nb, 2000).sum(1)
            ref = sums["ref"]
            print(f"\nS3  {name}: nascent RNA / reference along {len(sel)} active host genes (2-kb bins from the host start)")
            labs = [l for l in covs if l != "ref"]
            print("  kb      " + "".join(f"{l:>10}" for l in labs))
            for i in range(nb):
                vals = [sums[l][i] / max(ref[i], 1e-9) for l in labs]
                rows += [dict(cell_line=name, label=l, kb_from=2 * i, kb_to=2 * i + 2, ratio=v) for l, v in zip(labs, vals)]
                if i % 2 == 0 or i < 10:
                    print(f"  {2 * i:>3}-{2 * i + 2:<4}" + "".join(f"{v:>10.2f}" for v in vals))
        pd.DataFrame(rows).to_csv(f"{a.out_prefix}.S3_wave.tsv", sep="\t", index=False)

    # ---------------- S6 ----------------
    if a.context:
        rows = []
        print("\nS6  stable RNA by the rise measure: % of sites above the 95th percentile of random intronic positions")
        for name, path in kv(a.context).items():
            S = pd.read_csv(path, sep="\t", low_memory=False)
            col = "RNA_skok" if "RNA_skok" in S else [c for c in S.columns if "skok" in c.lower()][0]
            thr = np.nanpercentile(S[S.grupa == "kontrole"][col], 95)
            st3 = S[S.grupa.isin(["intron_akt", "intron_cichy", "miedzygen"])].copy()
            st3["third"] = pd.qcut(st3.cap_sum.rank(method="first"), 3, labels=["weak", "medium", "strong"])
            print(f"  {name} (threshold {thr:.2f}; column {col})")
            for th in ("weak", "medium", "strong"):
                cells = []
                for g in ("intron_akt", "intron_cichy", "miedzygen"):
                    x = st3[(st3.grupa == g) & (st3.third == th)][col].dropna()
                    pct = 100 * (x > thr).mean() if len(x) else np.nan
                    rows.append(dict(cell_line=name, third=th, group=g, n=len(x), pct=pct))
                    cells.append(f"{g}: {pct:5.1f}% (n={len(x)})")
                print(f"    {th:<7} " + "   ".join(cells))
            ip = S[S.grupa == "wewn_prom"][col].dropna()
            rows.append(dict(cell_line=name, third="all", group="wewn_prom", n=len(ip), pct=100 * (ip > thr).mean()))
            print(f"    internal promoters: {100 * (ip > thr).mean():5.1f}% (n={len(ip)})")
        pd.DataFrame(rows).to_csv(f"{a.out_prefix}.S6_rise.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

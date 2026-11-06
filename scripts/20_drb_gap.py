#!/usr/bin/env python3
"""
iSTART, module A: DRB gap test for the sites from step 11 (iTSS, internal promoters, random controls).

Idea: after DRB wash-out the host polymerase wave reaches about 12-13 kb from the host start site in 10 min.
Beyond the wave the DRB10 signal falls to background (about half the Bru10 level). An independent intronic
start site produces a "burst" of DRB10 signal downstream of itself and only background upstream. Measure per site:

  gap_index = log2[(DRB_down + 0.5) / (DRB_up + 0.5)] - log2[(Bru_down + 0.5) / (Bru_up + 0.5)]

i.e. how much more the signal rises downstream of the start site after DRB than at steady state (Bru10).
Normalising to Bru10 removes the effect of local expression and mappability. Significance: one-sided Fisher
test on read counts + BH correction separately for each site type; calibration on controls (criterion G2).
At low depth the per-site test has little power, so the script also computes an empirical calibration:
threshold = 95th percentile of gap_index of random controls; for each group the fraction of sites above the
threshold, the estimated share of sites with a true effect ((f - 0.05) / 0.95), the empirical FDR and a
Mann-Whitney test against controls. Flag drb_positive_emp = gap_index above the threshold.

Input:
  --evidence   *.evidence.tsv from 11_evidence.py
  --drb-plus/--drb-minus, --bru-plus/--bru-minus   coverage bigWigs by READ strand
Output (--out-prefix):
  .drb.tsv       measures for each site
  .summary.txt   orientation, G2 calibration, 2x2 classification, aggregate profiles
"""
import argparse
import sys

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import fisher_exact, mannwhitneyu


class Cov:
    def __init__(self, plus, minus):
        self.bw = {"plus": pyBigWig.open(plus), "minus": pyBigWig.open(minus)}
        self.flip = False

    def values(self, chrom, s, e, rna_strand):
        key = "plus" if rna_strand == "+" else "minus"
        if self.flip:
            key = "minus" if key == "plus" else "plus"
        bw = self.bw[key]
        L = bw.chroms(chrom)
        out = np.zeros(max(0, e - s))
        if L is None:
            return out
        s2, e2 = max(0, s), min(L, e)
        if e2 > s2:
            v = bw.values(chrom, s2, e2, numpy=True)
            if len(v) == e2 - s2:
                out[s2 - s:e2 - s] = np.abs(np.nan_to_num(v))
        return out


def along(pos, strand, a, b):
    """Interval [a, b) measured along the direction of transcription from position pos -> genomic coordinates."""
    return (pos + a, pos + b) if strand == "+" else (pos - b + 1, pos - a + 1)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    r = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--drb-plus", required=True)
    ap.add_argument("--drb-minus", required=True)
    ap.add_argument("--bru-plus", required=True)
    ap.add_argument("--bru-minus", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--read-len", type=float, default=51)
    ap.add_argument("--min-dist", type=int, default=20000, help="min. distance from the host start site (active host)")
    ap.add_argument("--win", type=int, default=3000, help="window downstream of the start site (nt)")
    ap.add_argument("--gap", type=int, default=3000, help="window upstream of the start site (nt)")
    ap.add_argument("--margin", type=int, default=200, help="gap between the 'upstream' window and the start site (nt)")
    ap.add_argument("--min-reads", type=float, default=10, help="min. total reads in the 4 windows")
    ap.add_argument("--qmax", type=float, default=0.1)
    ap.add_argument("--flip", choices=["auto", "yes", "no"], default="auto")
    a = ap.parse_args()

    ev = pd.read_csv(a.evidence, sep="\t", low_memory=False)
    drb = Cov(a.drb_plus, a.drb_minus)
    bru = Cov(a.bru_plus, a.bru_minus)
    log = []

    # --- read orientation: first 5 kb of active hosts ---
    hosts = ev[ev.host_active.astype(str) == "True"].dropna(subset=["host_tss_used"]) \
        .drop_duplicates("host_gene_id")
    same = opp = 0.0
    for r in hosts.itertuples():
        s, e = along(int(r.host_tss_used), r.strand, 0, 5000)
        same += bru.values(r.chrom, s, e, r.strand).sum()
        opp += bru.values(r.chrom, s, e, "-" if r.strand == "+" else "+").sum()
    flip = (opp > same) if a.flip == "auto" else (a.flip == "yes")
    drb.flip = bru.flip = flip
    log.append(f"Orientation (Bru10, 5 kb downstream of the start of {len(hosts)} active hosts): RNA strand = "
               f"{same:.0f}, opposite strand = {opp:.0f} -> {'reads on the opposite strand (swapping)' if flip else 'reads on the RNA strand'}")

    # --- measures for each site ---
    rows = []
    for r in ev.itertuples(index=False):
        active = str(r.host_active) == "True"
        upstream = str(getattr(r, "upstream_of_used_host_tss", "False")) == "True"
        d = r.dist_used_tss
        elig = (not upstream) and ((not active) or (pd.notna(d) and d >= a.min_dist))
        ds, de = along(r.pos, r.strand, 0, a.win)
        us, ue = along(r.pos, r.strand, -(a.margin + a.gap), -a.margin)
        D_dn = drb.values(r.chrom, ds, de, r.strand).sum() / a.read_len
        D_up = drb.values(r.chrom, us, ue, r.strand).sum() / a.read_len
        B_dn = bru.values(r.chrom, ds, de, r.strand).sum() / a.read_len
        B_up = bru.values(r.chrom, us, ue, r.strand).sum() / a.read_len
        tot = D_dn + D_up + B_dn + B_up
        gi = np.log2((D_dn + .5) / (D_up + .5)) - np.log2((B_dn + .5) / (B_up + .5))
        testable = elig and tot >= a.min_reads
        p = fisher_exact([[round(D_dn), round(D_up)], [round(B_dn), round(B_up)]], alternative="greater")[1] \
            if testable else np.nan
        rows.append(dict(site_id=r.site_id, site_type=r.site_type, host_active=active,
                         tss_supported=str(r.tss_supported) == "True",
                         evidence_class=getattr(r, "evidence_class", "-"),
                         dist_used_tss=d, eligible=elig, testable=testable,
                         drb_down=D_dn, drb_up=D_up, bru_down=B_dn, bru_up=B_up, gap_index=gi, pval=p))
    df = pd.DataFrame(rows)
    df["qval"] = np.nan
    for t, idx in df[df.testable].groupby("site_type").groups.items():
        df.loc[idx, "qval"] = bh(df.loc[idx, "pval"].values)
    df["drb_positive"] = df.testable & (df.qval < a.qmax) & (df.gap_index > 0)
    ctrl = df[(df.site_type == "random_control") & df.host_active & df.testable].gap_index
    thr = float(np.percentile(ctrl, 95)) if len(ctrl) else np.inf
    df["drb_positive_emp"] = df.testable & (df.gap_index > thr)
    full = ev.merge(df[["site_id", "eligible", "testable", "drb_down", "drb_up", "bru_down", "bru_up",
                        "gap_index", "pval", "qval", "drb_positive", "drb_positive_emp"]], on="site_id", how="left")
    full.to_csv(f"{a.out_prefix}.drb.tsv", sep="\t", index=False)

    # --- calibration (G2) ---
    def line(name, m):
        sub = df[m]
        t = sub[sub.testable]
        med = t.gap_index.median() if len(t) else np.nan
        return (f"  {name:<52} eligible {int(sub.eligible.sum()):>5}, testable {len(t):>5}, "
                f"positive {int(t.drb_positive.sum()):>5} ({(t.drb_positive.mean()*100 if len(t) else 0):5.1f}%), "
                f"median gap_index {med:5.2f}")

    cand = df.site_type == "candidate"
    log.append("")
    log.append(f"Calibration (G2): fraction positive among testable sites (q < {a.qmax}, gap_index > 0)")
    log.append(line("random controls, active host", (df.site_type == "random_control") & df.host_active))
    log.append(line("capped internal promoters, active host",
                    (df.site_type == "internal_promoter") & df.host_active & df.tss_supported))
    log.append(line("capped candidates, inactive host", cand & ~df.host_active & df.tss_supported))
    log.append(line("capped candidates, active host", cand & df.host_active & df.tss_supported))
    log.append(line("uncapped candidates, active host", cand & df.host_active & ~df.tss_supported))

    groups = [
        ("random controls, active host", (df.site_type == "random_control") & df.host_active),
        ("capped internal promoters, active host",
         (df.site_type == "internal_promoter") & df.host_active & df.tss_supported),
        ("capped candidates, inactive host", cand & ~df.host_active & df.tss_supported),
        ("capped candidates, active host", cand & df.host_active & df.tss_supported),
        ("uncapped candidates, active host", cand & df.host_active & ~df.tss_supported),
    ]
    log.append("")
    log.append(f"Empirical calibration: threshold = 95th percentile of gap_index of random controls = {thr:.2f}")
    log.append("  group                                                 n   >thr    share with effect*  emp. FDR   p (M-W vs controls)")
    for name, m in groups:
        g = df[m & df.testable].gap_index
        if not len(g):
            continue
        f = float((g > thr).mean())
        pi1 = max(0.0, (f - 0.05) / 0.95)
        is_ctrl = name.startswith("random")
        fdr = min(1.0, 0.05 / f) if (f > 0 and not is_ctrl) else np.nan
        pmw = mannwhitneyu(g, ctrl, alternative="greater").pvalue if (len(ctrl) and not is_ctrl) else np.nan
        log.append(f"  {name:<52} {len(g):>5}  {f*100:5.1f}%   {pi1*100:5.1f}%             {fdr:5.2f}      {pmw:.1e}")
    log.append("  * estimate: (fraction above threshold - 5%) / 95%; emp. FDR = 5% / fraction above threshold")

    # --- 2x2 classification (candidates, active host, testable) ---
    t = df[cand & df.host_active & df.testable]
    log.append("")
    log.append("2x2 classification (candidates with an active host, testable; empirical flag):")
    for cap_ in (True, False):
        for drb_ in (True, False):
            n = int(((t.tss_supported == cap_) & (t.drb_positive_emp == drb_)).sum())
            log.append(f"  cap {'yes' if cap_ else 'no'} / signal ahead of the wave {'yes' if drb_ else 'no'}: {n}")

    # --- aggregate DRB10/Bru10 profile around sites (from -10 to +10 kb) ---
    def meta(sel):
        sub = ev[ev.site_id.isin(df[sel].site_id)]
        Dsum = np.zeros(20)
        Bsum = np.zeros(20)
        for r in sub.itertuples():
            s, e = along(r.pos, r.strand, -10000, 10000)
            dv = drb.values(r.chrom, s, e, r.strand)
            bv = bru.values(r.chrom, s, e, r.strand)
            if r.strand == "-":
                dv, bv = dv[::-1], bv[::-1]
            Dsum += dv.reshape(20, 1000).sum(1)
            Bsum += bv.reshape(20, 1000).sum(1)
        return Dsum / np.maximum(Bsum, 1e-9), len(sub)

    # group labels are data (read back by 50_figure_data.py): iTSS_akt / iTSS_nieakt = capped iTSS with an
    # active / inactive host, wewn_prom = internal promoters, kontrole = random controls
    prof = [("iTSS_akt", cand & df.host_active & df.tss_supported & df.eligible),
            ("iTSS_nieakt", cand & ~df.host_active & df.tss_supported & df.eligible),
            ("wewn_prom", (df.site_type == "internal_promoter") & df.host_active & df.tss_supported & df.eligible),
            ("kontrole", (df.site_type == "random_control") & df.host_active & df.eligible)]
    res = [(n,) + meta(m) for n, m in prof]
    log.append("")
    # 50_figure_data.py finds this block by the line starting with "Profil" and its "name n=N" entries
    log.append("Profile DRB10/Bru10 around sites (kb relative to the site): " +
               ", ".join(f"{n} n={k}" for n, _, k in res))
    log.append("  kb\t" + "\t".join(n for n, _, _ in res))
    for i in range(20):
        log.append(f"  {i - 10:+d}..{i - 9:+d}\t" + "\t".join(f"{v[i]:.2f}" for _, v, _ in res))
    log.append("  step (mean +1..+6 kb / mean -6..-1 kb):\t" +
               "\t".join(f"{v[11:16].mean() / max(v[4:9].mean(), 1e-9):.2f}" for _, v, _ in res))

    txt = "\n".join(log)
    open(f"{a.out_prefix}.summary.txt", "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()

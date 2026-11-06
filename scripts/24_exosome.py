#!/usr/bin/env python3
"""
iSTART: is RNA from intronic start sites degraded by the nuclear exosome?

For one knockout (KO) and its control (RNA-seq bigWigs, plus/minus strand, replicates already merged):
 1. Positive control of knockout efficiency: antisense RNAs just upstream of active host promoters
    (PROMPTs, window -2000..-200 nt, opposite strand) are classic exosome substrates.
    Measure: log2(KO/control) in the PROMPT window minus log2(KO/control) over the first 3 kb of the gene.
 2. Test: change in the RNA rise after the site, Δ = [log2(after/before)]_KO - [log2(after/before)]_control,
    after = 0..+3 kb, before = -3.2..-0.2 kb, for the site groups of the description table (step 21).
    Threshold = 95th percentile of Δ in random controls; Mann-Whitney test against the controls.
"""
import argparse

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import mannwhitneyu


class Cov:
    def __init__(self, plus, minus):
        self.bw = {"+": pyBigWig.open(plus), "-": pyBigWig.open(minus)}
        self.flip = False

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


def win(pos, strand, a, b):
    return (pos + a, pos + b) if strand == "+" else (pos - b + 1, pos - a + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True, help="itss.describe.tsv (step 21)")
    ap.add_argument("--ko-plus", required=True)
    ap.add_argument("--ko-minus", required=True)
    ap.add_argument("--ctrl-plus", required=True)
    ap.add_argument("--ctrl-minus", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--eps", type=float, default=0.01)
    a = ap.parse_args()
    eps = a.eps

    s = pd.read_csv(a.sites, sep="\t", low_memory=False)
    for col in ("host_active", "tss_supported"):
        s[col] = s[col].astype(str) == "True"
    ko, ct = Cov(a.ko_plus, a.ko_minus), Cov(a.ctrl_plus, a.ctrl_minus)

    # orientation (control): first 5 kb downstream of host start sites
    hp = s[s.site_type == "host_promoter"]
    same = sum(ct.mean(r.chrom, *win(int(r.pos), r.strand, 0, 5000), r.strand) for r in hp.head(2000).itertuples())
    opp = sum(ct.mean(r.chrom, *win(int(r.pos), r.strand, 0, 5000), "-" if r.strand == "+" else "+")
              for r in hp.head(2000).itertuples())
    ko.flip = ct.flip = opp > same
    log = [f"[{a.name}] orientation: RNA strand {same:.1f}, opposite {opp:.1f} -> {'swapping' if ko.flip else 'no swap'}"]

    # 1. positive control: PROMPTs
    pr = []
    for r in hp.itertuples():
        p, st = int(r.pos), r.strand
        anti = "-" if st == "+" else "+"
        ps, pe = win(p, st, -2000, -200)
        gs, ge = win(p, st, 0, 3000)
        pk, pc = ko.mean(r.chrom, ps, pe, anti), ct.mean(r.chrom, ps, pe, anti)
        gk, gc = ko.mean(r.chrom, gs, ge, st), ct.mean(r.chrom, gs, ge, st)
        pr.append((np.log2((pk + eps) / (pc + eps)), np.log2((gk + eps) / (gc + eps))))
    pr = np.array(pr)
    d_prompt = pr[:, 0] - pr[:, 1]
    # this line is parsed by 50_figure_data.py (prompt_check_txt): keep its wording in step with that regex
    log.append(f"[{a.name}] positive control (PROMPT, {len(pr)} promoters): log2FC PROMPT median {np.median(pr[:, 0]):+.2f}, "
               f"gene start {np.median(pr[:, 1]):+.2f}, difference {np.median(d_prompt):+.2f}; "
               f"PROMPT higher than gene in {np.mean(d_prompt > 0) * 100:.0f}% of promoters")

    # 2. change in the RNA rise after the site
    rows = []
    sel = s.site_type.isin(["candidate", "internal_promoter", "random_control"])
    for r in s[sel].itertuples():
        p, st = int(r.pos), r.strand
        ds, de = win(p, st, 0, 3000)
        us, ue = win(p, st, -3200, -200)
        kd, ku = ko.mean(r.chrom, ds, de, st), ko.mean(r.chrom, us, ue, st)
        cd, cu = ct.mean(r.chrom, ds, de, st), ct.mean(r.chrom, us, ue, st)
        rows.append(dict(site_id=r.site_id, delta=np.log2((kd + eps) / (ku + eps)) - np.log2((cd + eps) / (cu + eps)),
                         lfc_down=np.log2((kd + eps) / (cd + eps)), covered=(kd + ku + cd + cu) > 4 * eps))
    d = s[sel].merge(pd.DataFrame(rows), on="site_id")
    d.to_csv(f"{a.out_prefix}.{a.name}.tsv", sep="\t", index=False)
    d = d[d.covered]
    cand = d.site_type == "candidate"
    ctrl = d[(d.site_type == "random_control") & d.host_active].delta
    thr = np.percentile(ctrl, 95)
    groups = [("random controls", (d.site_type == "random_control") & d.host_active),
              ("capped internal promoters", (d.site_type == "internal_promoter") & d.host_active & d.tss_supported),
              ("capped iTSS, active host", cand & d.host_active & d.tss_supported),
              ("  of which promoter-like", cand & d.host_active & d.tss_supported & (d.typ == "promotor-podobny")),
              ("  of which enhancer-like", cand & d.host_active & d.tss_supported & (d.typ == "enhancer-podobny")),
              ("capped iTSS, inactive host", cand & ~d.host_active & d.tss_supported),
              ("uncapped iTSS", cand & ~d.tss_supported)]
    log.append(f"[{a.name}] change in the RNA rise after the site upon knockout "
               f"(Δ; threshold = 95th percentile of controls = {thr:.2f}):")
    log.append(f"  {'group':<42}{'n':>6}{'median Δ':>11}{'Δ>thr':>9}{'median log2FC down':>19}{'p (M-W)':>10}")
    for name, m in groups:
        g = d[m]
        if not len(g):
            continue
        p = mannwhitneyu(g.delta, ctrl, alternative="greater").pvalue if not name.startswith("random") else np.nan
        log.append(f"  {name:<42}{len(g):>6}{g.delta.median():>11.2f}{(g.delta > thr).mean() * 100:>8.0f}%"
                   f"{g.lfc_down.median():>19.2f}{p:>10.1e}")
    txt = "\n".join(log)
    open(f"{a.out_prefix}.{a.name}.txt", "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Step 39: Fig. 5a with the 3-kb rule, from the tables results/context_{K562,HCT116}.tsv.
Value per initiation event = RNA_nadw (29_context.py already divides by cap_sum). First a comparison with the original table
results/final.fig5_robust.tsv (as a check), then the result after the rule: sites in introns of active
hosts and internal promoters only >= 3 kb upstream of the next host exon.
Run from the project directory:
    python scripts/39_fig5a.py > results/exon_dist/fig5a.txt 2>&1
"""
import numpy as np
import pandas as pd

FAR = 3000
rng = np.random.default_rng(1)
pd.set_option("display.width", 220)


def wilson(k, n):
    if n == 0:
        return np.nan, np.nan
    z, p = 1.96, k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return 100 * (c - h), 100 * (c + h)


def boot_med(v, nb=1000):
    v = np.asarray(v, float)
    if len(v) < 3:
        return np.nan, np.nan
    m = np.median(rng.choice(v, (nb, len(v)), replace=True), axis=1)
    return np.percentile(m, 2.5), np.percentile(m, 97.5)


def summ(v):
    k, n = int((v > 0).sum()), len(v)
    lo, hi = wilson(k, n)
    mlo, mhi = boot_med(v)
    return f"{n:>6} {100 * k / max(n, 1):5.1f}% [{lo:4.1f}-{hi:4.1f}] {1000 * np.median(v) if n else np.nan:8.3f} [{1000 * mlo:7.3f},{1000 * mhi:7.3f}]"


print("======== Original table for Fig. 5 (stable RNA) ========")
o = pd.read_csv("results/final.fig5_robust.tsv", sep="\t")
print("thirds:", sorted(o.third.astype(str).unique()))
print(o[o.measure.astype(str).str.contains("stable", case=False)].to_string(index=False))

for line, base in (("K562", "results/v2"), ("HCT116", "results/hct116")):
    print(f"\n======== {line} ========")
    c = pd.read_csv(f"results/context_{line}.tsv", sep="\t", low_memory=False)
    print("groups:", dict(c.grupa.astype(str).value_counts()))
    d = pd.read_csv(f"{base}/itss.describe.tsv", sep="\t", low_memory=False,
                    usecols=["site_id", "strand", "pos", "intron_start", "intron_end"])
    for col in ("pos", "intron_start", "intron_end"):
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d["dd"] = np.where(d.strand == "+", d.intron_end - d.pos, d.pos - d.intron_start)
    c = c.merge(d[["site_id", "dd"]], on="site_id", how="left")
    cap = pd.to_numeric(c.cap_sum, errors="coerce")
    c["v"] = np.where(cap > 0, pd.to_numeric(c.RNA_nadw, errors="coerce"), np.nan)
    print(f"distance to the exon matched for {c.dd.notna().sum()} of {len(c)} rows")
    print(f"\n  {'group':<16} {'n':>6} {'% > 0 [95% CI]':>20} {'med x1000 [95% CI]':>26} | {'n (3 kb)':>6} {'% > 0 [95% CI]':>20} {'med x1000 [95% CI]':>26}")
    for g, x in c.groupby(c.grupa.astype(str)):
        x = x[np.isfinite(x.v)]
        gl = g.lower()
        rule = (("akt" in gl and "nieakt" not in gl) or "wewn" in gl or "prom" in gl or "kontrol" in gl) and x.dd.notna().any()
        y = x[x.dd >= FAR] if rule else x
        print(f"  {g:<16} {summ(x.v.values)} | {summ(y.v.values)}{'' if rule else '   (rule not applied)'}")
print("\n==== end ====")

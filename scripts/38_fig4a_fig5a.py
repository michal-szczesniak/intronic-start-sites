#!/usr/bin/env python
"""Step 38: Fig. 4a (MTREX) and Fig. 5a (stable RNA per initiation event) with the 3-kb rule.
First reproduces the original numbers (as a check), then gives the values after the rule: sites in active
hosts (intronic start sites, internal promoters, random positions) only >= 3 kb upstream of the next exon.
Sites in silent hosts and between genes are not subject to the rule.
Run from the project directory:
    python scripts/38_fig4a_fig5a.py > results/exon_dist/fig4a_fig5a.txt 2>&1
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, norm

FAR = 3000
rng = np.random.default_rng(1)
pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 40)


def b(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


def sec(t):
    print(f"\n{'=' * 8} {t} {'=' * 8}", flush=True)


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


def dist_down(d):
    for c in ("pos", "intron_start", "intron_end"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return np.where(d.strand == "+", d.intron_end - d.pos, d.pos - d.intron_start)


# ---------- Fig. 4a: MTREX ----------
sec("Fig. 4a: MTREX (sites with coverage in all four windows, as in the original)")
e = pd.read_csv("results/v2/exosome/exo.MTREX.tsv", sep="\t", low_memory=False)
e = e[b(e.covered)] if "covered" in e else e
e["dd"] = dist_down(e)
e["v"] = pd.to_numeric(e["delta"], errors="coerce")
st, act, sup = e.site_type.astype(str), b(e.host_active), b(e.tss_supported)
typ = e["typ"].astype(str) if "typ" in e else pd.Series("", index=e.index)
for lab, far_rule in (("original (all distances)", False), ("3-kb rule", True)):
    far = (e.dd >= FAR) if far_rule else pd.Series(True, index=e.index)
    rnd = e[(st == "random_control") & act & far].v.dropna()
    thr = np.percentile(rnd, 95)
    print(f"\n  {lab}: threshold = {thr:.3f}")
    groups = [("Random intronic positions", (st == "random_control") & act & far),
              ("Internal promoters", (st == "internal_promoter") & sup & act & far),
              ("Intronic start sites, active host", (st == "candidate") & sup & act & far),
              ("… promoter-like", (st == "candidate") & sup & act & far & typ.str.startswith("promotor")),
              ("… enhancer-like", (st == "candidate") & sup & act & far & typ.str.startswith("enhancer")),
              ("Intronic start sites, silent host", (st == "candidate") & sup & ~act)]
    for g, m in groups:
        v = e[m].v.dropna()
        if g.startswith("Random"):
            print(f"    {g:<36} n {len(v):>5}  {100 * (v > thr).mean():5.1f}%")
            continue
        p2 = mannwhitneyu(v, rnd).pvalue if len(v) > 2 else np.nan
        p1 = mannwhitneyu(v, rnd, alternative="greater").pvalue if len(v) > 2 else np.nan
        print(f"    {g:<36} n {len(v):>5}  {100 * (v > thr).mean():5.1f}%   p two-sided = {p2:.2g}, one-sided = {p1:.2g}")

# ---------- Fig. 5a: context tables ----------
for line, base in (("K562", "results/v2"), ("HCT116", "results/hct116")):
    sec(f"Fig. 5a: {line}")
    path = f"results/context_{line}.tsv"
    if not os.path.exists(path):
        print(f"[missing file] {path}")
        continue
    c = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"{path}: {len(c)} rows; columns: {', '.join(map(str, c.columns))}")
    gcol = next((x for x in c.columns if c[x].astype(str).str.contains("between|intronic|internal|międzygen|wewn", case=False).mean() > 0.3), None)
    cand_v = [x for x in c.columns if any(k in str(x).lower() for k in ("per_init", "per_cap", "stable_per", "rna_per", "per_read", "na_odczyt"))]
    print(f"group column: {gcol}; candidate columns for the per-initiation-event value: {cand_v}")
    if gcol:
        print(f"groups: {dict(c[gcol].astype(str).value_counts())}")
    if not (gcol and cand_v):
        print(c.head(3).to_string())
        continue
    vcol = cand_v[0]
    if {"intron_start", "intron_end", "pos", "strand"} <= set(c.columns):
        c["dd"] = dist_down(c)
    else:
        d = pd.read_csv(f"{base}/itss.describe.tsv", sep="\t", low_memory=False,
                        usecols=["site_id", "chrom", "strand", "pos", "intron_start", "intron_end"])
        d["dd"] = dist_down(d)
        key = "site_id" if "site_id" in c else None
        if key:
            c = c.merge(d[["site_id", "dd"]], on="site_id", how="left")
        else:
            pc = next((x for x in ("pos", "position", "summit", "mid") if x in c), None)
            c = c.merge(d[["chrom", "strand", "pos", "dd"]].rename(columns={"pos": pc}), on=["chrom", "strand", pc], how="left")
    c["v"] = pd.to_numeric(c[vcol], errors="coerce")
    print(f"value: {vcol}; distance matched for {c.dd.notna().sum()} of {len(c)} rows")
    print(f"\n  {'group':<34} {'n':>6} {'% > 0 [95% CI]':>22} {'median x1000 [95% CI]':>30} | {'n (3 kb)':>8} {'% > 0 [95% CI]':>22} {'median x1000 [95% CI]':>30}")
    for g, x in c.groupby(c[gcol].astype(str)):
        x = x[x.v.notna()]
        active_like = ("active" in g.lower() or "aktyw" in g.lower() or "internal" in g.lower() or "wewn" in g.lower()) \
            and not ("silent" in g.lower() or "nieaktyw" in g.lower() or "wycisz" in g.lower())
        y = x[x.dd >= FAR] if active_like and x.dd.notna().any() else x
        out = []
        for z in (x, y):
            k, n = int((z.v > 0).sum()), len(z)
            lo, hi = wilson(k, n)
            mlo, mhi = boot_med(z.v)
            out.append(f"{n:>6} {100 * k / max(n, 1):>6.1f}% [{lo:4.1f}-{hi:4.1f}] {1000 * z.v.median():>10.3f} [{1000 * mlo:7.3f},{1000 * mhi:7.3f}]")
        print(f"  {g[:34]:<34} {out[0]} | {out[1]}{'' if active_like else '   (rule not applied)'}")
print("\n==== end ====")

#!/usr/bin/env python
"""Step 40: inclusion (PSI) of the next host exon and intronic start sites, in both cell lines.
Step 36 missed the HCT116 junction file data/hct116/RNA.junctions.tsv (built in jobs/hct116.sh)."""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

JUNCTIONS = {"K562": "data/rnaseq/K562_total.junctions.tsv", "HCT116": "data/hct116/RNA.junctions.tsv"}
sites = pd.read_csv("results/exon_dist/sites_exon_dist.tsv", sep="\t", low_memory=False)


def read_junctions(path):
    j = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False)
    if not str(j.iloc[0, 1]).replace(".", "").isdigit():
        j = j.iloc[1:].reset_index(drop=True)
    count_col = 3 if j.shape[1] == 4 else 6
    out = pd.DataFrame({"chrom": j.iloc[:, 0].astype(str), "js": pd.to_numeric(j.iloc[:, 1], errors="coerce"),
                        "je": pd.to_numeric(j.iloc[:, 2], errors="coerce"),
                        "n": pd.to_numeric(j.iloc[:, count_col], errors="coerce")})
    return out.dropna()


for line, path in JUNCTIONS.items():
    j = read_junctions(path)
    x = sites[(sites.line == line) & (sites.ex_last == False) & sites.ex_s.notna() & (sites.dist_down < 10000)].copy()
    hb = j[(j.chrom == "chr11") & (j.js > 5_240_000) & (j.je < 5_260_000)].n.sum()
    print(f"\n==== {line}: {path}, {len(j)} junctions, reads in the HBG1/2 locus {hb:.0f} (high in K562 only)")
    ends, starts = set(zip(j.chrom, j.je.astype(np.int64))), set(zip(j.chrom, j.js.astype(np.int64)))
    be = max((sum((c, int(s) - 1 + o) in ends for c, s in zip(x.chrom, x.ex_s)), o) for o in (-1, 0, 1))
    bs = max((sum((c, int(e) + 1 + o) in starts for c, e in zip(x.chrom, x.ex_e)), o) for o in (-1, 0, 1))
    off_e, off_s = be[1], bs[1]
    print(f"coordinate offsets: ends {off_e:+d}, starts {off_s:+d} (matches {be[0] + bs[0]})")
    idx = {c: (g.js.values, g.je.values, g.n.values) for c, g in j.groupby("chrom")}

    def psi(c, s, e):
        if c not in idx:
            return np.nan
        js, je, n = idx[c]
        inc = (n[je == s - 1 + off_e].sum() + n[js == e + 1 + off_s].sum()) / 2
        skip = n[(js - off_s < s) & (je - off_e > e)].sum()
        return inc / (inc + skip) if inc + skip >= 10 else np.nan

    x["psi"] = [psi(c, int(s), int(e)) for c, s, e in zip(x.chrom, x.ex_s, x.ex_e)]
    x["cls"] = np.select([x.psi >= 0.9, x.psi < 0.9], ["constitutive", "alternative"], "no data")
    print(f"PSI computed for {x.psi.notna().sum()} of {len(x)} sites < 10 kb from a non-last exon")
    for k, g in x.groupby("cls"):
        a, r = g[g.grp == "start"], g[g.grp == "losowa"]
        da, dr = a.drb.dropna().astype(bool), r.drb.dropna().astype(bool)
        print(f"  {k:<13} starts {len(a):>4}, random {len(r):>5}; median capped reads (starts) {a.cap_sum.median():>5.0f}; "
              f"DRB response {100 * da.mean() if len(da) else float('nan'):5.1f}% / {100 * dr.mean() if len(dr) else float('nan'):5.1f}%")
    alt = x[(x.grp == "start") & (x.cls == "alternative")].cap_sum
    con = x[(x.grp == "start") & (x.cls == "constitutive")].cap_sum
    if len(alt) > 5 and len(con) > 5:
        print(f"  capped reads of starts: alternative {alt.median():.0f} (n = {len(alt)}) vs constitutive {con.median():.0f} "
              f"(n = {len(con)}), Mann-Whitney p = {mannwhitneyu(alt, con).pvalue:.2g}")
print("\n==== done ====")

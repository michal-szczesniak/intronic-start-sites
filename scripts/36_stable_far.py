#!/usr/bin/env python
"""Step 36: stable RNA without the host exon artefact.

Step 35 showed that the 0-3 kb window downstream of the site overlaps the next host exon when that exon lies
closer than 3 kb, and mRNA reads then look like stable RNA (random positions 1-2 kb from the exon:
approx. 27% "stable"). Here we recompute the stable-RNA-based results only for sites at least
3 kb upstream of the next exon of the host reference transcript, with the threshold taken from random positions
that satisfy the same rule. Numbers for: Fig. 3a, Fig. 4b, Fig. 2d (stable RNA), Fig. 5a.
At the end: inclusion of the downstream exon (PSI) from the junction files, with automatic format detection.

Run from the project directory (without srun if you are already on a node):
    python scripts/36_stable_far.py > results/exon_dist/stable_far.txt 2>&1
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact, mannwhitneyu

LINES = {"K562": "results/v2", "HCT116": "results/hct116"}
FAR = 3000
pd.set_option("display.width", 200)


def b(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


def sec(t):
    print(f"\n{'=' * 8} {t} {'=' * 8}", flush=True)


def pc(x):
    return f"{100 * np.mean(x):5.1f}% (n = {len(x)})" if len(x) else "  -  (n = 0)"


def mh(tabs):
    """Mantel-Haenszel odds ratio and test (with continuity correction) for a list of 2x2 tables."""
    num = den = sa = sea = sva = 0.0
    for a, bb, c, d in tabs:
        n = a + bb + c + d
        if n < 2:
            continue
        num += a * d / n
        den += bb * c / n
        sa += a
        sea += (a + bb) * (a + c) / n
        sva += (a + bb) * (c + d) * (a + c) * (bb + d) / (n * n * (n - 1))
    orr = num / den if den else np.nan
    p = chi2.sf((abs(sa - sea) - 0.5) ** 2 / sva, 1) if sva else np.nan
    return orr, p


supp = pd.read_csv("results/final.supp_table1_intronic_starts.tsv", sep="\t", low_memory=False)

for line, base in LINES.items():
    d = pd.read_csv(f"results/repeats.{line}.tsv", sep="\t", low_memory=False)
    u = pd.read_csv(f"{base}/u1pas.tsv", sep="\t", usecols=["site_id", "protected"], low_memory=False)
    d = d.merge(u.drop_duplicates("site_id"), on="site_id", how="left")
    st, act, sup = d.site_type.astype(str), b(d.host_active), b(d.tss_supported)
    for c in ("pos", "intron_start", "intron_end", "cap_sum", "rna_step"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["dist_down"] = np.where(d.strand == "+", d.intron_end - d.pos, d.pos - d.intron_start)
    d["far"] = d.dist_down >= FAR
    itss = d[(st == "candidate") & sup & act].copy()
    ip = d[(st == "internal_promoter") & sup & act].copy()
    rnd = d[(st == "random_control") & act].copy()
    thr_all = np.nanpercentile(rnd.rna_step, 95)
    thr_far = np.nanpercentile(rnd.loc[rnd.far, "rna_step"], 95)
    sec(f"{line}: stable RNA thresholds (95th percentile of rna_step at random positions)")
    agree = (itss.rna_step > thr_all) == b(itss.stabilne_RNA)
    print(f"all distances: {thr_all:.3f} (agreement with the stabilne_RNA column: {100 * agree.mean():.1f}%); "
          f">= 3 kb from the exon: {thr_far:.3f}")
    for name, x in (("intronic start sites", itss), ("internal promoters", ip), ("random positions", rnd)):
        print(f"  {name:<22} closer than 3 kb to the exon: {100 * (~x.far).mean():5.1f}% of {len(x)}")
    for x in (itss, ip, rnd):
        x["stab_all"] = x.rna_step > thr_all
        x["stab_far"] = x.rna_step > thr_far

    # Fig. 3a: thirds of initiation strength with shared boundaries
    e1, e2 = np.nanpercentile(pd.concat([itss.cap_sum, ip.cap_sum]), [100 / 3, 200 / 3])
    for x in (itss, ip):
        x["third"] = pd.cut(x.cap_sum, [-np.inf, e1, e2, np.inf], labels=["weak", "medium", "strong"])
    sec(f"{line}: Fig. 3a - stable RNA by initiation-strength third (boundaries {e1:.0f}, {e2:.0f})")
    print(f"  {'third':<8} {'start sites, all':>24} {'start sites, >= 3 kb':>24} {'int. prom., all':>24} {'int. prom., >= 3 kb':>24}")
    for t in ["weak", "medium", "strong"]:
        a, p = itss[itss.third == t], ip[ip.third == t]
        print(f"  {t:<8} {pc(a.stab_all):>24} {pc(a.loc[a.far, 'stab_far']):>24} {pc(p.stab_all):>24} {pc(p.loc[p.far, 'stab_far']):>24}")
    print(f"  random positions >= 3 kb: {pc(rnd.loc[rnd.far, 'stab_far'])} (approx. 5% by definition)")

    # Fig. 4b: donor site before the first PAS
    sec(f"{line}: Fig. 4b - protected sites vs. stable RNA (start sites >= 3 kb from the exon)")
    f = itss[itss.far]
    pr = b(f.protected)
    for lab, x, col in (("all distances", itss, "stab_all"), (">= 3 kb", f, "stab_far")):
        prx = b(x.protected)
        t = [[int((x[col] & prx).sum()), int((~x[col] & prx).sum())], [int((x[col] & ~prx).sum()), int((~x[col] & ~prx).sum())]]
        print(f"  {lab:<20} protected {pc(x.loc[prx, col])}  unprotected {pc(x.loc[~prx, col])}  Fisher p = {fisher_exact(t)[1]:.2g}")
    print(f"  fraction protected: all {100 * b(itss.protected).mean():.1f}%, >= 3 kb {100 * pr.mean():.1f}%, "
          f"< 3 kb {100 * b(itss.loc[~itss.far, 'protected']).mean():.1f}%")

    # Fig. 2d: odds ratio for stable RNA, shared vs. line-specific, within strength thirds (CMH)
    s = supp[supp.cell_line == line][["site_id", "shared_with_other_cell_line", "stable_rna_per_initiation"]]
    j = itss.merge(s, on="site_id", how="left")
    if j.shared_with_other_cell_line.isna().mean() > 0.5:
        print("  [note] poor site_id match to Supp. Table 1; skipping Fig. 2d and 5a")
    else:
        j["shared"] = b(j.shared_with_other_cell_line)
        sec(f"{line}: Fig. 2d - stable RNA, shared vs. line-specific start sites (CMH within strength thirds)")
        for lab, x, col in (("all distances", j, "stab_all"), (">= 3 kb", j[j.far], "stab_far")):
            tabs = []
            for t in ["weak", "medium", "strong"]:
                g = x[x.third == t]
                tabs.append((int((g.shared & g[col]).sum()), int((g.shared & ~g[col]).sum()),
                             int((~g.shared & g[col]).sum()), int((~g.shared & ~g[col]).sum())))
            orr, p = mh(tabs)
            print(f"  {lab:<20} OR = {orr:.2f}, p = {p:.2g}; stable: shared {pc(x.loc[x.shared, col])}, specific {pc(x.loc[~x.shared, col])}")
        # Fig. 5a: fraction with a positive excess and median per read (start sites in active hosts)
        sec(f"{line}: Fig. 5a - stable RNA per initiation event (start sites in active hosts)")
        v = pd.to_numeric(j.stable_rna_per_initiation, errors="coerce")
        for lab, m in (("all distances", v.notna()), (">= 3 kb", v.notna() & j.far), ("< 3 kb", v.notna() & ~j.far)):
            print(f"  {lab:<20} positive excess {pc(v[m] > 0)}; median x1000 = {1000 * v[m].median():.3f}")

    # Sentences in Results: two numbers from the text
    sec(f"{line}: numbers for the text")
    print(f"  stable RNA, start sites (>= 3 kb): {pc(itss.loc[itss.far, 'stab_far'])}; internal promoters (>= 3 kb): {pc(ip.loc[ip.far, 'stab_far'])}")
    print(f"  in the strongest third: start sites {pc(itss.loc[itss.far & (itss.third == 'strong'), 'stab_far'])}, "
          f"internal promoters {pc(ip.loc[ip.far & (ip.third == 'strong'), 'stab_far'])}")


# ---------- downstream exon inclusion (PSI): detecting the junction file format ----------
sec("Junction files: format")
files = sorted(os.path.join("data/rnaseq", f) for f in os.listdir("data/rnaseq") if f.endswith("junctions.tsv"))
for p in files:
    with open(p) as fh:
        head = [next(fh, "").rstrip("\n") for _ in range(3)]
    print(f"--- {p}")
    for ln in head:
        print("   " + ln[:200])
S = "results/exon_dist/sites_exon_dist.tsv"
if os.path.exists(S) and files:
    sites = pd.read_csv(S, sep="\t", low_memory=False)

    def read_any(p):
        j = pd.read_csv(p, sep="\t", header=None, comment="#", low_memory=False)
        if not str(j.iloc[0, 1]).replace(".", "").isdigit():  # header row
            j = j.iloc[1:].reset_index(drop=True)
        num = [c for c in j.columns[1:] if pd.to_numeric(j[c], errors="coerce").notna().mean() > 0.95]
        if len(num) < 3:
            raise ValueError(f"too few numeric columns: {list(j.columns)}")
        out = pd.DataFrame({"chrom": j.iloc[:, 0].astype(str), "js": pd.to_numeric(j[num[0]]), "je": pd.to_numeric(j[num[1]]),
                            "uniq": pd.to_numeric(j[num[-1]] if j.shape[1] == 4 else j[6])})
        return out.dropna()

    sec("Downstream exon inclusion (PSI): sites < 10 kb from the exon, exon not the last one")
    hb = {}
    J = {}
    for p in files:
        try:
            J[p] = read_any(p)
            hb[p] = J[p][(J[p].chrom == "chr11") & (J[p].js > 5_240_000) & (J[p].je < 5_260_000)].uniq.sum()
            print(f"  {p}: {len(J[p])} junctions, reads in HBG1/2 {hb[p]:.0f}")
        except Exception as e:
            print(f"  {p}: could not read ({e})")
    if len(J) >= 2:
        k562 = next((p for p in J if "K562" in os.path.basename(p)), max(J, key=hb.get))
        hct = min(J, key=hb.get)
        print(f"  K562 <- {k562}; HCT116 <- {hct}")
        for line, p in (("K562", k562), ("HCT116", hct)):
            j = J[p]
            x = sites[(sites.line == line) & (sites.ex_last == False) & sites.ex_s.notna() & (sites.dist_down < 10000)].copy()
            # calibrate coordinate offsets (junction ends and starts separately) at exon boundaries
            ends = set(zip(j.chrom, j.je.astype(np.int64)))
            starts = set(zip(j.chrom, j.js.astype(np.int64)))
            be = max((sum((c, int(s0) - 1 + o) in ends for c, s0 in zip(x.chrom, x.ex_s)), o) for o in (-1, 0, 1))
            bs = max((sum((c, int(e0) + 1 + o) in starts for c, e0 in zip(x.chrom, x.ex_e)), o) for o in (-1, 0, 1))
            off_e, off_s = be[1], bs[1]
            best = (be[0] + bs[0], f"ends {off_e:+d}, starts {off_s:+d}")
            idx = {ch: (g.js.values, g.je.values, g.uniq.values) for ch, g in j.groupby("chrom")}

            def psi(c, s0, e0):
                if c not in idx:
                    return np.nan
                js_, je_, u_ = idx[c]
                left = u_[je_ == s0 - 1 + off_e].sum()
                right = u_[js_ == e0 + 1 + off_s].sum()
                skip = u_[(js_ - off_s < s0) & (je_ - off_e > e0)].sum()
                inc = (left + right) / 2
                return inc / (inc + skip) if inc + skip >= 10 else np.nan

            x["psi2"] = [psi(c, int(s0), int(e0)) for c, s0, e0 in zip(x.chrom, x.ex_s, x.ex_e)]
            print(f"\n  {line}: offsets {best[1]} (matches {best[0]}); PSI for {x.psi2.notna().sum()} of {len(x)} sites")
            x["klasa"] = np.select([x.psi2 >= 0.9, x.psi2 < 0.9], ["constitutive", "alternative"], "no data")
            for k, g in x.groupby("klasa"):
                a, r = g[g.grp == "start"], g[g.grp == "losowa"]
                da, dr = a.drb.dropna().astype(bool), r.drb.dropna().astype(bool)
                print(f"   {k:<14} start sites {len(a):>4}, random {len(r):>5}; start-site cap_sum {a.cap_sum.median():>5.0f}; "
                      f"DRB {100 * da.mean() if len(da) else float('nan'):5.1f}% / {100 * dr.mean() if len(dr) else float('nan'):5.1f}%")
            a1 = x[(x.grp == "start") & (x.klasa == "alternative")].cap_sum
            a2 = x[(x.grp == "start") & (x.klasa == "constitutive")].cap_sum
            if len(a1) > 5 and len(a2) > 5:
                print(f"   start-site cap_sum: alternative {a1.median():.0f} vs. constitutive {a2.median():.0f}, "
                      f"Mann-Whitney p = {mannwhitneyu(a1, a2).pvalue:.2g}")

print("\n==== end ====")

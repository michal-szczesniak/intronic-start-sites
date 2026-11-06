#!/usr/bin/env python
"""Step 37: the 3-kb rule (step 36) for the remaining stable-RNA-based results.
  1. Fig. 5a: stable RNA per initiation event, all groups (the Fig. 5a source file is detected automatically).
  2. Fig. 4a: MTREX test (change in the RNA step downstream of the site after MTREX depletion).
  3. Supp. Table 1 v2: column distance_to_next_host_exon; stable RNA measures left empty for sites < 3 kb.
Run from the project directory:
    python scripts/37_stable_far_rest.py > results/exon_dist/stable_far_rest.txt 2>&1
"""
import glob, os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

FAR = 3000
LINES = {"K562": "results/v2", "HCT116": "results/hct116"}
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 30)


def sec(t):
    print(f"\n{'=' * 8} {t} {'=' * 8}", flush=True)


def b(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


def find_col(d, *subs):
    for s in subs:
        for c in d.columns:
            if s in str(c).lower():
                return c
    return None


# distance to the next host exon for every intronic site (start sites, promoters, random positions)
dist = []
for line, base in LINES.items():
    d = pd.read_csv(f"{base}/itss.describe.tsv", sep="\t", low_memory=False,
                    usecols=["site_id", "site_type", "strand", "pos", "intron_start", "intron_end", "host_active", "tss_supported"])
    for c in ("pos", "intron_start", "intron_end"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["dist_down"] = np.where(d.strand == "+", d.intron_end - d.pos, d.pos - d.intron_start)
    d["line"] = line
    dist.append(d)
D = pd.concat(dist, ignore_index=True)

# ---------- 1. Fig. 5a ----------
sec("1. Fig. 5a: source file")
cands = sorted(set(glob.glob("results/*fig5*") + glob.glob("results/**/*fig5*", recursive=True)))
print("candidates:", cands)
for path in [c for c in cands if c.endswith(".tsv")]:
    f5 = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"--- {path}: {len(f5)} rows; columns: {', '.join(map(str, f5.columns))}")
    print(f5.head(3).to_string())
    sid = "site_id" if "site_id" in f5 else None
    grp = find_col(f5, "grupa", "group", "grp", "klasa", "kind")
    lin = find_col(f5, "cell_line", "line", "linia")
    val = find_col(f5, "per_init", "per_read", "per_cap", "na_odczyt", "excess", "nadwyzka", "nadwyżka")
    print(f"detected: site_id={sid}, group={grp}, line={lin}, value={val}")
    if not (sid and grp and val):
        print("could not recognise the columns; check the header above")
        continue
    if lin is None:
        f5["_line"] = np.where(f5[sid].isin(D[D.line == "K562"].site_id), "K562", "HCT116")
        lin = "_line"
    m = f5.merge(D[["site_id", "line", "dist_down"]], left_on=[sid, lin], right_on=["site_id", "line"], how="left") \
        if lin != "_line" else f5.merge(D[["site_id", "dist_down"]], left_on=sid, right_on="site_id", how="left")
    m["v"] = pd.to_numeric(m[val], errors="coerce")
    m["intron"] = m.dist_down.notna()
    keep = ~m.intron | (m.dist_down >= FAR)
    print(f"\n{'line':<7} {'group':<40} {'n all':>11} {'% > 0':>7} {'med x1000':>10} | {'n (3-kb rule)':>15} {'% > 0':>7} {'med x1000':>10}")
    for (L, g), x in m.groupby([lin, grp]):
        y = x[keep.loc[x.index]]
        print(f"{str(L):<7} {str(g)[:40]:<40} {x.v.notna().sum():>11} {100 * (x.v > 0).mean():>6.1f}% {1000 * x.v.median():>10.3f} | "
              f"{y.v.notna().sum():>15} {100 * (y.v > 0).mean():>6.1f}% {1000 * y.v.median():>10.3f}")
    print("(sites without a host intron, e.g. between genes, are not subject to the rule and remain unchanged)")

# ---------- 2. Fig. 4a: MTREX ----------
sec("2. Fig. 4a: MTREX with the 3-kb rule")
mt = [p for p in glob.glob("results/v2/exosome/*MTREX*") if p.endswith(".tsv")]
print("files:", mt)
for path in mt:
    e = pd.read_csv(path, sep="\t", low_memory=False)
    print(f"--- {path}: {len(e)} rows; columns: {', '.join(map(str, e.columns))}")
    dcol = next((c for c in e.columns if str(c).lower() in ("delta", "d", "δ", "step_change", "diff")), None) \
        or find_col(e, "delta", "δ", "diff", "zmiana")
    if "site_id" not in e or dcol is None:
        print(e.head(3).to_string())
        continue
    x = e.merge(D[D.line == "K562"][["site_id", "site_type", "dist_down", "host_active", "tss_supported"]], on="site_id", how="left",
                suffixes=("", "_d"))
    st = x["site_type"].astype(str)
    act, sup = b(x["host_active"]), b(x["tss_supported"])
    x["v"] = pd.to_numeric(x[dcol], errors="coerce")
    for lab, far in (("all distances", pd.Series(True, index=x.index)), ("3-kb rule", x.dist_down >= FAR)):
        rnd = x[(st == "random_control") & act & far].v
        thr = np.nanpercentile(rnd, 95)
        print(f"\n  {lab}: threshold (95th percentile of random positions) = {thr:.3f}")
        groups = {"random positions": (st == "random_control") & act,
                  "internal promoters": (st == "internal_promoter") & sup & act,
                  "intronic start sites, active host": (st == "candidate") & sup & act,
                  "intronic start sites, silent host": (st == "candidate") & sup & ~act}
        for g, msk in groups.items():
            msk = msk & (far if g != "intronic start sites, silent host" else True)
            v = x[msk].v.dropna()
            p = mannwhitneyu(v, rnd.dropna(), alternative="greater").pvalue if g != "random positions" and len(v) > 2 else float("nan")
            print(f"    {g:<42} n {len(v):>5}  above threshold {100 * (v > thr).mean():5.1f}%  p = {p:.2g}")

# ---------- 3. Supp. Table 1 v2 ----------
sec("3. Supp. Table 1 v2")
s = pd.read_csv("results/final.supp_table1_intronic_starts.tsv", sep="\t", low_memory=False)
s = s.merge(D[["site_id", "line", "dist_down"]].rename(columns={"line": "cell_line", "dist_down": "distance_to_next_host_exon"}),
            on=["site_id", "cell_line"], how="left")
near = s.distance_to_next_host_exon < FAR
for c in [c for c in s.columns if c.startswith("stable_rna")]:
    s.loc[near, c] = np.nan
cols = list(s.columns)
cols.insert(cols.index("position") + 1 if "position" in cols else 4, cols.pop(cols.index("distance_to_next_host_exon")))
out = "results/final.supp_table1_intronic_starts.v2.tsv"
s[cols].to_csv(out, sep="\t", index=False)
print(f"wrote {out}: {len(s)} rows; distance matched for {s.distance_to_next_host_exon.notna().sum()}; "
      f"stable RNA measures cleared for {int(near.sum())} sites < 3 kb from the exon")
print("\n==== end ====")

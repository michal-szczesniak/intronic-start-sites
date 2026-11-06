#!/usr/bin/env python
"""Step 45: depth-matched deficit of the sense arm of divergent pairs, compared across definitions
(as in step 42 and with other matching rules). Input: results/{v2,hct116}/directionality.tsv (step 26).
DI = log2((S+1)/(A+1)); deficit = 1 - 2^x, x = a difference in DI between intronic and matched intergenic pairs."""
import numpy as np
import pandas as pd

LINES = {"K562": "results/v2/directionality.tsv", "HCT116": "results/hct116/directionality.tsv"}
TOL, R, B = 0.10, 200, 1000


def match(tot_i, tv, rng):
    """each intronic pair -> one unused intergenic pair with total reads within +-TOL (random among candidates)"""
    used = np.zeros(len(tv), bool)
    out = np.full(len(tot_i), -1)
    for k in rng.permutation(len(tot_i)):
        t = tot_i[k]
        i0, i1 = np.searchsorted(tv, t * (1 - TOL)), np.searchsorted(tv, t * (1 + TOL), side="right")
        free = np.flatnonzero(~used[i0:i1])
        if len(free):
            j = i0 + free[rng.integers(len(free))]
            used[j] = True
            out[k] = j
    return out


def measures(a, m):
    return {"diff of medians": 1 - 2 ** (np.median(a) - np.median(m)),
            "median of paired diffs": 1 - 2 ** np.median(a - m),
            "diff of means": 1 - 2 ** (np.mean(a) - np.mean(m))}


def run(a_di, a_tot, iv_di, iv_tot, rng):
    idx = match(a_tot, iv_tot, rng)
    ok = idx >= 0
    return measures(a_di[ok], iv_di[idx[ok]]), ok.mean()


for line, path in LINES.items():
    d = pd.read_csv(path, sep="\t", low_memory=False)
    for c in ("sens", "anty", "DI", "host_signal"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["tot"] = d.sens + d.anty
    inter = d[d.klasa == "międzygenowa"].dropna(subset=["DI", "tot"]).sort_values("tot")
    intr = d[d.klasa != "międzygenowa"].dropna(subset=["DI", "tot"])
    act = intr[intr.host_signal >= 5]
    print(f"\n==== {line} ====")
    print(f"intronic pairs: all {len(intr)}, active (host_signal >= 5) {len(act)}, host_signal > 5 {int((intr.host_signal > 5).sum())}, "
          f"genes with active pairs {act.host.nunique()}; intergenic pairs {len(inter)}")
    a_di, a_tot, hosts = act.DI.values, act.tot.values, act.host.values
    iv_di, iv_tot = inter.DI.values, inter.tot.values
    # point estimates: many random matchings
    reps = [run(a_di, a_tot, iv_di, iv_tot, np.random.default_rng(s))[0] for s in range(R)]
    pt = {k: np.median([r[k] for r in reps]) for k in reps[0]}
    rg = {k: (np.percentile([r[k] for r in reps], 2.5), np.percentile([r[k] for r in reps], 97.5)) for k in reps[0]}
    # bootstrap CIs: pairs, or whole genes (re-matching in every resample)
    rng = np.random.default_rng(1)
    ug = pd.unique(hosts)
    by_gene = {g: np.flatnonzero(hosts == g) for g in ug}
    bs_pairs, bs_genes = [], []
    for b in range(B):
        k = rng.integers(0, len(a_di), len(a_di))
        bs_pairs.append(run(a_di[k], a_tot[k], iv_di, iv_tot, rng)[0])
        gk = np.concatenate([by_gene[g] for g in rng.choice(ug, len(ug))])
        bs_genes.append(run(a_di[gk], a_tot[gk], iv_di, iv_tot, rng)[0])
    print(f"{'measure':24s} {'estimate (200 matchings)':>26s} {'range over matchings':>22s} {'95% CI, pairs':>16s} {'95% CI, genes':>16s}")
    for k in pt:
        cp = np.percentile([r[k] for r in bs_pairs], [2.5, 97.5])
        cg = np.percentile([r[k] for r in bs_genes], [2.5, 97.5])
        print(f"{k:24s} {100 * pt[k]:25.1f}% {100 * rg[k][0]:9.1f} to {100 * rg[k][1]:5.1f}% "
              f"{100 * cp[0]:7.1f}-{100 * cp[1]:5.1f}% {100 * cg[0]:7.1f}-{100 * cg[1]:5.1f}%")
print("\nStep 42 used one random matching and pair-level bootstrap with 'diff of medians'.")

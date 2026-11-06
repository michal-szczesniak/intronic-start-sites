#!/usr/bin/env python
"""Step 45b: depth-matched deficit of the sense arm, with 1:1 nearest-neighbour matching, for active and silent genes.
Matching: 1:1 greedy nearest neighbour on depth (sens + anty), without replacement, intronic pairs in random order
(seed 1), ties -> control with larger or equal depth, tolerance |d_ctrl - d_intr| <= 0.1 * max(d_intr, 1).
Effect: d = median DI (intronic) - median DI (matched controls); deficit = 1 - 2^d.
CI: percentile bootstrap over whole host genes (10,000; seed 1001), each intronic pair resampled with its control.
Also: median of paired differences, Wilcoxon signed-rank on paired differences."""
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

LINES = {"K562": "results/v2/directionality.tsv", "HCT116": "results/hct116/directionality.tsv"}


def nn_match(tot_i, tv, order):
    used = np.zeros(len(tv), bool)
    out = np.full(len(tot_i), -1)
    for k in order:
        t = tot_i[k]
        tol = 0.1 * max(t, 1)
        r = np.searchsorted(tv, t)          # first control with depth >= t
        l = r - 1
        while r < len(tv) and used[r]:
            r += 1
        while l >= 0 and used[l]:
            l -= 1
        dr = tv[r] - t if r < len(tv) else np.inf
        dl = t - tv[l] if l >= 0 else np.inf
        j, dj = (r, dr) if dr <= dl else (l, dl)
        if dj <= tol:
            used[j] = True
            out[k] = j
    return out


def analyse(name, grp, ctrl, seed_order=1, seed_bs=1001, nb=10000):
    ctrl = ctrl.sort_values("tot", kind="stable")
    tv, cdi = ctrl.tot.values, ctrl.DI.values
    a_tot, a_di, hosts = grp.tot.values, grp.DI.values, grp.host.values
    idx = nn_match(a_tot, tv, np.random.default_rng(seed_order).permutation(len(a_tot)))
    ok = idx >= 0
    a, m, h = a_di[ok], cdi[idx[ok]], hosts[ok]
    rel = np.abs(tv[idx[ok]] - a_tot[ok]) / np.maximum(a_tot[ok], 1)
    d = np.median(a) - np.median(m)
    rng = np.random.default_rng(seed_bs)
    ug = pd.unique(h)
    by = {g: np.flatnonzero(h == g) for g in ug}
    bs = []
    for _ in range(nb):
        k = np.concatenate([by[g] for g in rng.choice(ug, len(ug))])
        bs.append(1 - 2 ** (np.median(a[k]) - np.median(m[k])))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    p = wilcoxon(a - m).pvalue
    print(f"  {name:8s} matched {ok.sum()} of {len(ok)} pairs in {len(ug)} genes; identical depth {100 * (rel == 0).mean():.1f}%, "
          f"max rel. diff {100 * rel.max():.1f}%")
    print(f"           d = {d:.3f} ({np.median(a):.3f} vs {np.median(m):.3f}): deficit {100 * (1 - 2 ** d):.1f}% "
          f"(95% CI {100 * lo:.1f}-{100 * hi:.1f}%, genes); Wilcoxon p = {p:.2g}; "
          f"median of paired differences {np.median(a - m):.3f} ({100 * (1 - 2 ** np.median(a - m)):.1f}%)")


for line, path in LINES.items():
    d = pd.read_csv(path, sep="\t", low_memory=False)
    for c in ("sens", "anty", "DI", "host_signal"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["tot"] = d.sens + d.anty
    ctrl = d[d.klasa == "międzygenowa"].dropna(subset=["DI", "tot"])
    intr = d[d.klasa != "międzygenowa"].dropna(subset=["DI", "tot"])
    print(f"\n==== {line} ==== controls {len(ctrl)}")
    analyse("active", intr[intr.host_signal >= 5], ctrl)
    analyse("silent", intr[intr.host_signal < 5], ctrl)

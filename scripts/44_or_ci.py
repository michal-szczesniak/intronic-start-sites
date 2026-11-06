#!/usr/bin/env python
"""Step 44: 95% CIs of the Mantel-Haenszel odds ratios in Fig. 2d (transposon with TE classes only; stable RNA)."""
import numpy as np
import pandas as pd
from scipy.stats import chi2

SUPP = "results/final.supp_table1_intronic_starts.v2.tsv"
TE = ("LINE", "SINE", "LTR", "DNA", "Retroposon", "RC")


def b(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


def mh_or(tables):
    R = S = PR = PSQR = QS = num = den = 0.0
    for a, b_, c, d in tables:
        n = a + b_ + c + d
        if n < 2:
            continue
        r, s = a * d / n, b_ * c / n
        p, q = (a + d) / n, (b_ + c) / n
        R += r; S += s; PR += p * r; PSQR += p * s + q * r; QS += q * s
        num += a - (a + b_) * (a + c) / n
        den += (a + b_) * (c + d) * (a + c) * (b_ + d) / (n * n * (n - 1))
    o = R / S
    se = np.sqrt(PR / (2 * R ** 2) + PSQR / (2 * R * S) + QS / (2 * S ** 2))
    return o, np.exp(np.log(o) - 1.96 * se), np.exp(np.log(o) + 1.96 * se), chi2.sf((abs(num) - 0.5) ** 2 / den, 1)


t = pd.read_csv(SUPP, sep="\t", low_memory=False)
print("repeat_class values (confirmed sites):")
print(t[b(t.capped_start_confirmed)].repeat_class.astype(str).value_counts().head(25).to_string())
stab = [c for c in t.columns if "stable" in c.lower()]
print("stable-RNA columns:", stab)
for line in ("K562", "HCT116"):
    x = t[(t.cell_line == line) & b(t.capped_start_confirmed) & b(t.host_active)].copy()
    x["third"] = pd.qcut(pd.to_numeric(x.procap_reads).rank(method="first"), 3, labels=False)
    sh = b(x.shared_with_other_cell_line)
    rc = x.repeat_class.astype(str)
    feats = {"conserved (phyloP > 1)": (pd.to_numeric(x.phylop_mean, errors="coerce") > 1, pd.Series(True, index=x.index)),
             "in a transposon (TE classes)": (rc.str.split("/").str[0].isin(TE), pd.Series(True, index=x.index))}
    for c in stab:
        v = x[c].astype(str).str.strip().str.lower()
        ok = v.isin(["true", "false", "1", "0", "1.0", "0.0", "yes", "no"])
        if ok.sum() > 50:
            feats[f"stable RNA ({c})"] = (v.isin(["true", "1", "1.0", "yes"]), ok)
    for f, (y, ok) in feats.items():
        tabs = []
        for th in range(3):
            m = (x.third == th) & ok
            tabs.append(((sh & y & m).sum(), (sh & ~y & m).sum(), (~sh & y & m).sum(), (~sh & ~y & m).sum()))
        o, lo, hi, p = mh_or(tabs)
        print(f"{line:6s} {f:40s} MH OR {o:.2f} (95% CI {lo:.2f}-{hi:.2f}), CMH p = {p:.2g}, n = {int(ok.sum())}")

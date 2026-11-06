#!/usr/bin/env python
"""Step 42: additional controls reported in the text.
Each section is independent; a failing section prints its error and the script goes on."""
import traceback

import numpy as np
import pandas as pd
from scipy.stats import chi2, mannwhitneyu, norm, wilcoxon

LINES = {"K562": "results/v2", "HCT116": "results/hct116"}
SUPP = "results/final.supp_table1_intronic_starts.v2.tsv"
RNG = np.random.default_rng(42)


def b(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


def col(df, *names):
    return next((n for n in names if n in df.columns), None)


def section(title):
    def wrap(f):
        print(f"\n======== {title} ========", flush=True)
        try:
            f()
        except Exception:
            traceback.print_exc()
        return f
    return wrap


def boot_median(v, n=1000):
    v = np.asarray(v, float)
    m = [np.median(RNG.choice(v, len(v))) for _ in range(n)]
    return np.percentile(m, 2.5), np.percentile(m, 97.5)


def mh_or(tables):
    """Mantel-Haenszel odds ratio with Robins-Breslow-Greenland 95% CI and CMH p. tables: list of (a, b, c, d)."""
    R = S = PR = PS_QR = QS = 0.0
    num = den = var = 0.0
    for a, b_, c, d in tables:
        n = a + b_ + c + d
        if n < 2:
            continue
        r, s = a * d / n, b_ * c / n
        p, q = (a + d) / n, (b_ + c) / n
        R += r; S += s; PR += p * r; PS_QR += p * s + q * r; QS += q * s
        e = (a + b_) * (a + c) / n
        v = (a + b_) * (c + d) * (a + c) * (b_ + d) / (n * n * (n - 1))
        num += a - e; den += v
    orr = R / S
    se = np.sqrt(PR / (2 * R ** 2) + PS_QR / (2 * R * S) + QS / (2 * S ** 2))
    stat = (abs(num) - 0.5) ** 2 / den
    return orr, np.exp(np.log(orr) - 1.96 * se), np.exp(np.log(orr) + 1.96 * se), chi2.sf(stat, 1)


@section("1. PAS density by DRB class, confirmed intronic start sites only (HCT116)")
def pas_drb():
    t = pd.read_csv(SUPP, sep="\t", low_memory=False)
    x = t[(t.cell_line == "HCT116") & b(t.capped_start_confirmed)]
    cls = x.drb_class_10min.astype(str).str.lower()
    pas = pd.to_numeric(x.pas_per_kb_transcribed_strand, errors="coerce")
    for lab, m in (("active hosts", b(x.host_active)), ("all hosts", pd.Series(True, index=x.index))):
        g = {c: pas[m & cls.str.startswith(c)].dropna() for c in ("productive", "short", "no")}
        resp = pd.concat([g["productive"], g["short"]])
        print(f"  {lab}: " + ", ".join(f"{c} {v.mean():.2f} (n = {len(v)})" for c, v in g.items())
              + f"; productive vs short p = {mannwhitneyu(g['productive'], g['short']).pvalue:.2g}"
              + f"; responsive {resp.mean():.2f} vs no response p = {mannwhitneyu(resp, g['no']).pvalue:.2g}")
        prot = b(x.splice_donor_before_first_pas)
        tab = pd.crosstab(cls[m & cls.str.match('productive|short|no')].str.split().str[0], prot[m])
        print("  protected arrangement by class (share protected): "
              + ", ".join(f"{k} {100 * r[True] / r.sum():.1f}%" for k, r in tab.iterrows() if True in r))


@section("2. Odds ratios with 95% CI, shared vs cell-type-specific confirmed sites in active hosts (Fig. 2d)")
def shared_or():
    t = pd.read_csv(SUPP, sep="\t", low_memory=False)
    for line in ("K562", "HCT116"):
        x = t[(t.cell_line == line) & b(t.capped_start_confirmed) & b(t.host_active)].copy()
        x["third"] = pd.qcut(pd.to_numeric(x.procap_reads).rank(method="first"), 3, labels=False)
        sh = b(x.shared_with_other_cell_line)
        feats = {"conserved (phyloP > 1)": pd.to_numeric(x.phylop_mean, errors="coerce") > 1,
                 "in a transposon": x.repeat_class.notna() & ~x.repeat_class.astype(str).str.lower().isin(["", "nan", "none", "."])}
        for f, y in feats.items():
            tabs = []
            for th in range(3):
                m = x.third == th
                tabs.append(((sh & y & m).sum(), (sh & ~y & m).sum(), (~sh & y & m).sum(), (~sh & ~y & m).sum()))
            o, lo, hi, p = mh_or(tabs)
            print(f"  {line:6s} {f:24s} MH OR {o:.2f} (95% CI {lo:.2f}-{hi:.2f}), CMH p = {p:.2g}; n shared {int(sh.sum())}, specific {int((~sh).sum())}")
        try:
            import statsmodels.formula.api as smf
            x["shared"] = sh.astype(int)
            x["cons"] = feats["conserved (phyloP > 1)"].astype(int)
            x["te"] = feats["in a transposon"].astype(int)
            x["lp"] = np.log10(pd.to_numeric(x.procap_reads) + 1)
            r = smf.logit("shared ~ cons + te + lp", data=x).fit(disp=0)
            ci = r.conf_int()
            print(f"  {line:6s} logistic regression shared ~ conserved + transposon + log10 PRO-cap: " + "; ".join(
                f"{k} OR {np.exp(r.params[k]):.2f} ({np.exp(ci.loc[k, 0]):.2f}-{np.exp(ci.loc[k, 1]):.2f}), p = {r.pvalues[k]:.2g}"
                for k in ("cons", "te", "lp")))
        except ImportError:
            print("  (statsmodels not installed: logistic regression skipped)")


@section("3. Divergent pairs: bootstrap CI of medians, depth matching with explicit tests, silent genes (Supp. Fig. S8)")
def pairs():
    pooled_p = []
    for line, base in LINES.items():
        d = pd.read_csv(f"{base}/directionality.tsv", sep="\t", low_memory=False)
        d["tot"] = pd.to_numeric(d.sens, errors="coerce") + pd.to_numeric(d.anty, errors="coerce")
        inter = d[d.klasa == "międzygenowa"].dropna(subset=["DI", "tot"])
        intr = d[d.klasa != "międzygenowa"].dropna(subset=["DI", "tot"]).copy()
        hs = pd.to_numeric(intr.host_signal, errors="coerce")
        groups = {"between genes": inter, "intron, active": intr[hs >= 5], "intron, silent": intr[hs < 5],
                  "intron, silent (no capped read at host start)": intr[hs == 0]}
        k36 = b(groups["intron, active"].k36)
        groups["active, H3K36me3"] = groups["intron, active"][k36]
        groups["active, no H3K36me3"] = groups["intron, active"][~k36]
        for g, x in groups.items():
            lo, hi = boot_median(x.DI)
            print(f"  {line:6s} {g:45s} n = {len(x):5d}  median ratio {2 ** x.DI.median():.2f} (95% CI {2 ** lo:.2f}-{2 ** hi:.2f}), "
                  f"Wilcoxon vs 0 p = {wilcoxon(x.DI).pvalue:.2g}")
        # depth matching: each intronic pair to one pair between genes with total reads within 10%
        pool = inter.sort_values("tot")
        tv, dv = pool.tot.values, pool.DI.values
        for g in ("intron, active", "intron, silent"):
            x = groups[g].sample(frac=1, random_state=1)
            used = np.zeros(len(tv), bool)
            a, m = [], []
            for tot, di in zip(x.tot.values, x.DI.values):
                i0, i1 = np.searchsorted(tv, tot * 0.9), np.searchsorted(tv, tot * 1.1, side="right")
                free = np.flatnonzero(~used[i0:i1])
                if len(free):
                    j = i0 + free[RNG.integers(len(free))]
                    used[j] = True
                    a.append(di); m.append(dv[j])
            a, m = np.array(a), np.array(m)
            deficit = 1 - 2 ** (np.median(a) - np.median(m))
            bs = []
            for _ in range(1000):
                k = RNG.integers(0, len(a), len(a))
                bs.append(1 - 2 ** (np.median(a[k]) - np.median(m[k])))
            p_w = wilcoxon(a - m).pvalue
            p_mw = mannwhitneyu(a, m).pvalue
            print(f"  {line:6s} {g:15s} matched {len(a)} of {len(x)}: deficit {100 * deficit:.1f}% "
                  f"(95% CI {100 * np.percentile(bs, 2.5):.1f} to {100 * np.percentile(bs, 97.5):.1f}%), "
                  f"Wilcoxon signed-rank on matched differences p = {p_w:.2g}, Mann-Whitney p = {p_mw:.2g}")
        pooled_p.append(wilcoxon(groups["intron, silent"].DI, alternative="less").pvalue)
    stat = -2 * np.sum(np.log(pooled_p))
    print(f"  silent genes, one-sided Wilcoxon (sense arm weaker) per line p = "
          + ", ".join(f"{p:.2g}" for p in pooled_p) + f"; combined (Fisher) p = {chi2.sf(stat, 2 * len(pooled_p)):.2g}")


@section("4. Stable-RNA excess at random intronic positions and other groups (baseline for Fig. 5a)")
def baseline():
    for line in LINES:
        cx = pd.read_csv(f"results/context_{line}.tsv", sep="\t", low_memory=False)
        gc = col(cx, "grupa", "group", "grp")
        num = [c for c in cx.columns if c not in (gc,) and pd.api.types.is_numeric_dtype(cx[c])
               and any(k in c.lower() for k in ("rna", "dal", "skok", "excess", "po", "przed", "after", "before"))]
        print(f"  {line}: numeric RNA-type columns: {', '.join(num)}")
        for c in num:
            g = cx.groupby(gc)[c]
            print(f"    {c:22s} " + "; ".join(f"{k} median {v.median():.3g}, >0 {100 * (v > 0).mean():.0f}% (n = {v.notna().sum()})"
                                             for k, v in g))


@section("5. DRB group test without sites that have an active internal promoter up to 20 kb upstream (HCT116, 10 min)")
def drb_ip():
    d = pd.read_csv("results/hct116/ttdrb_10.drb.tsv", sep="\t", low_memory=False)
    print("  columns: " + ", ".join(d.columns[:40]))
    gi = col(d, "gap_index", "gap")
    desc = pd.read_csv("results/hct116/itss.describe.tsv", sep="\t", low_memory=False)
    ip = desc[(desc.site_type == "internal_promoter") & b(desc.tss_supported)]
    ipk = {}
    for c, s, p in zip(ip.chrom, ip.strand, ip.pos):
        ipk.setdefault((c, s), []).append(int(p))
    ipk = {k: np.sort(v) for k, v in ipk.items()}

    def has_ip(c, s, p):
        v = ipk.get((c, s))
        if v is None:
            return False
        lo, hi = (p - 20000, p - 1) if s == "+" else (p + 1, p + 20000)
        i = np.searchsorted(v, lo)
        return i < len(v) and v[i] <= hi
    ctrl = d[(d.site_type == "random_control") & b(d.host_active) & b(d.testable)]
    thr = np.nanpercentile(ctrl[gi], 95)
    x = d[(d.site_type == "candidate") & b(d.tss_supported) & b(d.host_active) & b(d.testable)].copy()
    x["ip_up"] = [has_ip(c, s, int(p)) for c, s, p in zip(x.chrom, x.strand, x.pos)]
    for lab, y in (("all", x), ("with active internal promoter <= 20 kb upstream", x[x.ip_up]),
                   ("without", x[~x.ip_up])):
        v = y[gi].dropna()
        print(f"  {lab:48s} n = {len(v):5d}  above threshold {100 * (v > thr).mean():.1f}%  "
              f"Mann-Whitney vs random p = {mannwhitneyu(v, ctrl[gi].dropna(), alternative='greater').pvalue:.2g}")


print("\n==== done ====")

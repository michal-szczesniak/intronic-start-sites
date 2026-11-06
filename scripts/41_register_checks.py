#!/usr/bin/env python
"""Step 41: values reported in the text that no earlier step printed.
Each section is independent; a failing section prints its error and the script goes on."""
import bisect
import glob
import importlib.util
import subprocess
import traceback

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kruskal, mannwhitneyu, spearmanr, wilcoxon

LINES = {"K562": "results/v2", "HCT116": "results/hct116"}


def b(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


def col(df, *names):
    """First of the given column names present in df (None if none is)."""
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


@section("1. Divergent pairs: medians, Wilcoxon against 0, active vs silent, H3K36me3 (Supp. Fig. S8)")
def pairs():
    for line, base in LINES.items():
        d = pd.read_csv(f"{base}/directionality.tsv", sep="\t", low_memory=False)
        intr = d[d.klasa != "międzygenowa"].copy()
        intr["active"] = pd.to_numeric(intr.host_signal, errors="coerce") >= 5
        groups = {"between genes": d[d.klasa == "międzygenowa"], "intron, active": intr[intr.active],
                  "intron, silent": intr[~intr.active]}
        for g, x in groups.items():
            di = x.DI.dropna()
            print(f"  {line:6s} {g:15s} n = {len(di):5d}  median ratio {2 ** di.median():.2f}  "
                  f"Wilcoxon vs 0 p = {wilcoxon(di).pvalue:.2g}")
        a, s = groups["intron, active"].DI.dropna(), groups["intron, silent"].DI.dropna()
        print(f"  {line:6s} active vs silent: Mann-Whitney p = {mannwhitneyu(a, s).pvalue:.2g}")
        act = groups["intron, active"]
        k = b(act.k36)
        print(f"  {line:6s} active, with vs without H3K36me3: n = {k.sum()} / {(~k).sum()}, medians "
              f"{2 ** act.DI[k].median():.2f} / {2 ** act.DI[~k].median():.2f}, "
              f"Mann-Whitney p = {mannwhitneyu(act.DI[k].dropna(), act.DI[~k].dropna()).pvalue:.2g}")
        r = spearmanr(pd.to_numeric(act.host_signal, errors="coerce"), act.DI, nan_policy="omit")
        print(f"  {line:6s} active: Spearman DI vs host promoter signal r = {r.statistic:.3f}, p = {r.pvalue:.2g}")


@section("2. Splice-donor strength of the next host exon vs initiation and DRB response (sites < 10 kb)")
def donor():
    s = pd.read_csv("results/exon_dist/sites_exon_dist.tsv", sep="\t", low_memory=False)
    print("  columns: " + ", ".join(s.columns))
    dc = col(s, "donor_score", "donor_pwm", "pwm_donor", "donor", "ss_score", "donor_down")
    gc = col(s, "grp", "group", "kind", "site_type", "typ")
    dr = col(s, "drb", "drb_pos", "drb_positive", "drb_above", "drb_hit")
    if gc:
        print(f"  group column {gc}: " + ", ".join(f"{k} {v}" for k, v in s[gc].astype(str).value_counts().items()))
    if dc is None:
        print("  no donor-score column found: please send the column list above")
        return
    start = next((v for v in ("start", "itss", "candidate", "iTSS") if gc and (s[gc].astype(str) == v).any()), None)
    for line in ("K562", "HCT116"):
        m = (s.line == line) & (pd.to_numeric(s.dist_down, errors="coerce") < 10000) & ~b(s.ex_last) & s[dc].notna()
        if start:
            m &= s[gc].astype(str) == start
        x = s[m].copy()
        x["donor_score"] = pd.to_numeric(x[dc], errors="coerce")
        x["drb"] = x[dr] if dr else np.nan
        r = spearmanr(x.donor_score, x.cap_sum, nan_policy="omit")
        x["t"] = pd.qcut(x.donor_score.rank(method="first"), 3, labels=["weak", "medium", "strong"])
        kw = kruskal(*[g.cap_sum.dropna() for _, g in x.groupby("t", observed=True)])
        print(f"  {line}: n = {len(x)}, Spearman donor score vs capped reads r = {r.statistic:.3f}, p = {r.pvalue:.2g}; "
              f"Kruskal-Wallis over thirds p = {kw.pvalue:.2g}; medians "
              + ", ".join(f"{t} {g.cap_sum.median():.0f}" for t, g in x.groupby("t", observed=True)))
        y = x[x.drb.notna()].copy()
        if len(y) > 30:
            y["hit"] = b(y.drb) if y.drb.dtype == object or y.drb.dtype == bool else pd.to_numeric(y.drb) > 0
            tab = pd.crosstab(y.t, y.hit)
            print(f"         DRB response ({dr}) by third: " + ", ".join(
                f"{t} {100 * y[y.t == t].hit.mean():.1f}% (n = {(y.t == t).sum()})" for t in tab.index)
                  + f"; chi-square p = {chi2_contingency(tab).pvalue:.2g}")


@section("3. Intronic start sites at any annotated transcript start (±500 nt, same strand), current sets")
def annotated():
    spec = importlib.util.spec_from_file_location("c10", "scripts/10_candidates.py")
    c10 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(c10)
    tx, _ = c10.parse_gtf("ref/gencode.v47.annotation.gtf.gz")
    tss = {}
    for t in tx.values():
        tss.setdefault((t["chrom"], t["strand"]), []).append(c10.tss_of(t))
    tss = {k: sorted(v) for k, v in tss.items()}
    for line, base in LINES.items():
        d = pd.read_csv(f"{base}/itss.describe.tsv", sep="\t", low_memory=False)
        d = d[d.site_type == "candidate"]

        def near(c, s, p):
            tl = tss.get((c, s), [])
            j = bisect.bisect_left(tl, p - 500)
            return j < len(tl) and tl[j] <= p + 500
        hit = np.array([near(c, s, int(p)) for c, s, p in zip(d.chrom, d.strand, d.pos)])
        conf = b(d.tss_supported).values
        print(f"  {line}: all {len(d)} sites, at an annotated start {hit.sum()} ({100 * hit.mean():.1f}%), "
              f"not at one {100 * (1 - hit.mean()):.1f}%; confirmed {conf.sum()}: not at one {100 * (1 - hit[conf].mean()):.1f}%")


@section("4. Polyadenylation-signal ratio, transcribed / opposite strand, per group (Fig. 4c)")
def pas():
    for line, base in LINES.items():
        u = pd.read_csv(f"{base}/u1pas.tsv", sep="\t", low_memory=False)
        act, sup = b(u.host_active), b(u.tss_supported)
        groups = {"host promoters": u.site_type == "host_promoter",
                  "internal promoters in use": (u.site_type == "internal_promoter") & sup & act,
                  "all internal promoters": u.site_type == "internal_promoter",
                  "intronic, confirmed, active host": (u.site_type == "candidate") & sup & act,
                  "random, active host": (u.site_type == "random_control") & act,
                  "random, all": u.site_type == "random_control"}
        for g, m in groups.items():
            x = u[m]
            print(f"  {line:6s} {g:34s} n = {len(x):5d}  PAS per kb {x.pas_sens.mean():.2f} / {x.pas_anty.mean():.2f}"
                  f"  ratio {x.pas_sens.mean() / x.pas_anty.mean():.2f}  donor before first PAS "
                  f"{100 * b(x.protected).mean():.2f}% ({int(b(x.protected).sum())} of {len(x)})")
        cx = pd.read_csv(f"results/context_{line}.tsv", sep="\t", low_memory=False)
        gc = col(cx, "grupa", "group", "grp", "site_type")
        if gc is None or "pas_sens" not in cx or "pas_anty" not in cx:
            print(f"  {line:6s} context file columns: " + ", ".join(cx.columns))
            continue
        for g, x in cx.groupby(gc):
            print(f"  {line:6s} context: {g:28s} n = {len(x):5d}  PAS per kb {x.pas_sens.mean():.2f} / "
                  f"{x.pas_anty.mean():.2f}  ratio {x.pas_sens.mean() / x.pas_anty.mean():.2f}")


@section("5. PRO-cap depth: summed signal of the capped bigWigs (millions of reads)")
def depth():
    import pyBigWig
    for f in ["data/cap_plus.bw", "data/cap_minus.bw"] + sorted(glob.glob("data/hct116/*line_*.bw")):
        bw = pyBigWig.open(f)
        tot = sum(sum(abs(v) * (e - s) for s, e, v in (bw.intervals(c) or [])) for c in bw.chroms())
        print(f"  {f}: {tot / 1e6:.2f} million")


@section("6. Nascent RNA per initiation event, medians (Fig. 5b)")
def fig5b():
    for f in sorted(glob.glob("results/final.fig5*.tsv")):
        t = pd.read_csv(f, sep="\t")
        m = t.apply(lambda r: r.astype(str).str.contains("nascent|DRB", case=False).any(), axis=1)
        print(f"--- {f}")
        print(t[m].to_string(index=False) if m.any() else t.head(20).to_string(index=False))


@section("7. Version of pyGenomeTracks and the wave-front margin in step 27")
def misc():
    for e in sorted(glob.glob("envs/*")):
        r = subprocess.run(f"conda list -p {e} 2>/dev/null | grep -i -E '^pygenometracks '", shell=True,
                           capture_output=True, text=True)
        if r.stdout.strip():
            print(f"  {e}: {r.stdout.strip()}")
    for cmd in ("grep -n -i 'front' scripts/27_classify.py | head -15",
                "ls scripts | grep -E '^29'; grep -n -i -E 'front|22200|40000|min.?dist|drb' scripts/29_*.py | head -30",
                "grep -n -i -E 'front|dist' jobs/context.sh | head -10",
                "grep -n -E 'cap_uncap|uncap_sum' scripts/*.py | head -20"):
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        print(f"$ {cmd}\n{r.stdout}{r.stderr}")


@section("8. Supplementary Table 1 v2: columns and rows")
def supp():
    f = "results/final.supp_table1_intronic_starts.v2.tsv"
    t = pd.read_csv(f, sep="\t", low_memory=False)
    print(f"  {f}: {len(t)} rows, {t.shape[1]} columns")
    print("  " + ", ".join(t.columns))
    lc = col(t, "cell_line", "line")
    if lc:
        print(t.groupby(lc).size().to_string())
    if "distance_to_next_host_exon" in t and "stable_rna_rise_log2" in t:
        near = pd.to_numeric(t.distance_to_next_host_exon, errors="coerce") < 3000
        print(f"  sites < 3 kb from the next exon with a stable RNA value (should be 0): "
              f"{int((near & t.stable_rna_rise_log2.notna()).sum())}")

@section("9. Additional checks: Fig. 5b distance to the host start, PAS by DRB class, bidirectional peaks")
def additional_checks():
    t = pd.read_csv("results/final.supp_table1_intronic_starts.v2.tsv", sep="\t", low_memory=False)
    lc = col(t, "cell_line", "line")
    for line, base in LINES.items():
        d = pd.read_csv(f"{base}/itss.describe.tsv", sep="\t", low_memory=False)
        d = d[d.site_type == "candidate"]
        conf = d[b(d.tss_supported)]
        if "peak_kind" in d:
            print(f"  {line}: peak kind of the {len(conf)} confirmed sites: " + ", ".join(
                f"{k} {v} ({100 * v / len(conf):.1f}%)" for k, v in conf.peak_kind.astype(str).value_counts().items()))
        if line != "HCT116":
            continue
        x = t[t[lc] == line].merge(d[["site_id", "dist_used_tss"]], on="site_id", how="left")
        print(f"  {line}: Supp. Table rows {len(x)}, matched to a host distance {x.dist_used_tss.notna().sum()}")
        act = b(x.host_active)
        far = pd.to_numeric(x.drb_far_per_initiation, errors="coerce")
        dist = pd.to_numeric(x.dist_used_tss, errors="coerce").abs()
        for lab, m in (("all", act), (">= 22.2 kb", act & (dist >= 22200)), (">= 40 kb", act & (dist >= 40000))):
            v = far[m].dropna()
            if len(v):
                bs = [np.median(np.random.default_rng(i).choice(v, len(v))) for i in range(1000)]
                print(f"  Fig. 5b intronic, active host, {lab:10s}: n = {len(v)}, median {v.median() * 1e3:.3f} x 10^-3 "
                      f"(95% CI {np.percentile(bs, 2.5) * 1e3:.3f}-{np.percentile(bs, 97.5) * 1e3:.3f}), "
                      f"distance to host start median {dist[m & far.notna()].median() / 1000:.1f} kb")
        cls = x.drb_class_10min.astype(str).str.lower()
        pas = pd.to_numeric(x.pas_per_kb_transcribed_strand, errors="coerce")
        for lab, m in (("all hosts", x.drb_class_10min.notna()), ("active hosts", x.drb_class_10min.notna() & act)):
            g = {c: pas[m & cls.str.startswith(c)].dropna() for c in ("productive", "short", "no")}
            resp = pd.concat([g["productive"], g["short"]])
            print(f"  PAS per kb by DRB class, {lab}: " + ", ".join(f"{c} {v.mean():.2f} (n = {len(v)})" for c, v in g.items())
                  + f"; productive vs short p = {mannwhitneyu(g['productive'], g['short']).pvalue:.2g}; responsive vs no "
                  f"response p = {mannwhitneyu(resp, g['no']).pvalue:.2g}; Kruskal-Wallis p = {kruskal(*g.values()).pvalue:.2g}")


print("\n==== done ====")

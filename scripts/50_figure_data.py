#!/usr/bin/env python3
"""Step 50: the plotted values of every data figure panel, computed from the result tables.

Writes one tidy TSV per panel to results/figures/data/<fig>_<panel>.tsv, plus
  results/figures/data/_status.tsv    one row per panel: ok / missing input / error, rows written, inputs used
  results/figures/data/_sources.tsv   every input file read, with size and modification time
and prints a short report.

Each panel function names the analysis step whose definition it re-uses (same filters, thresholds, windows,
confidence intervals).  Where a value exists only in a text file written by an analysis step (*.summary.txt,
exo.*.txt) or can be recomputed from the signal files only, this is said in the panel's docstring; by default the
signal files are used when pyBigWig and the files are available, and the text file is read as a fallback and as a
cross-check.

Run from the project directory:
    python scripts/50_figure_data.py                      # all panels
    python scripts/50_figure_data.py --panels fig2_* figS8 # some panels (shell-style patterns)
    python scripts/50_figure_data.py --no-bigwig          # profiles from the *.summary.txt files only
    python scripts/50_figure_data.py --skip-gtf           # skip Fig. 1b (genomic location of all PRO-cap starts)
Only pandas, numpy and scipy are needed (pyBigWig optionally, for the profile panels).
"""
import argparse
import fnmatch
import importlib.util
import math
import os
import re
import sys
import time
import traceback
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import chi2, mannwhitneyu, spearmanr

LINES = {"K562": "results/v2", "HCT116": "results/hct116"}
SUPP = "results/final.supp_table1_intronic_starts.v2.tsv"
OUT = "results/figures/data"
FAR = 3000  # 3-kb rule (steps 36-39): stable RNA only >= 3 kb upstream of the next exon of the host
GTF = "ref/gencode.v47.annotation.gtf.gz"
PINTS = {"uni": "data/encode/ENCSR261KBX/ENCFF636REV.bed.gz", "bi": "data/encode/ENCSR261KBX/ENCFF819CPA.bed.gz"}
SOURCES = []  # (panel, path)
ARGS = None


# ----------------------------------------------------------------------------------------------- helpers
class MissingInput(Exception):
    pass


def b(s):
    """Boolean columns as written by pandas ('True'/'False'), as in steps 35-44."""
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


def num(s):
    return pd.to_numeric(s, errors="coerce")


CUR = {"panel": ""}


def need(path):
    if not os.path.exists(path):
        raise MissingInput(path)
    SOURCES.append((CUR["panel"], path))
    return path


def read(path, **kw):
    return pd.read_csv(need(path), sep="\t", low_memory=False, **kw)


def wilson(k, n):
    """Wilson 95% interval in percent (as in steps 38-39)."""
    if n == 0:
        return np.nan, np.nan
    z, p = 1.96, k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return 100 * (c - h), 100 * (c + h)


def dist_down(d):
    """Distance to the next exon of the host's reference transcript, in the direction of transcription (step 36)."""
    pos, s, e = num(d["pos"]), num(d["intron_start"]), num(d["intron_end"])
    return pd.Series(np.where(d["strand"].astype(str) == "+", e - pos, pos - s), index=d.index)


def load_script(name):
    path = need(os.path.join("scripts", name))
    spec = importlib.util.spec_from_file_location("m" + name.split("_")[0], path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def mwu_greater(a, b_):
    a, b_ = np.asarray(a, float), np.asarray(b_, float)
    if len(a) < 1 or len(b_) < 1:
        return np.nan
    return float(mannwhitneyu(a, b_, alternative="greater").pvalue)


PANELS = []


def panel(name, slow=False):
    def wrap(f):
        PANELS.append((name, f, slow))
        return f
    return wrap


def report(msg):
    print(f"    {msg}", flush=True)


# ------------------------------------------------------------------------------------------- Figure 1
def _read_uni(path):
    """PINTS unidirectional peaks, as in 12_itss.read_uni (column 15 = summit)."""
    import gzip
    out = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            s, e = int(f[1]), int(f[2])
            summit = int(f[14]) if len(f) > 14 and f[14].lstrip("-").isdigit() else (s + e) // 2
            if not (s <= summit < e):
                summit = (s + e) // 2
            out.append((f[0], f[5], summit, "uni", np.nan))
    return out


def _read_bi(path):
    """PINTS bidirectional peaks, as in 12_itss.read_bi (columns 5 and 6 = the two starts)."""
    import gzip
    out = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[4]), int(f[5]), f[3]))
    return out


def pints_starts():
    """All K562 PRO-cap starts: unidirectional summits + both arms of bidirectional peaks (helpers/starts_breakdown.py)."""
    uni, bi = _read_uni(need(PINTS["uni"])), _read_bi(need(PINTS["bi"]))
    starts = [(c, s, p, "single") for c, s, p, _, _ in uni]
    for c, p5, p6, _ in bi:
        starts += [(c, "+", p5, "pair"), (c, "-", p6, "pair")]
    return starts


@panel("fig1_b_waffle", slow=True)
def fig1_b_waffle():
    """Fig. 1b, left: genomic location of all K562 PRO-cap starts.
    Replicates scripts/helpers/starts_breakdown.py (7 categories; the figure merges 2+4+5 into 'elsewhere in genes,
    on their strand').  Inputs: GENCODE GTF and the two PINTS peak files (this is the only panel that needs them).
    Squares: 200 (0.5% each), largest-remainder rounding."""
    c10 = load_script("10_candidates.py")
    tx, exons = c10.parse_gtf(need(GTF))
    tss, ex, pc, oth = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
    for t in tx.values():
        k = (t["chrom"], t["strand"])
        tss[k].append(c10.tss_of(t))
        ex[k] += exons[t["transcript_id"]]
        (pc if t["gene_type"] == "protein_coding" else oth)[k].append((t["start"], t["end"]))
    tss = {k: sorted(v) for k, v in tss.items()}
    EX = {k: c10.IntervalSet(v) for k, v in ex.items()}
    PC = {k: c10.IntervalSet(v) for k, v in pc.items()}
    OT = {k: c10.IntervalSet(v) for k, v in oth.items()}
    E = c10.IntervalSet([])
    import bisect

    def cat(c, s, p):
        k, o = (c, s), (c, "-" if s == "+" else "+")
        tl = tss.get(k, [])
        j = bisect.bisect_left(tl, p - 500)
        if j < len(tl) and tl[j] <= p + 500:
            return 1
        if EX.get(k, E).dist(p) == 0:
            return 2
        if PC.get(k, E).overlaps(p, p + 1):
            return 3 if EX.get(k, E).dist(p) >= 1000 else 4
        if OT.get(k, E).overlaps(p, p + 1):
            return 5
        if PC.get(o, E).overlaps(p, p + 1) or OT.get(o, E).overlaps(p, p + 1):
            return 6
        return 7

    starts = pints_starts()
    cnt = defaultdict(int)
    for c, s, p, _ in starts:
        cnt[cat(c, s, p)] += 1
    tot = len(starts)
    groups = [("tss", "at annotated transcript starts", [1]), ("anti", "antisense to genes", [6]),
              ("inter", "between genes", [7]), ("other", "elsewhere in genes, on their strand", [2, 4, 5]),
              ("itss", "in introns, host strand, >=1 kb from exons", [3])]
    rows = [dict(category=g, label=lab, categories="+".join(map(str, cs)), n_starts=sum(cnt[c] for c in cs),
                 total=tot) for g, lab, cs in groups]
    d = pd.DataFrame(rows)
    d["pct"] = 100 * d.n_starts / tot
    raw = d.n_starts / tot * 200
    sq = np.floor(raw).astype(int)
    for i in np.argsort(-(raw - sq).values)[: 200 - sq.sum()]:
        sq.iloc[i] += 1
    d["squares"] = sq
    report(f"{tot} PRO-cap starts; " + ", ".join(f"{r.category} {r.pct:.1f}%" for r in d.itertuples()))
    detail = pd.DataFrame([dict(category_no=k, n_starts=v, pct=100 * v / tot) for k, v in sorted(cnt.items())])
    return {"fig1_b_waffle.tsv": d, "fig1_b_waffle_7categories.tsv": detail}


@panel("fig1_b_flow")
def fig1_b_flow():
    """Fig. 1b, right: intronic start sites through the analysis (step 11 evidence tables; step 12 linked_annot).
    'All PRO-cap start sites' (K562) = unidirectional summits + 2 x bidirectional peaks (PINTS files)."""
    rows = []
    for line, base in LINES.items():
        e = read(f"{base}/itss.evidence.tsv")
        c = e[e.site_type == "candidate"]
        sup, act = b(c.tss_supported), b(c.host_active)
        vals = dict(intronic=len(c), capped=int(sup.sum()), capped_active=int((sup & act).sum()),
                    capped_silent=int((sup & ~act).sum()), near_annot=int(b(c.linked_annot).sum()))
        rows += [dict(line=line, step=k, n=v) for k, v in vals.items()]
    try:
        n = len(pints_starts())
    except MissingInput as m:
        report(f"PINTS peak files not found ({m}); 'All PRO-cap start sites' left empty")
        n = np.nan
    rows.append(dict(line="K562", step="all_starts", n=n))
    rows.append(dict(line="HCT116", step="all_starts", n=np.nan))
    d = pd.DataFrame(rows)
    report("; ".join(f"{r.line} {r.step} {r.n}" for r in d.itertuples() if r.step != "all_starts"))
    return {"fig1_b_flow.tsv": d}


@panel("fig1_c")
def fig1_c():
    """Fig. 1c: nested protein-coding genes (either strand, any host transcript) vs stand-alone protein-coding genes,
    % with a detected start by expression bin (step 30 with --orient both --introns any: known_genes_wide.tsv)."""
    k = read("results/known_genes_wide.tsv")
    rows, summ = [], []
    pcm = k["type"] == "protein_coding"
    for line in LINES:
        det, binc = b(k[f"det_{line}"]), num(k[f"bin_{line}"])
        for kind, m in (("nested", pcm & (k.kind == "nested")), ("alone", pcm & (k.kind == "standalone"))):
            for bb in range(6):
                mm = m & (binc == bb)
                rows.append(dict(line=line, kind=kind, bin=bb, n=int(mm.sum()), k=int(det[mm].sum()),
                                 pct=100 * det[mm].mean() if mm.sum() else np.nan))
        nest = pcm & (k.kind == "nested")
        clear = nest & b(k[f"clear_{line}"])
        summ.append(dict(line=line, n_nested_pc=int(nest.sum()), n_nested_pc_same=int((nest & (k.orient == "same")).sum()),
                         n_nested_pc_anti=int((nest & (k.orient == "anti")).sum()), n_clear=int(clear.sum()),
                         n_clear_detected=int(det[clear].sum())))
        report(f"{line}: {summ[-1]['n_clear_detected']} of {summ[-1]['n_clear']} clearly expressed nested genes detected; "
               f"{summ[-1]['n_nested_pc']} nested protein-coding genes")
    return {"fig1_c.tsv": pd.DataFrame(rows), "fig1_c_summary.tsv": pd.DataFrame(summ)}


# ------------------------------------------------------------------------------------------- Figure 2
@panel("fig2_a")
def fig2_a():
    """Fig. 2a: % of sites with a histone peak within 1 kb (step 21, itss.describe.tsv; groups as in 21_describe)."""
    rows = []
    for line, base in LINES.items():
        d = read(f"{base}/itss.describe.tsv")
        st, act, sup = d.site_type.astype(str), b(d.host_active), b(d.tss_supported)
        groups = {"host": st == "host_promoter", "int": (st == "internal_promoter") & act & sup,
                  "itss": (st == "candidate") & act & sup, "rand": (st == "random_control") & act}
        for mark in ("H3K27ac", "H3K4me3", "H3K4me1"):
            v = b(d[mark])
            for g, m in groups.items():
                rows.append(dict(line=line, mark=mark, group=g, n=int(m.sum()), k=int(v[m].sum()),
                                 pct=100 * v[m].mean() if m.sum() else np.nan))
    out = pd.DataFrame(rows)
    report("n per group: " + "; ".join(f"{r.line} {r.group} {r.n}" for r in out[out.mark == "H3K27ac"].itertuples()))
    return {"fig2_a.tsv": out}


@panel("fig2_b")
def fig2_b():
    """Fig. 2b: transposable-element class at the start nucleotide, log2(share of intronic start sites / share of
    random intronic positions); confirmed sites in active hosts vs random positions in active hosts
    (step 28, repeats.<line>.tsv, column grupa_rep)."""
    cls = {"none": "brak", "alu": "SINE/Alu", "l1": "LINE/L1", "dna": "DNA", "ltr": "LTR/ERV"}
    rows = []
    from scipy.stats import fisher_exact
    for line in LINES:
        d = read(f"results/repeats.{line}.tsv")
        st, act, sup = d.site_type.astype(str), b(d.host_active), b(d.tss_supported)
        it, ct = d[(st == "candidate") & act & sup], d[(st == "random_control") & act]
        for c, g in cls.items():
            a1, a2 = int((it.grupa_rep == g).sum()), int((ct.grupa_rep == g).sum())
            p = fisher_exact([[a1, len(it) - a1], [a2, len(ct) - a2]])[1]
            rows.append(dict(line=line, cls=c, rep_group=g, n_itss=len(it), k_itss=a1, share_itss=100 * a1 / len(it),
                             n_rand=len(ct), k_rand=a2, share_rand=100 * a2 / len(ct),
                             log2_ratio=np.log2((a1 / len(it)) / (a2 / len(ct))) if a1 and a2 else np.nan, p_fisher=p))
        te_i = ~it.grupa_rep.isin(["brak", "proste/niska złożoność"])
        te_c = ~ct.grupa_rep.isin(["brak", "proste/niska złożoność"])
        rows.append(dict(line=line, cls="any_te", rep_group="not brak / simple", n_itss=len(it), k_itss=int(te_i.sum()),
                         share_itss=100 * te_i.mean(), n_rand=len(ct), k_rand=int(te_c.sum()),
                         share_rand=100 * te_c.mean(), log2_ratio=np.log2(te_i.mean() / te_c.mean()), p_fisher=np.nan))
        report(f"{line}: in a transposable element {100 * te_i.mean():.0f}% vs random {100 * te_c.mean():.0f}%")
    return {"fig2_b.tsv": pd.DataFrame(rows)}


AP1 = ("FOS", "FRA", "JUN", "ATF3", "BATF", "AP-1")
MARE = ("BACH", "NRF2", "NF-E2", "NFE2", "MAF")


def motif_name(full):
    """Short motif name from a HOMER known-motif name, e.g. 'Gata6(Zf)/HUVEC-Gata6-ChIP-Seq(GSE51936)/Homer' -> GATA6."""
    n = full.split("/")[0]
    n = n.split("(")[0].strip()
    n = re.sub(r"-fusion$", "", n)
    return n.upper()


def motif_family(full, name):
    fam = re.search(r"\(([^)]*)\)", full.split("/")[0])
    fam = fam.group(1) if fam else ""
    if name.startswith("GATA") or name == "TRPS1":
        return "GATA"
    if name.startswith(MARE):
        return "MARE"
    if name.startswith(AP1):
        return "AP-1"
    if "ETS" in fam:
        return "ETS"
    if name.startswith("KLF"):
        return "KLF"
    return "other"


@panel("fig2_c")
def fig2_c():
    """Fig. 2c: HOMER known motifs, top 25 per cell line (jobs/homer.sh: results/motifs/<line>/knownResults.txt,
    intronic start sites, confirmed, active host, 300 nt, vs random intronic positions in active hosts).
    -log10 p from the 'Log P-value' column (natural log); short names and families are derived from the HOMER
    motif names by rule (motif_name / motif_family)."""
    rows = []
    for line in LINES:
        p = need(f"results/motifs/{line}/knownResults.txt")
        k = pd.read_csv(p, sep="\t")
        cols = list(k.columns)
        logp = next(c for c in cols if c.lower().startswith("log p"))
        tgt = next(c for c in cols if c.startswith("% of Target"))
        bg = next(c for c in cols if c.startswith("% of Background"))
        k = k.head(25).reset_index(drop=True)
        for i, r in k.iterrows():
            full = str(r[cols[0]])
            nm = motif_name(full)
            tg, bgv = float(str(r[tgt]).rstrip("%")), float(str(r[bg]).rstrip("%"))
            rows.append(dict(line=line, rank=i + 1, motif=full, name=nm, family=motif_family(full, nm),
                             log10p=-float(r[logp]) / math.log(10), pct_target=tg, pct_background=bgv,
                             log2_enrichment=np.log2(tg / bgv) if bgv > 0 else np.nan))
    d = pd.DataFrame(rows)
    d["top_in_family"] = False
    for (line, fam), g in d.groupby(["line", "family"]):
        if fam != "other":
            d.loc[g.log10p.idxmax(), "top_in_family"] = True
    report("top motif: " + "; ".join(f"{l} {g.iloc[0]['name']} (-log10 p {g.iloc[0].log10p:.0f})" for l, g in d.groupby("line")))
    return {"fig2_c.tsv": d}


def mh_or(tables):
    """Mantel-Haenszel OR, Robins-Breslow-Greenland 95% CI and CMH p (scripts/44_or_ci.py, verbatim)."""
    R = S = PR = PSQR = QS = num_ = den = 0.0
    for a, b_, c, d in tables:
        n = a + b_ + c + d
        if n < 2:
            continue
        r, s = a * d / n, b_ * c / n
        p, q = (a + d) / n, (b_ + c) / n
        R += r; S += s; PR += p * r; PSQR += p * s + q * r; QS += q * s
        num_ += a - (a + b_) * (a + c) / n
        den += (a + b_) * (c + d) * (a + c) * (b_ + d) / (n * n * (n - 1))
    o = R / S
    se = np.sqrt(PR / (2 * R ** 2) + PSQR / (2 * R * S) + QS / (2 * S ** 2))
    return o, np.exp(np.log(o) - 1.96 * se), np.exp(np.log(o) + 1.96 * se), chi2.sf((abs(num_) - 0.5) ** 2 / den, 1)


@panel("fig2_d")
def fig2_d():
    """Fig. 2d: overlap of confirmed sites (all hosts) and CMH odds ratios, shared vs cell-type-specific sites in
    active hosts within thirds of PRO-cap reads (Mantel-Haenszel OR with Robins-Breslow-Greenland CI as in
    scripts/44_or_ci.py, on Supplementary Table 1 v2).  'transposon' follows 31_final_tables.py (any element except
    none and simple/low complexity); 44_or_ci.py's 'TE classes' variant is
    written as feature transposon_te_classes."""
    TE = ("LINE", "SINE", "LTR", "DNA", "Retroposon", "RC")
    t = read(SUPP)
    venn, ors = [], []
    for line in LINES:
        x = t[t.cell_line == line]
        conf, sh = b(x.capped_start_confirmed), b(x.shared_with_other_cell_line)
        venn.append(dict(line=line, n_confirmed=int(conf.sum()), n_shared=int((conf & sh).sum()),
                         n_specific=int((conf & ~sh).sum())))
        x = x[conf & b(x.host_active)].copy()
        x["third"] = pd.qcut(num(x.procap_reads).rank(method="first"), 3, labels=False)
        sh = b(x.shared_with_other_cell_line)
        rc = x.repeat_class.astype(str)
        # transposon as in 31_final_tables.py and Fig. 2b (any RepeatMasker element except simple/low complexity);
        # 44_or_ci.py also printed a variant restricted to TE classes (LINE/SINE/LTR/DNA/...): transposon_te_classes
        feats = {"conserved": num(x.phylop_mean) > 1,
                 "transposon": ~rc.isin(["none", "simple/low complexity"]),
                 "transposon_te_classes": rc.str.split("/").str[0].isin(TE)}
        for f, y in feats.items():
            tabs = []
            for th in range(3):
                m = x.third == th
                tabs.append((int((sh & y & m).sum()), int((sh & ~y & m).sum()), int((~sh & y & m).sum()), int((~sh & ~y & m).sum())))
            o, lo, hi, p = mh_or(tabs)
            ors.append(dict(line=line, feature=f, **{"or": o}, ci_low=lo, ci_high=hi, p=p, n_shared=int(sh.sum()),
                            n_specific=int((~sh).sum())))
            report(f"{line} {f}: OR {o:.2f} ({lo:.2f}-{hi:.2f}), p = {p:.2g}")
    return {"fig2_d_venn.tsv": pd.DataFrame(venn), "fig2_d_or.tsv": pd.DataFrame(ors)}


# ------------------------------------------------------------------------------------------- Figure 3
@panel("fig3_a")
def fig3_a():
    """Fig. 3a: % with stable RNA by third of PRO-cap reads, intronic start sites and internal promoters (confirmed,
    active host), only sites >= 3 kb upstream of the next host exon, threshold = 95th percentile of rna_step of random
    positions meeting the same rule; thirds with common boundaries over both groups (all distances)
    (scripts/36_stable_far.py, section 'Fig. 3a', column '>= 3 kb')."""
    rows = []
    for line in LINES:
        d = read(f"results/repeats.{line}.tsv")
        st, act, sup = d.site_type.astype(str), b(d.host_active), b(d.tss_supported)
        for c in ("pos", "intron_start", "intron_end", "cap_sum", "rna_step"):
            d[c] = num(d[c])
        d["dist_down"] = dist_down(d)
        d["far"] = d.dist_down >= FAR
        itss = d[(st == "candidate") & sup & act].copy()
        ip = d[(st == "internal_promoter") & sup & act].copy()
        rnd = d[(st == "random_control") & act].copy()
        thr_far = np.nanpercentile(rnd.loc[rnd.far, "rna_step"], 95)
        e1, e2 = np.nanpercentile(pd.concat([itss.cap_sum, ip.cap_sum]), [100 / 3, 200 / 3])
        for g, x in (("itss", itss), ("int", ip)):
            x["third"] = pd.cut(x.cap_sum, [-np.inf, e1, e2, np.inf], labels=["weak", "medium", "strong"]).astype(str)
            x["stab_far"] = x.rna_step > thr_far
            for th in ("weak", "medium", "strong"):
                y = x[(x.third == th) & x.far]
                rows.append(dict(line=line, group=g, third=th, n=len(y), k=int(y.stab_far.sum()),
                                 pct=100 * y.stab_far.mean() if len(y) else np.nan, n_all_distances=int((x.third == th).sum()),
                                 edge_low=e1, edge_high=e2, threshold=thr_far))
        y = rnd[rnd.far]
        rows.append(dict(line=line, group="rand", third="all", n=len(y), k=int((y.rna_step > thr_far).sum()),
                         pct=100 * (y.rna_step > thr_far).mean(), n_all_distances=len(rnd), edge_low=e1, edge_high=e2,
                         threshold=thr_far))
        report(f"{line}: threshold {thr_far:.2f}, third boundaries {e1:.0f} / {e2:.0f}; strongest third: "
               + ", ".join(f"{r['group']} {r['pct']:.1f}%" for r in rows[-7:] if r["third"] == "strong"))
    return {"fig3_a.tsv": pd.DataFrame(rows)}


@panel("fig3_b")
def fig3_b():
    """Fig. 3b: single-site classes 10 min vs 40 min after DRB release, both replicates (step 27, classes_10.tsv);
    percentages among resolved sites."""
    d = read("results/hct116/classes_10.tsv")
    st, act, sup = d.site_type.astype(str), b(d.host_active), b(d.tss_supported)
    groups = {"rand": (st == "random_control") & act, "int": (st == "internal_promoter") & act & sup,
              "itss": (st == "candidate") & act & sup, "itss_silent": (st == "candidate") & ~act & sup}
    rows = []
    for g, m in groups.items():
        k = d.loc[m, "klasa"].astype(str)
        n = {c: int((k == c).sum()) for c in ("produktywny", "krótki", "bez odpowiedzi", "nierozstrzygnięty")}
        res = n["produktywny"] + n["krótki"] + n["bez odpowiedzi"]
        rows.append(dict(group=g, n_all=int(m.sum()), n_resolved=res, n_productive=n["produktywny"], n_short=n["krótki"],
                         n_none=n["bez odpowiedzi"], n_unresolved=n["nierozstrzygnięty"],
                         pct_productive=100 * n["produktywny"] / res if res else np.nan,
                         pct_short=100 * n["krótki"] / res if res else np.nan,
                         pct_none=100 * n["bez odpowiedzi"] / res if res else np.nan))
        report(f"{g}: productive {rows[-1]['pct_productive']:.1f}%, short {rows[-1]['pct_short']:.1f}% (n = {res})")
    return {"fig3_b.tsv": pd.DataFrame(rows)}


class Cov:
    """bigWig pair by READ strand, as 20_drb_gap.Cov (values on the RNA strand, optional swap)."""

    def __init__(self, plus, minus):
        import pyBigWig
        self.bw = {"plus": pyBigWig.open(plus), "minus": pyBigWig.open(minus)}
        self.flip = False

    def values(self, chrom, s, e, rna_strand):
        key = "plus" if rna_strand == "+" else "minus"
        if self.flip:
            key = "minus" if key == "plus" else "plus"
        bw = self.bw[key]
        L = bw.chroms(chrom)
        out = np.zeros(max(0, e - s))
        if L is None:
            return out
        s2, e2 = max(0, s), min(L, e)
        if e2 > s2:
            v = bw.values(chrom, s2, e2, numpy=True)
            if len(v) == e2 - s2:
                out[s2 - s:e2 - s] = np.abs(np.nan_to_num(v))
        return out


def along(pos, strand, a, b_):
    return (pos + a, pos + b_) if strand == "+" else (pos - b_ + 1, pos - a + 1)


def drb_profile_bigwig(drb_tsv, test, ref):
    """Meta-profile test/ref in 1-kb bins from -10 to +10 kb, summed over sites (20_drb_gap.py, 'Profil')."""
    ev = read(drb_tsv)
    for p in test + ref:
        need(p)
    T, R = Cov(*test), Cov(*ref)
    hosts = ev[ev.host_active.astype(str) == "True"].dropna(subset=["host_tss_used"]).drop_duplicates("host_gene_id")
    same = opp = 0.0
    for r in hosts.itertuples():
        s, e = along(int(r.host_tss_used), r.strand, 0, 5000)
        same += R.values(r.chrom, s, e, r.strand).sum()
        opp += R.values(r.chrom, s, e, "-" if r.strand == "+" else "+").sum()
    T.flip = R.flip = opp > same
    act, sup, el = b(ev.host_active), b(ev.tss_supported), b(ev.eligible)
    st = ev.site_type.astype(str)
    groups = {"itss": (st == "candidate") & act & sup & el, "itss_silent": (st == "candidate") & ~act & sup & el,
              "int": (st == "internal_promoter") & act & sup & el, "rand": (st == "random_control") & act & el}
    rows = []
    for g, m in groups.items():
        sub = ev[m]
        D, B = np.zeros(20), np.zeros(20)
        for r in sub.itertuples():
            s, e = along(int(r.pos), r.strand, -10000, 10000)
            dv, bv = T.values(r.chrom, s, e, r.strand), R.values(r.chrom, s, e, r.strand)
            if r.strand == "-":
                dv, bv = dv[::-1], bv[::-1]
            D += dv.reshape(20, 1000).sum(1)
            B += bv.reshape(20, 1000).sum(1)
        ratio = D / np.maximum(B, 1e-9)
        for i in range(20):
            rows.append(dict(group=g, kb_from=i - 10, kb_to=i - 9, kb_mid=i - 9.5, ratio=ratio[i], n_sites=len(sub),
                             source="bigwig"))
    return pd.DataFrame(rows), f"flip={T.flip}"


def drb_profile_summary(path):
    """Same profile read from the text file written by 20_drb_gap.py (<out-prefix>.summary.txt)."""
    gm = {"iTSS_akt": "itss", "iTSS_nieakt": "itss_silent", "wewn_prom": "int", "kontrole": "rand"}
    lines = open(need(path), encoding="utf-8").read().splitlines()
    i = next(j for j, l in enumerate(lines) if l.startswith("Profil"))
    ns = dict(re.findall(r"(\w+) n=(\d+)", lines[i]))
    hdr = lines[i + 1].split()[1:]
    rows = []
    for l in lines[i + 2:i + 22]:
        f = l.split()
        a = int(f[0].split("..")[0])
        for name, v in zip(hdr, f[1:]):
            rows.append(dict(group=gm.get(name, name), kb_from=a, kb_to=a + 1, kb_mid=a + 0.5, ratio=float(v),
                             n_sites=int(ns.get(name, -1)), source="summary.txt"))
    return pd.DataFrame(rows)


def profile_panel(drb_tsv, summary, test, ref):
    out, note = None, ""
    bw_ok = not ARGS.no_bigwig
    if bw_ok:
        try:
            import pyBigWig  # noqa: F401
        except ImportError:
            bw_ok, note = False, "pyBigWig not importable"
    if bw_ok and all(os.path.exists(p) for p in test + ref):
        out, note = drb_profile_bigwig(drb_tsv, test, ref)
    elif bw_ok:
        note = "bigWig files not found: " + ", ".join(p for p in test + ref if not os.path.exists(p))
    try:
        sm = drb_profile_summary(summary)
    except MissingInput:
        sm = None
    if out is None:
        if sm is None:
            raise MissingInput(f"{summary} (and no bigWig: {note})")
        report(f"profile from {summary} ({note or 'no bigWig'})")
        return sm
    if sm is not None:
        m = out.merge(sm, on=["group", "kb_from"], suffixes=("", "_txt"))
        dmax = (m.ratio.round(2) - m.ratio_txt).abs().max()
        nok = (m.n_sites == m.n_sites_txt).all()
        report(f"profile recomputed from bigWig ({note}); vs {os.path.basename(summary)}: max |diff| {dmax:.3f}, "
               f"n identical: {nok}")
        out["ratio_summary_txt"] = m.ratio_txt.values
    return out


@panel("fig3_c")
def fig3_c():
    """Fig. 3c: nascent RNA 10 or 20 min / 40 min around sites beyond the host's wave front, summed over sites, 1-kb
    bins (20_drb_gap.py meta-profile; jobs/hct116_drb.sh: --min-dist 40000 / 70000, groups use 'eligible').
    The per-bin values are not in any table: recomputed from data/hct116_drb/TT_*min.read_*.bw with the eligible flags
    of ttdrb_<m>.drb.tsv, or read from results/hct116/ttdrb_<m>.summary.txt (written by 20_drb_gap.py)."""
    D = "data/hct116_drb"
    outs = []
    for m in (10, 20):
        p = profile_panel(f"results/hct116/ttdrb_{m}.drb.tsv", f"results/hct116/ttdrb_{m}.summary.txt",
                          [f"{D}/TT_{m}min.read_plus.bw", f"{D}/TT_{m}min.read_minus.bw"],
                          [f"{D}/TT_40min.read_plus.bw", f"{D}/TT_40min.read_minus.bw"])
        p.insert(0, "minutes", m)
        outs.append(p)
    d = pd.concat(outs, ignore_index=True)
    report("n sites: " + "; ".join(f"{m} min {g} {n}" for (m, g), n in d.groupby(["minutes", "group"]).n_sites.first().items()))
    return {"fig3_c.tsv": d}


# ------------------------------------------------------------------------------------------- Figure 4
@panel("fig4_a")
def fig4_a():
    """Fig. 4a: % of sites whose stable-RNA rise increases in MTREX knockout cells (delta) above the 95th percentile
    of random positions; covered sites; 3-kb rule for active hosts (incl. the threshold); silent hosts without the rule;
    one-sided Mann-Whitney vs random positions (scripts/38_fig4a_fig5a.py, 3-kb rule section)."""
    e = read("results/v2/exosome/exo.MTREX.tsv")
    e = e[b(e.covered)].copy() if "covered" in e else e
    e["dd"] = dist_down(e)
    e["v"] = num(e["delta"])
    st, act, sup = e.site_type.astype(str), b(e.host_active), b(e.tss_supported)
    typ = e["typ"].astype(str) if "typ" in e else pd.Series("", index=e.index)
    far = e.dd >= FAR
    rnd = e[(st == "random_control") & act & far].v.dropna()
    thr = np.percentile(rnd, 95)
    groups = [("rand", (st == "random_control") & act & far), ("int", (st == "internal_promoter") & sup & act & far),
              ("itss", (st == "candidate") & sup & act & far),
              ("prom", (st == "candidate") & sup & act & far & typ.str.startswith("promotor")),
              ("enh", (st == "candidate") & sup & act & far & typ.str.startswith("enhancer")),
              ("silent", (st == "candidate") & sup & ~act)]
    rows = []
    for g, m in groups:
        v = e[m].v.dropna()
        rows.append(dict(group=g, n=len(v), k=int((v > thr).sum()), pct=100 * (v > thr).mean(),
                         p=np.nan if g == "rand" else (mwu_greater(v, rnd) if len(v) > 2 else np.nan), threshold=thr))
    d = pd.DataFrame(rows)
    report(f"threshold {thr:.3f}; " + "; ".join(f"{r.group} {r.pct:.1f}% (n = {r.n})" for r in d.itertuples()))
    return {"fig4_a.tsv": d}


def u1pas_groups(u):
    st, act, sup = u.site_type.astype(str), b(u.host_active), b(u.tss_supported)
    return {"host": st == "host_promoter", "int": (st == "internal_promoter") & sup & act,
            "intall": st == "internal_promoter", "itss": (st == "candidate") & sup & act,
            "rand": (st == "random_control") & act}


@panel("fig4_b")
def fig4_b():
    """Fig. 4b: % of sites with a splice-donor motif before the first PAS in the first kb (or no PAS)
    (step 25, u1pas.tsv, column protected; groups as in 41_register_checks.py section 4)."""
    rows = []
    for line, base in LINES.items():
        u = read(f"{base}/u1pas.tsv")
        pr = b(u.protected)
        for g, m in u1pas_groups(u).items():
            rows.append(dict(line=line, group=g, n=int(m.sum()), k=int(pr[m].sum()), pct=100 * pr[m].mean()))
    d = pd.DataFrame(rows)
    report("; ".join(f"{r.line} {r.group} {r.pct:.1f}%" for r in d.itertuples()))
    return {"fig4_b.tsv": d}


@panel("fig4_c")
def fig4_c():
    """Fig. 4c: PAS per kb in the first kb, transcribed vs opposite strand, and their ratio of mean densities
    (u1pas.tsv as in 41 section 4; start sites between genes from context_<line>.tsv, group miedzygen, step 29)."""
    rows = []
    for line, base in LINES.items():
        u = read(f"{base}/u1pas.tsv")
        for g, m in u1pas_groups(u).items():
            if g == "intall":
                continue
            x = u[m]
            rows.append(dict(line=line, group=g, n=len(x), pas_transcribed=num(x.pas_sens).mean(),
                             pas_opposite=num(x.pas_anty).mean()))
        c = read(f"results/context_{line}.tsv")
        x = c[c.grupa == "miedzygen"]
        rows.append(dict(line=line, group="inter", n=len(x), pas_transcribed=num(x.pas_sens).mean(),
                         pas_opposite=num(x.pas_anty).mean()))
    d = pd.DataFrame(rows)
    d["ratio"] = d.pas_transcribed / d.pas_opposite
    report("ratio: " + "; ".join(f"{r.line} {r.group} {r.ratio:.2f}" for r in d.itertuples()))
    return {"fig4_c.tsv": d}


# ------------------------------------------------------------------------------------------- Figure 5
@panel("fig5_a")
def fig5_a():
    """Fig. 5a: share of sites with a positive stable-RNA excess (Wilson 95% CI) and median excess per PRO-cap read,
    unannotated sites passing the cap test (step 29, context_<line>.tsv), 3-kb rule for intronic sites in active hosts,
    internal promoters and random positions (distance from itss.describe.tsv) (scripts/39_fig5a.py).
    v = RNA_nadw, the excess of RNA-seq coverage after the site per PRO-cap read (step 29).  Bootstrap CIs of the
    medians: 1,000 resamples, default_rng(1), consumed in the order of 39_fig5a.py."""
    rng = np.random.default_rng(1)

    def boot_med(v, nb=1000):
        v = np.asarray(v, float)
        if len(v) < 3:
            return np.nan, np.nan
        m = np.median(rng.choice(v, (nb, len(v)), replace=True), axis=1)
        return np.percentile(m, 2.5), np.percentile(m, 97.5)

    slug = {"intron_akt": "act", "intron_cichy": "sil", "miedzygen": "inter", "wewn_prom": "prom", "kontrole": "rand"}
    rows = []
    for line, base in LINES.items():
        c = read(f"results/context_{line}.tsv")
        d = read(f"{base}/itss.describe.tsv", usecols=["site_id", "strand", "pos", "intron_start", "intron_end"])
        d["dd"] = dist_down(d)
        c = c.merge(d[["site_id", "dd"]], on="site_id", how="left")
        cap = num(c.cap_sum)
        c["v"] = np.where(cap > 0, num(c.RNA_nadw), np.nan)
        for g, x in c.groupby(c.grupa.astype(str)):
            x = x[np.isfinite(x.v)]
            gl = g.lower()
            rule = (("akt" in gl and "nieakt" not in gl) or "wewn" in gl or "prom" in gl or "kontrol" in gl) and x.dd.notna().any()
            y = x[x.dd >= FAR] if rule else x
            boot_med(x.v.values)  # consumes the generator as 39_fig5a.py does (all distances first)
            k, n = int((y.v > 0).sum()), len(y)
            lo, hi = wilson(k, n)
            mlo, mhi = boot_med(y.v.values)
            rows.append(dict(line=line, group=slug.get(g, g), grupa=g, three_kb_rule=rule, n=n, k_positive=k,
                             pct_positive=100 * k / n if n else np.nan, pct_lo=lo, pct_hi=hi,
                             median_x1e6=1e6 * np.median(y.v) if n else np.nan, median_lo_x1e6=1e6 * mlo,
                             median_hi_x1e6=1e6 * mhi, median_rna_nadw_x1e6=1e6 * np.nanmedian(num(y.RNA_nadw)),
                             n_all_distances=len(x)))
        report(f"{line}: " + "; ".join(f"{r['group']} {r['pct_positive']:.1f}% (n = {r['n']}), median {r['median_x1e6']:.2f}"
                                       for r in rows if r["line"] == line))
    return {"fig5_a.tsv": pd.DataFrame(rows)}


@panel("fig5_b")
def fig5_b():
    """Fig. 5b: nascent RNA 2-8 kb after the site per PRO-cap read, 10 min after DRB release (step 29,
    context_HCT116.tsv, DRB_dal; set only for sites > 22.2 kb from the host start in use or in silent hosts),
    median x 1000 with a 95% bootstrap CI (1,000 resamples, default_rng(1), groups in the order drawn).
    If results/final.fig5_robust.tsv exists, its values are copied to the columns robust_* for comparison."""
    c = read("results/context_HCT116.tsv")
    rng = np.random.default_rng(1)
    rows = []
    for g, slug in (("intron_akt", "act"), ("intron_cichy", "sil"), ("miedzygen", "inter"), ("wewn_prom", "prom")):
        v = num(c.loc[c.grupa == g, "DRB_dal"]).dropna().values
        m = np.median(rng.choice(v, (1000, len(v)), replace=True), axis=1) if len(v) > 2 else np.array([np.nan])
        rows.append(dict(group=slug, grupa=g, n=len(v), median_x1e3=1e3 * np.median(v) if len(v) else np.nan,
                         lo_x1e3=1e3 * np.percentile(m, 2.5), hi_x1e3=1e3 * np.percentile(m, 97.5)))
    d = pd.DataFrame(rows)
    rob = "results/final.fig5_robust.tsv"
    if os.path.exists(rob):
        r = read(rob)
        r = r[r.measure.astype(str).str.contains("nascent|DRB", case=False) & (r.third.astype(str) == "all")
              & (r.cell_line == "HCT116")]
        gm = {"intronic, active host": "act", "intronic, silent host": "sil", "between genes": "inter",
              "internal promoters": "prom"}
        r = r.assign(group=r.group.map(gm))[["group", "n", "median", "median_lo", "median_hi"]]
        r.columns = ["group", "robust_n", "robust_median", "robust_lo", "robust_hi"]
        d = d.merge(r, on="group", how="left")
        report("final.fig5_robust.tsv nascent rows: " + "; ".join(
            f"{x.group} n {x.robust_n} median {x.robust_median:.3f} [{x.robust_lo:.3f}, {x.robust_hi:.3f}]"
            for x in d.itertuples() if pd.notna(x.robust_n)))
    report("; ".join(f"{x.group} {x.median_x1e3:.3f} [{x.lo_x1e3:.3f}, {x.hi_x1e3:.3f}] n = {x.n}" for x in d.itertuples()))
    return {"fig5_b.tsv": d}


# ------------------------------------------------------------------------------------- Supplementary
@panel("figS1")
def figS1():
    """Supp. Fig. S1: cap test, % passing in active hosts with Wilson 95% CI (step 34, section S1, from the step 11
    evidence tables; annotated intronic RNAs from annot.evidence.tsv)."""
    rows = []
    for line, base in LINES.items():
        e = read(f"{base}/itss.evidence.tsv")
        a = read(f"{base}/annot.evidence.tsv")
        sets = {"rand": e[e.site_type == "random_control"], "annot": a[a.site_type == "candidate"],
                "int": e[e.site_type == "internal_promoter"], "itss": e[e.site_type == "candidate"]}
        for g, x in sets.items():
            for host, m in (("active", b(x.host_active)), ("silent", ~b(x.host_active))):
                y = b(x.loc[m, "tss_supported"])
                k, n = int(y.sum()), len(y)
                lo, hi = wilson(k, n)
                rows.append(dict(line=line, group=g if host == "active" else g + "_silent", host=host, n=n, k=k,
                                 pct=100 * k / n if n else np.nan, lo=lo, hi=hi))
    d = pd.DataFrame(rows)
    report("; ".join(f"{r.line} {r.group} {r.pct:.2f}% ({r.k}/{r.n})" for r in d[d.host == "active"].itertuples()))
    return {"figS1.tsv": d}


@panel("figS2")
def figS2():
    """Supp. Fig. S2: confirmed HCT116 sites (CTCF line, active host) passing the cap test in the RAD21 line, by
    third of CTCF-line reads, Wilson 95% CI; Spearman correlation of reads at sites passing in both
    (step 34, section S2: results/supp/supp.S2_sites.tsv)."""
    x = read("results/supp/supp.S2_sites.tsv")
    x["pass"] = b(x.rad21_pass)
    rows = []
    for th in ("weak", "medium", "strong", "all"):
        g = x if th == "all" else x[x.third.astype(str) == th]
        k, n = int(g["pass"].sum()), len(g)
        lo, hi = wilson(k, n)
        rows.append(dict(third=th, n=n, k=k, pct=100 * k / n, lo=lo, hi=hi, median_reads=num(g.cap_sum).median()))
    sh = x[x["pass"]]
    r, p = spearmanr(num(sh.cap_sum), num(sh.rad21_reads))
    report(f"all {rows[-1]['pct']:.1f}% ({rows[-1]['k']}/{rows[-1]['n']}); Spearman r = {r:.2f} (n = {len(sh)})")
    return {"figS2.tsv": pd.DataFrame(rows),
            "figS2_cor.tsv": pd.DataFrame([dict(set="passing_both", spearman_r=r, p=p, n=len(sh))])}


@panel("figS3")
def figS3():
    """Supp. Fig. S3: the DRB wave along 1,000 active host genes, nascent RNA / reference in 2-kb bins
    (step 34, section S3: results/supp/supp.S3_wave.tsv).  'shown' marks the bins drawn (and printed by step 34):
    every bin up to 20 kb, then every other bin."""
    d = read("results/supp/supp.S3_wave.tsv")
    d["kb_mid"] = (num(d.kb_from) + num(d.kb_to)) / 2
    i = (num(d.kb_from) // 2).astype(int)
    d["shown"] = (i < 10) | (i % 2 == 0)
    report("labels: " + ", ".join(f"{l}/{s}" for l, s in d[["cell_line", "label"]].drop_duplicates().itertuples(index=False)))
    return {"figS3.tsv": d.rename(columns={"cell_line": "line"})}


@panel("figS4")
def figS4():
    """Supp. Fig. S4 (K562, BruDRB-seq 10 min / Bru-seq; results/v2/itss_w8.drb.tsv, 20_drb_gap.py, --win 8000 --gap 5000,
    --min-dist 20000).  (a) meta-profile as in Fig. 3c (bigWig recompute or results/v2/itss_w8.summary.txt);
    (b) share of testable sites above the 95th percentile of random positions (active host, testable),
    Wilson 95% CI, one-sided Mann-Whitney vs random positions (empirical calibration section of 20_drb_gap.py)."""
    B = "results/brudrb"
    prof = profile_panel("results/v2/itss_w8.drb.tsv", "results/v2/itss_w8.summary.txt",
                         [f"{B}/K562_DRB10.read_plus.bw", f"{B}/K562_DRB10.read_minus.bw"],
                         [f"{B}/K562_Bru10.read_plus.bw", f"{B}/K562_Bru10.read_minus.bw"])
    d = read("results/v2/itss_w8.drb.tsv")
    st, act, sup, te = d.site_type.astype(str), b(d.host_active), b(d.tss_supported), b(d.testable)
    gi = num(d.gap_index)
    ctrl = gi[(st == "random_control") & act & te]
    thr = float(np.percentile(ctrl, 95))
    groups = {"rand": (st == "random_control") & act, "int": (st == "internal_promoter") & act & sup,
              "itss": (st == "candidate") & act & sup, "sil": (st == "candidate") & ~act & sup}
    rows = []
    for g, m in groups.items():
        v = gi[m & te]
        k, n = int((v > thr).sum()), len(v)
        lo, hi = wilson(k, n)
        rows.append(dict(group=g, n=n, k=k, pct=100 * k / n, lo=lo, hi=hi,
                         p=np.nan if g == "rand" else mwu_greater(v, ctrl), threshold=thr))
    out = pd.DataFrame(rows)
    report(f"threshold {thr:.2f}; " + "; ".join(f"{r.group} {r.pct:.1f}% (n = {r.n})" for r in out.itertuples()))
    return {"figS4_a.tsv": prof, "figS4_b.tsv": out}


def prompt_check_bigwig(ko):
    """PROMPT check of 24_exosome.py (positive control), recomputed from the knockout/control bigWigs."""
    import pyBigWig

    class C:
        def __init__(self, plus, minus):
            self.bw = {"+": pyBigWig.open(plus), "-": pyBigWig.open(minus)}
            self.flip = False

        def mean(self, chrom, s, e, strand):
            if self.flip:
                strand = "-" if strand == "+" else "+"
            bw = self.bw[strand]
            L = bw.chroms(chrom)
            if L is None:
                return 0.0
            s, e = max(0, s), min(L, e)
            if e <= s:
                return 0.0
            v = bw.values(chrom, s, e, numpy=True)
            return float(np.abs(np.nan_to_num(v)).mean()) if len(v) else 0.0

    X = "data/exosome"
    files = [f"{X}/{ko}_KO_plus.bw", f"{X}/{ko}_KO_minus.bw", f"{X}/{ko}_CTRL_plus.bw", f"{X}/{ko}_CTRL_minus.bw"]
    for f in files:
        need(f)
    s = read("results/v2/itss.describe.tsv")
    k_, c_ = C(*files[:2]), C(*files[2:])
    hp = s[s.site_type == "host_promoter"]
    win = along
    same = sum(c_.mean(r.chrom, *win(int(r.pos), r.strand, 0, 5000), r.strand) for r in hp.head(2000).itertuples())
    opp = sum(c_.mean(r.chrom, *win(int(r.pos), r.strand, 0, 5000), "-" if r.strand == "+" else "+")
              for r in hp.head(2000).itertuples())
    k_.flip = c_.flip = opp > same
    eps, pr = 0.01, []
    for r in hp.itertuples():
        p, st = int(r.pos), r.strand
        anti = "-" if st == "+" else "+"
        ps, pe = win(p, st, -2000, -200)
        gs, ge = win(p, st, 0, 3000)
        pk, pc = k_.mean(r.chrom, ps, pe, anti), c_.mean(r.chrom, ps, pe, anti)
        gk, gc = k_.mean(r.chrom, gs, ge, st), c_.mean(r.chrom, gs, ge, st)
        pr.append((np.log2((pk + eps) / (pc + eps)), np.log2((gk + eps) / (gc + eps))))
    pr = np.array(pr)
    dp = pr[:, 0] - pr[:, 1]
    return dict(ko=ko, n_promoters=len(pr), prompt_median=np.median(pr[:, 0]), gene_median=np.median(pr[:, 1]),
                diff_median=np.median(dp), pct_prompt_higher=100 * np.mean(dp > 0), source="bigwig")


def prompt_check_txt(ko):
    """The same four numbers from results/v2/exosome/exo.<KO>.txt (written by 24_exosome.py, 2 decimals)."""
    txt = open(need(f"results/v2/exosome/exo.{ko}.txt"), encoding="utf-8").read()
    # matches the line in both wordings of 24_exosome.py: the original Polish one (files in the data archive)
    # and the English one (re-runs)
    m = re.search(r"PROMPT, (\d+) (?:promotor\w*|promoters?)\): log2FC PROMPT (?:mediana|median) ([-+\d.]+), "
                  r"(?:początek genu|gene start) ([-+\d.]+), (?:różnica|difference) ([-+\d.]+); "
                  r"PROMPT (?:wyżej niż gen w|higher than gene in) (\d+)%", txt)
    if not m:
        raise MissingInput(f"PROMPT line not found in results/v2/exosome/exo.{ko}.txt")
    return dict(ko=ko, n_promoters=int(m.group(1)), prompt_median=float(m.group(2)), gene_median=float(m.group(3)),
                diff_median=float(m.group(4)), pct_prompt_higher=float(m.group(5)), source="exo.txt")


@panel("figS5")
def figS5():
    """Supp. Fig. S5 (K562 exosome knockouts, step 24).  (a) PROMPT check: median log2(KO/control) 200-2,000 nt
    upstream of host promoters (opposite strand) and over the first 3 kb of the gene; the values exist only in
    results/v2/exosome/exo.<KO>.txt, so they are recomputed from data/exosome/<KO>_{KO,CTRL}_{plus,minus}.bw when
    possible.  (b) % of covered sites with delta above the 95th percentile of random positions (active host), all
    distances, one-sided Mann-Whitney vs random positions (exo.<KO>.tsv)."""
    a_rows, b_rows = [], []
    for ko in ("MTREX", "EXOSC10", "EXOSC5"):
        r = None
        if not ARGS.no_bigwig:
            try:
                r = prompt_check_bigwig(ko)
            except (MissingInput, ImportError) as e:
                report(f"{ko}: PROMPT check from bigWig not possible ({e}); reading exo.{ko}.txt")
        try:
            t = prompt_check_txt(ko)
        except MissingInput:
            t = None
        if r is None:
            if t is None:
                raise MissingInput(f"bigWigs and results/v2/exosome/exo.{ko}.txt")
            r = t
        elif t is not None:
            report(f"{ko}: bigWig PROMPT {r['prompt_median']:+.3f} / gene {r['gene_median']:+.3f} / {r['pct_prompt_higher']:.1f}%"
                   f"; exo.{ko}.txt {t['prompt_median']:+.2f} / {t['gene_median']:+.2f} / {t['pct_prompt_higher']:.0f}%")
        a_rows.append(r)
        d = read(f"results/v2/exosome/exo.{ko}.tsv")
        d = d[b(d.covered)]
        st, act, sup = d.site_type.astype(str), b(d.host_active), b(d.tss_supported)
        typ = d.typ.astype(str)
        delta = num(d.delta)
        ctrl = delta[(st == "random_control") & act]
        thr = np.percentile(ctrl, 95)
        cand = st == "candidate"
        groups = {"rand": (st == "random_control") & act, "int": (st == "internal_promoter") & act & sup,
                  "itss": cand & act & sup, "prom": cand & act & sup & (typ == "promotor-podobny"),
                  "enh": cand & act & sup & (typ == "enhancer-podobny"), "sil": cand & ~act & sup}
        for g, m in groups.items():
            v = delta[m]
            b_rows.append(dict(ko=ko, group=g, n=len(v), k=int((v > thr).sum()), pct=100 * (v > thr).mean(),
                               p=np.nan if g == "rand" else mwu_greater(v, ctrl), threshold=thr))
    a = pd.DataFrame(a_rows)
    report("; ".join(f"{r.ko} PROMPT {r.prompt_median:+.2f}, gene {r.gene_median:+.2f}, {r.pct_prompt_higher:.0f}% ({r.source})"
                     for r in a.itertuples()))
    return {"figS5_a.tsv": a, "figS5_b.tsv": pd.DataFrame(b_rows)}


@panel("figS6")
def figS6():
    """Supp. Fig. S6: % of unannotated confirmed sites with RNA_skok above the 95th percentile of all random positions
    (no 3-kb rule), thirds of PRO-cap reads over intronic and intergenic sites (step 34, section S6, context tables)."""
    rows = []
    for line in LINES:
        S = read(f"results/context_{line}.tsv")
        col = "RNA_skok"
        thr = np.nanpercentile(num(S[S.grupa == "kontrole"][col]), 95)
        st3 = S[S.grupa.isin(["intron_akt", "intron_cichy", "miedzygen"])].copy()
        st3["third"] = pd.qcut(num(st3.cap_sum).rank(method="first"), 3, labels=["weak", "medium", "strong"]).astype(str)
        for th in ("weak", "medium", "strong"):
            for g, slug in (("intron_akt", "act"), ("intron_cichy", "sil"), ("miedzygen", "inter")):
                x = num(st3[(st3.grupa == g) & (st3.third == th)][col]).dropna()
                rows.append(dict(line=line, group=slug, third=th, n=len(x), k=int((x > thr).sum()),
                                 pct=100 * (x > thr).mean() if len(x) else np.nan, threshold=thr))
        ip = num(S[S.grupa == "wewn_prom"][col]).dropna()
        rows.append(dict(line=line, group="prom", third="all", n=len(ip), k=int((ip > thr).sum()),
                         pct=100 * (ip > thr).mean(), threshold=thr))
        report(f"{line}: threshold {thr:.2f}; internal promoters {rows[-1]['pct']:.1f}%")
    return {"figS6.tsv": pd.DataFrame(rows)}


@panel("figS7")
def figS7():
    """Supp. Fig. S7: % of genes with a detected start by expression bin, stand-alone vs nested on the host's strand
    vs opposite strand, lncRNA and protein-coding; o/e over all bins with the expectation from stand-alone genes of the
    same type and bin (step 30, known_genes_wide.tsv)."""
    k = read("results/known_genes_wide.tsv")
    rows, oe = [], []
    for line in LINES:
        det, binc = b(k[f"det_{line}"]), num(k[f"bin_{line}"])
        for ty in ("lncRNA", "protein_coding"):
            sa = (k["type"] == ty) & (k.kind == "standalone")
            rate = det[sa].groupby(binc[sa]).mean()
            for g, m in (("stand", sa), ("same", (k["type"] == ty) & (k.kind == "nested") & (k.orient == "same")),
                         ("anti", (k["type"] == ty) & (k.kind == "nested") & (k.orient == "anti"))):
                for bb in range(6):
                    mm = m & (binc == bb)
                    rows.append(dict(line=line, type=ty, group=g, bin=bb, n=int(mm.sum()), k=int(det[mm].sum()),
                                     pct=100 * det[mm].mean() if mm.sum() else np.nan))
                if g != "stand":
                    exp_ = binc[m].map(rate).fillna(0).sum()
                    obs = int(det[m].sum())
                    oe.append(dict(line=line, type=ty, group=g, n=int(m.sum()), observed=obs, expected=exp_,
                                   oe=obs / exp_ if exp_ else np.nan))
    o = pd.DataFrame(oe)
    report("o/e: " + "; ".join(f"{r.line} {r.type} {r.group} {r.oe:.2f}" for r in o.itertuples()))
    return {"figS7.tsv": pd.DataFrame(rows), "figS7_oe.tsv": o}


@panel("figS8")
def figS8():
    """Supp. Fig. S8: median ratio of the two arms of divergent pairs (2^median DI) with 95% bootstrap CI of the median
    (step 26, directionality.tsv), as in scripts/42_extra_checks.py section 3, with the same random number
    sequence (default_rng(42), consumed group by group and by the depth matching).  Depth-matched deficits are
    written to figS8_matched.tsv."""
    RNG = np.random.default_rng(42)

    def boot_median(v, n=1000):
        v = np.asarray(v, float)
        m = [np.median(RNG.choice(v, len(v))) for _ in range(n)]
        return np.percentile(m, 2.5), np.percentile(m, 97.5)

    from scipy.stats import wilcoxon
    rows, match = [], []
    slug = {"between genes": "inter", "intron, active": "active", "intron, silent": "silent",
            "intron, silent (no capped read at host start)": "silent_zero", "active, H3K36me3": "k36",
            "active, no H3K36me3": "nok36"}
    for line, base in LINES.items():
        d = read(f"{base}/directionality.tsv")
        d["tot"] = num(d.sens) + num(d.anty)
        inter = d[d.klasa == "międzygenowa"].dropna(subset=["DI", "tot"])
        intr = d[d.klasa != "międzygenowa"].dropna(subset=["DI", "tot"]).copy()
        hs = num(intr.host_signal)
        groups = {"between genes": inter, "intron, active": intr[hs >= 5], "intron, silent": intr[hs < 5],
                  "intron, silent (no capped read at host start)": intr[hs == 0]}
        k36 = b(groups["intron, active"].k36)
        groups["active, H3K36me3"] = groups["intron, active"][k36]
        groups["active, no H3K36me3"] = groups["intron, active"][~k36]
        for g, x in groups.items():
            lo, hi = boot_median(x.DI)
            rows.append(dict(line=line, group=slug[g], n=len(x), median_DI=x.DI.median(), median_ratio=2 ** x.DI.median(),
                             ci_low=2 ** lo, ci_high=2 ** hi, p_wilcoxon=wilcoxon(x.DI).pvalue))
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
            bs = []
            for _ in range(1000):
                kk = RNG.integers(0, len(a), len(a))
                bs.append(1 - 2 ** (np.median(a[kk]) - np.median(m[kk])))
            match.append(dict(line=line, group=slug[g], n=len(x), matched=len(a),
                              deficit_pct=100 * (1 - 2 ** (np.median(a) - np.median(m))),
                              ci_low_pct=100 * np.percentile(bs, 2.5), ci_high_pct=100 * np.percentile(bs, 97.5),
                              p_wilcoxon=wilcoxon(a - m).pvalue))
    d = pd.DataFrame(rows)
    report("; ".join(f"{r.line} {r.group} {r.median_ratio:.2f} ({r.ci_low:.2f}-{r.ci_high:.2f}) n = {r.n}" for r in d.itertuples()
                     if r.group in ("inter", "silent", "active")))
    return {"figS8.tsv": d, "figS8_matched.tsv": pd.DataFrame(match)}


@panel("figS9")
def figS9():
    """Supp. Fig. S9 by distance to the next host exon (step 35, results/exon_dist/sites_exon_dist.tsv: confirmed
    intronic start sites and random positions in active hosts).  (a) % with stable RNA under the original threshold
    (column stable = stabilne_RNA, all random positions); (b) median PRO-cap reads of the start sites; (c) % above the
    random 95th percentile of the DRB test (column drb: K562 itss_w8, HCT116 ttdrb_10)."""
    BINS = [1000, 2000, 3000, 5000, 10000, 20000, np.inf]
    LABS = ["1-2 kb", "2-3 kb", "3-5 kb", "5-10 kb", "10-20 kb", ">20 kb"]
    S = read("results/exon_dist/sites_exon_dist.tsv")
    S["b6"] = pd.cut(num(S.dist_down), BINS, labels=range(6), right=False)
    rows = []
    for line in LINES:
        d = S[S.line == line]
        for gname, g in (("start", "itss"), ("losowa", "rand")):
            for bb in range(6):
                x = d[(d.grp == gname) & (d.b6 == bb)]
                st_ = b(x.stable)
                rows.append(dict(line=line, measure="stable", group=g, bin=bb, bin_label=LABS[bb], n=len(x),
                                 value=100 * st_.mean() if len(x) else np.nan))
                dr = num(x.drb).dropna().astype(bool)
                rows.append(dict(line=line, measure="drb", group=g, bin=bb, bin_label=LABS[bb], n=len(dr),
                                 value=100 * dr.mean() if len(dr) else np.nan))
                if g == "itss":
                    rows.append(dict(line=line, measure="capped_median", group=g, bin=bb, bin_label=LABS[bb], n=len(x),
                                     value=num(x.cap_sum).median()))
    d = pd.DataFrame(rows)
    report("stable 1-2 kb: " + "; ".join(f"{r.line} {r.group} {r.value:.1f}%" for r in d[(d.measure == "stable") & (d.bin == 0)].itertuples()))
    return {"figS9.tsv": d}


# ---------------------------------------------------------------------------------------- cross-checks
@panel("crosscheck_st1")
def crosscheck_st1():
    """Consistency of the result tables with Supplementary Table 1 v2 (counts per cell line); not a figure panel."""
    t = read(SUPP)
    rows = []
    for line, base in LINES.items():
        x = t[t.cell_line == line]
        e = read(f"{base}/itss.evidence.tsv")
        c = e[e.site_type == "candidate"]
        pairs = [("sites", len(x), len(c)),
                 ("confirmed", int(b(x.capped_start_confirmed).sum()), int(b(c.tss_supported).sum())),
                 ("confirmed, active", int((b(x.capped_start_confirmed) & b(x.host_active)).sum()),
                  int((b(c.tss_supported) & b(c.host_active)).sum()))]
        for item, a, bb in pairs:
            rows.append(dict(line=line, item=item, supp_table1=a, evidence_table=bb, ok=a == bb))
    d = pd.DataFrame(rows)
    report("; ".join(f"{r.line} {r.item}: {r.supp_table1} / {r.evidence_table} {'OK' if r.ok else 'DIFF'}" for r in d.itertuples()))
    return {"_crosscheck_supp_table1.tsv": d}


# ------------------------------------------------------------------------------------------------ main
def stale(name):
    """Remove the tables an earlier run wrote for a panel that now failed, so that 51 cannot draw old values."""
    old = [f for f in sorted(os.listdir(ARGS.out)) if f == f"{name}.tsv" or (f.startswith(f"{name}_") and f.endswith(".tsv"))]
    for f in old:
        os.remove(os.path.join(ARGS.out, f))
    return f"; removed the earlier table(s) {', '.join(old)}" if old else ""


def main():
    global ARGS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panels", nargs="*", help="panel names or shell patterns (default: all)")
    ap.add_argument("--no-bigwig", action="store_true", help="never read bigWig files (use the text outputs)")
    ap.add_argument("--skip-gtf", action="store_true", help="skip fig1_b_waffle (parses the GENCODE GTF, a few minutes)")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--list", action="store_true", help="list the panels and their sources, then exit")
    ARGS = ap.parse_args()
    if ARGS.list:
        for name, f, slow in PANELS:
            print(f"{name}{' (slow)' if slow else ''}: {' '.join((f.__doc__ or '').split())}\n")
        return
    os.makedirs(ARGS.out, exist_ok=True)
    print(f"50_figure_data.py  project {os.getcwd()}  python {sys.version.split()[0]}  pandas {pd.__version__}  "
          f"numpy {np.__version__}")
    status = []
    for name, f, slow in PANELS:
        if ARGS.panels and not any(fnmatch.fnmatch(name, p) for p in ARGS.panels):
            continue
        if slow and ARGS.skip_gtf:
            status.append(dict(panel=name, status="skipped", rows=0, files="", detail="--skip-gtf"))
            print(f"[{name}]\n  SKIP {name}: --skip-gtf (a table from an earlier run, if any, is kept)")
            continue
        CUR["panel"] = name
        t0 = time.time()
        print(f"[{name}]", flush=True)
        try:
            out = f()
            files = []
            for fn, df in out.items():
                df.to_csv(os.path.join(ARGS.out, fn), sep="\t", index=False)
                files.append(f"{fn}({len(df)})")
            status.append(dict(panel=name, status="ok", rows=sum(len(x) for x in out.values()), files=" ".join(files),
                               detail=f"{time.time() - t0:.0f} s"))
            print(f"  ok: {' '.join(files)}  [{time.time() - t0:.0f} s]")
        except MissingInput as m:
            status.append(dict(panel=name, status="skipped: missing input", rows=0, files="", detail=str(m)))
            print(f"  SKIP {name}: input not found: {m} (the other panels go on; 51 leaves this panel as a frame)"
                  + stale(name))
        except Exception as e:  # one panel must not stop the others
            tb = traceback.format_exc().strip().splitlines()
            status.append(dict(panel=name, status="error", rows=0, files="", detail=f"{type(e).__name__}: {e}"))
            print(f"  ERROR in {name} {type(e).__name__}: {e} (the other panels go on){stale(name)}\n    "
                  + "\n    ".join(tb[-6:]))
    st = pd.DataFrame(status)
    st.to_csv(os.path.join(ARGS.out, "_status.tsv"), sep="\t", index=False)
    src = []
    for pnl, p in SOURCES:
        try:
            s = os.stat(p)
            src.append(dict(panel=pnl, path=p, bytes=s.st_size, modified=time.strftime("%Y-%m-%d %H:%M", time.localtime(s.st_mtime))))
        except OSError:
            pass
    pd.DataFrame(src).drop_duplicates().to_csv(os.path.join(ARGS.out, "_sources.tsv"), sep="\t", index=False)
    n_ok = int((st.status == "ok").sum()) if len(st) else 0
    print(f"\n{n_ok} of {len(st)} panels written to {ARGS.out}/; "
          + ", ".join(f"{r.panel}: {r.status}" for r in st.itertuples() if r.status != "ok"))


if __name__ == "__main__":
    main()

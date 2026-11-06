#!/usr/bin/env python
"""Step 35: distance from the intronic start site to the next host exon vs. initiation and transcript fate.

Questions:
  EMATS (Fiszbein et al. 2019, 2023): is proximity to an efficiently spliced downstream exon
      associated with stronger initiation (cap_sum)?
  iSTART: is it associated with more productive transcription (stable RNA, response after DRB)?
Controls: random intronic positions in the same distance bins (the 0-3 kb stable RNA window
may overlap the exon), distance to the upstream exon as a comparison, intron length
and host start signal in the regression.

Run from the project directory (interactive node, approx. 5 min, 16 GB):
    python scripts/35_exon_distance.py > results/exon_dist/summary.txt 2>&1
"""
import gzip, os, re, subprocess
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

GTF = "ref/gencode.v47.annotation.gtf.gz"
GENOME = "ref/GRCh38.primary_assembly.genome.fa"
OUT = "results/exon_dist"
LINES = {"K562": "results/v2", "HCT116": "results/hct116"}
DRB = {"HCT116": "results/hct116/ttdrb_10.drb.tsv", "K562": "results/v2/itss_w8.drb.tsv"}
BINS = [1000, 2000, 3000, 5000, 10000, 20000, np.inf]
LABS = ["1-2 kb", "2-3 kb", "3-5 kb", "5-10 kb", "10-20 kb", ">20 kb"]
os.makedirs(OUT, exist_ok=True)
pd.set_option("display.width", 200)


def b(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


def sec(t):
    print(f"\n{'=' * 8} {t} {'=' * 8}", flush=True)


def pct(x):
    return f"{100 * x.mean():5.1f}%" if len(x) else "    -"


# ---------- 1. exons of the host reference transcripts (GTF) ----------
def load_sites():
    out = {}
    for line, base in LINES.items():
        d = pd.read_csv(f"results/repeats.{line}.tsv", sep="\t", low_memory=False)
        u = pd.read_csv(f"{base}/u1pas.tsv", sep="\t", usecols=["site_id", "protected"], low_memory=False)
        d = d.merge(u.drop_duplicates("site_id"), on="site_id", how="left")
        d["line"] = line
        out[line] = d
    return out


sites = load_sites()
ref_tx = set(pd.concat(sites.values())["host_ref_tx"].dropna().astype(str))
sec(f"Exons of {len(ref_tx)} reference transcripts from {GTF}")
rows = []
pat_tx = re.compile(r'transcript_id "([^"]+)"')
with gzip.open(GTF, "rt") as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.split("\t", 8)
        if f[2] != "exon":
            continue
        m = pat_tx.search(f[8])
        if m and m.group(1) in ref_tx:
            rows.append((m.group(1), f[0], int(f[3]), int(f[4]), f[6]))
ex = pd.DataFrame(rows, columns=["tx", "chrom", "start", "end", "strand"]).sort_values(["tx", "start"])
print(f"exons: {len(ex)}")
EX = {tx: g[["start", "end"]].values for tx, g in ex.groupby("tx")}
STR = ex.groupby("tx")["strand"].first().to_dict()
CHR = ex.groupby("tx")["chrom"].first().to_dict()


def downstream_exon(tx, pos, strand):
    """Next exon in the direction of transcription: (start, end, is_last) or None."""
    e = EX.get(tx)
    if e is None:
        return None
    if strand == "+":
        idx = np.where(e[:, 0] > pos)[0]
        if not len(idx):
            return None
        i = idx[0]
        return e[i, 0], e[i, 1], i == len(e) - 1
    idx = np.where(e[:, 1] < pos)[0]
    if not len(idx):
        return None
    i = idx[-1]
    return e[i, 0], e[i, 1], i == 0


# ---------- 2. donor-site sequences and PWM score ----------
def donor_region(chrom, s, e, strand):
    # 9 nt: 3 exonic + 6 intronic (transcript orientation)
    return f"{chrom}:{e - 2}-{e + 6}" if strand == "+" else f"{chrom}:{s - 6}-{s + 2}"


def fetch(regions):
    reg = sorted(set(regions))
    with open(f"{OUT}/regions.txt", "w") as fo:
        fo.write("\n".join(reg) + "\n")
    fa = subprocess.run(["samtools", "faidx", GENOME, "-r", f"{OUT}/regions.txt"],
                        capture_output=True, text=True, check=True).stdout
    seqs, name = {}, None
    for ln in fa.splitlines():
        if ln.startswith(">"):
            name = ln[1:].strip()
            seqs[name] = ""
        elif name:
            seqs[name] += ln.strip().upper()
    return seqs


RC = str.maketrans("ACGTN", "TGCAN")
train = []
for tx, e in EX.items():
    st, ch = STR[tx], CHR[tx]
    inner = e[:-1] if st == "+" else e[1:]  # exons that have a donor site
    for s, en in inner:
        train.append((donor_region(ch, s, en, st), st))
seqs = fetch([r for r, _ in train])


def oriented(region, strand):
    s = seqs.get(region, "")
    return s.translate(RC)[::-1] if strand == "-" else s


tr = [oriented(r, st) for r, st in train]
tr = [s for s in tr if len(s) == 9 and s[3:5] == "GT"]
cnt = np.full((9, 4), 0.5)
for s in tr:
    for i, c in enumerate(s):
        if c in "ACGT":
            cnt[i, "ACGT".index(c)] += 1
PWM = np.log2(cnt / cnt.sum(1, keepdims=True) / 0.25)
print(f"Donor-site PWM from {len(tr)} internal exons; consensus: "
      + "".join("ACGT"[i] for i in PWM.argmax(1)))


def pwm_score(s):
    if len(s) != 9 or "N" in s:
        return np.nan
    return float(sum(PWM[i, "ACGT".index(c)] for i, c in enumerate(s)))


# ---------- 3. splice junctions: assigning files to cell lines, and exon inclusion ----------
def read_sj(path):
    j = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False)
    j = j.iloc[:, :7]
    j.columns = ["chrom", "js", "je", "jstr", "motif", "annot", "uniq"]
    for c in ("js", "je", "uniq"):
        j[c] = pd.to_numeric(j[c], errors="coerce")
    return j.dropna(subset=["js", "je", "uniq"])


sec("Splice junction files (assigned to cell lines by reads in the HBG1/HBG2 locus)")
SJ = {}
for p in sorted(os.path.join("data/rnaseq", f) for f in os.listdir("data/rnaseq") if f.endswith("junctions.tsv")):
    try:
        j = read_sj(p)
        hb = j[(j.chrom == "chr11") & (j.js > 5_240_000) & (j.je < 5_260_000)].uniq.sum()
        print(f"  {p}: {len(j)} junctions, reads in HBG1/2: {hb}; first row: {j.iloc[0].tolist()}")
        SJ[p] = (j, hb)
    except Exception as e:
        print(f"  {p}: could not read ({e})")
if len(SJ) >= 2:
    ranked = sorted(SJ, key=lambda p: SJ[p][1], reverse=True)
    SJ_LINE = {"K562": SJ[ranked[0]][0], "HCT116": SJ[ranked[1]][0]}
    print(f"  K562 <- {ranked[0]}, HCT116 <- {ranked[1]}")
else:
    SJ_LINE = {}
    print("  too few junction files; skipping exon inclusion")


def index_sj(j):
    idx = {}
    for ch, g in j.groupby("chrom"):
        js, je, u = g.js.values.astype(np.int64), g.je.values.astype(np.int64), g.uniq.values.astype(float)
        idx[ch] = (js, je, u, g.groupby("je").uniq.sum().to_dict(), g.groupby("js").uniq.sum().to_dict())
    return idx


def psi(ix, chrom, s, e):
    if chrom not in ix:
        return np.nan, 0.0
    js, je, u, by_end, by_start = ix[chrom]
    left, right = by_end.get(s - 1, 0.0), by_start.get(e + 1, 0.0)
    skip = u[(js < s) & (je > e)].sum()
    inc = (left + right) / 2
    return (inc / (inc + skip) if inc + skip >= 10 else np.nan), inc + skip


SJ_LINE = {k: index_sj(v) for k, v in SJ_LINE.items()}


# ---------- 4. site table ----------
recs = []
for line, d in sites.items():
    st = d.site_type.astype(str)
    keep = ((st == "candidate") & b(d.tss_supported) | (st == "random_control")) & b(d.host_active)
    d = d[keep & (pd.to_numeric(d.min_dist_host_exon, errors="coerce") >= 1000)].copy()
    d["grp"] = np.where(d.site_type == "candidate", "start", "losowa")
    plus = d.strand == "+"
    d["dist_down"] = np.where(plus, d.intron_end - d.pos, d.pos - d.intron_start)
    d["dist_up"] = np.where(plus, d.pos - d.intron_start, d.intron_end - d.pos)
    d["intron_len"] = d.intron_end - d.intron_start
    d["bin"] = pd.cut(d.dist_down, BINS, labels=LABS, right=False)
    d["stable"] = b(d.stabilne_RNA)
    dex = [downstream_exon(t, p, s) for t, p, s in zip(d.host_ref_tx.astype(str), d.pos, d.strand)]
    d["ex_s"] = [x[0] if x else np.nan for x in dex]
    d["ex_e"] = [x[1] if x else np.nan for x in dex]
    d["ex_last"] = [x[2] if x else np.nan for x in dex]
    recs.append(d)
S = pd.concat(recs, ignore_index=True)
need = S[(S.ex_last == False) & S.ex_s.notna()]
regs = [donor_region(c, int(s), int(e), st) for c, s, e, st in zip(need.chrom, need.ex_s, need.ex_e, need.strand)]
seqs.update(fetch(regs))
S["donor_score"] = np.nan
S.loc[need.index, "donor_score"] = [pwm_score(oriented(r, st)) for r, st in zip(regs, need.strand)]
S["psi"], S["sj_reads"] = np.nan, np.nan
for line, j in SJ_LINE.items():
    m = (S.line == line) & S.ex_s.notna()
    vals = [psi(j, c, int(s), int(e)) for c, s, e in zip(S.chrom[m], S.ex_s[m], S.ex_e[m])]
    S.loc[m, "psi"] = [v[0] for v in vals]
    S.loc[m, "sj_reads"] = [v[1] for v in vals]

# DRB response: above the 95th percentile of random positions in active hosts
S["drb"] = np.nan
for line, path in DRB.items():
    if not os.path.exists(path):
        continue
    try:
        g = pd.read_csv(path, sep="\t", low_memory=False)
        g = g[b(g.testable)] if "testable" in g else g
        thr = np.nanpercentile(g.loc[(g.site_type == "random_control") & b(g.host_active), "gap_index"], 95)
        resp = dict(zip(g.site_id, g.gap_index > thr))
        m = S.line == line
        S.loc[m, "drb"] = S.loc[m, "site_id"].map(resp).astype(float)
        print(f"DRB {line}: {path}, threshold {thr:.2f}, sites with a value {int(S.loc[m, 'drb'].notna().sum())}")
    except Exception as e:
        print(f"DRB {line}: skipped ({type(e).__name__}: {e})")
S.to_csv(f"{OUT}/sites_exon_dist.tsv", sep="\t", index=False)


def logit_irls(X, y, names, it=50):
    """Logistic regression (Newton-Raphson); returns coefficients and Wald p-values."""
    from scipy.stats import norm
    X = np.column_stack([np.ones(len(X)), X])
    beta = np.zeros(X.shape[1])
    for _ in range(it):
        p = 1 / (1 + np.exp(-X @ beta))
        W = p * (1 - p)
        H = X.T @ (X * W[:, None])
        step = np.linalg.solve(H, X.T @ (y - p))
        beta += step
        if np.abs(step).max() < 1e-8:
            break
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    return {n: (beta[i + 1], 2 * norm.sf(abs(beta[i + 1] / se[i + 1]))) for i, n in enumerate(names)}


# ---------- 5. results ----------
def table(d, title):
    print(f"\n{title}")
    print(f"  {'distance':<9} {'n start':>7} {'n rnd.':>7} {'cap_sum':>8} {'stab. start':>11} {'stab. rnd.':>10}"
          f" {'DRB start':>9} {'DRB rnd.':>8} {'prot.':>7}")
    for lab in LABS:
        a, r = d[(d.grp == "start") & (d.bin == lab)], d[(d.grp == "losowa") & (d.bin == lab)]
        da, dr = a.drb.dropna().astype(bool), r.drb.dropna().astype(bool)
        print(f"  {lab:<9} {len(a):>7} {len(r):>7} {a.cap_sum.median():>8.0f} {pct(a.stable):>11} {pct(r.stable):>10}"
              f" {pct(da):>9} {pct(dr):>8} {pct(b(a.protected)) if 'protected' in a else '':>7}")


for line in LINES:
    d = S[S.line == line]
    a = d[d.grp == "start"]
    sec(f"{line}: {len(a)} start sites (supported, active host), {int((d.grp == 'losowa').sum())} random positions")
    print(f"distance to the downstream exon, median: start sites {a.dist_down.median():.0f} nt, random {d[d.grp == 'losowa'].dist_down.median():.0f} nt")
    table(d, "A. By distance to the next exon (stab. = stable RNA, DRB = above the random-position threshold, prot. = donor before PAS)")
    # EMATS: initiation vs. distance
    r1 = spearmanr(a.dist_down, a.cap_sum, nan_policy="omit")
    r2 = spearmanr(a.dist_up, a.cap_sum, nan_policy="omit")
    print(f"\nB. Initiation (cap_sum) vs. distance: to the downstream exon rho = {r1.statistic:.3f} (p = {r1.pvalue:.2g}); "
          f"to the upstream exon rho = {r2.statistic:.3f} (p = {r2.pvalue:.2g})")
    a2 = a.assign(t=pd.qcut(a.intron_len.rank(method="first"), 3, labels=["short", "medium", "long"]))
    for t, g in a2.groupby("t", observed=True):
        r = spearmanr(g.dist_down, g.cap_sum, nan_policy="omit")
        print(f"   introns {t:<8} (n = {len(g)}): rho = {r.statistic:.3f}, p = {r.pvalue:.2g}")
    # without start sites closer than 3 kb (the stable RNA window does not touch the exon)
    far = d[d.dist_down >= 3000]
    print(f"\nC. Only sites >= 3 kb from the exon: stable RNA, start sites {pct(far[far.grp == 'start'].stable)}, random {pct(far[far.grp == 'losowa'].stable)}")
    # logistic regression (variables in log10)
    for grp in ("start", "losowa"):
        try:
            x = d[d.grp == grp].copy()
            names = ["dist_down", "dist_up", "host_tss_signal"] + (["cap_sum"] if grp == "start" else [])
            X = np.column_stack([np.log10(pd.to_numeric(x[c], errors="coerce").clip(lower=1)) for c in names])
            ok = np.isfinite(X).all(1)
            X = (X[ok] - X[ok].mean(0)) / X[ok].std(0)
            res = logit_irls(X, x.stable.values[ok].astype(float), names)
            print(f"\nD. Logistic regression of stable RNA ({grp}, n = {int(ok.sum())}): coefficient per 1 SD of log10, and p")
            for k, (c, p) in res.items():
                print(f"   {k:<16} {c:+.3f}  p = {p:.2g}")
        except Exception as e:
            print(f"\nD. regression ({grp}) skipped ({type(e).__name__}: {e})")
    # donor strength and exon inclusion (sites up to 10 kb from the exon, exon not the last one)
    near = d[(d.dist_down < 10000) & (d.ex_last == False)].copy()
    if near.donor_score.notna().sum() > 30:
        near["donor"] = pd.qcut(near.donor_score.rank(method="first"), 3, labels=["weak", "medium", "strong"])
        print("\nE. Sites < 10 kb from the exon: donor-site strength of the downstream exon (PWM thirds)")
        for t, g in near.groupby("donor", observed=True):
            ga, gr = g[g.grp == "start"], g[g.grp == "losowa"]
            print(f"   donor {t:<6} n {len(ga):>4}/{len(gr):>5}  cap_sum {ga.cap_sum.median():>5.0f}  stable {pct(ga.stable)} / {pct(gr.stable)}"
                  f"  DRB {pct(ga.drb.dropna().astype(bool))} / {pct(gr.drb.dropna().astype(bool))}")
    if near.psi.notna().sum() > 30:
        near["incl"] = np.select([near.psi >= 0.9, near.psi < 0.9], ["constitutive (PSI >= 0.9)", "alternative (PSI < 0.9)"], "too few reads")
        print("\nF. Sites < 10 kb from the exon: inclusion of the downstream exon (RNA-seq junctions)")
        for t, g in near.groupby("incl"):
            ga, gr = g[g.grp == "start"], g[g.grp == "losowa"]
            print(f"   {t:<27} n {len(ga):>4}/{len(gr):>5}  cap_sum {ga.cap_sum.median():>5.0f}  stable {pct(ga.stable)} / {pct(gr.stable)}"
                  f"  DRB {pct(ga.drb.dropna().astype(bool))} / {pct(gr.drb.dropna().astype(bool))}")
        alt = near[(near.grp == "start") & (near.psi < 0.9)].cap_sum
        con = near[(near.grp == "start") & (near.psi >= 0.9)].cap_sum
        if len(alt) > 5 and len(con) > 5:
            print(f"   start-site cap_sum: alternative {alt.median():.0f} vs. constitutive {con.median():.0f}, "
                  f"Mann-Whitney p = {mannwhitneyu(alt, con).pvalue:.2g}")
    # nearest vs. farthest: Fisher test for stable RNA (start sites minus random)
    a1, a3 = a[a.dist_down < 3000], a[a.dist_down >= 10000]
    t = [[int(a1.stable.sum()), int((~a1.stable).sum())], [int(a3.stable.sum()), int((~a3.stable).sum())]]
    if min(len(a1), len(a3)) > 0:
        print(f"\nG. Start sites < 3 kb vs. >= 10 kb from the exon: stable RNA {pct(a1.stable)} vs. {pct(a3.stable)}, "
              f"OR {fisher_exact(t)[0]:.2f}, p = {fisher_exact(t)[1]:.2g}")

print(f"\nSite table: {OUT}/sites_exon_dist.tsv\n==== end ====")

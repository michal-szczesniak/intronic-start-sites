#!/usr/bin/env python3
"""
iSTART: final tables for the manuscript.

1. Supplementary Table 1: all intronic start sites of both cell lines with every measure (English column names).
2. Figure 5: RNA per initiation event with 95% bootstrap confidence intervals, for intronic start sites in
   active and silent host genes, start sites between genes and internal promoters, within thirds of strength.
   Pooled estimate = sum(excess signal after the start over the level before it) / sum(PRO-cap reads).
3. Start sites shared by both cell lines versus specific to one, compared within thirds of strength
   (Fisher test per third, Cochran-Mantel-Haenszel test across thirds).
Input file names follow the repository layout (results/...); run from the project directory.
"""
import argparse

import numpy as np
import pandas as pd
from scipy.stats import chi2, fisher_exact

LINES = {
    "K562": dict(rep="results/repeats.K562.tsv", u1="results/v2/u1pas.tsv", exo="results/v2/exosome/exo.MTREX.tsv",
                 ctx="results/context_K562.tsv", cls=None),
    "HCT116": dict(rep="results/repeats.HCT116.tsv", u1="results/hct116/u1pas.tsv", exo=None,
                   ctx="results/context_HCT116.tsv", cls="results/hct116/classes_10.tsv"),
}
TYP = {"promotor-podobny": "promoter-like", "enhancer-podobny": "enhancer-like", "bez znaczników": "unmarked"}
REP = {"brak": "none", "proste/niska złożoność": "simple/low complexity", "inne": "other", "SINE/inne": "SINE/other",
       "LINE/inne": "LINE/other"}
ORI = {"zgodna": "same", "przeciwna": "opposite"}
CLS = {"produktywny": "productive", "krótki": "short", "bez odpowiedzi": "no response", "nierozstrzygnięty": "unresolved"}


def shared_flag(a, b, w=200):
    idx = {}
    for r in b.itertuples():
        idx.setdefault((r.chrom, r.strand), []).append(int(r.pos))
    idx = {k: np.sort(np.array(v)) for k, v in idx.items()}
    out = []
    for r in a.itertuples():
        v = idx.get((r.chrom, r.strand))
        ok = False
        if v is not None:
            j = np.searchsorted(v, int(r.pos))
            ok = any(0 <= jj < len(v) and abs(v[jj] - int(r.pos)) <= w for jj in (j - 1, j))
        out.append(ok)
    return np.array(out)


def cmh(tables):
    """Cochran-Mantel-Haenszel test for a list of 2x2 tables [[a, b], [c, d]]."""
    num, var, or_n, or_d = 0.0, 0.0, 0.0, 0.0
    for (a, b), (c, d) in tables:
        n = a + b + c + d
        if n < 2:
            continue
        e = (a + b) * (a + c) / n
        num += a - e
        var += (a + b) * (c + d) * (a + c) * (b + d) / (n * n * (n - 1))
        or_n += a * d / n
        or_d += b * c / n
    stat = (abs(num) - 0.5) ** 2 / var if var > 0 else np.nan
    return (or_n / or_d if or_d else np.nan), (chi2.sf(stat, 1) if var > 0 else np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--out-prefix", default="results/final")
    a = ap.parse_args()
    rng = np.random.default_rng(1)

    # ---------- 1. Supplementary Table 1 ----------
    tabs = {}
    for name, f in LINES.items():
        d = pd.read_csv(f["rep"], sep="\t", low_memory=False)
        for c in ("host_active", "tss_supported"):
            d[c] = d[c].astype(str) == "True"
        d = d[d.site_type == "candidate"].copy()
        u = pd.read_csv(f["u1"], sep="\t", low_memory=False)[["site_id", "pas_sens", "pas_anty", "prot_sens"]]
        d = d.merge(u, on="site_id", how="left")
        if f["exo"]:
            x = pd.read_csv(f["exo"], sep="\t", low_memory=False)[["site_id", "delta"]]
            d = d.merge(x, on="site_id", how="left")
        ctx = pd.read_csv(f["ctx"], sep="\t", low_memory=False)
        keep = [c for c in ("site_id", "RNA_nadw", "DRB_dal") if c in ctx]
        d = d.merge(ctx[keep], on="site_id", how="left")
        if f["cls"] and "klasa" not in d:
            k = pd.read_csv(f["cls"], sep="\t", low_memory=False)[["site_id", "klasa"]]
            d = d.merge(k, on="site_id", how="left")
        tabs[name] = d
    tabs["K562"]["shared"] = shared_flag(tabs["K562"], tabs["HCT116"][tabs["HCT116"].tss_supported])
    tabs["HCT116"]["shared"] = shared_flag(tabs["HCT116"], tabs["K562"][tabs["K562"].tss_supported])
    out = []
    for name, d in tabs.items():
        t = pd.DataFrame({
            "cell_line": name, "site_id": d.site_id, "chrom": d.chrom, "strand": d.strand, "position": d.pos.astype(int),
            "host_gene": d.get("host_gene_name", d.get("host_name", pd.Series("", index=d.index))),
            "host_active": d.host_active, "capped_start_confirmed": d.tss_supported, "procap_reads": d.cap_sum,
            "chromatin_type": d.typ.map(TYP).fillna(d.typ) if "typ" in d else np.nan})
        for c in d.columns:
            if c.startswith("H3K"):
                t[c] = d[c]
        t["phylop_mean"] = d.phylop
        t["stable_rna_rise_log2"] = d.rna_step
        t["stable_rna_per_initiation"] = d.RNA_nadw if "RNA_nadw" in d else np.nan
        t["pas_per_kb_transcribed_strand"] = d.pas_sens
        t["pas_per_kb_opposite_strand"] = d.pas_anty
        t["splice_donor_before_first_pas"] = d.prot_sens
        t["repeat_class"] = d.grupa_rep.map(REP).fillna(d.grupa_rep)
        t["repeat_element"] = d.element
        t["repeat_orientation"] = d.orientacja.map(ORI).fillna(d.orientacja)
        t["drb_class_10min"] = d.klasa.map(CLS) if "klasa" in d else np.nan
        t["drb_far_per_initiation"] = d.DRB_dal if "DRB_dal" in d else np.nan
        t["mtrex_knockout_change"] = d.delta if "delta" in d else np.nan
        t["shared_with_other_cell_line"] = d.shared
        out.append(t)
    S1 = pd.concat(out, ignore_index=True)
    S1.to_csv(f"{a.out_prefix}.supp_table1_intronic_starts.tsv", sep="\t", index=False)
    print(f"Supplementary Table 1: {len(S1)} intronic start sites "
          f"({', '.join(f'{n} {int((S1.cell_line == n).sum())}' for n in tabs)}), {S1.shape[1]} columns -> "
          f"{a.out_prefix}.supp_table1_intronic_starts.tsv")

    # ---------- 2. Figure 5: RNA per initiation with bootstrap CIs ----------
    groups = {"intron_akt": "intronic, active host", "intron_cichy": "intronic, silent host",
              "miedzygen": "between genes", "wewn_prom": "internal promoters"}
    rows = []
    print("\nFigure 5: signal per initiation event (pooled: sum of excess / sum of PRO-cap reads, x1000), 95% bootstrap CI")
    for name, f in LINES.items():
        S = pd.read_csv(f["ctx"], sep="\t", low_memory=False)
        st3 = S[S.grupa.isin(["intron_akt", "intron_cichy", "miedzygen"])].copy()
        st3["third"] = pd.qcut(st3.cap_sum.rank(method="first"), 3, labels=["weak", "medium", "strong"])
        S = S.merge(st3[["site_id", "third"]], on="site_id", how="left")
        S["third"] = S.third.astype(object).fillna("all")
        for col, lab in (("RNA_nadw", "stable RNA"), ("DRB_dal", "nascent RNA 2-8 kb (DRB 10 min)")):
            if col not in S or S[col].isna().all():
                continue
            print(f"\n  {name} – {lab}")
            print(f"  {'third':<8}" + "".join(f"{groups[g]:>36}" for g in groups))
            for th in ("weak", "medium", "strong", "all"):
                cells = []
                for g in groups:
                    x = S[(S.grupa == g) & ((S.third == th) if th != "all" else True)].dropna(subset=[col])
                    if g == "wewn_prom" and th != "all":
                        cells.append(f"{'':>36}")
                        continue
                    if len(x) < 5:
                        cells.append(f"{'–':>36}")
                        continue
                    cap = np.maximum(x.cap_sum.values, 1)
                    exc = x[col].values * cap
                    est = exc.sum() / cap.sum() * 1000
                    bs = []
                    for _ in range(a.boot):
                        i = rng.integers(0, len(x), len(x))
                        bs.append(exc[i].sum() / cap[i].sum() * 1000)
                    lo, hi = np.percentile(bs, [2.5, 97.5])
                    rows.append(dict(cell_line=name, measure=lab, third=th, group=groups[g], n=len(x), estimate=est, ci_low=lo, ci_high=hi))
                    cells.append(f"{est:>12.3f} [{lo:.3f}, {hi:.3f}] n={len(x):<5}")
                print(f"  {th:<8}" + "".join(cells))
    pd.DataFrame(rows).to_csv(f"{a.out_prefix}.fig5_per_initiation.tsv", sep="\t", index=False)

    # ---------- 3. shared vs specific within strength thirds ----------
    print("\nShared vs cell-type-specific intronic start sites (capped, active host), within thirds of strength")
    for name, d in tabs.items():
        x = d[d.tss_supported & d.host_active].copy()
        x["third"] = pd.qcut(x.cap_sum.rank(method="first"), 3, labels=["weak", "medium", "strong"])
        x["conserved"] = x.phylop > 1
        x["transposon"] = ~x.grupa_rep.isin(["brak", "proste/niska złożoność"])
        x["stable"] = x.stabilne_RNA.astype(str) == "True" if "stabilne_RNA" in x else np.nan
        print(f"\n  {name}: shared {int(x.shared.sum())}, specific {int((~x.shared).sum())}")
        feats = [("conserved", "phyloP > 1"), ("transposon", "in a transposon"), ("stable", "stable RNA")]
        if "klasa" in x:
            x["short"] = x.klasa == "krótki"
            feats.append(("short", "short DRB burst"))
        print(f"  {'feature':<18}" + "".join(f"{t + ' (shared / specific)':>34}" for t in ("weak", "medium", "strong"))
              + f"{'CMH odds ratio, p':>22}")
        for col, lab in feats:
            cells, tabs2 = [], []
            for th in ("weak", "medium", "strong"):
                y = x[x.third == th]
                s1, s0 = y[y.shared][col].astype(bool), y[~y.shared][col].astype(bool)
                tab = [[int(s1.sum()), len(s1) - int(s1.sum())], [int(s0.sum()), len(s0) - int(s0.sum())]]
                tabs2.append(tab)
                p = fisher_exact(tab)[1] if len(s1) and len(s0) else np.nan
                cells.append(f"{s1.mean() * 100:>13.0f}% / {s0.mean() * 100:3.0f}% (p={p:.1e})")
            orr, p = cmh(tabs2)
            print(f"  {lab:<18}" + "".join(f"{c:>34}" for c in cells) + f"{orr:>12.2f}, {p:.1e}")


if __name__ == "__main__":
    main()

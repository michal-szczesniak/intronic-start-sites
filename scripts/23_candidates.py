#!/usr/bin/env python3
"""
iSTART: list of the strongest candidates among capped intronic start sites.

Each start site gets 1 point for each of: (1) a DRB effect above the 95th percentile of random controls,
(2) a stable-RNA rise above the 95th percentile of controls, (3) conservation (phyloP > 1),
(4) promoter-like chromatin. Sorted by points, then by start-site strength (PRO-cap).
Output: the full table (--out) and the top candidates on screen, with a link to the UCSC browser (±5 kb).
"""
import argparse

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--describe", default="results/modA/itss.describe.tsv")
    ap.add_argument("--phylop", default="results/modA/sites_phylop.tab")
    ap.add_argument("--expr", default="results/modA/expr.itss_expression.tsv")
    ap.add_argument("--out", default="results/modA/candidates_ranked.tsv")
    ap.add_argument("--top", type=int, default=25)
    a = ap.parse_args()

    d = pd.read_csv(a.describe, sep="\t", low_memory=False).reset_index()
    t = pd.read_csv(a.phylop, sep="\t", header=None, names=["n", "size", "covered", "sum", "mean0", "mean"])
    t["index"] = t.n.str[1:].astype(int)
    d = d.merge(t[["index", "mean"]].rename(columns={"mean": "phylop"}), on="index", how="left")
    e = pd.read_csv(a.expr, sep="\t", low_memory=False)
    d = d.merge(e[["site_id", "rna_step", "junc_down"]], on="site_id", how="left")
    for col in ("host_active", "tss_supported", "testable"):
        d[col] = d[col].astype(str) == "True"

    ctrl = d[(d.site_type == "random_control") & d.host_active]
    thr_drb = np.percentile(ctrl[ctrl.testable].gap_index, 95)
    thr_rna = np.percentile(ctrl.rna_step.dropna(), 95)
    c = d[(d.site_type == "candidate") & d.tss_supported].copy()
    if "at_other_pc_tss" in c:
        nested = c.at_other_pc_tss.astype(str) == "True"
        print(f"Skipping {int(nested.sum())} start sites lying at the start of another protein-coding gene (±500 nt)")
        c = c[~nested]
    c["drb_ok"] = c.testable & (c.gap_index > thr_drb)
    c["rna_ok"] = c.rna_step > thr_rna
    c["cons_ok"] = c.phylop > 1
    c["prom_ok"] = c.typ == "promotor-podobny"
    c["score"] = c[["drb_ok", "rna_ok", "cons_ok", "prom_ok"]].sum(axis=1)
    c = c.sort_values(["score", "cap_sum"], ascending=[False, False])
    c["ucsc"] = [f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={r.chrom}:{int(r.pos) - 5000}-{int(r.pos) + 5000}"
                 for r in c.itertuples()]
    c.to_csv(a.out, sep="\t", index=False)

    print(f"Thresholds: DRB gap_index > {thr_drb:.2f}, RNA rise > {thr_rna:.2f}, phyloP > 1")
    print(f"Capped start sites: {len(c)}; number of criteria met (0-4):")
    for s, n in c.score.value_counts().sort_index(ascending=False).items():
        print(f"  {s} criteria: {n}")
    print(f"\nTop candidates ({a.top}):")
    print(f"  {'site':<28}{'host':<14}{'act.':<6}{'type':<18}{'PRO-cap':>8}{'DRB':>6}{'RNA':>6}{'phyloP':>8}{'cl.':>4}  annotation")
    for r in c.head(a.top).itertuples():
        print(f"  {r.chrom}:{int(r.pos)}:{r.strand:<8}{str(r.host_gene_name)[:13]:<14}{'yes' if r.host_active else 'no':<6}"
              f"{str(r.typ)[:17]:<18}{r.cap_sum:>8.0f}{r.gap_index:>6.2f}{r.rna_step:>6.2f}{r.phylop:>8.2f}"
              f"{str(r.evidence_class):>4}  {str(r.cand_gene_name)[:20]}")
    print(f"\nFull table with UCSC links: {a.out}")


if __name__ == "__main__":
    main()

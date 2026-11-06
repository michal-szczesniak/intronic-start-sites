#!/usr/bin/env python3
"""
iSTART: do intronic RNAs exist in the cell as stable RNAs (total RNA-seq, K562)?

A. Annotated RNAs (candidates from step 10, gene level):
   - multi-exon: reads on splice junctions that exactly match the candidate's introns
     (the host pre-mRNA cannot produce such junctions) -> "confirmed by splicing" if >= --min-junc;
   - all: coverage enrichment (same strand) over the candidate's exons relative to 2-kb flanks.
   Results split by: capped start yes/no, active host yes/no.
B. Intronic start sites (iTSS) and controls: rise in stable-RNA coverage after the start site,
   log2[(coverage 0..+3 kb) / (coverage -3.2..-0.2 kb)], and splice junctions in the host intron
   starting up to 20 kb after the start site; comparison with random controls.

Input: --gtf, --cand-tx (candidates_tx.tsv), --cand-genes (annot.genes.tsv from step 11),
       --sites (itss.describe.tsv from step 21), --junctions (chrom, start, end, count; 0-based),
       --plus/--minus (RNA-seq bigWigs; orientation is checked automatically).
"""
import argparse
import bisect
import importlib.util
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("c10", os.path.join(HERE, "10_candidates.py"))
c10 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c10)


class Cov:
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--cand-tx", required=True)
    ap.add_argument("--cand-genes", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--plus", required=True)
    ap.add_argument("--minus", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--min-junc", type=int, default=2)
    ap.add_argument("--flank", type=int, default=2000)
    a = ap.parse_args()
    log = []

    # --- splice junctions ---
    J = {}
    jidx = defaultdict(list)
    with open(a.junctions) as fh:
        for line in fh:
            c, s, e, n = line.split()[:4]
            s, e, n = int(s), int(e), int(n)
            J[(c, s, e)] = n
            jidx[c].append((s, e, n))
    for c in jidx:
        jidx[c].sort()
    jstart = {c: [x[0] for x in v] for c, v in jidx.items()}
    log.append(f"Splice junctions: {len(J)}")

    cov = Cov(a.plus, a.minus)

    # --- orientation: first 5 kb downstream of active host start sites ---
    sites = pd.read_csv(a.sites, sep="\t", low_memory=False)
    for col in ("host_active", "tss_supported", "testable"):
        if col in sites:
            sites[col] = sites[col].astype(str) == "True"
    hp = sites[sites.site_type == "host_promoter"].head(2000)
    same = sum(cov.mean(r.chrom, int(r.pos), int(r.pos) + 5000, r.strand) if r.strand == "+" else
               cov.mean(r.chrom, int(r.pos) - 5000, int(r.pos) + 1, r.strand) for r in hp.itertuples())
    cov.flip = True
    opp = sum(cov.mean(r.chrom, int(r.pos), int(r.pos) + 5000, r.strand) if r.strand == "+" else
              cov.mean(r.chrom, int(r.pos) - 5000, int(r.pos) + 1, r.strand) for r in hp.itertuples())
    cov.flip = opp > same
    log.append(f"RNA-seq orientation (5 kb after {len(hp)} host start sites): files 'plus/minus' = RNA strand: {same:.1f}, "
               f"swapped: {opp:.1f} -> {'swapping' if cov.flip else 'no swap'}")

    # --- A. annotated RNAs ---
    ct = pd.read_csv(a.cand_tx, sep="\t", low_memory=False)
    tx, exons = c10.parse_gtf(a.gtf)
    rows = []
    for r in ct.itertuples():
        ex = exons.get(r.cand_tx_id)
        if not ex:
            continue
        introns = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
        jmax = max([J.get((r.chrom, s, e), 0) for s, e in introns], default=0)
        exon_cov = np.mean([cov.mean(r.chrom, s, e, r.strand) for s, e in ex])
        fl = (cov.mean(r.chrom, r.start - a.flank, r.start, r.strand) + cov.mean(r.chrom, r.end, r.end + a.flank, r.strand)) / 2
        rows.append(dict(cand_tx_id=r.cand_tx_id, cand_gene_id=r.cand_gene_id, n_exons=len(ex),
                         junc_max=jmax, exon_cov=exon_cov, flank_cov=fl,
                         enrich=np.log2((exon_cov + 1e-3) / (fl + 1e-3))))
    t = pd.DataFrame(rows)
    gt = t.groupby("cand_gene_id").agg(multi=("n_exons", lambda x: (x > 1).any()), junc_max=("junc_max", "max"),
                                       enrich=("enrich", "max"), exon_cov=("exon_cov", "max")).reset_index()
    gn = pd.read_csv(a.cand_genes, sep="\t", low_memory=False)
    for col in ("host_active", "tss_supported"):
        gn[col] = gn[col].astype(str) == "True"
    gt = gt.merge(gn[["cand_gene_id", "tss_supported", "host_active", "evidence_class"]], on="cand_gene_id", how="left")
    gt["junction_ok"] = gt.multi & (gt.junc_max >= a.min_junc)
    gt["visible"] = gt.enrich >= 1
    gt.to_csv(f"{a.out_prefix}.annot_expression.tsv", sep="\t", index=False)

    log.append("")
    log.append(f"A. Annotated RNAs (genes). Splicing: >= {a.min_junc} reads on a junction matching a candidate intron; "
               f"coverage: exons >= 2x flanks.")
    log.append(f"  {'group':<46}{'genes':>6}{'multi-ex.':>11}{'splicing':>10}{'% multi-ex.':>13}{'coverage ≥2x':>14}")

    def arow(name, g):
        m = g[g.multi]
        return (f"  {name:<46}{len(g):>6}{len(m):>11}{int(m.junction_ok.sum()):>10}"
                f"{(m.junction_ok.mean() * 100 if len(m) else 0):>12.0f}%{g.visible.mean() * 100:>13.0f}%")

    log.append(arow("all", gt))
    log.append(arow("active host", gt[gt.host_active == True]))
    log.append(arow("active host, capped start", gt[(gt.host_active == True) & (gt.tss_supported == True)]))
    log.append(arow("active host, uncapped", gt[(gt.host_active == True) & (gt.tss_supported != True)]))
    log.append(arow("inactive host", gt[gt.host_active != True]))

    # --- B. intronic start sites and controls ---
    def step(r):
        p = int(r.pos)
        if r.strand == "+":
            dn, up = cov.mean(r.chrom, p, p + 3000, "+"), cov.mean(r.chrom, p - 3200, p - 200, "+")
        else:
            dn, up = cov.mean(r.chrom, p - 2999, p + 1, "-"), cov.mean(r.chrom, p + 201, p + 3201, "-")
        return np.log2((dn + 1e-3) / (up + 1e-3))

    def junc_down(r):
        p, c = int(r.pos), r.chrom
        if c not in jidx or pd.isna(r.intron_start):
            return 0
        ia, ib = int(r.intron_start), int(r.intron_end)
        st, v = jstart[c], jidx[c]
        n = 0
        if r.strand == "+":
            for k in range(bisect.bisect_left(st, p), bisect.bisect_left(st, min(p + 20000, ib))):
                s, e, cnt = v[k]
                if e <= ib:
                    n += cnt
        else:
            for k in range(bisect.bisect_left(st, ia), bisect.bisect_left(st, p)):
                s, e, cnt = v[k]
                if p - 20000 <= e <= p:
                    n += cnt
        return n

    sel = sites.site_type.isin(["candidate", "random_control", "internal_promoter"])
    s2 = sites[sel].copy()
    s2["rna_step"] = [step(r) for r in s2.itertuples()]
    s2["junc_down"] = [junc_down(r) for r in s2.itertuples()]
    s2.to_csv(f"{a.out_prefix}.itss_expression.tsv", sep="\t", index=False)
    ctrl = s2[(s2.site_type == "random_control") & s2.host_active]
    thr = np.percentile(ctrl.rna_step, 95)
    jthr = max(1, np.percentile(ctrl.junc_down, 95))
    cand = s2.site_type == "candidate"
    groups = [("random controls", (s2.site_type == "random_control") & s2.host_active),
              ("capped internal promoters", (s2.site_type == "internal_promoter") & s2.host_active & s2.tss_supported),
              ("capped iTSS, active host", cand & s2.host_active & s2.tss_supported),
              ("  of which promoter-like", cand & s2.host_active & s2.tss_supported & (s2.typ == "promotor-podobny")),
              ("  of which enhancer-like", cand & s2.host_active & s2.tss_supported & (s2.typ == "enhancer-podobny")),
              ("capped iTSS, inactive host", cand & ~s2.host_active & s2.tss_supported),
              ("uncapped iTSS", cand & ~s2.tss_supported)]
    log.append("")
    log.append(f"B. Stable RNA after the start site (coverage rise 0..+3 kb vs -3.2..-0.2 kb; "
               f"threshold = 95th percentile of controls = {thr:.2f}) and splice junctions in the intron "
               f"up to 20 kb after the start site (threshold = {jthr:.0f} reads):")
    log.append(f"  {'group':<42}{'n':>6}{'median rise':>15}{'rise>thr':>11}{'splicing>thr':>15}")
    for name, m in groups:
        g = s2[m]
        if not len(g):
            continue
        log.append(f"  {name:<42}{len(g):>6}{g.rna_step.median():>15.2f}{(g.rna_step > thr).mean() * 100:>10.0f}%"
                   f"{(g.junc_down > jthr).mean() * 100:>14.0f}%")
    txt = "\n".join(log)
    open(f"{a.out_prefix}.expression.txt", "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()

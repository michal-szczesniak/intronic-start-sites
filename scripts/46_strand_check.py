#!/usr/bin/env python3
"""
iSTART step 46: strand decomposition of the detection of nested genes (Supplementary Fig. 8 and Supplementary Table 1), checked against
the way expression is measured and against clustering by host gene.

Groups: nested protein-coding and lncRNA genes on the host's strand (same) and on the opposite strand (anti), each
compared with stand-alone genes of the same type in the same expression bin (expected = sum of stand-alone detection
rates; O/E over bins 1-5). Strand effect = O/E(anti) / O/E(same); excess = strand effect lncRNA / strand effect pc.

Expression estimators:
  cov        step 30 as published: mean exon coverage minus mean coverage 2-5 kb on both sides (same strand)
  cov_multi  the same, restricted to genes with at least one intron of their own (the gene set of 'junc')
  cov_clear  the same, restricted to loci whose exon signal is >= 3x their surroundings (any bin >= 1)
  junc       host-robust: maximum number of RNA-seq reads on a splice junction that exactly matches one of the gene's
             own introns (introns shared with any other GENCODE gene are ignored); host pre-mRNA cannot produce
             these junctions; expressed = >= --min-junc reads; bins = fifths of expressed stand-alone genes
Confidence intervals (95%, percentile): binomial bootstrap within bins (genes independent) and Poisson cluster
bootstrap over host genes (all nested genes of one host resampled together; each stand-alone gene its own cluster).
Also: Mantel-Haenszel risk ratio host-strand lncRNA vs host-strand protein-coding (strata = bins 1-5) and the
expected number of host-strand lncRNAs with a de novo intronic start site within 500 nt of an annotated start
(column itss_<line> of step 30) if each had its own start detected like stand-alone lncRNAs of the same expression.

Input : --genes results/known_genes_wide.tsv  --gtf ref/gencode.v47.annotation.gtf.gz
        --line "NAME|junctions.tsv|rna_plus.bw,rna_minus.bw" (repeatable; bigWigs optional, for an orientation check)
Output: <out-prefix>.counts.tsv, <out-prefix>.summary.tsv; report on stdout
"""
import argparse
import gzip
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

ATTR = re.compile(r'(\w+) "([^"]*)"')
GROUPS = ["pc_stand", "lnc_stand", "pc_same", "pc_anti", "lnc_same", "lnc_anti"]
REF = {"pc_same": "pc_stand", "pc_anti": "pc_stand", "lnc_same": "lnc_stand", "lnc_anti": "lnc_stand"}
GI = {g: i for i, g in enumerate(GROUPS)}
B15 = [1, 2, 3, 4, 5]


def tf(s):
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "1.0"])


def read_gtf_introns(path, keep):
    """introns (chrom, s, e; 0-based, as in the junction tables) of every transcript; owners over all genes"""
    tx_gene, tx_ex = {}, defaultdict(list)
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] != "exon":
                continue
            a = dict(ATTR.findall(f[8]))
            t = a["transcript_id"]
            tx_gene[t] = (a["gene_id"].split(".")[0], f[0])
            tx_ex[t].append((int(f[3]) - 1, int(f[4])))
    owners = defaultdict(set)
    gene_introns = defaultdict(set)
    for t, ex in tx_ex.items():
        gid, chrom = tx_gene[t]
        ex.sort()
        for (s1, e1), (s2, e2) in zip(ex[:-1], ex[1:]):
            if s2 > e1:
                key = (chrom, e1, s2)
                owners[key].add(gid)
                if gid in keep:
                    gene_introns[gid].add(key)
    own = {g: sorted(k for k in v if len(owners[k]) == 1) for g, v in gene_introns.items()}
    return own


def read_junctions(path, needed):
    out = {}
    with open(path) as fh:
        for line in fh:
            c, s, e, n = line.split()[:4]
            k = (c, int(s), int(e))
            if k in needed:
                out[k] = out.get(k, 0) + int(n)
    return out


def tables(cell, det, w=None):
    N = np.bincount(cell, weights=w, minlength=36).reshape(6, 6)
    K = np.bincount(cell, weights=det if w is None else w * det, minlength=36).reshape(6, 6)
    return N, K


def stats(N, K):
    out = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        for g, r in REF.items():
            ref = K[GI[r]] / N[GI[r]]
            o = K[GI[g], 1:6].sum()
            e = np.nansum(N[GI[g], 1:6] * ref[1:6])
            out["OE_" + g] = o / e if e > 0 else np.nan
        out["strand_lnc"] = out["OE_lnc_anti"] / out["OE_lnc_same"]
        out["strand_pc"] = out["OE_pc_anti"] / out["OE_pc_same"]
        out["excess"] = out["strand_lnc"] / out["strand_pc"]
    return out


def mh_rr(N, K, a, c):
    """Mantel-Haenszel risk ratio of group a vs group c over bins 1-5, Greenland-Robins 95% CI"""
    num = den = var = 0.0
    for b in B15:
        n1, k1, n0, k0 = N[GI[a], b], K[GI[a], b], N[GI[c], b], K[GI[c], b]
        T = n1 + n0
        if n1 == 0 or n0 == 0:
            continue
        num += k1 * n0 / T
        den += k0 * n1 / T
        var += (n1 * n0 * (k1 + k0) - k1 * k0 * T) / T ** 2
    if num == 0 or den == 0:
        return np.nan, np.nan, np.nan
    rr = num / den
    se = np.sqrt(var / (num * den))
    return rr, rr * np.exp(-1.96 * se), rr * np.exp(1.96 * se)


def analyse(d, line, label, boot, rng, out_counts, out_sum):
    d = d[d.bin.notna()].copy()
    d["bin"] = d.bin.astype(int)
    cell = (d.grp.map(GI).values * 6 + d.bin.values).astype(int)
    det = d.det.values.astype(float)
    N, K = tables(cell, det)
    for g in GROUPS:
        for b in range(6):
            out_counts.append(dict(line=line, estimator=label, group=g, bin=b, n=int(N[GI[g], b]), k=int(K[GI[g], b])))
    pt = stats(N, K)
    rr = mh_rr(N, K, "lnc_same", "pc_same")
    # binomial bootstrap within cells
    p = np.divide(K, N, out=np.zeros_like(K), where=N > 0)
    bb = defaultdict(list)
    for _ in range(boot):
        Kb = rng.binomial(N.astype(int), p).astype(float)
        for k, v in stats(N, Kb).items():
            bb[k].append(v)
    # Poisson cluster bootstrap over hosts
    cl = d.cluster.astype("category").cat.codes.values
    ncl = cl.max() + 1
    cb = defaultdict(list)
    for _ in range(boot):
        w = rng.poisson(1.0, ncl)[cl].astype(float)
        Nb, Kb = tables(cell, det, w)
        for k, v in stats(Nb, Kb).items():
            cb[k].append(v)
    nest = d[d.kind == "nested"]
    hc = nest.groupby("cluster").agg(n=("gene_id", "size"), types=("type", "nunique"))
    print(f"\n--- {line} | {label} ---")
    print(f"  hosts with nested genes: {len(hc)}; with >= 2 nested genes: {int((hc.n >= 2).sum())} "
          f"(holding {int(hc.n[hc.n >= 2].sum())} of {len(nest)} nested genes); with both lncRNA and protein-coding: "
          f"{int((hc.types == 2).sum())}")
    print("  " + "  ".join(f"{g}: n(1-5)={int(N[GI[g], 1:6].sum())} k={int(K[GI[g], 1:6].sum())} n(1-3)={int(N[GI[g], 1:4].sum())}"
                           for g in ("pc_anti", "pc_same", "lnc_anti", "lnc_same")))
    for k in ["OE_pc_anti", "OE_pc_same", "OE_lnc_anti", "OE_lnc_same", "strand_lnc", "strand_pc", "excess"]:
        lb, hb = np.nanpercentile(bb[k], [2.5, 97.5])
        lc, hc_ = np.nanpercentile(cb[k], [2.5, 97.5])
        print(f"  {k:<12} {pt[k]:.3f}   binomial {lb:.2f}-{hb:.2f}   host-cluster {lc:.2f}-{hc_:.2f}")
        out_sum.append(dict(line=line, estimator=label, stat=k, value=pt[k], binom_lo=lb, binom_hi=hb, clust_lo=lc, clust_hi=hc_))
    print(f"  MH risk ratio lnc_same vs pc_same (bins 1-5): {rr[0]:.3f} ({rr[1]:.2f}-{rr[2]:.2f})")
    out_sum.append(dict(line=line, estimator=label, stat="MH_RR_lnc_same_vs_pc_same", value=rr[0], binom_lo=rr[1], binom_hi=rr[2],
                        clust_lo=np.nan, clust_hi=np.nan))


def orientation(G, line, bws):
    try:
        import pyBigWig
    except ImportError:
        print(f"  orientation check skipped ({line}): pyBigWig not available")
        return
    p, m = bws.split(",")
    bw = {"+": pyBigWig.open(p), "-": pyBigWig.open(m)}
    top = G[(G.kind == "standalone") & (G.type == "protein_coding")].nlargest(300, f"expr_{line}")
    r = []
    for x in top.itertuples():
        o = "-" if x.strand == "+" else "+"
        try:
            a = abs(bw[x.strand].stats(x.chrom, int(x.start), int(x.end), type="mean")[0] or 0)
            b = abs(bw[o].stats(x.chrom, int(x.start), int(x.end), type="mean")[0] or 0)
        except RuntimeError:
            continue
        r.append((a + 1e-3) / (b + 1e-3))
    med = float(np.median(r)) if r else float("nan")
    print(f"  orientation check ({line}): 300 most expressed stand-alone protein-coding genes, coverage on own strand / "
          f"opposite strand, median {med:.1f} -> {'OK' if med > 1 else 'STRANDS SWAPPED'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--line", action="append", required=True)
    ap.add_argument("--min-junc", type=int, default=2)
    ap.add_argument("--boot", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=46)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    G = pd.read_csv(a.genes, sep="\t", low_memory=False)
    G = G[G["type"].isin(["protein_coding", "lncRNA"]) & G.kind.isin(["standalone", "nested"])].copy()
    G["orient"] = G.orient.fillna("").astype(str)
    G = G[(G.kind == "standalone") | G.orient.isin(["same", "anti"])]
    pre = np.where(G["type"] == "protein_coding", "pc", "lnc")
    G["grp"] = pd.Series(pre, index=G.index) + "_" + np.where(G.kind == "standalone", "stand", G.orient)
    G["cluster"] = np.where(G.kind == "nested", "H:" + G.host.astype(str), "G:" + G.gene_id.astype(str))
    print(f"genes: {len(G)}; " + ", ".join(f"{g} {int((G.grp == g).sum())}" for g in GROUPS))

    own = read_gtf_introns(a.gtf, set(G.gene_id))
    G["n_own_introns"] = G.gene_id.map(lambda g: len(own.get(g, [])))
    needed = {k for v in own.values() for k in v}
    print(f"genes with >= 1 intron of their own: {int((G.n_own_introns > 0).sum())}; introns: {len(needed)}")

    counts, summ = [], []
    for spec in a.line:
        line, jf, bws = (spec.split("|") + ["", ""])[:3]
        print(f"\n===================== {line} =====================")
        if bws:
            orientation(G, line, bws)
        det = tf(G[f"det_{line}"])
        J = read_junctions(jf, needed)
        jmax = G.gene_id.map(lambda g: max([J.get(k, 0) for k in own.get(g, [])], default=0)).astype(float)
        multi = G.n_own_introns > 0
        ref = jmax[(G.kind == "standalone") & multi & (jmax >= a.min_junc)]
        edges = np.quantile(ref, [0.2, 0.4, 0.6, 0.8])
        jbin = pd.Series(np.where(jmax >= a.min_junc, np.searchsorted(edges, jmax, side="right") + 1, 0), index=G.index)
        print(f"  junction bins: edges {edges.tolist()} (reads on the best own junction); expressed stand-alone multi-exon genes "
              f"{len(ref)}; sizes " + ", ".join(f"{b}:{int(((jbin == b) & (G.kind == 'standalone') & multi).sum())}" for b in range(6)))
        cbin = pd.to_numeric(G[f"bin_{line}"], errors="coerce")
        clear = (pd.to_numeric(G[f"exon_{line}"], errors="coerce") >= 3 * pd.to_numeric(G[f"flank_{line}"], errors="coerce"))
        base = G[["gene_id", "type", "kind", "grp", "cluster"]].copy()
        base["det"] = det.values
        for label, b in (("cov", cbin), ("cov_multi", cbin.where(multi)), ("cov_clear", cbin.where(clear | (cbin == 0))),
                         ("junc", jbin.where(multi).astype(float))):
            d = base.copy()
            d["bin"] = b
            analyse(d, line, label, a.boot, rng, counts, summ)

        # expected de novo start sites near host-strand lncRNAs if each had its own start
        it = tf(G[f"itss_{line}"]) if f"itss_{line}" in G else None
        if it is not None and it.any():
            ns = G.kind.eq("nested") & G.orient.eq("same")
            p1 = it[ns & det].mean()
            p0 = it[ns & ~det].mean()
            sa = (G.grp == "lnc_stand")
            rate = det[sa].groupby(cbin[sa]).mean()
            ls = G.grp == "lnc_same"
            r = cbin[ls].map(rate).fillna(0)
            pcs = G.grp == "pc_same"
            pc_rate = det[G.grp == "pc_stand"].groupby(cbin[G.grp == "pc_stand"]).mean()
            oe_pc_same = det[pcs].sum() / cbin[pcs].map(pc_rate).fillna(0).sum()
            e1 = float((r * p1 + (1 - r) * p0).sum())
            r2 = (r * oe_pc_same).clip(upper=1)
            e2 = float((r2 * p1 + (1 - r2) * p0).sum())
            print(f"\n  de novo start sites within 500 nt of the annotated start, host-strand lncRNAs ({int(ls.sum())} genes): "
                  f"observed {int(it[ls].sum())}; P(de novo | cap test passed) = {p1:.3f}, P(de novo | not passed) = {p0:.3f} "
                  f"(host-strand nested genes); expected if each had its own start: {e1:.1f} (as stand-alone lncRNAs), "
                  f"{e2:.1f} (with the host-strand shortfall of protein-coding genes, O/E {oe_pc_same:.3f})")
            summ.append(dict(line=line, estimator="cov", stat="denovo_observed", value=int(it[ls].sum())))
            summ.append(dict(line=line, estimator="cov", stat="denovo_expected_standalone", value=e1))
            summ.append(dict(line=line, estimator="cov", stat="denovo_expected_pc_same_shortfall", value=e2))

    pd.DataFrame(counts).to_csv(f"{a.out_prefix}.counts.tsv", sep="\t", index=False)
    pd.DataFrame(summ).to_csv(f"{a.out_prefix}.summary.tsv", sep="\t", index=False, float_format="%.4f")
    print(f"\nwritten: {a.out_prefix}.counts.tsv, {a.out_prefix}.summary.tsv")


if __name__ == "__main__":
    main()

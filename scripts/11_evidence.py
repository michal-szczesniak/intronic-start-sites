#!/usr/bin/env python3
"""
iSTART, Stage 0, step 2: evidence for the start (PRO-cap) and the end (PolyASite) of the sites from step 1.

Input:
  --sites        sites.tsv from 10_candidates.py
  --host-tss     host_tss.tsv from 10_candidates.py
  --cap-plus/--cap-minus        bigWig of 5' ends of capped RNA (e.g. CoPRO capped, ENCSR261KBX)
  --uncap-plus/--uncap-minus    (optional) bigWig of uncapped RNA (e.g. ENCSR114HGS) - background control
  --pas          (optional) PolyASite 2.0 clusters in BED (chromosome names with or without "chr")
Output:
  --out-prefix.evidence.tsv     all sites with measures and flags
  --out-prefix.genes.tsv        candidates at the gene level (best start unit)
  --out-prefix.summary.txt      orientation check, calibration on controls, class counts, G1 criterion

Evidence classes (candidates): B = start + end; C = start only; '-' = no start.
Class A (start + end + long read) requires optional long reads - column long_read = NA.
"""
import argparse
import bisect
import gzip
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pyBigWig
from scipy.stats import poisson


def open_any(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


class Track:
    """Pair of bigWigs (plus/minus); values are taken as |x| (some pipelines store the minus strand as negative)."""

    def __init__(self, plus, minus, name):
        self.name = name
        self.bw = {"+": pyBigWig.open(plus), "-": pyBigWig.open(minus)}
        self.flip = False
        self.total = {k: abs(v.header().get("sumData", 0)) or 1.0 for k, v in self.bw.items()}

    def _bw(self, strand):
        if self.flip:
            strand = "-" if strand == "+" else "+"
        return self.bw[strand]

    def values(self, chrom, s, e, strand):
        bw = self._bw(strand)
        L = bw.chroms(chrom)
        if L is None:
            return np.zeros(max(0, e - s))
        s2, e2 = max(0, s), min(L, e)
        out = np.zeros(e - s)
        if e2 > s2:
            v = bw.values(chrom, s2, e2, numpy=True)
            out[s2 - s:e2 - s] = np.abs(np.nan_to_num(v))
        return out

    def depth(self):
        return self.total["+"] + self.total["-"]


def orientation_check(track, ref_tss, w=50):
    """Signal summed at host reference start sites: same strand vs opposite strand; meta-profile."""
    same = opp = 0.0
    prof = np.zeros(2 * w + 1)
    for r in ref_tss.itertuples():
        v = track.values(r.chrom, r.pos - w, r.pos + w + 1, r.strand)
        o = track.values(r.chrom, r.pos - w, r.pos + w + 1, "-" if r.strand == "+" else "+")
        same += v.sum()
        opp += o.sum()
        prof += v if r.strand == "+" else v[::-1]
    offset = int(np.argmax(prof)) - w if prof.sum() > 0 else 0
    return same, opp, offset


def load_pas(path):
    idx = defaultdict(list)
    with open_any(path) as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            chrom = f[0]
            if not chrom.startswith("chr"):
                chrom = "chrM" if chrom in ("MT", "M") else "chr" + chrom
            idx[(chrom, f[5])].append((int(f[1]), int(f[2])))
    out = {}
    for k, v in idx.items():
        v.sort()
        out[k] = (v, [s for s, _ in v], max(e - s for s, e in v))
    print(f"[PAS] clusters: {sum(len(v[0]) for v in out.values())}", file=sys.stderr)
    return out


def pas_hit(pas, chrom, strand, pos, w):
    """Nearest PAS cluster within pos±w on the same strand; returns its distance or None."""
    if (chrom, strand) not in pas:
        return None
    ivs, starts, maxlen = pas[(chrom, strand)]
    j0 = bisect.bisect_left(starts, pos - w - maxlen)
    j1 = bisect.bisect_right(starts, pos + w)
    best = None
    for s, e in ivs[j0:j1]:
        if e > pos - w and s <= pos + w:
            d = 0 if s <= pos < e else min(abs(pos - s), abs(pos - (e - 1)))
            best = d if best is None else min(best, d)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--host-tss", required=True)
    ap.add_argument("--cap-plus", required=True)
    ap.add_argument("--cap-minus", required=True)
    ap.add_argument("--uncap-plus")
    ap.add_argument("--uncap-minus")
    ap.add_argument("--pas")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--window", type=int, default=100, help="window around the start site (± nt)")
    ap.add_argument("--flank-in", type=int, default=500)
    ap.add_argument("--flank-out", type=int, default=2000)
    ap.add_argument("--bg-floor", type=float, default=1e-3, help="background floor (reads/nt)")
    ap.add_argument("--min-reads", type=float, default=5)
    ap.add_argument("--pval", type=float, default=1e-3)
    ap.add_argument("--min-ratio", type=float, default=2.0, help="min. capped/uncapped ratio (after normalisation)")
    ap.add_argument("--min-exon-dist", type=int, default=1000)
    ap.add_argument("--pas-window", type=int, default=200)
    ap.add_argument("--g1-min", type=int, default=500, help="G1 threshold: number of class B genes")
    ap.add_argument("--no-auto-flip", action="store_true")
    a = ap.parse_args()

    sites = pd.read_csv(a.sites, sep="\t", dtype={"intron_no": str, "intron_start": str, "intron_end": str, "tes_list": str})
    htss = pd.read_csv(a.host_tss, sep="\t")
    pas = load_pas(a.pas) if a.pas else None
    cap = Track(a.cap_plus, a.cap_minus, "capped")
    uncap = Track(a.uncap_plus, a.uncap_minus, "uncapped") if a.uncap_plus and a.uncap_minus else None
    log = []

    # --- 1. strand orientation and peak offset at known start sites ---
    ref = htss[htss.is_ref == 1].drop_duplicates(["chrom", "strand", "pos"])
    same, opp, offset = orientation_check(cap, ref)
    log.append(f"Orientation (capped, ±50 nt around {len(ref)} MANE start sites): same strand = {same:.0f}, "
               f"opposite = {opp:.0f}; meta-profile peak shifted by {offset} nt")
    if opp > same and not a.no_auto_flip:
        cap.flip = True
        if uncap:
            uncap.flip = True
        log.append("WARNING: more signal on the opposite strand - swapping plus/minus (check the ENCODE file descriptions)")
    if abs(offset) > 10:
        log.append("WARNING: meta-profile peak > 10 nt from the annotated start site - check that the bigWig holds 5' ends")

    # --- 2. host start site used in K562 ---
    W = a.window
    hu = {}
    for g, df in htss.drop_duplicates(["host_gene_id", "pos"]).groupby("host_gene_id"):
        best_pos, best_sig = None, -1.0
        for r in df.itertuples():
            sig = cap.values(r.chrom, r.pos - W, r.pos + W + 1, r.strand).sum()
            if sig > best_sig:
                best_pos, best_sig = r.pos, sig
        hu[g] = (best_pos, best_sig)

    # --- 3. start-site measures for each site ---
    rows = []
    cache = {}
    for r in sites.itertuples(index=False):
        key = (r.chrom, r.strand, r.pos)
        if key not in cache:
            v = cap.values(r.chrom, r.pos - W, r.pos + W + 1, r.strand)
            csum = float(v.sum())
            peak = int(r.pos - W + np.argmax(v)) if csum > 0 else r.pos
            fl = np.concatenate([
                cap.values(r.chrom, r.pos - a.flank_out, r.pos - a.flank_in, r.strand),
                cap.values(r.chrom, r.pos + a.flank_in, r.pos + a.flank_out, r.strand)])
            bg = max(float(fl.mean()) if fl.size else 0.0, a.bg_floor)
            lam = bg * (2 * W + 1)
            p = float(poisson.sf(csum - 1, lam)) if csum > 0 else 1.0
            if uncap:
                usum = float(uncap.values(r.chrom, r.pos - W, r.pos + W + 1, r.strand).sum())
                ratio = ((csum + 1) / cap.depth()) / ((usum + 1) / uncap.depth())
            else:
                usum, ratio = np.nan, np.nan
            cache[key] = (csum, peak, bg, p, usum, ratio)
        csum, peak, bg, p, usum, ratio = cache[key]
        hpos, hsig = hu.get(r.host_gene_id, (None, 0.0))
        d_used = (r.pos - hpos if r.strand == "+" else hpos - r.pos) if hpos is not None else np.nan
        exon_ok = r.min_dist_host_exon >= a.min_exon_dist
        tss_ok = (csum >= a.min_reads) and (p < a.pval) and exon_ok and \
                 (np.isnan(ratio) or ratio >= a.min_ratio)
        pas_d = None
        if pas and r.site_type == "candidate" and isinstance(r.tes_list, str) and r.tes_list != "NA":
            for t in r.tes_list.split(","):
                d = pas_hit(pas, r.chrom, r.strand, int(t), a.pas_window)
                if d is not None:
                    pas_d = d if pas_d is None else min(pas_d, d)
        rows.append(dict(
            **r._asdict(), cap_sum=csum, cap_peak_pos=peak, cap_peak_shift=peak - r.pos if r.strand == "+" else r.pos - peak,
            bg_per_nt=bg, pval=p, uncap_sum=usum, cap_uncap_ratio=ratio, exon_dist_ok=exon_ok,
            tss_supported=tss_ok, host_tss_used=hpos, host_tss_signal=hsig,
            host_active=hsig >= a.min_reads, dist_used_tss=d_used,
            upstream_of_used_host_tss=(d_used < 0) if not np.isnan(d_used) else False,
            pas_dist=pas_d if pas_d is not None else np.nan,
            pas_supported=pas_d is not None, long_read="NA"))
    ev = pd.DataFrame(rows)

    # BH correction separately for each site type (for information)
    ev["qval"] = np.nan
    for t, idx in ev.groupby("site_type").groups.items():
        pv = ev.loc[idx, "pval"].values
        order = np.argsort(pv)
        ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
        q = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty_like(q)
        out[order] = np.minimum(q, 1.0)
        ev.loc[idx, "qval"] = out

    def cls(r):
        if r.site_type != "candidate" or not r.tss_supported:
            return "-"
        return "B" if r.pas_supported else "C"

    ev["evidence_class"] = ev.apply(cls, axis=1)
    ev.to_csv(f"{a.out_prefix}.evidence.tsv", sep="\t", index=False)

    # --- 4. gene level: best start unit ---
    order = {"B": 0, "C": 1, "-": 2}
    cand = ev[ev.site_type == "candidate"].copy()
    cand["_o"] = cand.evidence_class.map(order)
    genes = cand.sort_values(["cand_gene_id", "_o", "cap_sum"], ascending=[True, True, False]) \
                .drop_duplicates("cand_gene_id").drop(columns="_o")
    genes.to_csv(f"{a.out_prefix}.genes.tsv", sep="\t", index=False)

    # --- 5. summary and calibration ---
    def frac(df):
        return f"{df.tss_supported.mean()*100:.1f}% ({int(df.tss_supported.sum())}/{len(df)})" if len(df) else "n/a"

    act = ev[ev.host_active]
    log.append("")
    log.append("Calibration of start-site calls (active hosts only):")
    log.append(f"  internal promoters (expected high, but not all of them are used in K562): "
               f"{frac(act[act.site_type == 'internal_promoter'])}")
    log.append(f"  random controls (≈ false-positive rate, expected low): "
               f"{frac(act[act.site_type == 'random_control'])}")
    log.append(f"  candidates (start units): {frac(act[act.site_type == 'candidate'])}")
    g = genes
    log.append("")
    log.append(f"Candidate genes: {len(g)}; with an active host: {int(g.host_active.sum())}")
    for c in ("B", "C", "-"):
        log.append(f"  class {c}: {int((g.evidence_class == c).sum())} "
                   f"(with an active host: {int(((g.evidence_class == c) & g.host_active).sum())})")
    nB = int(((g.evidence_class == "B") & g.host_active).sum())
    log.append("")
    log.append(f"G1 (number of iRNAs): {nB} class B genes with an active host; threshold {a.g1_min} → "
               f"{'MET' if nB >= a.g1_min else 'NOT MET - consider class C'}")
    up = int((g.upstream_of_used_host_tss & (g.evidence_class != '-')).sum())
    log.append(f"Note: {up} candidate genes with a start site lie upstream of the host start site used in K562 "
               f"(excluded from module A)")
    txt = "\n".join(log)
    open(f"{a.out_prefix}.summary.txt", "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()

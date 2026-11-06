#!/usr/bin/env python3
"""
iSTART, Stage 0, step 1b: intronic start sites (iTSS) detected directly from PRO-cap peaks.

Input:
  --gtf          GENCODE (the same as in step 10)
  --uni          unidirectional PRO-cap peaks (PINTS, ENCODE: ENCFF636REV.bed.gz)
  --bi           bidirectional PRO-cap peaks (PINTS, ENCODE: ENCFF819CPA.bed.gz)
  --cap-plus/--cap-minus   capped bigWigs (to check the orientation of bidirectional peaks)
  --pas          PolyASite 2.0 clusters
  --annot-sites  sites.tsv from step 10 (annotated candidates and internal promoters)
Output:
  --out          sites file in the step-10 format (input for 11_evidence.py), containing:
                 iTSS (site_type=candidate, source=procap_uni/procap_bi/procap_both),
                 internal promoters from step 10 and new random controls matched to the iTSS.

iTSS = summit of a PRO-cap peak on the host strand, in an intron of the reference transcript,
at least --min-exon-dist from host exons. For each iTSS:
  - end of the unit: first PolyASite cluster downstream of the start site in the same intron
    (at least --pas-min-dist from the start site and from the intron end) -> column tes_list,
  - link to the annotation: step-10 candidates within --link-window nt,
  - type: unidirectional peak, bidirectional peak (divergent, enhancer/promoter-like) or both.
"""
import argparse
import bisect
import gzip
import importlib.util
import os
import random
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def load(fname):
    spec = importlib.util.spec_from_file_location("m_" + fname.split("_")[0], os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


c10 = load("10_candidates.py")
c11 = load("11_evidence.py")


def open_any(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_uni(path):
    """PINTS unidirectional peaks: chr, start, end, name, p, strand, counts, padj, ..., summit (col. 15)."""
    out = []
    with open_any(path) as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            s, e = int(f[1]), int(f[2])
            summit = int(f[14]) if len(f) > 14 and f[14].lstrip("-").isdigit() else (s + e) // 2
            if not (s <= summit < e):
                summit = (s + e) // 2
            cnt = float(f[6]) if len(f) > 6 else np.nan
            out.append((f[0], f[5], summit, "uni", cnt))
    return out


def read_bi(path):
    """PINTS bidirectional peaks: chr, start, end, confidence, start on strand A (col. 5), start on strand B (col. 6)."""
    out = []
    with open_any(path) as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[4]), int(f[5]), f[3]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--uni", required=True)
    ap.add_argument("--bi", required=True)
    ap.add_argument("--cap-plus", required=True)
    ap.add_argument("--cap-minus", required=True)
    ap.add_argument("--pas", required=True)
    ap.add_argument("--annot-sites", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chrom")
    ap.add_argument("--host-gene-types", default="protein_coding")
    ap.add_argument("--host-exon-tx-types", default=c10.HOST_EXON_TX_TYPES_DEFAULT)
    ap.add_argument("--min-exon-dist", type=int, default=1000)
    ap.add_argument("--merge", type=int, default=50, help="distance for merging summits on the same strand (nt)")
    ap.add_argument("--link-window", type=int, default=500)
    ap.add_argument("--pas-min-dist", type=int, default=200)
    ap.add_argument("--ctrl-per-site", type=int, default=3)
    ap.add_argument("--ctrl-exclude", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261004)
    a = ap.parse_args()

    # --- hosts ---
    tx, exons = c10.parse_gtf(a.gtf)
    if a.chrom:
        tx = {k: v for k, v in tx.items() if v["chrom"] == a.chrom}
    hosts = c10.build_hosts(tx, exons, set(a.host_gene_types.split(",")), set(a.host_exon_tx_types.split(",")))
    hidx = c10.index_hosts(hosts)
    pc_tss = defaultdict(list)
    for t in tx.values():
        if t["gene_type"] == "protein_coding":
            pc_tss[(t["chrom"], t["strand"])].append(c10.tss_of(t))
    for k in pc_tss:
        pc_tss[k].sort()

    # --- peaks ---
    uni = read_uni(a.uni)
    bi = read_bi(a.bi)
    if a.chrom:
        uni = [x for x in uni if x[0] == a.chrom]
        bi = [x for x in bi if x[0] == a.chrom]
    cap = c11.Track(a.cap_plus, a.cap_minus, "capped")

    def sig(chrom, pos, strand, w=2):
        return cap.values(chrom, pos - w, pos + w + 1, strand).sum()

    # orientation of columns 5/6 in bidirectional peaks (A = plus?) - checked against the signal
    sA = sum(sig(c, p5, "+") + sig(c, p6, "-") for c, p5, p6, _ in bi[:5000])
    sB = sum(sig(c, p5, "-") + sig(c, p6, "+") for c, p5, p6, _ in bi[:5000])
    swap = sB > sA
    print(f"[BI] pairs: {len(bi)}; signal assuming col5=plus: {sA:.0f}, reversed: {sB:.0f} -> "
          f"{'col5 = minus' if swap else 'col5 = plus'}", file=sys.stderr)
    starts = list(uni)
    for c, p5, p6, conf in bi:
        if swap:
            p5, p6 = p6, p5
        starts.append((c, "+", p5, "bi", sig(c, p5, "+", 0)))
        starts.append((c, "-", p6, "bi", sig(c, p6, "-", 0)))
    print(f"[PEAKS] start sites (unidirectional: {len(uni)}, from bidirectional pairs: {2 * len(bi)})", file=sys.stderr)

    # --- merging summits on the same strand ---
    by = defaultdict(list)
    for c, s, p, kind, cnt in starts:
        by[(c, s)].append((p, kind, cnt))
    merged = []
    all_peak_pos = defaultdict(list)
    for (c, s), lst in by.items():
        lst.sort()
        cur = []
        for item in lst:
            if cur and item[0] - cur[-1][0] <= a.merge:
                cur.append(item)
            else:
                if cur:
                    merged.append((c, s, cur))
                cur = [item]
        if cur:
            merged.append((c, s, cur))
        all_peak_pos[(c, s)] = sorted(p for p, _, _ in lst)

    # --- start site in a host intron, on the host strand ---
    pas = c11.load_pas(a.pas)
    annot = pd.read_csv(a.annot_sites, sep="\t", dtype=str)
    cand_annot = annot[annot.site_type == "candidate"]
    ann_idx = defaultdict(list)
    for r in cand_annot.itertuples():
        ann_idx[(r.chrom, r.strand)].append((int(r.pos), r.cand_gene_id, r.cand_gene_name))
    for k in ann_idx:
        ann_idx[k].sort()

    sites = []
    n_in = 0
    for c, s, grp in merged:
        best = max(grp, key=lambda x: (sig(c, x[0], s), -abs(x[0] - grp[0][0])))
        pos = best[0]
        kinds = {g[1] for g in grp}
        hs = c10.containing_hosts(hidx, c, s, pos, pos + 1)
        hs = [h for h in hs if c10.ref_intron_of(h, pos, pos + 1) and h["exon_set"].dist(pos) >= a.min_exon_dist]
        if not hs:
            continue
        h = min(hs, key=lambda x: x["end"] - x["start"])
        n_in += 1
        _, ia, ib = c10.ref_intron_of(h, pos, pos + 1)
        # first PAS cluster downstream of the start site in the same intron
        tes, unit_len = "NA", np.nan
        key = (c, s)
        if key in pas:
            ivs, st, _ = pas[key]
            if s == "+":
                lo, hi = pos + a.pas_min_dist, ib - a.pas_min_dist
                j = bisect.bisect_left(st, lo)
                if j < len(ivs) and ivs[j][1] <= hi:
                    tes = (ivs[j][0] + ivs[j][1] - 1) // 2
            else:
                lo, hi = ia + a.pas_min_dist, pos - a.pas_min_dist
                j = bisect.bisect_right(st, hi) - 1
                while j >= 0 and ivs[j][1] > hi:
                    j -= 1
                if j >= 0 and ivs[j][0] >= lo:
                    tes = (ivs[j][0] + ivs[j][1] - 1) // 2
            if tes != "NA":
                unit_len = abs(tes - pos)
        # link to the annotation
        al = ann_idx.get(key, [])
        j0 = bisect.bisect_left(al, (pos - a.link_window,))
        links = [x for x in al[j0:] if x[0] <= pos + a.link_window]
        src = "procap_both" if kinds == {"uni", "bi"} else f"procap_{kinds.pop()}"
        sid = f"iTSS:{c}:{pos}:{s}"
        rec = c10.site_record(sid, "candidate", h, pos, dict(
            cand_gene_id=sid, cand_gene_name=",".join(x[2] for x in links) or "NA", source=src,
            cand_tx_ids=",".join(x[1] for x in links) or "NA", tes_list=str(tes), context="ref_intron"))
        pl = pc_tss.get(key, [])
        jj = bisect.bisect_left(pl, pos - 500)
        at_pc = jj < len(pl) and pl[jj] <= pos + 500
        rec.update(peak_kind=src, peak_count=best[2], unit_len=unit_len, linked_annot=len(links) > 0,
                   at_other_pc_tss=at_pc)
        sites.append(rec)
    print(f"[iTSS] start sites in host introns (sense strand, ≥{a.min_exon_dist} nt from exons): {n_in}; "
          f"with a PAS downstream in the same intron: {sum(r['tes_list'] != 'NA' for r in sites)}; "
          f"linked to the annotation: {sum(r['linked_annot'] for r in sites)}", file=sys.stderr)

    print(f"[iTSS] of these, at the start of another coding gene (±500 nt): {sum(r['at_other_pc_tss'] for r in sites)}",
          file=sys.stderr)

    # --- random controls matched to the iTSS (same host) ---
    rng = random.Random(a.seed)
    excl = defaultdict(list)
    for k, v in all_peak_pos.items():
        excl[k] += v
    for k, v in ann_idx.items():
        excl[k] += [x[0] for x in v]
    for r in annot[annot.site_type == "internal_promoter"].itertuples():
        excl[(r.chrom, r.strand)].append(int(r.pos))
    for k in excl:
        excl[k].sort()
    cache = {}

    def allowed(h):
        if h["gene_id"] in cache:
            return cache[h["gene_id"]]
        kt = excl[(h["chrom"], h["strand"])]
        ivs = []
        for s0, e0 in h["introns"]:
            s2, e2 = s0 + a.min_exon_dist, e0 - a.min_exon_dist
            if e2 <= s2:
                continue
            cut = [(max(s2, p - a.ctrl_exclude), min(e2, p + a.ctrl_exclude + 1))
                   for p in kt[bisect.bisect_left(kt, s2 - a.ctrl_exclude):bisect.bisect_right(kt, e2 + a.ctrl_exclude)]]
            cur = s2
            for cs, ce in c10.merge_intervals([x for x in cut if x[1] > x[0]]):
                if cs > cur:
                    ivs.append((cur, cs))
                cur = max(cur, ce)
            if cur < e2:
                ivs.append((cur, e2))
        cache[h["gene_id"]] = ivs
        return ivs

    ctrls = []
    for srec in sites:
        h = hosts[srec["host_gene_id"]]
        ivs = allowed(h)
        tot = sum(e - b for b, e in ivs)
        if tot <= 0:
            continue
        for i in range(a.ctrl_per_site):
            for _ in range(20):
                x = rng.randrange(tot)
                for b, e in ivs:
                    if x < e - b:
                        p = b + x
                        break
                    x -= e - b
                if h["exon_set"].dist(p) >= a.min_exon_dist:
                    break
            else:
                continue
            rec = c10.site_record(f"{srec['site_id']}:ctrl{i+1}", "random_control", h, p, dict(
                cand_gene_id=srec["cand_gene_id"], cand_gene_name="NA", source="random",
                cand_tx_ids="NA", tes_list="NA", context="ref_intron"))
            rec.update(peak_kind="NA", peak_count=np.nan, unit_len=np.nan, linked_annot=False, at_other_pc_tss=False)
            ctrls.append(rec)
    print(f"[CTRL-] random controls: {len(ctrls)}", file=sys.stderr)

    # --- output: iTSS + controls + internal promoters from step 10 ---
    df = pd.DataFrame(sites + ctrls)
    ip = annot[annot.site_type == "internal_promoter"].copy()
    for col in ("peak_kind", "peak_count", "unit_len", "linked_annot", "at_other_pc_tss"):
        ip[col] = "NA" if col == "peak_kind" else (False if col in ("linked_annot", "at_other_pc_tss") else np.nan)
    if a.chrom:
        ip = ip[ip.chrom == a.chrom]
    out = pd.concat([df, ip[df.columns]], ignore_index=True)
    out["pos"] = out["pos"].astype(int)
    out = out.sort_values(["chrom", "pos"])
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)
    print(f"[OUT] {a.out}: iTSS {len(sites)}, controls {len(ctrls)}, internal promoters {len(ip)}", file=sys.stderr)


if __name__ == "__main__":
    main()
